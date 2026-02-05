module fsps_ssp
    !> @brief Core module for Simple Stellar Population (SSP) generation.
    !>
    !> @details
    !> This module calculates the time-evolution of a Single Stellar Population (SSP),
    !> which is the fundamental building block of stellar population synthesis.
    !> It combines:
    !> 1.  **Isochrone Data:** Stellar evolution tracks (Mass, L, T, g vs Age).
    !> 2.  **Spectral Libraries:** Flux distributions for individual stars.
    !> 3.  **IMF Weights:** Initial Mass Function weighting.
    !> 4.  **Stellar Physics:** Horizontal Branch, Blue Stragglers, AGB modifications.
    !>
    !> **Key Architecture:**
    !> This module replaces the legacy `SSP_GEN` subroutine. It uses a 
    !> "Recyclable Buffer" pattern (`isochrone_buffer_t`) to process one time-step 
    !> at a time using contiguous 1D arrays. This significantly improves L1/L2 
    !> cache locality compared to the legacy strided access on `(NT, NM)` global arrays.
    !>
    !> **Optimization Strategy:**
    !> - **Contiguous Memory:** Data for the current age is gathered into a compact buffer.
    !> - **Cache Blocking:** Spectral summation uses a block-loop (Batch Size 32) 
    !>   to keep the accumulator in L1 cache, maximizing vectorization throughput.
    !> - **Vectorization:** Critical loops are marked `contiguous` to enable SIMD.

    use fsps_precision, only: WP
    use fsps_constants, only: NM, M_SOL, L_SOL, VERBOSE, BHB_SBS_TIME, TIME_RES_INCR, SAFE_FLOOR
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t
    
    !> Physics Modules
    use fsps_imf, only: compute_imf_weights
    use fsps_stellar_modifications, only: modify_horizontal_branch, &
                                          apply_blue_stragglers, &
                                          modify_giant_branch, &
                                          add_remnant_mass, &
                                          add_xray_binaries
    use fsps_gas, only: apply_nebular_emission
    use fsps_smoothing, only: apply_smoothing
    use fsps_spectral_library, only: get_stellar_spectrum

    implicit none
    private

    public :: generate_ssp_grid
    public :: isochrone_buffer_t
    public :: init_isochrone_buffer, reset_buffer, free_isochrone_buffer
    public :: load_timestep_data
    public :: apply_isochrone_physics
    public :: compute_integrated_properties, accumulate_spectrum
    public :: interpolate_time_grid
    public :: configure_imf_parameters

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------
    real(WP), parameter :: MIN_WEIGHT_CUTOFF = 1.0e-8_wp  ! Skip stars with negligible weight
    real(WP), parameter :: SMOOTHING_SIGMA_FLAG = 99.0_wp ! Flag to use variable sigma from LSF
    real(WP), parameter :: LN10 = log(10.0_wp)            ! Natural log of 10
    integer, parameter :: BATCH_SIZE = 32                 ! Process stars in chunks for L1 cache locality

    ! ------------------------------------------------------------------------
    ! DERIVED TYPES
    ! ------------------------------------------------------------------------
    
    !> @brief A recyclable buffer for processing a single isochrone time-step.
    !> @details
    !> Holds all properties for the stars alive at a specific age. 
    !> These are 1D CONTIGUOUS arrays of size NM (max stars).
    !> We allocate this ONCE and reuse it to avoid heap fragmentation.
    type :: isochrone_buffer_t
        integer :: n_stars
        
        ! Changed to POINTER to allow them to be targets for 2D views.
        ! We manually manage their memory in init/free routines.
        real(WP), pointer :: initial_mass(:) => null()
        real(WP), pointer :: current_mass(:) => null()
        real(WP), pointer :: log_lum(:)      => null()
        real(WP), pointer :: log_teff(:)     => null()
        real(WP), pointer :: log_g(:)        => null()
        real(WP), pointer :: phase(:)        => null()
        real(WP), pointer :: co_ratio(:)     => null()
        real(WP), pointer :: log_mdot(:)     => null()
        real(WP), pointer :: weights(:)      => null()
    end type isochrone_buffer_t

contains

    ! ------------------------------------------------------------------------
    ! MAIN DRIVER
    ! ------------------------------------------------------------------------

    !> @brief Main driver routine to generate the SSP grid.
    !>
    !> @details
    !> Orchestrates the extraction of isochrones, application of stellar physics,
    !> IMF weighting, spectral integration, and post-processing.
    !>
    !> **Pipeline:**
    !> 1. **Setup:** Initialize buffers and handle special cases (BPASS).
    !> 2. **Evolution Loop:** For each time step:
    !>    - Load isochrone data into contiguous buffer.
    !>    - Apply stellar physics (HB, BS, AGB).
    !>    - Sum integrated properties (Mass, Lbol).
    !>    - Accumulate Spectra.
    !> 3. **Post-Processing:**
    !>    - Interpolate time grid (if high-res).
    !>    - Add Nebular Emission.
    !>    - Add X-Ray Binaries.
    !>    - Apply Instrumental Smoothing (LSF).
    !>
    !> @param[in,out] ctx       The main FSPS context.
    !> @param[in]     pset      User parameter set.
    !> @param[out]    mass_grid Output mass history (time dependent).
    !> @param[out]    lbol_grid Output bolometric luminosity history.
    !> @param[out]    spec_grid Output spectral grid (wavelength, time).
    subroutine generate_ssp_grid(ctx, pset, mass_grid, lbol_grid, spec_grid)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        real(WP), intent(out), contiguous   :: mass_grid(:)
        real(WP), intent(out), contiguous   :: lbol_grid(:)
        real(WP), intent(out), contiguous   :: spec_grid(:,:)
        real(WP), allocatable               :: temp_spec_grid(:,:)

        ! Local variables
        type(isochrone_buffer_t) :: buf
        integer :: n_times, i_time, out_idx
        real(WP) :: time_log_yr

        ! Workspace for spectral generation (allocated on heap to prevent stack overflow)
        real(WP), allocatable :: work_spec(:)
        
        ! --------------------------------------------------------------------
        ! 1. INITIALIZATION
        ! --------------------------------------------------------------------
        
        ! Initialize outputs to zero
        mass_grid = 0.0_wp
        lbol_grid = 0.0_wp
        spec_grid = 0.0_wp

        ! Validation: Check for setup flag
        if (ctx%state%check_sps_setup == 0) then
            write(*,*) 'SSP_GEN ERROR: SPS_SETUP must be run before SSP_GEN.'
            stop
        end if

        ! Validation: Check metallicity index
        if (pset%zmet < 1 .or. pset%zmet > ctx%state%nz) then
            write(*,*) 'SSP_GEN ERROR: Metallicity index out of range:', pset%zmet
            stop
        end if

        ! --------------------------------------------------------------------
        ! 2. SPECIAL CASE: BPASS MODELS
        ! --------------------------------------------------------------------
        ! BPASS models are pre-computed and stored in bpass_* arrays.
        ! We do not synthesize them star-by-star.
        if (trim(ctx%state%isoc_type) == 'bpss') then
            
            ! Copy pre-computed data for the requested metallicity
            spec_grid = ctx%state%bpass_spec_ssp(:, :, pset%zmet)
            mass_grid = ctx%state%bpass_mass_ssp(:, pset%zmet)
            
            ! BPASS usually doesn't provide separate Lbol history in the same way,
            ! or it's handled differently, but we leave it 0 or computed elsewhere.
            ! (Legacy code implies simple copy for spec/mass)
            
            return ! Exit immediately
        end if

        ! --------------------------------------------------------------------
        ! 3. SETUP & CONFIGURATION
        ! --------------------------------------------------------------------

        ! Configure IMF parameters (maps pset -> context global vars)
        ! Also handles legacy file reading for user-defined IMFs.
        call configure_imf_parameters(ctx, pset)

        ! Initialize the reusable memory buffer
        call init_isochrone_buffer(buf)

        ! Allocation of workspace (Size of wavelength grid)
        allocate(work_spec(size(spec_grid, 1)))

        ! --------------------------------------------------------------------
        ! 4. EVOLUTION LOOP
        ! --------------------------------------------------------------------
        n_times = ctx%state%nt

        do i_time = 1, n_times
            
            ! Optimization: Skip ages not requested by the user
            if (pset%ssp_gen_age(i_time) == 0) cycle

            ! Map internal time index to output grid index
            ! (Supports legacy spacing logic via TIME_RES_INCR)
            out_idx = 1 + (i_time - 1) * TIME_RES_INCR

            ! Get current Age (log years)
            time_log_yr = ctx%state%timestep_isoc(pset%zmet, i_time)

            ! A. LOAD DATA
            !    Copy strided global data into our contiguous buffer
            call load_timestep_data(ctx, pset%zmet, i_time, buf)

            ! B. APPLY PHYSICS
            !    IMF, Horizontal Branch, Blue Stragglers, Giant Branch modifications
            call apply_isochrone_physics(ctx, pset, time_log_yr, buf)

            ! C. COMPUTE INTEGRATED PROPERTIES
            !    Mass and Bolometric Luminosity
            call compute_integrated_properties(ctx, buf, &
                                               mass_grid(out_idx), &
                                               lbol_grid(out_idx))

            ! D. ACCUMULATE SPECTRA (Pass workspace)
            !    Sum individual stellar spectra into the grid
            call accumulate_spectrum(ctx, pset, buf, spec_grid(:, out_idx), work_spec)

        end do

        ! Release buffer memory
        call free_isochrone_buffer(buf)
        if (allocated(work_spec)) deallocate(work_spec)

        ! --------------------------------------------------------------------
        ! 5. POST-PROCESSING
        ! --------------------------------------------------------------------

        ! A. TIME INTERPOLATION
        !    If the calculation was done on a coarse grid (TIME_RES_INCR > 1),
        !    interpolate to fill the gaps.
        if (TIME_RES_INCR > 1) then
            call interpolate_time_grid(ctx, pset%zmet, mass_grid, lbol_grid, spec_grid)
        end if

        ! B. NEBULAR EMISSION
        if (ctx%add_neb_emission_val == 2) then
            ! Allocate a temp grid to avoid aliasing arguments
            allocate(temp_spec_grid(size(spec_grid,1), size(spec_grid,2)))

            ! Pass spec_grid as input, temp_grid as output
            call apply_nebular_emission(ctx, pset, spec_grid, temp_spec_grid)

            ! Copy back
            spec_grid = temp_spec_grid
            deallocate(temp_spec_grid)
        end if

        ! C. X-RAY BINARIES
        if (ctx%add_xrb_emission_val == 1) then
            ! Note: Calls into fsps_stellar_modifications
            ! Use temp buffer to prevent Aliasing (In/Out same array)
            allocate(temp_spec_grid(size(spec_grid,1), size(spec_grid,2)))
            
            call add_xray_binaries(ctx, pset, spec_grid, temp_spec_grid)
            spec_grid = temp_spec_grid
            
            deallocate(temp_spec_grid)
        end if

        ! D. INSTRUMENTAL SMOOTHING
        if (ctx%smooth_lsf_val == 1) then
             call apply_lsf_smoothing_grid(ctx, spec_grid)
        end if

    end subroutine generate_ssp_grid

    ! ------------------------------------------------------------------------
    ! BUFFER MANAGEMENT
    ! ------------------------------------------------------------------------

    !> @brief Allocates the recyclable isochrone buffer.
    subroutine init_isochrone_buffer(buf)
        type(isochrone_buffer_t), intent(out) :: buf
        
        ! Allocate to maximum size NM defined in fsps_constants.
        ! We do this once per SSP generation run.
        allocate(buf%initial_mass(NM))
        allocate(buf%current_mass(NM))
        allocate(buf%log_lum(NM))
        allocate(buf%log_teff(NM))
        allocate(buf%log_g(NM))
        allocate(buf%phase(NM))
        allocate(buf%co_ratio(NM))
        allocate(buf%log_mdot(NM))
        allocate(buf%weights(NM))
        
        call reset_buffer(buf)
    end subroutine init_isochrone_buffer

    !> @brief Resets scalars and zeroes arrays (soft clear).
    subroutine reset_buffer(buf)
        type(isochrone_buffer_t), intent(inout) :: buf
        buf%n_stars = 0
        buf%weights = 0.0_wp ! Essential to zero this as it's often accumulated
    end subroutine reset_buffer

    !> @brief Deallocates the isochrone buffer.
    subroutine free_isochrone_buffer(buf)
        type(isochrone_buffer_t), intent(inout) :: buf
        
        ! Pointers use 'associated', not 'allocated'
        if (associated(buf%initial_mass)) deallocate(buf%initial_mass)
        if (associated(buf%current_mass)) deallocate(buf%current_mass)
        if (associated(buf%log_lum))      deallocate(buf%log_lum)
        if (associated(buf%log_teff))     deallocate(buf%log_teff)
        if (associated(buf%log_g))        deallocate(buf%log_g)
        if (associated(buf%phase))        deallocate(buf%phase)
        if (associated(buf%co_ratio))     deallocate(buf%co_ratio)
        if (associated(buf%log_mdot))     deallocate(buf%log_mdot)
        if (associated(buf%weights))      deallocate(buf%weights)
        
        ! Good practice to nullify after deallocation
        buf%initial_mass => null()
        buf%current_mass => null()
        buf%log_lum      => null()
        buf%log_teff     => null()
        buf%log_g        => null()
        buf%phase        => null()
        buf%co_ratio     => null()
        buf%log_mdot     => null()
        buf%weights      => null()
    end subroutine free_isochrone_buffer


    !> @brief Copies strided global data into the contiguous local buffer.
    !> @details 
    !> Optimizes memory access by gathering data from the strided `ctx` structure
    !> into the hot L1/L2 cache-friendly `buf` once per timestep.
    subroutine load_timestep_data(ctx, z_idx, t_idx, buf)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)              :: z_idx, t_idx
        type(isochrone_buffer_t), intent(inout) :: buf
        
        integer :: n
        
        ! Get number of mass points, clamped to max NM
        ! Inside load_timestep_data
        if (ctx%state%nmass_isoc(z_idx, t_idx) > NM) then
            if (VERBOSE == 1) &
                write(*,*) 'SSP_GEN WARNING: Isochrone truncated at NM stars. Increase NM.'
        endif
        n = min(ctx%state%nmass_isoc(z_idx, t_idx), NM)
        buf%n_stars = n

        ! Vectorized Copy: Global(z, t, 1:n) -> Buffer(1:n)
        ! This pays the "stride tax" exactly once per timestep.
        buf%initial_mass(1:n) = ctx%state%mini_isoc(z_idx, t_idx, 1:n)
        buf%current_mass(1:n) = ctx%state%mact_isoc(z_idx, t_idx, 1:n)
        buf%log_lum(1:n)      = ctx%state%logl_isoc(z_idx, t_idx, 1:n)
        buf%log_teff(1:n)     = ctx%state%logt_isoc(z_idx, t_idx, 1:n)
        buf%log_g(1:n)        = ctx%state%logg_isoc(z_idx, t_idx, 1:n)
        buf%phase(1:n)        = ctx%state%phase_isoc(z_idx, t_idx, 1:n)
        buf%co_ratio(1:n)     = ctx%state%ffco_isoc(z_idx, t_idx, 1:n)
        buf%log_mdot(1:n)     = ctx%state%lmdot_isoc(z_idx, t_idx, 1:n)
        
        ! Weights are calculated fresh every step, but zeroing is safe
        buf%weights(1:n)      = 0.0_wp

    end subroutine load_timestep_data

    ! ------------------------------------------------------------------------
    ! PHYSICS ORCHESTRATION
    ! ------------------------------------------------------------------------

    !> @brief Applies all stellar population modifications to the buffer.
    !>
    !> @details
    !> Transforms the raw isochrone data in the buffer by:
    !> 1. Computing IMF weights.
    !> 2. Redistributing Horizontal Branch stars (Blue HB).
    !> 3. Adding Blue Straggler stars.
    !> 4. Shifting and scaling Giant Branch (AGB/RGB) properties.
    !>
    !> @note
    !> Uses rank-remapping (get_2d_view) to interface the 1D buffer with 
    !> legacy physics routines that expect 2D arrays.
    !>
    !> @param[in,out] ctx          Simulation context.
    !> @param[in]     pset         User settings.
    !> @param[in]     time_log_yr  Current age in log(years).
    !> @param[in,out] buf          The working buffer (modified in place).
    subroutine apply_isochrone_physics(ctx, pset, time_log_yr, buf)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        real(WP), intent(in)                :: time_log_yr
        type(isochrone_buffer_t), intent(inout) :: buf

        ! Local variables
        real(WP) :: hb_weight
        
        ! Wrapper for n_stars to satisfy legacy array interfaces
        integer :: n_stars_wrap(1)
        
        ! Pointers for 2D views (Rank-2)
        real(WP), pointer :: p_mini(:,:), p_mact(:,:)
        real(WP), pointer :: p_logl(:,:), p_logt(:,:), p_logg(:,:)
        real(WP), pointer :: p_phase(:,:)

        ! Initialize wrapper with current count
        n_stars_wrap(1) = buf%n_stars
        hb_weight = 0.0_wp

        ! 1. Calculate IMF Weights
        call compute_imf_weights(ctx, buf%initial_mass(1:buf%n_stars), &
                                 buf%weights(1:buf%n_stars), buf%n_stars)

        ! --------------------------------------------------------------------
        ! PREPARE 2D VIEWS (Inline Rank Remapping)
        ! --------------------------------------------------------------------
        ! We map the 1D buffer (1:NM) to a 2D pointer (1:1, 1:NM).
        ! This is safe because buf%components are pointers (implicit targets).
        p_mini(1:1, 1:NM)  => buf%initial_mass
        p_mact(1:1, 1:NM)  => buf%current_mass
        p_logl(1:1, 1:NM)  => buf%log_lum
        p_logt(1:1, 1:NM)  => buf%log_teff
        p_logg(1:1, 1:NM)  => buf%log_g
        p_phase(1:1, 1:NM) => buf%phase

        ! --------------------------------------------------------------------
        ! 2. MODIFY HORIZONTAL BRANCH (HB)
        ! --------------------------------------------------------------------
        if (pset%fbhb > 0.0_wp .or. pset%sbss > 1.0e-3_wp) then
            
            ! We pass '1' as the time index because our view is (1, NM).
            ! n_stars_wrap is modified in-place if new stars are added.
            call modify_horizontal_branch( &
                ctx, 1, pset%fbhb, time_log_yr, hb_weight, &
                n_stars_wrap, &
                p_mini, p_mact, p_logl, p_logt, p_logg, p_phase, &
                buf%weights) ! weights is 1D in the interface
        end if

        ! --------------------------------------------------------------------
        ! 3. ADD BLUE STRAGGLERS (BS)
        ! --------------------------------------------------------------------
        if (time_log_yr >= BHB_SBS_TIME .and. pset%sbss > 1.0e-3_wp) then
            
            call apply_blue_stragglers( &
                ctx, 1, pset%zmet, pset%sbss, hb_weight, &
                n_stars_wrap, &
                p_mini, p_mact, p_logl, p_logt, p_logg, p_phase, &
                buf%weights)
        end if

        ! Sync buffer count with the wrapper (in case stars were added)
        buf%n_stars = n_stars_wrap(1)

        ! --------------------------------------------------------------------
        ! 4. MODIFY GIANT BRANCH (AGB / RGB / PAGB)
        ! --------------------------------------------------------------------
        ! This routine does NOT add stars, so we pass buf%n_stars directly.
        ! We must pass pset%zmet (metallicity index) from settings.
        
        call modify_giant_branch( &
            ctx, 1, pset%zmet, time_log_yr, buf%n_stars, &
            pset%delt, pset%dell, pset%pagb, pset%redgb, pset%agb, &
            p_logl, p_logt, p_phase, buf%weights)

    end subroutine apply_isochrone_physics

    ! ------------------------------------------------------------------------
    ! INTEGRATION KERNELS
    ! ------------------------------------------------------------------------

    !> @brief Computes integrated scalar properties (Mass, Lbol) for the population.
    !>
    !> @details
    !> 1. Sums the mass of living stars.
    !> 2. Adds the mass of stellar remnants (BH, NS, WD) based on the turn-off mass.
    !> 3. Sums the bolometric luminosity of all stars.
    !>
    !> @param[in,out] ctx      Simulation context.
    !> @param[in]     buf      The isochrone buffer (contains weights and properties).
    !> @param[out]    tot_mass Total Mass of the SSP [M_sol].
    !> @param[out]    tot_lbol Total Bolometric Luminosity [log10(L_sol)].
    subroutine compute_integrated_properties(ctx, buf, tot_mass, tot_lbol)
        type(fsps_context_t), intent(inout) :: ctx
        type(isochrone_buffer_t), intent(in) :: buf
        real(WP), intent(out) :: tot_mass, tot_lbol

        integer :: n
        real(WP) :: max_living_mass, linear_lum

        n = buf%n_stars
        
        ! Handle empty buffer edge case
        if (n <= 0) then
            tot_mass = 0.0_wp
            tot_lbol = -99.0_wp ! Arbitrary low value
            return
        end if

        ! 1. Calculate Mass of Living Stars
        !    Vectorized dot product: sum(weights * current_mass)
        tot_mass = sum(buf%weights(1:n) * buf%current_mass(1:n))
        
        ! 2. Add Remnant Mass (Black Holes, Neutron Stars, White Dwarfs)
        if (ctx%add_stellar_remnants_val == 1) then
             
             ! The turn-off mass is approximately the maximum initial mass 
             ! of stars still present in the isochrone.
             max_living_mass = maxval(buf%initial_mass(1:n))
             
             ! This routine integrates the IMF for dead stars and adds to tot_mass
             call add_remnant_mass(ctx, tot_mass, max_living_mass)
        end if

        ! 3. Calculate Total Bolometric Luminosity
        !    L_tot = sum( weight * 10^logL )
        !    We compute the linear sum first, then take log10.
        
        linear_lum = sum(buf%weights(1:n) * exp(buf%log_lum(1:n) * LN10))
        
        ! Prevent log(0)
        tot_lbol = log10(max(linear_lum, SAFE_FLOOR))

    end subroutine compute_integrated_properties

    !> @brief Sums the spectra of all stars in the buffer into the SSP spectrum.
    !>
    !> @details
    !> Iterates over every star in the isochrone buffer:
    !> 1. Checks filters (Weight > cutoff, EVTYPE, MASSCUT).
    !> 2. Generates the specific stellar spectrum for that star.
    !> 3. Accumulates the weighted spectrum into the total grid.
    !>
    !> @param[in,out] ctx       Simulation context.
    !> @param[in]     pset      User settings (for EVTYPE/MASSCUT).
    !> @param[in]     buf       The isochrone buffer.
    !> @param[in,out] spec_out  The output spectrum accumulator (L_sol/Hz).
    !> @param[in,out] work_spec Workspace array for single star spectrum.
    subroutine accumulate_spectrum(ctx, pset, buf, spec_out, work_spec)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        type(isochrone_buffer_t), intent(in):: buf
        real(WP), intent(inout), contiguous :: spec_out(:)
        real(WP), intent(inout), contiguous :: work_spec(:)

        ! Local variables
        integer :: j, k, batch_count
        real(WP) :: linear_lbol, current_weight

        ! L1 Cache Accumulator (Small enough to stay hot)
        real(WP) :: batch_sum(size(spec_out))

        batch_count = 0
        batch_sum   = 0.0_wp

        ! 1. Loop over all stars in the buffer
        do j = 1, buf%n_stars
            
            current_weight = buf%weights(j)

            ! Optimization: Skip stars with negligible contribution.
            ! CAUTION: Stars with weight ~0.0 must be genuinely dead/empty, 
            ! not just low-mass/rare.
            if (current_weight <= 0.0_wp) cycle

            ! ----------------------------------------------------------------
            ! FILTERS
            ! ----------------------------------------------------------------
            
            ! Filter by Evolutionary Phase (pset%evtype)
            ! -1 means "All Phases". Otherwise, match integer phase.
            if (pset%evtype /= -1) then
                if (int(buf%phase(j)) /= pset%evtype) cycle
            end if

            ! Filter by Initial Mass (pset%masscut)
            ! Only include stars below the mass cut (Original FSPS logic)
            if (buf%initial_mass(j) >= pset%masscut) cycle

            ! ----------------------------------------------------------------
            ! GENERATION
            ! ----------------------------------------------------------------
            
            ! Convert LogL -> Linear L (Required by get_stellar_spectrum)
            linear_lbol = exp(buf%log_lum(j) * LN10)

            ! Retrieve spectrum for this specific star
            call get_stellar_spectrum( &
                ctx, &
                pset, &
                buf%current_mass(j), &
                buf%log_teff(j), &
                linear_lbol, &
                buf%log_g(j), &
                buf%phase(j), &
                buf%co_ratio(j), &
                buf%log_mdot(j), &
                work_spec) ! Output

            ! ----------------------------------------------------------------
            ! BATCH ACCUMULATION (L1 Cache)
            ! ----------------------------------------------------------------
            ! Accumulate into local stack array (fastest access)
            ! Compiler will likely SIMDize this loop
            do k = 1, size(spec_out)
                batch_sum(k) = batch_sum(k) + (current_weight * work_spec(k))
            end do
            
            batch_count = batch_count + 1

            ! ----------------------------------------------------------------
            ! FLUSH TO MAIN MEMORY
            ! ----------------------------------------------------------------
            if (batch_count >= BATCH_SIZE) then
                spec_out = spec_out + batch_sum
                batch_sum = 0.0_wp
                batch_count = 0
            end if

        end do

        ! Flush remaining stars
        if (batch_count > 0) then
            spec_out = spec_out + batch_sum
        end if

    end subroutine accumulate_spectrum

    ! ------------------------------------------------------------------------
    ! POST-PROCESSING HELPERS
    ! ------------------------------------------------------------------------

    !> @brief Interpolates the sparse calculated time grid onto the full high-res grid.
    !>
    !> @details
    !> When `TIME_RES_INCR > 1`, the SSP generation only computes every Nth time step.
    !> This routine fills in the intermediate steps using interpolation.
    !>
    !> **Interpolation Logic:**
    !> - **Spectrum:** Log-linear interpolation (power law).
    !>   \f$ S_{new} = 10^{( (1-dt) \log S_1 + dt \log S_2 )} \f$
    !> - **Mass/Lbol:** Linear interpolation.
    !>   \f$ M_{new} = (1-dt) M_1 + dt M_2 \f$
    !>
    !> @param[in]     ctx       Simulation context (contains time arrays).
    !> @param[in]     z_idx     Metallicity index (to access correct time array).
    !> @param[in,out] mass_grid Integrated mass history (modified in place).
    !> @param[in,out] lbol_grid Integrated luminosity history (modified in place).
    !> @param[in,out] spec_grid Spectral grid (modified in place).
    subroutine interpolate_time_grid(ctx, z_idx, mass_grid, lbol_grid, spec_grid)
        type(fsps_context_t), intent(in)    :: ctx
        integer, intent(in)                 :: z_idx
        real(WP), intent(inout), contiguous :: mass_grid(:)
        real(WP), intent(inout), contiguous :: lbol_grid(:)
        real(WP), intent(inout), contiguous :: spec_grid(:,:)

        ! Local variables
        integer :: j
        integer :: idx_low, idx_high
        integer :: raw_idx_low
        real(WP) :: dt_weight
        real(WP) :: t_target, t1, t2
        integer :: nt, nt_full

        real(WP), pointer :: time_computed(:)
        real(WP), pointer :: time_full(:)

        time_computed => ctx%state%timestep_isoc(z_idx, :)
        time_full     => ctx%state%time_full

        nt      = ctx%state%nt
        nt_full = ctx%state%ntfull
        
        ! Optimization: Initialize "Hunter" index
        raw_idx_low = 1

        ! Iterate over the FULL high-res time grid
        do j = 1, nt_full
            
            if (mod(j - 1, TIME_RES_INCR) == 0) cycle

            t_target = time_full(j)
            
            ! Optimization: Sequential Hunt
            ! Since t_target increases monotonically, we only need to move forward.
            ! Check if we need to advance the index
            do while (raw_idx_low < nt - 1)
                if (time_computed(raw_idx_low + 1) > t_target) exit
                raw_idx_low = raw_idx_low + 1
            end do
            
            ! 2. Calculate Interpolation Weight (dt)
            t1 = time_computed(raw_idx_low)
            t2 = time_computed(raw_idx_low + 1)
            
            if (abs(t2 - t1) > tiny(0.0_wp)) then
                dt_weight = (t_target - t1) / (t2 - t1)
            else
                dt_weight = 0.0_wp
            end if

            ! 3. Map Sparse Indices
            idx_low  = 1 + (raw_idx_low - 1) * TIME_RES_INCR
            idx_high = idx_low + TIME_RES_INCR

            ! 4. Interpolate Scalar Properties
            mass_grid(j) = (1.0_wp - dt_weight) * mass_grid(idx_low) + &
                           (dt_weight)          * mass_grid(idx_high)

            lbol_grid(j) = (1.0_wp - dt_weight) * lbol_grid(idx_low) + &
                           (dt_weight)          * lbol_grid(idx_high)

            ! 5. Interpolate Spectrum
            spec_grid(:, j) = exp( &
                LN10 * (1.0_wp - dt_weight) * log10(max(spec_grid(:, idx_low), SAFE_FLOOR)) + &
                LN10 * (dt_weight)          * log10(max(spec_grid(:, idx_high), SAFE_FLOOR)) )

        end do

    end subroutine interpolate_time_grid
    
    !> @brief Configures IMF parameters and handles user-defined IMF files.
    !>
    !> @details
    !> Maps the user parameters (pset) into the global context state.
    !> If a user-defined IMF is selected (imf_type=5), this routine handles
    !> opening and parsing the definition file.
    !>
    !> @param[in,out] ctx   Simulation context.
    !> @param[in]     pset  User settings.
    subroutine configure_imf_parameters(ctx, pset)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset

        ! Local variables
        integer :: imf_type
        integer :: unit_imf, io_status
        integer :: i
        character(len=1024) :: filepath

        imf_type = ctx%imf_type_val

        ! 1. Validate IMF Type
        !    Supported types: 0=Salpeter, 1=Chabrier, 2=Kroupa, 
        !                     3=van Dokkum, 4=Dave, 5=User-defined
        if (imf_type < 0 .or. imf_type > 5) then
            write(*,*) 'SSP_GEN ERROR: IMF type outside of range [0-5]:', imf_type
            stop
        end if

        ! 2. Transfer Scalar Parameters
        ctx%state%imf_alpha(1) = pset%imf1
        ctx%state%imf_alpha(2) = pset%imf2
        ctx%state%imf_alpha(3) = pset%imf3
        ctx%state%imf_vdmc     = pset%vdmc
        ctx%state%imf_mdave    = pset%mdave

        ! 3. Handle User-Defined IMF (Type 5)
        if (imf_type == 5) then
            
            ! Determine filename
            if (len_trim(pset%imf_filename) == 0) then
                filepath = trim(ctx%sps_home) // '/data/imf.dat'
            else
                filepath = trim(ctx%sps_home) // '/data/' // trim(pset%imf_filename)
            end if

            ! Open file
            open(newunit=unit_imf, file=trim(filepath), status='old', &
                 action='read', iostat=io_status)
            
            if (io_status /= 0) then
                write(*,*) 'SSP_GEN ERROR: Could not open IMF file: ', trim(filepath)
                stop
            end if

            ! Read Loop
            i = 0
            do
                ! Read into row i+1
                read(unit_imf, *, iostat=io_status) &
                     ctx%state%imf_user_alpha(1, i+1), & ! Mass Lower
                     ctx%state%imf_user_alpha(2, i+1), & ! Mass Upper
                     ctx%state%imf_user_alpha(3, i+1)    ! Slope
                
                if (io_status /= 0) exit ! EOF or Error
                
                i = i + 1
                
                ! Safety check for array bounds (100 is typical hard limit in legacy)
                if (i >= 100) then
                    write(*,*) 'SSP_GEN WARNING: Reached max segments (100) in user IMF.'
                    exit
                end if
            end do
            
            close(unit_imf)

            if (i == 0) then
                write(*,*) 'SSP_GEN ERROR: User IMF file was empty or invalid.'
                stop
            end if

            ! Store count and limits
            ctx%state%n_user_imf = i
            
            ! Define limits based on first lower bound and last upper bound
            ctx%state%imf_lower_limit = ctx%state%imf_user_alpha(1, 1)
            ctx%state%imf_upper_limit = ctx%state%imf_user_alpha(2, i)

        end if
        
        ! Verbose Logging (Optional, matching legacy behavior)
        if (VERBOSE == 1) then
            write(*,*) ''
            write(*,'("   Log(Z/Zsol): ",F6.3)') log10(ctx%state%zlegend(pset%zmet)/0.019_wp)
            write(*,'("   Fraction of blue HB stars: ",F6.3)') pset%fbhb
            write(*,'("   Ratio of BS to HB stars  : ",F6.3)') pset%sbss
            write(*,'("   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
                 pset%delt, pset%dell
            
            select case (imf_type)
            case (2)
                write(*,'("   IMF: ",I1,", slopes= ",3F4.1)') imf_type, ctx%state%imf_alpha
            case (3)
                write(*,'("   IMF: ",I1,", cut-off= ",F4.2)') imf_type, ctx%state%imf_vdmc
            case default
                write(*,'("   IMF: ",I1)') imf_type
            end select
        end if

    end subroutine configure_imf_parameters
    
    !> @brief Applies instrumental Line Spread Function (LSF) smoothing to the grid.
    !>
    !> @details
    !> Iterates over every time step in the spectral grid and convolves the 
    !> spectrum with the user-provided LSF (if loaded).
    !>
    !> @param[in]     ctx       Simulation context (contains LSF data).
    !> @param[in,out] spec_grid The full spectral grid (Lambda, Time).
    subroutine apply_lsf_smoothing_grid(ctx, spec_grid)
        type(fsps_context_t), intent(in)    :: ctx
        real(WP), intent(inout), contiguous :: spec_grid(:,:)

        integer :: j, n_times

        n_times = size(spec_grid, 2)

        ! Loop over time columns
        do j = 1, n_times
            
            ! Call the smoothing kernel for this specific time step.
            ! Note: The '99.0_wp' is a special flag in fsps_smoothing indicating
            ! that we should use the tabulated LSF kernel rather than a 
            ! simple Gaussian velocity dispersion.
            
            call apply_smoothing( &
                ctx, &
                ctx%state%spec_lambda, &      ! Wavelength grid
                spec_grid(:, j), &            ! Spectrum (in/out)
                SMOOTHING_SIGMA_FLAG, &       ! Sigma flag
                ctx%state%lsfinfo%minlam, &   ! LSF range min
                ctx%state%lsfinfo%maxlam, &   ! LSF range max
                ctx%state%lsfinfo%lsf)        ! LSF kernel array

        end do

    end subroutine apply_lsf_smoothing_grid

end module fsps_ssp