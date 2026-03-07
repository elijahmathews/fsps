module fsps_csp
    !> @brief Core module for Composite Stellar Population (CSP) synthesis.
    !>
    !> @details
    !> This module calculates the spectral energy distribution of a composite
    !> population by integrating Simple Stellar Populations (SSPs) over a 
    !> user-defined Star Formation History (SFH).
    !>
    !> **Key Architecture:**
    !> 1. **Driver (`compute_csp_scenario`):** Orchestrates the loop over requested 
    !>    output ages and manages high-level post-processing (IGM, AGN, Smoothing).
    !> 2. **Kernel (`integrate_csp_step`):** The vectorization engine. It computes 
    !>    SFH weights and performs the linear combination of SSP spectra.
    !>    It separates flux into "Young" and "Old" components to handle 
    !>    age-dependent dust attenuation efficiently.
    !> 3. **Helpers:** Dedicated routines for SFH weight calculation and 
    !>    complex physics corrections (Dust, Nebular).
    !>
    !> **Optimization Strategy:**
    !> - **Matrix Operations:** Where possible, spectral summation is treated as 
    !>   matrix multiplication (Weights x SSP_Grid).
    !> - **Pre-calculation:** SFH weights are computed into a contiguous buffer 
    !>   before spectral integration to maximize memory locality.

    use fsps_precision, only: WP
    use fsps_constants, only: NEMLINE, SAFE_FLOOR, C_LIGHT, VERBOSE
    use fsps_types, only: params, sfhparams, compspout
    use fsps_context_types, only: fsps_context_t

    !> Physics Modules
    use fsps_sfh, only: get_sfh_properties_at_age, compute_ssp_weights
    use fsps_dust, only: apply_dust_attenuation_and_emission, apply_agn_dust_emission
    use fsps_gas, only: apply_nebular_emission
    use fsps_smoothing, only: apply_smoothing
    use fsps_cosmology, only: compute_igm_transmission

    !> Math Modules
    use fsps_interpolation, only: find_interval, interpolate_linear

    !> Spectra Modules
    use fsps_photometry, only: compute_magnitudes
    use fsps_spectral_indices, only: compute_spectral_indices

    implicit none
    private

    !> Public API
    public :: compute_csp_scenario
    ! Expose internal kernels for unit testing
    public :: compute_sfh_weights
    public :: convert_sfhparams
    public :: integrate_csp_step
    public :: apply_dust_physics
    public :: apply_post_processing

contains

    ! ------------------------------------------------------------------------
    ! MAIN DRIVER
    ! ------------------------------------------------------------------------

    !> @brief Main driver to generate CSP properties for a set of parameters.
    !>
    !> @details
    !> This routine replaces the legacy `COMPSP`. It performs the following steps:
    !> 1. **Setup:** Creates local copies of the SSP grids and adds Nebular Emission 
    !>    (if enabled) to the SSPs *before* integration.
    !> 2. **Loop:** Iterates over the requested output ages.
    !>    - If `pset%tage > 0`: Computes a single snapshot at that specific age.
    !>    - If `pset%tage <= 0`: Computes the full time evolution history.
    !> 3. **Integration:** Calls `integrate_csp_step` to sum the SSPs.
    !> 4. **Physics:** Calls `apply_dust_physics` and `apply_post_processing`.
    !> 5. **Output:** Returns an allocatable array of `compspout` structures.
    !>
    !> @param[in,out] ctx       The main FSPS context.
    !> @param[in]     pset      User parameter set.
    !> @param[in]     nzin      Number of input metallicities to use (usually 1 or ctx%nz).
    !> @param[out]    results   Allocatable array of output structures.
    !> @param[out]    status    (Optional) Error status (0 = Success).
    subroutine compute_csp_scenario(ctx, pset, nzin, results, status)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        integer, intent(in)                 :: nzin
        type(compspout), allocatable, intent(out) :: results(:)
        integer, intent(out), optional      :: status

        ! Local variables
        real(WP) :: mass_csp, lbol_csp, mdust_total
        real(WP) :: target_age
        
        integer :: i, nt, nspec, n_outputs, start_idx
        
        if (present(status)) status = 0

        ! 1. SETUP
        nt    = ctx%state%ntfull
        nspec = ctx%state%nspec
        
        if (ctx%state%check_sps_setup == 0) then
            if (present(status)) status = 1; return
        end if

        if (.not. allocated(ctx%state%ssp_basis_spec) .or. &
            .not. allocated(ctx%state%ssp_basis_mass) .or. &
            .not. allocated(ctx%state%ssp_basis_lbol)) then
            if (present(status)) status = 3
            return
        end if
        
        ! Initialize host-side working grids from SSP inputs.
        ctx%state%csp_ssp_grid(:,:,1:nzin) = ctx%state%ssp_basis_spec(:,:,1:nzin)
        ctx%state%csp_emlin_grid(:,:,1:nzin) = 0.0_wp

        ! Explicitly push the initialized grids to the device
        !$acc update device(ctx%state%csp_ssp_grid, ctx%state%csp_emlin_grid)

        if (ctx%add_neb_emission_val == 1) then
            if (nzin > 1) then
                 if (present(status)) status = 2; return
            end if

            call apply_nebular_emission(ctx, pset, ctx%state%ssp_basis_spec(:,:,1), &
                                        ctx%state%csp_ssp_grid(:,:,1), ctx%state%csp_emlin_grid(:,:,1))
            
            ! Pull updated data to host for the integrator
            !$acc update host(ctx%state%csp_ssp_grid(:,:,1), ctx%state%csp_emlin_grid(:,:,1))
        end if

        ! 2. OPTIMIZATIONS (Pre-calculations)
        ! -----------------------------------
        
        ! A. Linearize Luminosity (Avoids 10**x inside hot loops)
        !$acc kernels present(ctx)
        ctx%state%csp_ssp_lum_linear(:,1:nzin) = 10.0_wp**ctx%state%ssp_basis_lbol(:,1:nzin)
        !$acc end kernels

        ! B. Pre-calculate IGM Transmission (Constant for this PSET)
        if (ctx%add_igm_absorption_val == 1 .and. pset%zred > SAFE_FLOOR) then
            call compute_igm_transmission(ctx%state%spec_lambda, &
                                          pset%zred, pset%igm_factor, &
                                          ctx%state%csp_igm_transmission)
        end if

        if (pset%tage > 0.0_wp) then
            n_outputs = 1; start_idx = 0 
        elseif (pset%tage == -99.0_wp .and. (pset%sfh == 2 .or. pset%sfh == 3)) then
            n_outputs = 1; start_idx = -99
        else
            n_outputs = nt; start_idx = 1
        end if
        
        allocate(results(n_outputs))

        ! 3. MAIN GENERATION LOOP
        ! -----------------------
        !$acc data present(ctx)
        do i = 1, n_outputs
            ! Ensure output arrays are allocated
            if (.not. allocated(results(i)%spec)) then
                allocate(results(i)%spec(nspec))
            end if
            if (.not. allocated(results(i)%emlines)) then
                allocate(results(i)%emlines(NEMLINE))
            end if
            
            ! Determine Target Age
            if (pset%tage > 0.0_wp) then
                target_age = pset%tage
            elseif (start_idx == -99) then
                target_age = maxval(ctx%state%sfh_tab(1, 1:ctx%state%ntabsfh)) / 1.0e9_wp
            else
                target_age = 10.0_wp**(ctx%state%time_full(i) - 9.0_wp)
            end if

            ! Integration Kernel (Passes pre-calculated Linear Lum)
            call integrate_csp_step(ctx, pset, target_age, nzin, &
                                    ctx%state%csp_ssp_grid, ctx%state%csp_emlin_grid, ctx%state%ssp_basis_mass, &
                                    ctx%state%csp_ssp_lum_linear, &
                                    mass_csp, lbol_csp)

            ! Dust Physics
            call apply_dust_physics(ctx, pset, ctx%state%csp_spec_final, ctx%state%csp_emlin_final, mdust_total)

            ! Post-Processing (Passes pre-calculated IGM)
            call apply_post_processing(ctx, pset, target_age, &
                                       mass_csp, lbol_csp, mdust_total, &
                                       ctx%state%csp_spec_final, ctx%state%csp_emlin_final, &
                                       ctx%state%csp_igm_transmission, &
                                       results(i))
        end do
        !$acc end data

    end subroutine compute_csp_scenario

    ! ------------------------------------------------------------------------
    ! CORE INTEGRATION KERNEL
    ! ------------------------------------------------------------------------

    !> @brief Integrates SSPs to produce raw Young/Old composite spectra.
    !>
    !> @details
    !> 1. Calculates SFH weights.
    !> 2. Determines the age index `i_tesc` that separates "Young" (birth cloud)
    !>    from "Old" (diffuse ISM) populations based on `pset%dust_tesc`.
    !> 3. Accumulates the weighted sum of SSP spectra and emission lines into 
    !>    the buffer, preserving the Young/Old separation.
    !> 4. Computes total Stellar Mass and Bolometric Luminosity.
    !>
    !> @param[inout]  ctx        Context.
    !> @param[in]     pset       User parameters.
    !> @param[in]     tage       Target age of the galaxy [Gyr].
    !> @param[in]     nzin       Number of metallicities.
    !> @param[in]     ssp_grid   Input SSP grid [Lambda, Age, Z].
    !> @param[in]     emlin_grid Input Emission Line grid [Line, Age, Z].
    !> @param[out]    mass_csp   Total stellar mass formed (normalized).
    !> @param[out]    lbol_csp   Total bolometric luminosity [log10(L_sol)].
    subroutine integrate_csp_step(ctx, pset, tage, nzin, ssp_grid, emlin_grid, &
                                  mass_ssp, ssp_lum_linear, mass_csp, lbol_csp)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        real(WP), intent(in)                :: tage
        integer, intent(in)                 :: nzin
        real(WP), intent(in), contiguous    :: ssp_grid(:,:,:)
        real(WP), intent(in), contiguous    :: emlin_grid(:,:,:)
        real(WP), intent(in), contiguous    :: mass_ssp(:,:)
        real(WP), intent(in), contiguous    :: ssp_lum_linear(:,:)
        real(WP), intent(out)               :: mass_csp, lbol_csp

        integer :: i, k, nt, i_tesc
        integer :: i_spec, i_em
        integer :: nspec, nem
        real(WP) :: dust_age_log
        real(WP) :: linear_lbol_sum
        real(WP) :: weight_ik
        real(WP) :: sum_young, sum_old
        real(WP) :: sum_em_young, sum_em_old
        real(WP) :: dev_scalars(2)
        
        nt = ctx%state%ntfull
        nspec = size(ssp_grid, 1)
        nem   = size(emlin_grid, 1)

        ! 1. Clear Accumulators
        ! Clear host arrays
        ctx%state%spec_young  = 0.0_wp
        ctx%state%spec_old    = 0.0_wp
        ctx%state%csp_emlin_young = 0.0_wp
        ctx%state%csp_emlin_old   = 0.0_wp

        ! Clear device arrays
        !$acc kernels present(ctx)
        ctx%state%spec_young  = 0.0_wp
        ctx%state%spec_old    = 0.0_wp
        ctx%state%csp_emlin_young = 0.0_wp
        ctx%state%csp_emlin_old   = 0.0_wp
        !$acc end kernels

        ! 2. Compute SFH Weights
        !$acc serial present(ctx)
        call compute_sfh_weights(ctx, pset, tage, nzin, ctx%state%csp_weights)
        !$acc end serial

        ! 3. Determine Dust Separation Index
        if (pset%dust_tesc > SAFE_FLOOR) then
            dust_age_log = pset%dust_tesc
        else
            dust_age_log = 7.0_wp
        end if
        i_tesc = max(1, min(find_interval(ctx%state%time_full, dust_age_log), nt))

        ! 4. Integration Loop for Scalars
        ! Executed sequentially on the device using persistent context memory
        ! to bypass NVHPC compiler bugs with dummy argument copyout and local register arrays.
        !$acc serial present(ctx, mass_ssp, ssp_lum_linear)
        ctx%state%scalar_reductions(1) = 0.0_wp
        ctx%state%scalar_reductions(2) = 0.0_wp

        do k = 1, nzin
            do i = 1, nt
                if (ctx%state%csp_weights(i, k) > SAFE_FLOOR) then
                    ctx%state%scalar_reductions(1) = ctx%state%scalar_reductions(1) + &
                                                     (ctx%state%csp_weights(i, k) * mass_ssp(i, k))
                    ctx%state%scalar_reductions(2) = ctx%state%scalar_reductions(2) + &
                                                     (ctx%state%csp_weights(i, k) * ssp_lum_linear(i, k))
                end if
            end do
        end do
        !$acc end serial

        ! Fetch results safely back to the host. Explicit array bounds bypass dope-vector bugs.
        !$acc update host(ctx%state%scalar_reductions(1:2))
        mass_csp        = ctx%state%scalar_reductions(1)
        linear_lbol_sum = ctx%state%scalar_reductions(2)

        if (linear_lbol_sum > SAFE_FLOOR) then
            lbol_csp = log10(linear_lbol_sum)
        else
            lbol_csp = 0.0_wp
        end if

        ! 5. Integration Loop for Spectra
#ifdef _OPENACC

        ! =====================================================================
        ! GPU OPTIMIZED PATH (Max Occupancy, No Atomics, Coalesced Memory)
        ! =====================================================================
        !$acc parallel loop gang vector present(ctx, ssp_grid) private(sum_young, sum_old, weight_ik)
        do i_spec = 1, nspec
            sum_young = 0.0_wp
            do k = 1, nzin
                do i = 1, i_tesc
                    weight_ik = ctx%state%csp_weights(i, k)
                    if (weight_ik > SAFE_FLOOR) then
                        sum_young = sum_young + (weight_ik * ssp_grid(i_spec, i, k))
                    end if
                end do
            end do

            sum_old = 0.0_wp
            do k = 1, nzin
                do i = i_tesc + 1, nt
                    weight_ik = ctx%state%csp_weights(i, k)
                    if (weight_ik > SAFE_FLOOR) then
                        sum_old = sum_old + (weight_ik * ssp_grid(i_spec, i, k))
                    end if
                end do
            end do

            ctx%state%spec_young(i_spec) = ctx%state%spec_young(i_spec) + sum_young
            ctx%state%spec_old(i_spec)   = ctx%state%spec_old(i_spec) + sum_old
        end do

#else

        ! =====================================================================
        ! CPU OPTIMIZED PATH (Strict Stride-1 Cache Locality, Vectorized)
        ! =====================================================================
        do k = 1, nzin
            ! Young stars
            do i = 1, i_tesc
                weight_ik = ctx%state%csp_weights(i, k)
                if (weight_ik > SAFE_FLOOR) then
                    ! Stride-1 inner loop for CPU vectorization/cache
                    do i_spec = 1, nspec
                        ctx%state%spec_young(i_spec) = ctx%state%spec_young(i_spec) + &
                            (weight_ik * ssp_grid(i_spec, i, k))
                    end do
                end if
            end do
            ! Old stars
            do i = i_tesc + 1, nt
                weight_ik = ctx%state%csp_weights(i, k)
                if (weight_ik > SAFE_FLOOR) then
                    do i_spec = 1, nspec
                        ctx%state%spec_old(i_spec) = ctx%state%spec_old(i_spec) + &
                            (weight_ik * ssp_grid(i_spec, i, k))
                    end do
                end if
            end do
        end do
#endif

        ! 6. Integration Loop for Emission Lines
#ifdef _OPENACC

        ! =====================================================================
        ! GPU OPTIMIZED PATH (Max Occupancy, No Atomics, Coalesced Memory)
        ! =====================================================================
        !$acc parallel loop gang vector present(ctx, emlin_grid) private(sum_em_young, sum_em_old, weight_ik)
        do i_em = 1, nem
            sum_em_young = 0.0_wp
            do k = 1, nzin
                do i = 1, i_tesc
                    weight_ik = ctx%state%csp_weights(i, k)
                    if (weight_ik > SAFE_FLOOR) then
                        sum_em_young = sum_em_young + (weight_ik * emlin_grid(i_em, i, k))
                    end if
                end do
            end do

            sum_em_old = 0.0_wp
            do k = 1, nzin
                do i = i_tesc + 1, nt
                    weight_ik = ctx%state%csp_weights(i, k)
                    if (weight_ik > SAFE_FLOOR) then
                        sum_em_old = sum_em_old + (weight_ik * emlin_grid(i_em, i, k))
                    end if
                end do
            end do

            ctx%state%csp_emlin_young(i_em) = ctx%state%csp_emlin_young(i_em) + sum_em_young
            ctx%state%csp_emlin_old(i_em)   = ctx%state%csp_emlin_old(i_em) + sum_em_old
        end do

#else

        ! =====================================================================
        ! CPU OPTIMIZED PATH (Strict Stride-1 Cache Locality, Vectorized)
        ! =====================================================================
        do k = 1, nzin
            ! Young stars
            do i = 1, i_tesc
                weight_ik = ctx%state%csp_weights(i, k)
                if (weight_ik > SAFE_FLOOR) then
                    ! Stride-1 inner loop for CPU vectorization/cache
                    do i_em = 1, nem
                        ctx%state%csp_emlin_young(i_em) = ctx%state%csp_emlin_young(i_em) + &
                            (weight_ik * emlin_grid(i_em, i, k))
                    end do
                end if
            end do
            ! Old stars
            do i = i_tesc + 1, nt
                weight_ik = ctx%state%csp_weights(i, k)
                if (weight_ik > SAFE_FLOOR) then
                    do i_em = 1, nem
                        ctx%state%csp_emlin_old(i_em) = ctx%state%csp_emlin_old(i_em) + &
                            (weight_ik * emlin_grid(i_em, i, k))
                    end do
                end if
            end do
        end do
        
#endif

        ! 7. Fetch results back to host for subsequent physics steps
        !$acc update host(ctx%state%spec_young, ctx%state%spec_old, ctx%state%csp_emlin_young, ctx%state%csp_emlin_old)

        if (linear_lbol_sum > 0.0_wp) then
            lbol_csp = log10(linear_lbol_sum)
        else
            lbol_csp = -99.0_wp
        end if

    end subroutine integrate_csp_step

    ! ------------------------------------------------------------------------
    ! PHYSICS LAYERS
    ! ------------------------------------------------------------------------

    !> @brief Applies dust attenuation and IR emission models.
    !>
    !> @details
    !> Wraps the `fsps_dust` module to process the separated Young/Old spectra.
    !> 1. Applies birth cloud attenuation to `spec_young`.
    !> 2. Applies diffuse ISM attenuation to both `spec_young` and `spec_old`.
    !> 3. Computes IR dust re-emission via energy balance.
    !>
    !> @param[inout]  ctx           Context.
    !> @param[in]     pset          User parameters.
    !> @param[out]    spec_total    Final attenuated spectrum (L_sol/Hz).
    !> @param[out]    emlin_total   Final attenuated emission lines (L_sol).
    !> @param[out]    mdust_total   Total dust mass (M_sol).
    subroutine apply_dust_physics(ctx, pset, spec_total, emlin_total, mdust_total)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        real(WP), intent(out)               :: spec_total(:)
        real(WP), intent(out)               :: emlin_total(:)
        real(WP), intent(out)               :: mdust_total

        ! The fsps_dust routine handles all the heavy lifting, 
        ! including the combination of Young + Old components.
        call apply_dust_attenuation_and_emission( &
            ctx, &
            pset, &
            ctx%state%spec_young, &
            ctx%state%spec_old, &
            ctx%state%csp_emlin_young, &
            ctx%state%csp_emlin_old, &
            spec_total, &    ! Output
            mdust_total, &   ! Output
            emlin_total)     ! Output

    end subroutine apply_dust_physics

    !> @brief Applies final physical corrections and populates output.
    !>
    !> @details
    !> Performs the following operations in order:
    !> 1. **Mass Renormalization:** Calculates the "Surviving Mass" fraction if needed.
    !> 2. **Smoothing:** Convolves with LSF or velocity dispersion.
    !> 3. **IGM:** Applies Madau absorption if redshift > 0.
    !> 4. **AGN:** Adds AGN dust torus emission.
    !> 5. **Output:** Packages everything into the `compspout` structure.
    !>
    !> @param[in]     ctx      Context.
    !> @param[in]     pset     User parameters.
    !> @param[in]     tage     Current age [Gyr].
    !> @param[in]     mass_csp Total formed stellar mass (from integration).
    !> @param[in]     lbol_csp Total bolometric luminosity (log L_sol).
    !> @param[in]     mdust    Total dust mass.
    !> @param[in,out] spec     The spectrum (Input: attenuated; Output: final).
    !> @param[in]     emlines  Emission lines.
    !> @param[out]    result   The output structure to populate.
    subroutine apply_post_processing(ctx, pset, tage, mass_csp, lbol_csp, mdust, &
                                     spec, emlines, igm_transmission, result)
        
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: pset
        real(WP), intent(in)             :: tage
        real(WP), intent(in)             :: mass_csp, lbol_csp, mdust
        real(WP), intent(inout)          :: spec(:)
        real(WP), intent(in)             :: emlines(:)
        real(WP), intent(in), optional   :: igm_transmission(:)
        type(compspout), intent(out)     :: result

        real(WP) :: mass_frac, sfr_norm, frac_linear
        real(WP) :: current_mass_surviving, lbol_final, z_effective
        logical :: do_full_post

        do_full_post = .not. ctx%fast_mode

        ! 0. Allocate optional output arrays only in full post-processing mode.
        if (do_full_post) then
            allocate(result%mags(ctx%state%nbands))
            allocate(result%indx(ctx%state%nindx))
        end if

        ! 1. Calculate Mass/SFR Properties (Renormalization)
        call get_sfh_properties_at_age(ctx, pset, tage, mass_frac, sfr_norm, frac_linear)
        
        current_mass_surviving = mass_csp * mass_frac
        
        ! 2. Handle History vs Snapshot Normalization
        if (pset%tage <= 0.0_wp) then
            ! History Mode: Output represents the total galaxy luminosity at time T.
            !$acc kernels present(spec)
            spec       = spec * mass_frac
            !$acc end kernels
            lbol_final = lbol_csp + log10(max(mass_frac, SAFE_FLOOR))
        else
            ! Snapshot Mode: Output is normalized to 1 M_sol formed *total*.
            lbol_final = lbol_csp
            if (mass_frac > SAFE_FLOOR) then
                sfr_norm = sfr_norm / mass_frac
            end if
            mass_frac = 1.0_wp 
        end if

        ! 3. Instrumental Smoothing
        if (do_full_post .and. pset%sigma_smooth > 0.0_wp) then
            !$acc update host(spec)
            call apply_smoothing(ctx, ctx%state%spec_lambda, spec, &
                                 pset%sigma_smooth, &
                                 pset%min_wave_smooth, pset%max_wave_smooth)
            !$acc update device(spec)
        end if

        ! 4. IGM Absorption
        if (ctx%add_igm_absorption_val == 1 .and. pset%zred > SAFE_FLOOR .and. present(igm_transmission)) then
            !$acc kernels present(spec, igm_transmission)
            spec = spec * igm_transmission
            !$acc end kernels
        end if

        ! 5. AGN Dust Emission
        if (ctx%add_agn_dust_val == 1 .and. pset%fagn > SAFE_FLOOR) then
            call apply_agn_dust_emission(ctx, pset, ctx%state%spec_lambda, &
                                         lbol_final, spec)
        end if

        ! 6. Calculate Magnitudes and Spectral Indices (full post-processing only)
        if (do_full_post) then
            !$acc update host(spec)

            ! Redshift for magnitudes calculation
            if (ctx%redshift_colors_val == 1) then
                 ! Inverse lookup from Age (tage) to Redshift using pre-computed spline.
                 ! cosmospl(:,2) is Age(Gyr), cosmospl(:,1) is Redshift.
                 z_effective = interpolate_linear(ctx%state%cosmospl(:,2), &
                                                  ctx%state%cosmospl(:,1), &
                                                  tage)
                 z_effective = min(max(z_effective, 0.0_wp), 20.0_wp)
            else
                 z_effective = pset%zred
            end if

            ! Compute Magnitudes
            if (pset%compute_mags == 1) then
                call compute_magnitudes(ctx, z_effective, spec, result%mags, pset%mag_compute)
            else
                result%mags = -99.0_wp
            end if

            ! Compute Spectral Indices
            if (pset%compute_indices == 1) then
                call compute_spectral_indices(ctx, ctx%state%spec_lambda, spec, result%indx)
            else
                result%indx = -99.0_wp
            end if
        end if

        ! 7. Populate Output Structure
        ! Update scalars/arrays modified on device
        !$acc update host(spec, emlines)

        result%age      = log10(tage * 1.0e9_wp)
        result%mass_csp = current_mass_surviving
        result%lbol_csp = lbol_final
        result%sfr      = sfr_norm
        result%mdust    = mdust * mass_frac
        result%mformed  = mass_frac
        result%spec     = max(spec, SAFE_FLOOR)
        result%emlines  = emlines * mass_frac
        ! result%mags and result%indx populated above

    end subroutine apply_post_processing

    ! ------------------------------------------------------------------------
    ! SFH LOGIC
    ! ------------------------------------------------------------------------

    !> @brief Computes the mass weights for each SSP bin for the given SFH.
    !>
    !> @details
    !> Orchestrates the calculation of weights for different SFH types:
    !> - **SSP (0):** Single burst.
    !> - **Tau/Delayed (1,4):** Exponential models + Constant + Burst.
    !> - **Simha (5):** Mixed Linear + Delayed-Tau.
    !> - **Tabular (2,3):** Sum of linear segments from a table.
    !>
    !> Ensures that the final weights sum to 1.0 M_sol formed (or appropriate 
    !> mass fraction) so that the resulting spectrum is per unit mass formed.
    !>
    !> @param[inout] ctx     Context.
    !> @param[in]    pset    User parameters.
    !> @param[in]    tage    Age of the galaxy [Gyr].
    !> @param[in]    nzin    Number of metallicity bins.
    !> @param[inout] weights Weights [ntfull, nzin].
    subroutine compute_sfh_weights(ctx, pset, tage, nzin, weights)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        real(WP), intent(in)                :: tage
        integer, intent(in)                 :: nzin
        ! Note: weights is intent(inout) rather than intent(out) to bypass
        ! NVHPC array descriptor reallocation bugs on the device.
        real(WP), intent(inout), contiguous :: weights(:,:)

        !$acc routine seq

        ! Local variables
        type(sfhparams) :: sfh
        real(WP) :: mass1, mass2
        real(WP) :: frac_linear, mass_frac, sfr, fburst_val
        real(WP) :: t1, t2, dt, zbin, dz
        integer :: i, k, imin, imax
        integer :: nt, j

        ! 1. Setup
        nt = ctx%state%ntfull
        do k = 1, nzin
            do i = 1, nt
                weights(i, k) = 0.0_wp
            end do
        end do
        do i = 1, nt
            ctx%state%sfh_w_tmp1(i) = 0.0_wp
            ctx%state%sfh_w_tmp2(i) = 0.0_wp
        end do

        ! Initialize SFH struct with unit conversions (Gyr -> Yr)
        call convert_sfhparams(pset, tage, sfh)

        ! Optimization: Only calculate up to the age of the universe + buffer
        ! (Legacy logic: +2 steps to be safe against interpolation edge cases)
        imax = min(max(find_interval(ctx%state%time_full, log10(sfh%tage)) + 2, 1), nt)

        ! --------------------------------------------------------------------
        ! CASE: SSP (Type 0)
        ! --------------------------------------------------------------------
        if (pset%sfh == 0) then
            sfh%type = -1        ! -1 indicates Burst in compute_ssp_weights
            sfh%tb   = sfh%tage  ! Burst at 'now' (lookback time = tage)
            
            imin = max(imax - 2, 1) ! Optimization: SSP is local
            call compute_ssp_weights(ctx, sfh, imin, imax, weights(:, 1))
            return
        end if

        ! --------------------------------------------------------------------
        ! CASE: Tau (1) or Delayed Tau (4)
        ! --------------------------------------------------------------------
        if (pset%sfh == 1 .or. pset%sfh == 4) then
            ! A. Main Exponential Component
            sfh%type = pset%sfh
            imin = 0
            call compute_ssp_weights(ctx, sfh, imin, imax, weights(:, 1))

            ! Normalize to 1.0 mass
            mass1 = sum(weights(1:imax, 1))
            if (mass1 < SAFE_FLOOR) mass1 = 1.0_wp
            weights(:, 1) = weights(:, 1) / mass1

            ! B. Add Constant and Burst (if requested)
            if (pset%const > 0.0_wp .or. pset%fburst > SAFE_FLOOR) then
                
                ! Constant Component
                sfh%type = 0 ! Constant
                call compute_ssp_weights(ctx, sfh, imin, imax, ctx%state%sfh_w_tmp1)
                mass1 = 0.0_wp
                do i = 1, imax
                    mass1 = mass1 + ctx%state%sfh_w_tmp1(i)
                end do
                if (mass1 < SAFE_FLOOR) mass1 = 1.0_wp

                ! Burst Component
                ctx%state%sfh_w_tmp2 = 0.0_wp
                fburst_val = 0.0_wp
                
                if (sfh%tb >= 0.0_wp) then
                    sfh%type = -1 ! Burst
                    call compute_ssp_weights(ctx, sfh, imin, imax, ctx%state%sfh_w_tmp2)
                    fburst_val = pset%fburst
                    
                    ! Extend imax to include burst if it happened earlier
                    imax = max(imax, min(max(find_interval(ctx%state%time_full, log10(sfh%tb)) + 2, 1), nt))
                end if

                ! Combine: (1 - C - B) * Tau + C * Const + B * Burst
                do i = 1, nt
                    weights(i, 1) = (1.0_wp - pset%const - fburst_val) * weights(i, 1) + &
                                    pset%const * (ctx%state%sfh_w_tmp1(i) / mass1) + &
                                    fburst_val * ctx%state%sfh_w_tmp2(i)
                end do
            end if
            return
        end if

        ! --------------------------------------------------------------------
        ! CASE: Simha (5)
        ! --------------------------------------------------------------------
        if (pset%sfh == 5) then
            ! A. Delayed Tau Portion
            sfh%type = 4
            imin = 0
            call compute_ssp_weights(ctx, sfh, imin, imax, ctx%state%sfh_w_tmp1)
            mass1 = 0.0_wp
            do i = 1, imax
                mass1 = mass1 + ctx%state%sfh_w_tmp1(i)
            end do

            ! B. Linear Cutoff Portion
            sfh%type = 5
            sfh%use_simha_limits = 1
            call compute_ssp_weights(ctx, sfh, imin, imax, ctx%state%sfh_w_tmp2)
            sfh%use_simha_limits = 0
            mass2 = 0.0_wp
            do i = 1, imax
                mass2 = mass2 + ctx%state%sfh_w_tmp2(i)
            end do

            ! Normalize
            if (mass1 < SAFE_FLOOR) mass1 = 1.0_wp
            if (mass2 < SAFE_FLOOR) mass2 = 1.0_wp

            ! C. Get Mixing Fraction
            call get_sfh_properties_at_age(ctx, pset, tage, mass_frac, sfr, frac_linear)

            ! Combine
            do i = 1, nt
                weights(i, 1) = (ctx%state%sfh_w_tmp1(i) / mass1) * (1.0_wp - frac_linear) + &
                                (ctx%state%sfh_w_tmp2(i) / mass2) * frac_linear
            end do
            return
        end if

        ! --------------------------------------------------------------------
        ! CASE: Tabular (2=File, 3=Array)
        ! --------------------------------------------------------------------
        if (pset%sfh == 2 .or. pset%sfh == 3) then
            ! We treat the table as a series of linear SFH segments.
            ! sfh_tab is in [Years, SFR, Z]
            
            ! Loop over table bins
            do j = 1, ctx%state%ntabsfh - 1
                
                ! Bin Edges in Lookback Time (Note: Table is in Forward Time)
                ! t1 = Older edge, t2 = Younger edge
                t1 = tage * 1.0e9_wp - ctx%state%sfh_tab(1, j+1)
                t2 = tage * 1.0e9_wp - ctx%state%sfh_tab(1, j)

                ! Skip bins entirely in the future
                if (t2 < 0.0_wp) cycle

                ! Average Z of bin
                zbin = (ctx%state%sfh_tab(3, j) + ctx%state%sfh_tab(3, j+1)) / 2.0_wp

                ! Calculate Slope for "Linear" Type 5 model
                ! Slope = - dSFR / dt / SFR_start (approx normalization for Simha model logic)
                ! Legacy code uses this specific definition compatible with Type 5 math.
                ! Note: Denominator is dt = (t2 - t1)
                dt = max(t2 - t1, SAFE_FLOOR)
                sfh%sf_slope = -(ctx%state%sfh_tab(2, j+1) - ctx%state%sfh_tab(2, j)) / &
                               dt / max(ctx%state%sfh_tab(2, j+1), SAFE_FLOOR)

                ! Set Integration Limits for this bin
                sfh%type = 5
                sfh%tq   = min(max(t1, 10.0_wp**ctx%tiny_logt_val), 10.0_wp**ctx%state%time_full(nt))
                sfh%tage = min(max(t2, 10.0_wp**ctx%tiny_logt_val), 10.0_wp**ctx%state%time_full(nt))
                sfh%sf_trunc = sfh%tage - sfh%tq ! Width of integration

                ! Calculate Mass formed in this bin (Analytic integral of linear segment)
                ! Mass = SFR_end * (1 + slope/2 * width) * width
                mass2 = ctx%state%sfh_tab(2, j+1) * &
                        (1.0_wp + sfh%sf_slope/2.0_wp * (sfh%tage + sfh%tq - 2.0_wp*t1)) * &
                        (sfh%tage - sfh%tq)

                ! Get Weights for this segment
                imin = min(max(find_interval(ctx%state%time_full, log10(max(t1, SAFE_FLOOR))) - 1, 0), nt)
                imax = min(max(find_interval(ctx%state%time_full, log10(max(t2, SAFE_FLOOR))) + 2, 0), nt)
                
                call compute_ssp_weights(ctx, sfh, imin, imax, ctx%state%sfh_w_tmp1)
                
                mass1 = 0.0_wp
                do i = 1, nt
                    mass1 = mass1 + ctx%state%sfh_w_tmp1(i)
                end do
                if (mass1 < SAFE_FLOOR) mass1 = 1.0_wp

                ! Distribute to Metallicities
                if (nzin > 1) then
                    ! Find Z interval
                    k = max(min(find_interval(ctx%state%zlegend, zbin), ctx%state%nz - 1), 1)
                    
                    ! Calc interpolation factor dz
                    dz = (log10(zbin) - log10(ctx%state%zlegend(k))) / &
                         (log10(ctx%state%zlegend(k+1)) - log10(ctx%state%zlegend(k)))
                    dz = max(min(dz, 1.0_wp), -1.0_wp) ! Clamp extrapolation
                    
                    ! Vectorized Add
                    do i = 1, nt
                        weights(i, k)   = weights(i, k)   + (1.0_wp - dz) * ctx%state%sfh_w_tmp1(i) * (mass2 / mass1)
                        weights(i, k+1) = weights(i, k+1) + dz            * ctx%state%sfh_w_tmp1(i) * (mass2 / mass1)
                    end do
                else
                    ! Single Z
                    do i = 1, nt
                        weights(i, 1) = weights(i, 1) + ctx%state%sfh_w_tmp1(i) * (mass2 / mass1)
                    end do
                end if
            end do
        end if

    end subroutine compute_sfh_weights

    !> @brief Converts user parameters (pset) into internal physical units (sfhparams).
    !>
    !> @details
    !> Performs necessary unit conversions (Gyr -> Years), calculates lookback times,
    !> and sets up derived parameters like truncation times and zero-crossings.
    !>
    !> @param[in]  pset Parameters from user.
    !> @param[in]  tage Age of the galaxy [Gyr].
    !> @param[out] sfh  Internal SFH structure.
    pure subroutine convert_sfhparams(pset, tage, sfh)
        type(params), intent(in)     :: pset
        real(WP), intent(in)         :: tage
        type(sfhparams), intent(out) :: sfh

        real(WP) :: start_time
        
        ! 1. Determine Start Time (Offset)
        if (any([1, 4, 5] == pset%sfh)) then
            start_time = pset%sf_start * 1.0e9_wp
        else
            start_time = 0.0_wp
        end if

        ! 2. Basic Conversions (Gyr -> Yr)
        sfh%tage     = tage * 1.0e9_wp - start_time
        sfh%tburst   = pset%tburst * 1.0e9_wp - start_time
        sfh%sf_trunc = pset%sf_trunc * 1.0e9_wp - start_time
        sfh%tau      = pset%tau * 1.0e9_wp
        
        ! Note sign flip: pset slope is (+ = increasing in forward time)
        ! sfh slope is (+ = increasing in lookback time)
        sfh%sf_slope = -pset%sf_slope / 1.0e9_wp

        ! 3. Derived Times (Lookback)
        sfh%tb = sfh%tage - sfh%tburst

        ! Truncation Time (Lookback)
        if (sfh%sf_trunc <= 0.0_wp .or. sfh%sf_trunc > sfh%tage) then
            sfh%tq = 0.0_wp
        else
            sfh%tq = sfh%tage - sfh%sf_trunc
        end if

        ! 4. Zero Crossing for Simha/Linear models
        ! t0 = time where SFR hits zero
        if (sfh%sf_slope > SAFE_FLOOR) then
            sfh%t0 = sfh%tq - (1.0_wp / sfh%sf_slope)
        else
            sfh%t0 = 0.0_wp
        end if

        ! Clamp t0 to be physical
        if (sfh%t0 > sfh%tq .or. sfh%t0 <= 0.0_wp) then
            sfh%t0 = 0.0_wp
        end if
        
        ! Initialize Simha flag
        sfh%use_simha_limits = 0

    end subroutine convert_sfhparams

end module fsps_csp
