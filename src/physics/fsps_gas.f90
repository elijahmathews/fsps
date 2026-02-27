module fsps_gas
    !> @brief
    !> Handles nebular emission physics (gas phase).
    !>
    !> @details
    !> This module processes the reprocessing of ionizing stellar radiation by
    !> the Interstellar Medium (ISM). It includes routines to:
    !> 1. Calculate the rate of ionizing photons (Q_H) produced by the stellar population.
    !> 2. Attenuate the stellar UV flux due to absorption by neutral hydrogen.
    !> 3. Interpolate pre-computed CLOUDY photoionization grid models for
    !>    continuum and line emission.
    !> 4. Add smoothed emission lines and nebular continuum to the output spectrum.
    !>
    !> AD_NOTE: This module contains table lookups and discrete index finding
    !> (`find_interval`), which create discontinuities. Special care is required
    !> for Automatic Differentiation (AD).

    use fsps_precision, only: WP
    use fsps_constants, only: NEMLINE, NEBNAGE, NEBNZ, NEBNIP, C_LIGHT, &
                              H_PLANCK, L_SOL, SAFE_FLOOR, PI
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: find_interval
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none
    private

    ! Public Interface
    public :: apply_nebular_emission
    public :: process_ionizing_radiation
    public :: interpolate_zu_slice_point
    public :: get_grid_indices_weights

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------
    real(WP), parameter :: SQRT_2_PI  = sqrt(2.0_wp * PI)
    real(WP), parameter :: KM_S_TO_ANGSTROM_FACTOR = 1.0e13_wp / C_LIGHT

contains

    !> @brief
    !> Main driver for adding nebular emission to the stellar spectrum.
    !>
    !> @details
    !> Iterates over the age of the stellar population. For each age step:
    !> 1. Calculates the number of ionizing photons (Q) that are absorbed by gas.
    !> 2. Interpolates the nebular grid (CLOUDY) for the given metallicity, 
    !>    ionization parameter, and age.
    !> 3. Adds the scaled nebular continuum and emission lines to the spectrum.
    !>
    !> @param[inout] ctx         FSPS context (contains CLOUDY grids and wavelength arrays).
    !> @param[in]    pset        FSPS parameter structure (gas_logz, gas_logu, etc.).
    !> @param[in]    sspi        Input Stellar Spectrum (L_sol/Hz).
    !> @param[inout] sspo        Output Spectrum (Stellar + Nebular).
    !> @param[out]   nebemline   (Optional) Output array for individual line luminosities.
    subroutine apply_nebular_emission(ctx, pset, sspi, sspo, nebemline)

        type(fsps_context_t), intent(inout)       :: ctx
        type(params), intent(in)                  :: pset
        real(WP), dimension(:,:), intent(in)      :: sspi
        real(WP), dimension(:,:), intent(inout)   :: sspo
        real(WP), dimension(:,:), intent(inout), optional :: nebemline

        ! Loop variables
        integer :: t, k, max_neb_time_idx
        integer :: nspec
        
        ! Grid Interpolation Indices & Weights
        integer :: idx_z, idx_u, idx_a
        real(WP) :: w_z, w_u, w_a
        
        ! Physics variables
        real(WP) :: q_ionizing
        logical :: use_xrb_grid
        logical :: calc_lines, calc_cont

        ! Temporary "Reduced" Grids
        ! We use allocatable arrays but manage them on device.
        real(WP), allocatable, dimension(:,:) :: neb_cont_grid_reduced ! (n_wave, n_age_grid)
        real(WP), allocatable, dimension(:,:) :: neb_line_grid_reduced ! (n_lines, n_age_grid)

        ! Buffers for the current time step. 
        ! We use !acc enter data create for these inside loop or allocate once.
        real(WP), allocatable, dimension(:) :: current_step_cont
        real(WP), allocatable, dimension(:) :: current_step_lines_log

        nspec = size(sspi, 1)

        !$acc data pcopyin(sspi) pcopy(sspo)
        !$acc update device(sspi)

        ! 0. Initialization & Validation
        ! ------------------------------
        !$acc kernels present(sspo, sspi)
        sspo = sspi
        !$acc end kernels
        
        if (present(nebemline)) then
            nebemline = 0.0_wp
        end if

        ! Determine what needs calculating
        calc_cont = (ctx%add_neb_continuum_val == 1)
        calc_lines = (ctx%nebemlineinspec_val == 1 .or. present(nebemline))
        
        ! Logic flag for XRB (BPSS models)
        use_xrb_grid = (ctx%state%isoc_type == 'bpss') .and. (ctx%add_xrb_emission_val == 1)

        ! Ensure Gaussian smoothing kernels are ready (Run on Host or Device?)
        ! compute_line_gaussians is likely expensive.
        if (ctx%setup_nebular_gaussians_val == 0 .and. ctx%nebemlineinspec_val == 1) then
            call compute_line_gaussians(ctx, pset)
        end if

        ! 1. Pre-calculate Interpolation Weights for Z and U
        ! --------------------------------------------------
        ! These are constant for the entire subroutine call.
        call get_grid_indices_weights(ctx%state%nebem_logz, pset%gas_logz, NEBNZ, idx_z, w_z)
        call get_grid_indices_weights(ctx%state%nebem_logu, pset%gas_logu, NEBNIP, idx_u, w_u)


        ! 2. Dimensionality Reduction
        ! ----------------------------------------------
        
        if (calc_cont) then
            allocate(neb_cont_grid_reduced(nspec, NEBNAGE))
            !$acc enter data create(neb_cont_grid_reduced)

            do k = 1, NEBNAGE
                do t = 1, nspec
                    ! Inline interpolate_zu_slice logic or use routine seq
                    ! interpolate_zu_slice needs to be routine seq if used here.
                    if (use_xrb_grid) then
                       neb_cont_grid_reduced(t, k) = interpolate_zu_slice_point(ctx%state%xnebem_cont, &
                            t, k, idx_z, idx_u, w_z, w_u)
                    else
                       neb_cont_grid_reduced(t, k) = interpolate_zu_slice_point(ctx%state%nebem_cont, &
                            t, k, idx_z, idx_u, w_z, w_u)
                    end if
                end do
            end do
        end if

        if (calc_lines) then
            allocate(neb_line_grid_reduced(NEMLINE, NEBNAGE))
            !$acc enter data create(neb_line_grid_reduced)

            do k = 1, NEBNAGE
                do t = 1, NEMLINE
                    if (use_xrb_grid) then
                       neb_line_grid_reduced(t, k) = interpolate_zu_slice_point(ctx%state%xnebem_line, &
                            t, k, idx_z, idx_u, w_z, w_u)
                    else
                       neb_line_grid_reduced(t, k) = interpolate_zu_slice_point(ctx%state%nebem_line, &
                            t, k, idx_z, idx_u, w_z, w_u)
                    end if
                end do
            end do
        end if

        ! 3. Main Time Loop
        ! -----------------
        max_neb_time_idx = find_interval(ctx%state%time_full, ctx%state%nebem_age(NEBNAGE))
        
        ! Scratch arrays
        allocate(current_step_cont(nspec))
        allocate(current_step_lines_log(NEMLINE))
        !$acc enter data create(current_step_cont, current_step_lines_log)

        ! NOTE: The loop over time steps 't' must be sequential because process_ionizing_radiation
        ! and integration might be heavy, and we are updating sspo(:, t).
        ! Parallelizing over 't' is possible if independent.
        ! But process_ionizing_radiation does integration.
        ! We will keep the loop sequential but run kernels inside.
        
        do t = 1, max_neb_time_idx
            
            ! A. Calculate Ionizing Photons
            call process_ionizing_radiation(ctx, pset, sspi(:,t), sspo(:,t), q_ionizing)

            ! Optimization: Skip expensive math if there is no ionizing radiation
            if (q_ionizing <= SAFE_FLOOR) cycle

            ! B. Interpolate Age (1D)
            ! -----------------------
            call get_grid_indices_weights(ctx%state%nebem_age, ctx%state%time_full(t), &
                                          NEBNAGE, idx_a, w_a)

            ! C. Add Continuum
            ! ----------------
            if (calc_cont) then
                !$acc parallel loop present(sspo, neb_cont_grid_reduced, current_step_cont) &
                !$acc               firstprivate(idx_a, w_a, q_ionizing, t)
                do k = 1, nspec
                    current_step_cont(k) = (1.0_wp - w_a) * neb_cont_grid_reduced(k, idx_a) + &
                                           (         w_a) * neb_cont_grid_reduced(k, idx_a + 1)
                    
                    sspo(k,t) = sspo(k,t) + (10.0_wp**current_step_cont(k)) * q_ionizing
                end do
            end if

            ! D. Add Lines
            ! ------------
            if (calc_lines) then
                !$acc parallel loop present(current_step_lines_log, neb_line_grid_reduced) &
                !$acc               firstprivate(idx_a, w_a)
                do k = 1, NEMLINE
                    current_step_lines_log(k) = (1.0_wp - w_a) * neb_line_grid_reduced(k, idx_a) + &
                                                (         w_a) * neb_line_grid_reduced(k, idx_a + 1)
                end do
                !$acc update host(current_step_lines_log)
                
                if (present(nebemline)) then
                    do k = 1, NEMLINE
                        nebemline(k,t) = (10.0_wp**current_step_lines_log(k)) * q_ionizing
                    end do
                end if

                if (ctx%nebemlineinspec_val == 1) then
                    call add_lines_to_spectrum(ctx, sspo(:,t), current_step_lines_log, q_ionizing)
                end if
            end if

        end do

        ! Cleanup
        !$acc exit data delete(current_step_cont, current_step_lines_log)
        deallocate(current_step_cont, current_step_lines_log)
        
        if (allocated(neb_cont_grid_reduced)) then
            !$acc exit data delete(neb_cont_grid_reduced)
            deallocate(neb_cont_grid_reduced)
        end if
        if (allocated(neb_line_grid_reduced)) then
            !$acc exit data delete(neb_line_grid_reduced)
            deallocate(neb_line_grid_reduced)
        end if

        !$acc end data

    end subroutine apply_nebular_emission

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> @brief
    !> Calculates the number of ionizing photons and attenuates the stellar UV.
    !>
    !> @details
    !> 1. Attenuates the output spectrum (`spec_out`) below the Lyman limit 
    !>    based on `frac_obrun` (fraction of photons that escape the nebula).
    !> 2. Integrates the *input* spectrum (`spec_in`) to find the total 
    !>    ionizing flux available.
    !> 3. Returns `q_val`, the number of photons *absorbed* by the gas.
    !>
    !> @param[in]  ctx       FSPS context (wavelength arrays).
    !> @param[in]  pset      Parameters (frac_obrun).
    !> @param[in]  spec_in   Intrinsic stellar spectrum.
    !> @param[out] spec_out  Attenuated stellar spectrum.
    !> @param[out] q_val     Number of ionizing photons absorbed.
    subroutine process_ionizing_radiation(ctx, pset, spec_in, spec_out, q_val)
        type(fsps_context_t), intent(in)    :: ctx
        type(params), intent(in)            :: pset
        real(WP), dimension(:), intent(in)  :: spec_in
        real(WP), dimension(:), intent(inout) :: spec_out
        real(WP), intent(out)               :: q_val
        
        integer :: whlylim
        real(WP) :: integral_flux
        integer :: i
        real(WP) :: y1, y2
        real(WP) :: frac_obrun_clamped
        
        whlylim = ctx%state%whlylim
        frac_obrun_clamped = max(0.0_wp, min(pset%frac_obrun, 1.0_wp))

        !$acc data pcopyin(spec_in) pcopy(spec_out)
        !$acc update device(spec_in)

        ! 1. Attenuate Output Spectrum (EUV < 912 A)
        if (whlylim > 0) then
            !$acc parallel loop present(ctx, spec_in, spec_out) firstprivate(frac_obrun_clamped)
            do i = 1, whlylim
                spec_out(i) = spec_in(i) * frac_obrun_clamped
            end do
        end if

        ! 2. Calculate Total Ionizing Photons (Q)
        if (whlylim < 2) then
            q_val = 0.0_wp
        else
            integral_flux = 0.0_wp
            ! Use specialized loop to avoid array temp (spec_in / spec_nu)
            !$acc parallel loop reduction(+:integral_flux) present(ctx, spec_in)
            do i = 1, whlylim - 1
                y1 = spec_in(i) / ctx%state%spec_nu(i)
                y2 = spec_in(i+1) / ctx%state%spec_nu(i+1)
                integral_flux = integral_flux + 0.5_wp * abs(ctx%state%spec_nu(i+1) - ctx%state%spec_nu(i)) * (y1 + y2)
            end do
            
            ! We check NaN on host? Or device?
            ! Can't check IEEE NaN easily on device across reduction.
            ! Assuming logic holds.
            q_val = abs((integral_flux / H_PLANCK * L_SOL) * (1.0_wp - pset%frac_obrun))
        end if

        !$acc end data

    end subroutine process_ionizing_radiation

    !> @brief
    !> Pre-computes Gaussian profiles for emission lines.
    !>
    !> @details
    !> Computes the normalized Gaussian profile `gaussnebarr` for each emission line,
    !> broadened by `sigma_smooth`.
    !>
    !> AD_NOTE: This modifies `ctx%state`. In a purely functional AD context, 
    !> this might need to return a local array instead.
    subroutine compute_line_gaussians(ctx, pset)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: pset
        
        integer :: i
        real(WP) :: sigma_angstroms, norm_factor
        real(WP), dimension(size(ctx%state%spec_lambda)) :: lambda_grid

        lambda_grid = ctx%state%spec_lambda

        do i = 1, NEMLINE
            ! Determine width in Angstroms
            if (ctx%smooth_velocity_val == 1) then
                ! velocity smoothing (sigma_smooth in km/s)
                sigma_angstroms = ctx%state%nebem_line_pos(i) * pset%sigma_smooth * KM_S_TO_ANGSTROM_FACTOR
            else
                ! wavelength smoothing (sigma_smooth in Angstroms)
                sigma_angstroms = pset%sigma_smooth
            end if

            ! Enforce minimum resolution (Nyquist sampling or instrumental limit)
            ! The factor of 2 ensures we don't alias on the grid.
            sigma_angstroms = max(sigma_angstroms, ctx%state%neb_res_min(i) * 2.0_wp)

            ! Compute Gaussian
            ! Profile = (1 / (sqrt(2pi)*sigma)) * exp( -0.5 * ((lam - lam_0)/sigma)^2 )
            ! We also convert luminosity units here if needed, but original code divides by C_LIGHT*lambda^2?
            ! Original: ... / C_LIGHT * nebem_line_pos(i)**2
            ! This converts from L_lambda to L_nu? 
            ! Yes: L_nu = L_lambda * lambda^2 / c. 
            norm_factor = (1.0_wp / (SQRT_2_PI * sigma_angstroms)) * &
                          (ctx%state%nebem_line_pos(i)**2 / C_LIGHT)

            ctx%state%gaussnebarr(:,i) = norm_factor * &
                exp( -0.5_wp * ((lambda_grid - ctx%state%nebem_line_pos(i)) / sigma_angstroms)**2 )
        end do
    end subroutine compute_line_gaussians

    !> @brief
    !> Adds emission lines to the spectrum.
    subroutine add_lines_to_spectrum(ctx, spectrum, line_lum_log, q_val)
        type(fsps_context_t), intent(in)    :: ctx
        real(WP), dimension(:), intent(inout) :: spectrum
        real(WP), dimension(:), intent(in)    :: line_lum_log
        real(WP), intent(in)                  :: q_val
        
        integer :: i, j
        real(WP) :: sum_val

        ! Manual Matmul
        ! spectrum(j) = sum(gauss(j, i) * flux(i))
        !$acc parallel loop gang vector present(ctx, spectrum, line_lum_log) private(sum_val)
        do j = 1, size(spectrum)
            sum_val = 0.0_wp
            do i = 1, NEMLINE
                sum_val = sum_val + ctx%state%gaussnebarr(j, i) * (10.0_wp**line_lum_log(i)) * q_val
            end do
            spectrum(j) = spectrum(j) + sum_val
        end do

    end subroutine add_lines_to_spectrum

    !> @brief
    !> Finds grid index and linear interpolation weight for a single value.
    !>
    !> @details
    !> AD_NOTE: `find_interval` introduces a discontinuity in the derivative 
    !> of the index w.r.t the input value. The weight `w` is differentiable.
    subroutine get_grid_indices_weights(grid, val, n_grid, idx, weight)
        real(WP), dimension(:), intent(in) :: grid
        real(WP), intent(in)               :: val
        integer, intent(in)                :: n_grid
        integer, intent(out)               :: idx
        real(WP), intent(out)              :: weight
        
        ! Find index in sorted array
        idx = find_interval(grid, val)
        
        ! Clamp to valid range [1, N-1] for interpolation safety
        idx = max(1, min(idx, n_grid - 1))
        
        ! Calculate weight (0.0 to 1.0)
        weight = (val - grid(idx)) / (grid(idx+1) - grid(idx))
        
        ! Clamp weight to avoid extrapolation
        weight = max(0.0_wp, min(weight, 1.0_wp))
    end subroutine get_grid_indices_weights

    !> @brief
    !> Scalar version of interpolate_zu_slice for use inside parallel loops.
    !> Returns single value at index `i_wave`.
    pure function interpolate_zu_slice_point(grid, i_wave, idx_age, idx_z, idx_u, w_z, w_u) result(val)
        !$acc routine seq
        real(WP), dimension(:,:,:,:), intent(in) :: grid
        integer, intent(in)  :: i_wave, idx_age, idx_z, idx_u
        real(WP), intent(in) :: w_z, w_u
        real(WP) :: val
        
        real(WP) :: c00, c01, c10, c11

        c00 = (1.0_wp - w_z) * (1.0_wp - w_u)
        c01 = (1.0_wp - w_z) * (         w_u)
        c10 = (         w_z) * (1.0_wp - w_u)
        c11 = (         w_z) * (         w_u)

        val = c00 * grid(i_wave, idx_z,   idx_age, idx_u  ) + &
              c01 * grid(i_wave, idx_z,   idx_age, idx_u+1) + &
              c10 * grid(i_wave, idx_z+1, idx_age, idx_u  ) + &
              c11 * grid(i_wave, idx_z+1, idx_age, idx_u+1)

    end function interpolate_zu_slice_point

end module fsps_gas