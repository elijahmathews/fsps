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

    use fsps_types, only: sp, params, nemline, nebnage, nebnz, nebnip, &
                          clight, mypi, hplank, lsun
    use fsps_context_types, only: fsps_context_t
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: find_interval
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan

    implicit none
    private

    ! Public Interface
    public :: apply_nebular_emission
    public :: process_ionizing_radiation
    public :: interpolate_zu_slice
    public :: get_grid_indices_weights

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------
    real(sp), parameter :: SAFE_FLOOR = tiny(0.0_sp)
    real(sp), parameter :: PI = acos(-1.0_sp)
    real(sp), parameter :: SQRT_2_PI  = sqrt(2.0_sp * PI)
    real(sp), parameter :: KM_S_TO_ANGSTROM_FACTOR = 1.0e13_sp / clight

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
        real(sp), dimension(:,:), intent(in)      :: sspi
        real(sp), dimension(:,:), intent(inout)   :: sspo
        real(sp), dimension(:,:), intent(inout), optional :: nebemline

        ! Loop variables
        integer :: t, k, max_neb_time_idx
        
        ! Grid Interpolation Indices & Weights
        integer :: idx_z, idx_u, idx_a
        real(sp) :: w_z, w_u, w_a
        
        ! Physics variables
        real(sp) :: q_ionizing
        logical :: use_xrb_grid
        logical :: calc_lines, calc_cont

        ! Temporary "Reduced" Grids
        ! We collapse the 4D grid (Wave, Z, Age, U) -> 2D (Wave, Age)
        ! These hold the spectrum/lines for the specific Z and U of this call,
        ! for every Age in the original grid.
        real(sp), allocatable, dimension(:,:) :: neb_cont_grid_reduced ! (n_wave, n_age_grid)
        real(sp), allocatable, dimension(:,:) :: neb_line_grid_reduced ! (n_lines, n_age_grid)

        ! Buffers for the current time step (interpolated from the reduced grids)
        ! If this winds up being a performance bottleneck, we can pre-allocate
        ! these as permanent arrays in the context.
        real(sp), dimension(size(sspi, 1)) :: current_step_cont
        real(sp), dimension(nemline)       :: current_step_lines_log

        ! 0. Initialization & Validation
        ! ------------------------------
        sspo = sspi
        if (present(nebemline)) nebemline = 0.0_sp

        ! Determine what needs calculating
        calc_cont = (ctx%add_neb_continuum_val == 1)
        calc_lines = (ctx%nebemlineinspec_val == 1 .or. present(nebemline))
        
        ! Logic flag for XRB (BPSS models)
        use_xrb_grid = (ctx%state%isoc_type == 'bpss') .and. (ctx%add_xrb_emission_val == 1)

        ! Ensure Gaussian smoothing kernels are ready
        if (ctx%setup_nebular_gaussians_val == 0 .and. ctx%nebemlineinspec_val == 1) then
            call compute_line_gaussians(ctx, pset)
        end if

        ! 1. Pre-calculate Interpolation Weights for Z and U
        ! --------------------------------------------------
        ! These are constant for the entire subroutine call.
        call get_grid_indices_weights(ctx%state%nebem_logz, pset%gas_logz, nebnz, idx_z, w_z)
        call get_grid_indices_weights(ctx%state%nebem_logu, pset%gas_logu, nebnip, idx_u, w_u)


        ! 2. Dimensionality Reduction
        ! ----------------------------------------------
        ! Instead of doing 3D interpolation (Z, Age, U) inside the time loop,
        ! we collapse Z and U first.
        
        if (calc_cont) then
            allocate(neb_cont_grid_reduced(size(sspi, 1), nebnage))
            do k = 1, nebnage
                ! Interpolates Z and U planes for the specific Age slice 'k'
                ! Implementation note: Helper function `interpolate_zu_slice` handles
                ! the 2D interpolation for a fixed age index.
                if (use_xrb_grid) then
                    call interpolate_zu_slice(ctx%state%xnebem_cont, &
                        neb_cont_grid_reduced(:, k), k, idx_z, idx_u, w_z, w_u)
                else
                    call interpolate_zu_slice(ctx%state%nebem_cont, &
                        neb_cont_grid_reduced(:, k), k, idx_z, idx_u, w_z, w_u)
                end if
            end do
        end if

        if (calc_lines) then
            allocate(neb_line_grid_reduced(nemline, nebnage))
            do k = 1, nebnage
                if (use_xrb_grid) then
                    call interpolate_zu_slice(ctx%state%xnebem_line, &
                        neb_line_grid_reduced(:, k), k, idx_z, idx_u, w_z, w_u)
                else
                    call interpolate_zu_slice(ctx%state%nebem_line, &
                        neb_line_grid_reduced(:, k), k, idx_z, idx_u, w_z, w_u)
                end if
            end do
        end if

        ! 3. Main Time Loop
        ! -----------------
        ! Only calculate up to the max age supported by the nebular grid
        max_neb_time_idx = find_interval(ctx%state%time_full, ctx%state%nebem_age(nebnage))

        do t = 1, max_neb_time_idx
            
            ! A. Calculate Ionizing Photons
            call process_ionizing_radiation(ctx, pset, sspi(:,t), sspo(:,t), q_ionizing)

            ! Optimization: Skip expensive math if there is no ionizing radiation
            if (q_ionizing <= SAFE_FLOOR) cycle

            ! B. Interpolate Age (1D)
            ! -----------------------
            ! Now we simply interpolate our pre-reduced grids along the Age axis.
            call get_grid_indices_weights(ctx%state%nebem_age, ctx%state%time_full(t), &
                                          nebnage, idx_a, w_a)

            ! C. Add Continuum
            ! ----------------
            if (calc_cont) then
                ! 1D Linear Interpolation: (1-w)*grid(idx) + w*grid(idx+1)
                current_step_cont = (1.0_sp - w_a) * neb_cont_grid_reduced(:, idx_a) + &
                                    (         w_a) * neb_cont_grid_reduced(:, idx_a + 1)
                
                ! Add to spectrum
                sspo(:,t) = sspo(:,t) + (10.0_sp**current_step_cont) * q_ionizing
            end if

            ! D. Add Lines
            ! ------------
            if (calc_lines) then
                ! 1D Linear Interpolation
                current_step_lines_log = (1.0_sp - w_a) * neb_line_grid_reduced(:, idx_a) + &
                                         (         w_a) * neb_line_grid_reduced(:, idx_a + 1)
                
                ! Store raw luminosities if requested
                if (present(nebemline)) then
                    nebemline(:,t) = (10.0_sp**current_step_lines_log) * q_ionizing
                end if

                ! Add to spectrum (Using MATMUL optimization)
                if (ctx%nebemlineinspec_val == 1) then
                    call add_lines_to_spectrum(ctx, sspo(:,t), current_step_lines_log, q_ionizing)
                end if
            end if

        end do

        ! Cleanup
        if (allocated(neb_cont_grid_reduced)) deallocate(neb_cont_grid_reduced)
        if (allocated(neb_line_grid_reduced)) deallocate(neb_line_grid_reduced)

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
        real(sp), dimension(:), intent(in)  :: spec_in
        real(sp), dimension(:), intent(inout) :: spec_out
        real(sp), intent(out)               :: q_val
        
        integer :: whlylim
        real(sp) :: integral_flux
        
        whlylim = ctx%state%whlylim

        ! 1. Attenuate Output Spectrum (EUV < 912 A)
        ! ------------------------------------------
        ! frac_obrun is the fraction of "runaway" stars or leakage.
        ! These photons escape the HII region without processing.
        ! Conversely, (1 - frac_obrun) are absorbed.
        if (whlylim > 0) then
            spec_out(1:whlylim) = spec_in(1:whlylim) * max(0.0_sp, min(pset%frac_obrun, 1.0_sp))
        end if

        ! 2. Calculate Total Ionizing Photons (Q)
        ! ---------------------------------------
        ! Q = Integral(L_nu / h*nu) d_nu
        if (whlylim < 2) then
            q_val = 0.0_sp
        else
            ! Note: spec_nu is frequency. spec_in is L_sol/Hz.
            ! Result is photons/sec scaled by lsun/hplank.
            integral_flux = integrate_trapezoid_array( &
                ctx%state%spec_nu(:whlylim), &
                spec_in(:whlylim) / ctx%state%spec_nu(:whlylim) &
            )
            
            if (ieee_is_nan(integral_flux)) then
                q_val = 0.0_sp
            else
                q_val = (integral_flux / hplank * lsun) * (1.0_sp - pset%frac_obrun)
            end if
        end if

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
        real(sp) :: sigma_angstroms, norm_factor
        real(sp), dimension(size(ctx%state%spec_lambda)) :: lambda_grid

        lambda_grid = ctx%state%spec_lambda

        do i = 1, nemline
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
            sigma_angstroms = max(sigma_angstroms, ctx%state%neb_res_min(i) * 2.0_sp)

            ! Compute Gaussian
            ! Profile = (1 / (sqrt(2pi)*sigma)) * exp( -0.5 * ((lam - lam_0)/sigma)^2 )
            ! We also convert luminosity units here if needed, but original code divides by clight*lambda^2?
            ! Original: ... / clight * nebem_line_pos(i)**2
            ! This converts from L_lambda to L_nu? 
            ! Yes: L_nu = L_lambda * lambda^2 / c. 
            norm_factor = (1.0_sp / (SQRT_2_PI * sigma_angstroms)) * &
                          (ctx%state%nebem_line_pos(i)**2 / clight)

            ctx%state%gaussnebarr(:,i) = norm_factor * &
                exp( -0.5_sp * ((lambda_grid - ctx%state%nebem_line_pos(i)) / sigma_angstroms)**2 )
        end do
    end subroutine compute_line_gaussians

    !> @brief
    !> Adds emission lines to the spectrum using matrix multiplication.
    !> This is significantly faster than looping over lines due to reduced memory I/O.
    subroutine add_lines_to_spectrum(ctx, spectrum, line_lum_log, q_val)
        type(fsps_context_t), intent(in)    :: ctx
        real(sp), dimension(:), intent(inout) :: spectrum
        real(sp), dimension(:), intent(in)    :: line_lum_log
        real(sp), intent(in)                  :: q_val
        
        real(sp), dimension(nemline) :: line_flux_linear

        ! Vectorize the log -> linear conversion
        line_flux_linear = (10.0_sp**line_lum_log) * q_val
        
        ! Perform Matrix-Vector multiplication: 
        ! [Lambda x Lines] * [Lines] = [Lambda]
        ! This sums all Gaussian profiles weighted by their flux in one pass.
        spectrum = spectrum + matmul(ctx%state%gaussnebarr, line_flux_linear)

    end subroutine add_lines_to_spectrum

    !> @brief
    !> Finds grid index and linear interpolation weight for a single value.
    !>
    !> @details
    !> AD_NOTE: `find_interval` introduces a discontinuity in the derivative 
    !> of the index w.r.t the input value. The weight `w` is differentiable.
    subroutine get_grid_indices_weights(grid, val, n_grid, idx, weight)
        real(sp), dimension(:), intent(in) :: grid
        real(sp), intent(in)               :: val
        integer, intent(in)                :: n_grid
        integer, intent(out)               :: idx
        real(sp), intent(out)              :: weight
        
        ! Find index in sorted array
        idx = find_interval(grid, val)
        
        ! Clamp to valid range [1, N-1] for interpolation safety
        idx = max(1, min(idx, n_grid - 1))
        
        ! Calculate weight (0.0 to 1.0)
        weight = (val - grid(idx)) / (grid(idx+1) - grid(idx))
        
        ! Clamp weight to avoid extrapolation
        weight = max(0.0_sp, min(weight, 1.0_sp))
    end subroutine get_grid_indices_weights

    !> @brief
    !> Interpolates a 2D slice (Metallicity Z, Ionization U) from the 4D grid
    !> for a fixed Age index.
    !>
    !> @details
    !> This reduces the 4D grid (Wave, Z, Age, U) down to a 1D array (Wave)
    !> for the specific Z, U, and Age parameters provided.
    !>
    !> Logic:
    !> Res = (1-wz)(1-wu)*V00 + (1-wz)(wu)*V01 + (wz)(1-wu)*V10 + (wz)(wu)*V11
    !>
    !> @param[in] grid     The 4D nebular grid (Wave/Line, Z, Age, U)
    !> @param[in] idx_age  The fixed Age index to slice at
    !> @param[in] idx_z    Lower index for Metallicity
    !> @param[in] idx_u    Lower index for Ionization Parameter
    !> @param[in] w_z      Interpolation weight for Z
    !> @param[in] w_u      Interpolation weight for U
    !> @return             1D array of interpolated values (wavelengths or lines)
    pure subroutine interpolate_zu_slice(grid, res, idx_age, idx_z, idx_u, w_z, w_u)
        real(sp), dimension(:,:,:,:), intent(in) :: grid
        real(sp), dimension(:), intent(out)    :: res
        integer, intent(in)  :: idx_age, idx_z, idx_u
        real(sp), intent(in) :: w_z, w_u
        
        ! Local coefficients for bilinear interpolation
        real(sp) :: c00, c01, c10, c11

        ! Pre-calculate coefficients
        ! This avoids re-calculating (1-w_z) etc. for every wavelength point
        c00 = (1.0_sp - w_z) * (1.0_sp - w_u)
        c01 = (1.0_sp - w_z) * (         w_u)
        c10 = (         w_z) * (1.0_sp - w_u)
        c11 = (         w_z) * (         w_u)

        ! Vectorized Array Operation
        ! Fortran arrays are column-major. Since the first dimension (:) is 
        ! contiguous, this operation is highly cache-efficient and vectorizable.
        res = c00 * grid(:, idx_z,   idx_age, idx_u  ) + &
              c01 * grid(:, idx_z,   idx_age, idx_u+1) + &
              c10 * grid(:, idx_z+1, idx_age, idx_u  ) + &
              c11 * grid(:, idx_z+1, idx_age, idx_u+1)

    end subroutine interpolate_zu_slice

end module fsps_gas