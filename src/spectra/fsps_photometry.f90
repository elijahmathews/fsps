
module fsps_photometry
    !> @brief
    !> Module for calculating photometric magnitudes (AB or Vega) from spectra.
    !>
    !> @details
    !> This module replaces the legacy `spec_mags.f90`. It handles:
    !> 1. Redshifting the rest-frame spectrum to the observed frame.
    !> 2. Applying cosmological distance dimming (Distance Modulus).
    !> 3. Integrating the spectrum over filter transmission curves.
    !> 4. Converting fluxes to AB or Vega magnitudes.

    use fsps_precision, only: WP
    use fsps_constants, only: SAFE_FLOOR, ABS_MAG_ZEROPOINT_LOG
    use fsps_context_types, only: fsps_context_t
    use fsps_interpolation, only: interpolate_linear
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_special_functions, only: mag_from_flux
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan

    implicit none
    private

    public :: compute_magnitudes

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------
    real(WP), parameter :: AB_ZEROPOINT   = 48.60_wp   !< CGS to AB Mag constant
    real(WP), parameter :: LOG10_FACTOR   = 2.5_wp     !< Magnitude scaling factor
    integer, parameter  :: IDX_V_BAND     = 1          !< Index of V-band (for Vega norm)

contains

    !> @brief
    !> Calculates magnitudes for a given spectrum and redshift.
    !>
    !> @details
    !> The routine first shifts the input spectrum (Lsun/Hz) to the observed frame
    !> if z > 0. It then integrates over the transmission curves defined in
    !> `ctx%state%bands`.
    !>
    !> @param[in]    ctx          The FSPS context (contains filter curves, cosmology).
    !> @param[in]    zred         Redshift of the galaxy.
    !> @param[in]    spec         Rest-frame spectrum (Lsun/Hz).
    !> @param[out]   mags         Output magnitudes array.
    !> @param[in]    mag_compute  (Optional) Integer mask (1=compute, 0=skip) for bands.
    subroutine compute_magnitudes(ctx, zred, spec, mags, mag_compute)
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: zred
        real(WP), dimension(:), intent(in), contiguous :: spec
        real(WP), dimension(:), intent(out), contiguous :: mags
        integer, dimension(:), intent(in), optional :: mag_compute

        ! Local variables
        integer :: i, n_spec, n_bands
        integer, allocatable :: work_flags(:)
        real(WP), allocatable :: obs_frame_spec(:)
        real(WP), allocatable :: flux_over_lambda(:)
        real(WP) :: dist_mod_term, integrated_flux

        ! Context shortcuts
        logical :: do_vega, do_light_ages

        ! 1. Setup and Validation
        n_spec = size(spec)
        n_bands = size(mags)

        if (n_spec < 2) then
            mags = get_quiet_nan()
            return
        end if

        ! Initialize output
        mags = get_quiet_nan()

        ! Parse Context Flags
        do_vega = (ctx%compute_vega_mags_val == 1)
        do_light_ages = (ctx%compute_light_ages_val == 1)

        ! 2. Determine which bands to compute
        allocate(work_flags(n_bands))
        if (present(mag_compute)) then
            work_flags = mag_compute
            ! If calculating Vega mags, we MUST calculate V-band (index 1) for normalization
            if (do_vega) work_flags(IDX_V_BAND) = 1
            
            ! Optimization: Early exit if no bands requested
            if (all(work_flags == 0)) return
        else
            work_flags = 1 ! Default: compute all
        end if

        ! 3. Prepare Spectrum (Redshift & Distance Modulus)
        allocate(obs_frame_spec(n_spec))
        
        call prepare_observed_spectrum(ctx, zred, spec, obs_frame_spec, dist_mod_term)

        ! 4. Pre-calculate terms for Integration
        ! The integral is int(F * T * dlam / lam).
        ! We pre-calculate (F / lam) here to save N_bands divisions later.
        allocate(flux_over_lambda(n_spec))
        where (ctx%state%spec_lambda > tiny(0.0_wp))
            flux_over_lambda = obs_frame_spec / ctx%state%spec_lambda
        elsewhere
            flux_over_lambda = 0.0_wp
        end where

        ! 5. Filter Integration Loop
        ! Note: We cannot vectorize the loop over 'i' easily because 'bands' is large
        ! and we want to use the optimized `integrate_trapezoid_array`.
        do i = 1, n_bands
            if (work_flags(i) == 0) cycle

            ! Integrate: Trapz(x=lambda, y = (Flux/lambda) * Transmission)
            ! Note: ctx%state%bands is shape (n_lambda, n_bands)
            integrated_flux = integrate_trapezoid_array( &
                                ctx%state%spec_lambda, &
                                flux_over_lambda * ctx%state%bands(:, i) &
                              )

            if (.not. do_light_ages) then
                ! Convert Flux to AB Magnitude
                ! Mag = -2.5*log10(F) - 48.60 - ZeroPointCorrection + DistanceModulus
                mags(i) = mag_from_flux(integrated_flux) - &
                            AB_ZEROPOINT - &
                            (LOG10_FACTOR * ABS_MAG_ZEROPOINT_LOG) + &
                            dist_mod_term
            else
                ! For light ages, we just want the raw weight (or NaN if invalid)
                if (integrated_flux > SAFE_FLOOR) then
                    mags(i) = integrated_flux
                end if
                ! else: remains NaN from initialization
            end if
        end do

        ! 6. Vega System Correction
        ! Applies offset relative to V-band if requested.
        ! Formula: M_i = m_i - v_i + v_1 (Derived from legacy logic)
        if (do_vega .and. .not. do_light_ages) then
            ! We only need to check the Anchor (V-Band)
            if (.not. ieee_is_nan(mags(IDX_V_BAND))) then
                
                ! Vectorized update!
                ! Valid bands become (Mag - Vega + Vega_Ref)
                ! NaN bands stay NaN (because NaN - Vega = NaN)
                mags(2:n_bands) = mags(2:n_bands) - &
                                  ctx%state%magvega(2:n_bands) + &
                                  ctx%state%magvega(IDX_V_BAND)
            else
                ! If V-band is invalid, the whole Vega normalization is impossible.
                mags = get_quiet_nan()
            end if
        end if

    end subroutine compute_magnitudes

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPERS
    ! ------------------------------------------------------------------------

    !> @brief
    !> Handles redshifting and cosmology calculations.
    !> 
    !> @details
    !> If z > 0, interpolates the spectrum to the observed frame and calculates
    !> the distance modulus term.
    subroutine prepare_observed_spectrum(ctx, z, spec_rest, spec_obs, dist_term)
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: z
        real(WP), dimension(:), intent(in) :: spec_rest
        real(WP), dimension(:), intent(out) :: spec_obs
        real(WP), intent(out) :: dist_term

        real(WP) :: dm, z_factor, target_lam_rest
        integer :: k, idx

        if (abs(z) > SAFE_FLOOR) then
            ! --- High Redshift Case ---

            ! OPTIMIZATION: 
            ! 1. Map observed wavelength BACK to rest frame: lam_rest = lam_obs / (1+z)
            !    This avoids allocating a temporary 'source_grid' array.
            ! 2. Use monotonic 'hunting' (tracking 'idx') instead of binary search.
            !    This reduces complexity from O(N log N) to O(N).
            
            idx = 1
            z_factor = 1.0_wp + z
            
            do k = 1, size(spec_obs)
                ! The rest-frame wavelength corresponding to this observed grid point
                target_lam_rest = ctx%state%spec_lambda(k) / z_factor

                ! Hunt for the interval: only move forward
                ! We stop when spec_lambda(idx+1) is just above our target
                do while (idx < size(spec_rest) - 1)
                    if (ctx%state%spec_lambda(idx+1) >= target_lam_rest) exit
                    idx = idx + 1
                end do

                ! Manual Linear Interpolation
                ! y = y1 + (x - x1) * slope
                spec_obs(k) = spec_rest(idx) + &
                              (target_lam_rest - ctx%state%spec_lambda(idx)) * &
                              (spec_rest(idx+1) - spec_rest(idx)) / &
                              (ctx%state%spec_lambda(idx+1) - ctx%state%spec_lambda(idx))
                
                ! Physical Constraints
                if (ieee_is_nan(spec_obs(k))) spec_obs(k) = 0.0_wp
                spec_obs(k) = max(spec_obs(k), 0.0_wp)
            end do

            ! 2. Cosmology / Distance Modulus
            dm = interpolate_linear(ctx%state%cosmospl(:,1), &
                                    ctx%state%cosmospl(:,3), &
                                    z)

            if (ieee_is_nan(dm) .or. dm <= SAFE_FLOOR) then
                dist_term = 0.0_wp
            else
                dist_term = (5.0_wp * log10(dm / 10.0_wp)) + & 
                            mag_from_flux(z_factor) ! Bandwidth stretching
            end if

        else
            ! --- Zero Redshift Case ---
            spec_obs  = spec_rest
            dist_term = 0.0_wp
        end if
        
    end subroutine prepare_observed_spectrum

    !> @brief Helper to generate a Quiet NaN
    pure function get_quiet_nan() result(val)
        real(WP) :: val
        val = ieee_value(0.0_wp, ieee_quiet_nan)
    end function get_quiet_nan

end module fsps_photometry