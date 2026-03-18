module fsps_imf
    !> @brief
    !> Handles Initial Mass Function (IMF) definitions and weighting.
    !>
    !> @details
    !> Contains routines to evaluate the IMF (dn/dM) for various standard
    !> parameterizations (Salpeter, Chabrier, Kroupa, etc.) and to compute
    !> the mass weights for SSP generation.
    
    use fsps_precision, only: WP
    use fsps_context_types, only: fsps_context_t
    use fsps_integration, only: integrate_romberg
    implicit none

    private

    public :: compute_imf_weights
    public :: get_imf_value
    public :: integrate_imf_interval_analytic
    public :: CHAB_MC, CHAB_SIGMA2, CHAB_IND

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------

    ! Chabrier 2003
    real(WP), parameter :: CHAB_MC     = 0.08_wp
    real(WP), parameter :: CHAB_SIGMA2 = 0.69_wp**2  ! Explicit squaring
    real(WP), parameter :: CHAB_IND    = 1.3_wp
    real(WP), parameter :: CHAB_LOG_MC = log10(CHAB_MC)
    real(WP), parameter :: CHAB_INV_2SIG2 = 1.0_wp / (2.0_wp * CHAB_SIGMA2)
    real(WP), parameter :: CHAB_HIGH_NORM = exp(-(CHAB_LOG_MC * CHAB_LOG_MC) * CHAB_INV_2SIG2)

    ! van Dokkum 2008
    real(WP), parameter :: VD_SIGMA2 = 0.69_wp**2
    real(WP), parameter :: VD_AH     = 0.0443_wp
    real(WP), parameter :: VD_IND    = 1.3_wp
    real(WP), parameter :: VD_AL     = 0.14_wp
    real(WP), parameter :: VD_NC     = 25.0_wp

contains

    !> @brief
    !> Computes the number of stars in each mass bin, normalized to 1 Msun total mass.
    !>
    !> @details
    !> 1. Calculates the number of stars in each mass bin i (integrating dn/dm).
    !> 2. Calculates the total mass of the system (integrating m*dn/dm).
    !> 3. Normalizes the weights so the total mass of the population is 1.0.
    !>
    !> @param[inout] ctx    Simulation context (contains IMF params).
    !> @param[in]    mini   Initial masses of the isochrone points.
    !> @param[out]   weights Output weights (Number of stars per unit total mass).
    !> @param[in]    nmass  Number of mass points.
    subroutine compute_imf_weights(ctx, mini, weights, nmass)
        type(fsps_context_t), intent(in) :: ctx
        real(WP), dimension(:), intent(in) :: mini
        real(WP), dimension(:), intent(out) :: weights
        integer, intent(in) :: nmass

        real(WP) :: m1, m2, total_mass
        integer :: i, imf_type_base
        logical :: use_analytic

        ! Access limits from state
        real(WP) :: lower_limit, upper_limit, lower_bound

        lower_limit = ctx%state%imf_lower_limit
        upper_limit = ctx%state%imf_upper_limit
        lower_bound = ctx%state%imf_lower_bound
        imf_type_base = mod(ctx%imf_type_val, 10)

        select case (imf_type_base)
        case (0, 2, 4, 5)
            use_analytic = .true.
        case default
            use_analytic = .false.
        end select

        weights = 0.0_wp

        ! 1. Calculate Number of Stars in each bin
        do i = 1, nmass
            
            ! Skip bins outside the requested IMF integration range
            if (mini(i) < lower_limit .or. mini(i) > upper_limit) cycle

            ! Define Bin Edges
            if (i == 1) then
                ! Special case for first bin (Geneva models)
                m1 = lower_bound
            else
                m1 = mini(i) - 0.5_wp * (mini(i) - mini(i-1))
            end if

            if (i == nmass) then
                m2 = mini(i)
            else
                m2 = mini(i) + 0.5_wp * (mini(i+1) - mini(i))
            end if

            ! Validate Monotonicity
            if (m2 < m1) then
                write(*, '("[FSPS-IMF] Warning: Non-monotonic mass bin at index ", I0, &
                         & " (m1=", ES10.3, ", m2=", ES10.3, ")")') i, m1, m2
                cycle
            end if
            
            ! Skip empty bins
            if (m2 == m1) cycle

            ! Integrate dn/dM over [m1, m2]
            ! Uses the wrapper_imf_count to integrate number density
            if (use_analytic) then
                weights(i) = integrate_imf_interval_analytic(ctx, m1, m2, mass_weighted=.false.)
            else
                weights(i) = integrate_romberg(ctx, wrapper_imf_count, m1, m2)
            end if
        end do

        ! 2. Calculate Total Mass (Normalization Factor)
        ! Integrate M * dn/dM over the full range [lower_limit, upper_limit]
        if (use_analytic) then
            total_mass = integrate_imf_interval_analytic(ctx, lower_limit, upper_limit, mass_weighted=.true.)
        else
            total_mass = integrate_romberg(ctx, wrapper_imf_mass, lower_limit, upper_limit)
        end if

        ! 3. Normalize
        if (total_mass > 0.0_wp) then
            weights = weights / total_mass
        end if

    end subroutine compute_imf_weights

    ! ------------------------------------------------------------------------
    ! INTEGRATION WRAPPERS
    ! These match the interface required by fsps_integration
    ! ------------------------------------------------------------------------

    !> @brief Wrapper to integrate Number Density (dn/dM)
    pure function wrapper_imf_count(ctx, x) result(res)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: x
        real(WP) :: res
        
        res = get_imf_value(ctx, x, mass_weighted=.false.)
    end function wrapper_imf_count

    !> @brief Wrapper to integrate Mass Density (M * dn/dM)
    pure function wrapper_imf_mass(ctx, x) result(res)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: x
        real(WP) :: res
        
        res = get_imf_value(ctx, x, mass_weighted=.true.)
    end function wrapper_imf_mass

    !> @brief Analytic integral of IMF over [m1, m2] for power-law IMF families.
    !>
    !> @details
    !> Supports IMF base types:
    !> - 0: Salpeter
    !> - 2: Kroupa (piecewise power law)
    !> - 4: Dave (two-slope power law)
    !> - 5: User-defined piecewise power law
    !>
    !> Falls back to zero for unsupported types (caller guards usage).
    pure function integrate_imf_interval_analytic(ctx, m1, m2, mass_weighted) result(int_val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m1, m2
        logical, intent(in)  :: mass_weighted
        real(WP) :: int_val

        integer :: imf_type_base
        real(WP) :: lo, hi, a1, a2, a3, mdave, c2
        real(WP) :: seg_lo, seg_hi, imfcu
        integer  :: n, n_user

        int_val = 0.0_wp
        lo = min(m1, m2)
        hi = max(m1, m2)
        if (hi <= lo) return

        imf_type_base = mod(ctx%imf_type_val, 10)

        select case (imf_type_base)
        case (0)
            int_val = integrate_powerlaw(lo, hi, ctx%state%salp_ind, mass_weighted, 1.0_wp)

        case (2)
            a1 = ctx%state%imf_alpha(1)
            a2 = ctx%state%imf_alpha(2)
            a3 = ctx%state%imf_alpha(3)
            c2 = 0.5_wp**(-a1 + a2)

            seg_lo = max(lo, 0.08_wp)
            seg_hi = min(hi, 0.5_wp)
            if (seg_hi > seg_lo) int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, a1, mass_weighted, 1.0_wp)

            seg_lo = max(lo, 0.5_wp)
            seg_hi = min(hi, 1.0_wp)
            if (seg_hi > seg_lo) int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, a2, mass_weighted, c2)

            seg_lo = max(lo, 1.0_wp)
            seg_hi = hi
            if (seg_hi > seg_lo) int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, a3, mass_weighted, c2)

        case (4)
            a1 = ctx%state%imf_alpha(1)
            a2 = ctx%state%imf_alpha(2)
            mdave = ctx%state%imf_mdave

            seg_lo = max(lo, 0.08_wp)
            seg_hi = min(hi, mdave)
            if (seg_hi > seg_lo) int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, a1, mass_weighted, 1.0_wp)

            seg_lo = max(lo, mdave)
            seg_hi = hi
            if (seg_hi > seg_lo) then
                int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, a2, mass_weighted, mdave**(-a1 + a2))
            end if

        case (5)
            n_user = ctx%state%n_user_imf
            imfcu = 1.0_wp

            ! First segment
            if (n_user >= 1) then
                seg_lo = max(lo, ctx%state%imf_user_alpha(1,1))
                seg_hi = min(hi, ctx%state%imf_user_alpha(2,1))
                if (seg_hi > seg_lo) then
                    int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, ctx%state%imf_user_alpha(3,1), mass_weighted, 1.0_wp)
                end if
            end if

            do n = 2, n_user
                imfcu = imfcu * ctx%state%imf_user_alpha(1,n)**(-ctx%state%imf_user_alpha(3,n-1) + &
                        ctx%state%imf_user_alpha(3,n))

                seg_lo = max(lo, ctx%state%imf_user_alpha(1,n))
                seg_hi = min(hi, ctx%state%imf_user_alpha(2,n))
                if (seg_hi > seg_lo) then
                    int_val = int_val + integrate_powerlaw(seg_lo, seg_hi, ctx%state%imf_user_alpha(3,n), mass_weighted, imfcu)
                end if
            end do

        case default
            int_val = 0.0_wp
        end select
    end function integrate_imf_interval_analytic

    !> @brief Integral of coeff * m^(-alpha) or coeff * m^(1-alpha) over [lo, hi].
    pure function integrate_powerlaw(lo, hi, alpha, mass_weighted, coeff) result(val)
        !$omp declare target
        real(WP), intent(in) :: lo, hi, alpha, coeff
        logical, intent(in)  :: mass_weighted
        real(WP) :: val
        real(WP) :: denom

        if (hi <= lo) then
            val = 0.0_wp
            return
        end if

        if (mass_weighted) then
            denom = 2.0_wp - alpha
        else
            denom = 1.0_wp - alpha
        end if

        if (abs(denom) > 1.0e-12_wp) then
            if (mass_weighted) then
                val = coeff * (hi**(2.0_wp - alpha) - lo**(2.0_wp - alpha)) / denom
            else
                val = coeff * (hi**(1.0_wp - alpha) - lo**(1.0_wp - alpha)) / denom
            end if
        else
            ! log integral limit when exponent is -1
            val = coeff * log(hi / lo)
        end if
    end function integrate_powerlaw

    ! ------------------------------------------------------------------------
    ! CORE IMF LOGIC
    ! ------------------------------------------------------------------------

    !> @brief
    !> Evaluates the Initial Mass Function (IMF) at specific mass points.
    !>
    !> @details
    !> Calculates the unnormalized IMF value for the algorithm specified in `ctx%imf_type_val`.
    !> Supported IMF types (determined by `mod(ctx%imf_type_val, 10)`):
    !> - 0: Salpeter ([1955](https://ui.adsabs.harvard.edu/abs/1955ApJ...121..161S))
    !> - 1: Chabrier ([2003](https://ui.adsabs.harvard.edu/abs/2003PASP..115..763C))
    !> - 2: Kroupa ([2001](https://ui.adsabs.harvard.edu/abs/2001MNRAS.322..231K))
    !> - 3: van Dokkum ([2008](https://ui.adsabs.harvard.edu/abs/2008ApJ...674...29V))
    !> - 4: Dave ([2008](https://ui.adsabs.harvard.edu/abs/2008MNRAS.385..147D))
    !> - 5: User-defined (piecewise power law)
    !>
    !> @param[in] ctx            The simulation context containing IMF parameters.
    !> @param[in] mass           Array of stellar masses (in M_sun).
    !> @param[in] mass_weighted  Control flag for the output form:
    !>                           - .false.: Returns Number Density (dN/dM).
    !>                           - .true.:  Returns Mass Density (M * dN/dM approx dN/dlnM).
    !> @return    imf_val        The calculated IMF values.
    pure function get_imf_value(ctx, mass, mass_weighted) result(imf_val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: mass
        logical, intent(in) :: mass_weighted
        real(WP) :: imf_val
        
        integer :: imf_type_base
        real(WP) :: alpha1, alpha2, alpha3
        real(WP) :: log_m, vdmc, vd_break, vd_term, mdave
        real(WP) :: imfcu
        integer  :: n, n_user

        ! The original code used (type + 10) to signal mass-weighting.
        ! We strip that here to get the base algorithm type.
        imf_type_base = mod(ctx%imf_type_val, 10)

        select case (imf_type_base)
        case (0) 
            ! Salpeter (1955)
            imf_val = mass**(-ctx%state%salp_ind)
        case (1) 
            ! Chabrier (2003)
            if (mass < 1.0_wp) then
                log_m = log10(mass)
                imf_val = exp(-((log_m - CHAB_LOG_MC)**2) * CHAB_INV_2SIG2) / mass
            else
                imf_val = CHAB_HIGH_NORM * mass**(-(CHAB_IND + 1.0_wp))
            end if
        case (2) 
            ! Kroupa (2001)
            alpha1 = ctx%state%imf_alpha(1)
            alpha2 = ctx%state%imf_alpha(2)
            alpha3 = ctx%state%imf_alpha(3)

            if (mass >= 0.08_wp .and. mass < 0.5_wp) then
                imf_val = mass**(-alpha1)
            else if (mass >= 0.5_wp .and. mass < 1.0_wp) then
                imf_val = 0.5_wp**(-alpha1 + alpha2) * mass**(-alpha2)
            else if (mass >= 1.0_wp) then
                imf_val = 0.5_wp**(-alpha1 + alpha2) * mass**(-alpha3)
            else
                imf_val = 0.0_wp
            end if
        case (3) 
            ! van Dokkum (2008)
            vdmc = ctx%state%imf_vdmc
            vd_break = VD_NC * vdmc

            if (mass <= vd_break) then
                log_m = log10(mass)
                vd_term = ((log_m - log10(vdmc))**2) / (2.0_wp * VD_SIGMA2)
                imf_val = VD_AL * (0.5_wp * vd_break)**(-VD_IND) * exp(-vd_term)
            else
                imf_val = VD_AH * mass**(-VD_IND)
            end if

            ! Convert from dn/dlnM to dn/dM
            imf_val = imf_val / mass
        case (4) 
            ! Dave (2008)
            alpha1 = ctx%state%imf_alpha(1)
            alpha2 = ctx%state%imf_alpha(2)
            mdave  = ctx%state%imf_mdave

            if (mass >= 0.08_wp .and. mass < mdave) then
                imf_val = mass**(-alpha1)
            else if (mass >= mdave) then
                imf_val = mdave**(-alpha1 + alpha2) * mass**(-alpha2)
            else
                imf_val = 0.0_wp
            end if
        case (5) 
            ! User-defined
            n_user = ctx%state%n_user_imf
            imf_val = 0.0_wp

            if (mass >= ctx%state%imf_user_alpha(1,1) .and. &
                mass <  ctx%state%imf_user_alpha(2,1)) then
                imf_val = mass**(-ctx%state%imf_user_alpha(3,1))
            end if

            imfcu = 1.0_wp
            do n = 2, n_user
                if (mass >= ctx%state%imf_user_alpha(1,n) .and. &
                    mass <  ctx%state%imf_user_alpha(2,n)) then

                    imf_val = mass**(-ctx%state%imf_user_alpha(3,n)) * &
                              ctx%state%imf_user_alpha(1,n)**(-ctx%state%imf_user_alpha(3,n-1) + &
                              ctx%state%imf_user_alpha(3,n)) * imfcu
                end if

                imfcu = imfcu * ctx%state%imf_user_alpha(1,n)**(-ctx%state%imf_user_alpha(3,n-1) + &
                        ctx%state%imf_user_alpha(3,n))
            end do
        case default
            imf_val = 0.0_wp
        end select

        ! Apply Mass Weighting if requested
        if (mass_weighted) then
            imf_val = imf_val * mass
        end if

    end function get_imf_value

    ! ------------------------------------------------------------------------
    ! SPECIFIC IMF IMPLEMENTATIONS
    ! ------------------------------------------------------------------------

    pure subroutine imf_salpeter(ctx, m, val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        
        ! Power law: m^(-alpha)
        val = m**(-ctx%state%salp_ind)
    end subroutine imf_salpeter

    pure subroutine imf_chabrier(m, val)
        !$omp declare target
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        real(WP) :: log_m, log_mc, term

        log_mc = log10(CHAB_MC)

        if (m < 1.0_wp) then
            ! Log-normal part
            log_m = log10(m)
            term = (log_m - log_mc)**2 / (2.0_wp * CHAB_SIGMA2)
            val = exp(-term) / m
        else
            ! Power law part
            term = log_mc**2 / (2.0_wp * CHAB_SIGMA2)
            val = exp(-term) * m**(-CHAB_IND) / m
        end if
    end subroutine imf_chabrier

    pure subroutine imf_kroupa(ctx, m, val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        real(WP), dimension(3) :: alpha

        alpha = ctx%state%imf_alpha

        if (m >= 0.08_wp .and. m < 0.5_wp) then
            val = m**(-alpha(1))
        else if (m >= 0.5_wp .and. m < 1.0_wp) then
            val = 0.5_wp**(-alpha(1) + alpha(2)) * m**(-alpha(2))
        else if (m >= 1.0_wp) then
            val = 0.5_wp**(-alpha(1) + alpha(2)) * m**(-alpha(3))
        else
            val = 0.0_wp
        end if
    end subroutine imf_kroupa

    pure subroutine imf_vandokkum(ctx, m, val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        real(WP) :: breakpoint, term, log_m, log_mc

        breakpoint = VD_NC * ctx%state%imf_vdmc
        log_mc = log10(ctx%state%imf_vdmc)

        if (m <= breakpoint) then
            ! Lognormal-ish
            log_m = log10(m)
            term = (log_m - log_mc)**2 / (2.0_wp * VD_SIGMA2)
            val = VD_AL * (0.5_wp * breakpoint)**(-VD_IND) * exp(-term)
        else
            ! Power law
            val = VD_AH * m**(-VD_IND)
        end if

        ! Convert from dn/dlnM to dn/dM
        val = val / m
    end subroutine imf_vandokkum

    pure subroutine imf_dave(ctx, m, val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        real(WP), dimension(3) :: alpha
        real(WP) :: mdave

        alpha = ctx%state%imf_alpha
        mdave = ctx%state%imf_mdave

        if (m >= 0.08_wp .and. m < mdave) then
            val = m**(-alpha(1))
        else if (m >= mdave) then
            val = mdave**(-alpha(1) + alpha(2)) * m**(-alpha(2))
        else
            val = 0.0_wp
        end if
    end subroutine imf_dave

    pure subroutine imf_user_defined(ctx, m, val)
        !$omp declare target
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        integer :: n
        real(WP) :: imfcu
        integer :: n_user

        ! Alias for cleaner syntax
        n_user = ctx%state%n_user_imf

        val = 0.0_wp

        ! First segment
        if (m >= ctx%state%imf_user_alpha(1,1) .and. &
            m <  ctx%state%imf_user_alpha(2,1)) then
            val = m**(-ctx%state%imf_user_alpha(3,1))
        end if

        ! Subsequent segments
        imfcu = 1.0_wp
        do n = 2, n_user
            if (m >= ctx%state%imf_user_alpha(1,n) .and. &
                m <  ctx%state%imf_user_alpha(2,n)) then
                
                val = m**(-ctx%state%imf_user_alpha(3,n)) * &
                         ctx%state%imf_user_alpha(1,n)**(-ctx%state%imf_user_alpha(3,n-1) + &
                         ctx%state%imf_user_alpha(3,n)) * imfcu
            end if
            
            ! Update cumulative factor
            imfcu = imfcu * ctx%state%imf_user_alpha(1,n)**(-ctx%state%imf_user_alpha(3,n-1) + &
                    ctx%state%imf_user_alpha(3,n))
        end do
    end subroutine imf_user_defined

end module fsps_imf