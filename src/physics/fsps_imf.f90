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
    public :: CHAB_MC, CHAB_SIGMA2, CHAB_IND

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------

    ! Chabrier 2003
    real(WP), parameter :: CHAB_MC     = 0.08_wp
    real(WP), parameter :: CHAB_SIGMA2 = 0.69_wp**2  ! Explicit squaring
    real(WP), parameter :: CHAB_IND    = 1.3_wp

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
        integer :: i

        ! Access limits from state
        real(WP) :: lower_limit, upper_limit, lower_bound

        lower_limit = ctx%state%imf_lower_limit
        upper_limit = ctx%state%imf_upper_limit
        lower_bound = ctx%state%imf_lower_bound

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
            weights(i) = integrate_romberg(ctx, wrapper_imf_count, m1, m2)
        end do

        ! 2. Calculate Total Mass (Normalization Factor)
        ! Integrate M * dn/dM over the full range [lower_limit, upper_limit]
        total_mass = integrate_romberg(ctx, wrapper_imf_mass, lower_limit, upper_limit)

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
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: x
        real(WP) :: res
        
        res = get_imf_value(ctx, x, mass_weighted=.false.)
    end function wrapper_imf_count

    !> @brief Wrapper to integrate Mass Density (M * dn/dM)
    pure function wrapper_imf_mass(ctx, x) result(res)
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: x
        real(WP) :: res
        
        res = get_imf_value(ctx, x, mass_weighted=.true.)
    end function wrapper_imf_mass

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
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: mass
        logical, intent(in) :: mass_weighted
        real(WP) :: imf_val
        
        integer :: imf_type_base

        ! The original code used (type + 10) to signal mass-weighting.
        ! We strip that here to get the base algorithm type.
        imf_type_base = mod(ctx%imf_type_val, 10)

        select case (imf_type_base)
        case (0) 
            ! Salpeter (1955)
            call imf_salpeter(ctx, mass, imf_val)
        case (1) 
            ! Chabrier (2003)
            call imf_chabrier(mass, imf_val)
        case (2) 
            ! Kroupa (2001)
            call imf_kroupa(ctx, mass, imf_val)
        case (3) 
            ! van Dokkum (2008)
            call imf_vandokkum(ctx, mass, imf_val)
        case (4) 
            ! Dave (2008)
            call imf_dave(ctx, mass, imf_val)
        case (5) 
            ! User-defined
            call imf_user_defined(ctx, mass, imf_val)
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
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in) :: m
        real(WP), intent(out) :: val
        
        ! Power law: m^(-alpha)
        val = m**(-ctx%state%salp_ind)
    end subroutine imf_salpeter

    pure subroutine imf_chabrier(m, val)
        !$acc routine seq
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
        !$acc routine seq
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
        !$acc routine seq
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
        !$acc routine seq
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
        !$acc routine seq
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