module fsps_special_functions
    !> @brief
    !> Provides mathematical special functions needed for FSPS physics.
    !>
    !> @details
    !> Currently contains the Exponential Integral Ei(x).
    
    use fsps_precision, only: WP
    implicit none

    private

    public :: exponential_integral

contains

    !> @brief
    !> Computes the Exponential Integral Ei(x).
    !>
    !> @details
    !> Based on the implementation by Shanjie Zhang and Jianming Jin.
    !> Uses a power series expansion for small x (|x| <= 40) and an
    !> asymptotic expansion for large x.
    !>
    !> @param[in] x  The argument (must be > 0 for real result).
    !> @return    res The value of Ei(x). 
    !>                Returns -Infinity at x=0.
    !>                Returns NaN for x < 0.
    pure function exponential_integral(x) result(res)
        real(WP), intent(in) :: x
        real(WP) :: res

        integer :: k
        integer, parameter :: MAXIT = 1000
        real(WP) :: r
        
        ! Euler-Mascheroni constant
        real(WP), parameter :: GAMMA = 0.5772156649015328_wp
        real(WP), parameter :: EPS = 1.0e-20_wp
        
        if (x == 0.0_wp) then
            ! Return negative infinity
            res = get_neg_infinity()
            return
        else if (x < 0.0_wp) then
            ! Ei(x) is complex for x < 0; return NaN in this real implementation
            res = get_quiet_nan()
            return
        end if

        if (abs(x) <= 40.0_wp) then
            ! Power series expansion around x=0
            ! Ei(x) = Gamma + ln(x) + sum(x^k / (k * k!))
            res = 1.0_wp
            r = 1.0_wp
            
            do k = 1, MAXIT
                r = r * k * x / (real(k, WP) + 1.0_wp)**2
                res = res + r
                
                ! Convergence check
                if (abs(r/res) <= EPS) exit
                
                ! Check for non-convergence
                if (k == MAXIT) then
                    res = get_quiet_nan()
                    return
                end if
            end do
            
            res = GAMMA + log(x) + x * res

        else
            ! Asymptotic expansion for large x (divergent series, limited terms)
            ! Ei(x) ~ exp(x)/x * (1 + 1!/x + 2!/x^2 + ...)
            res = 1.0_wp
            r = 1.0_wp
            
            do k = 1, 20
                r = r * real(k, WP) / x
                res = res + r
            end do
            
            res = exp(x) / x * res
        end if

    end function exponential_integral

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> @brief Helper to generate a Quiet NaN
    pure function get_quiet_nan() result(val)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real(WP) :: val
        val = ieee_value(0.0_wp, ieee_quiet_nan)
    end function get_quiet_nan

    !> @brief Helper to generate Negative Infinity
    pure function get_neg_infinity() result(val)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf
        real(WP) :: val
        val = ieee_value(0.0_wp, ieee_negative_inf)
    end function get_neg_infinity

end module fsps_special_functions