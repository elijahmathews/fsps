module fsps_integration
    !> @brief
    !> Provides robust numerical integration routines.
    !>
    !> @details
    !> This module consolidates discrete array integration (Trapezoidal rule)
    !> and function integration (Romberg method with Richardson extrapolation).
    !>
    !> It supports both standard mathematical functions and context-aware functions
    !> used in FSPS physics calculations.

    use fsps_types, only: sp
    use fsps_context_types, only: fsps_context_t
    implicit none

    private

    ! Publicly expose the integration routines
    public :: integrate_trapezoid_array
    public :: integrate_romberg

    !> @brief
    !> Abstract interface for a standard vectorized function f(x).
    abstract interface
        pure function func_interface(x) result(res)
            import :: sp
            real(sp), dimension(:), intent(in) :: x
            real(sp), dimension(size(x)) :: res
        end function func_interface
    end interface

    !> @brief
    !> Abstract interface for a context-aware vectorized function f(ctx, x).
    abstract interface
        pure function func_ctx_interface(ctx, x) result(res)
            import :: sp, fsps_context_t
            type(fsps_context_t), intent(in) :: ctx
            real(sp), dimension(:), intent(in) :: x
            real(sp), dimension(size(x)) :: res
        end function func_ctx_interface
    end interface

    !> @brief
    !> Generic interface for Romberg Integration.
    !> Routes to either the simple or context-aware implementation.
    interface integrate_romberg
        module procedure integrate_romberg_simple
        module procedure integrate_romberg_context
    end interface

contains

    !> @brief
    !> Computes the integral of y(x) using the composite trapezoidal rule.
    !>
    !> @details
    !> Calculates the area under the curve defined by discrete points (x, y).
    !>
    !> @note
    !> This routine calculates the **unsigned** area. It uses the absolute difference
    !> in x (`|x_{i+1} - x_i|`). This ensures the calculated area is positive
    !> even if the x-axis is sorted in descending order.
    !>
    !> @param[in] x     The x-coordinates (abscissa).
    !> @param[in] y     The y-coordinates (ordinate).
    !>
    !> @return    area  The integrated area. Returns NaN on error.
    pure function integrate_trapezoid_array(x, y) result(area)
        real(sp), dimension(:), intent(in), contiguous :: x, y
        real(sp) :: area

        integer :: n

        n = size(x)

        ! Require at least 2 points and matching dimensions
        if (n < 2 .or. (n /= size(y))) then
            area = get_quiet_nan()
            return
        end if

        ! Vectorized trapezoidal calculation.
        ! Math: Area = 0.5 * Sum( |dx| * (y_i + y_{i+1}) )
        ! Optimizations:
        ! 1. Factor 0.5 out of the sum to replace N divisions with 1 multiplication.
        ! 2. Use array slicing for vectorization.
        area = 0.5_sp * sum( abs(x(2:n) - x(1:n-1)) * (y(2:n) + y(1:n-1)) )

    end function integrate_trapezoid_array

    !> @brief
    !> Integrates a function f(x) from a to b using Romberg (1995) integration.
    !>
    !> @details
    !> Uses the trapezoidal rule with iterative refinement and Neville's algorithm
    !> for Richardson extrapolation to step size h=0.
    !>
    !> @param[in] func  The function to integrate (must accept/return arrays).
    !> @param[in] a     Lower integration limit.
    !> @param[in] b     Upper integration limit.
    !> @return    res   The integral value. Returns NaN on non-convergence.
    pure function integrate_romberg_simple(func, a, b) result(res)
        procedure(func_interface) :: func
        real(sp), intent(in) :: a, b
        real(sp) :: res

        integer, parameter :: MAX_STEPS = 20
        integer, parameter :: K_ORDER = 5
        real(sp), parameter :: EPS = 1.0e-7_sp
        
        real(sp), dimension(MAX_STEPS + 1) :: h, s
        real(sp) :: dqromb, zero_h
        integer :: j

        zero_h = 0.0_sp
        h(1) = 1.0_sp

        ! Initial coarse step
        call refine_trapezoid_simple(func, a, b, s(1), 1)
        res = s(1)

        do j = 1, MAX_STEPS
            ! Copy previous estimate to current slot before refining
            s(j+1) = s(j)

            ! Refine the trapezoidal sum (reduce step size by half)
            call refine_trapezoid_simple(func, a, b, s(j+1), j+1)
            res = s(j+1)
            
            ! Record the relative step size for extrapolation
            h(j+1) = 0.25_sp * h(j)

            ! Perform Richardson Extrapolation after K steps
            if (j >= K_ORDER) then
                call polynomial_extrapolation(h(j-K_ORDER+1:j+1), s(j-K_ORDER+1:j+1), &
                                     zero_h, res, dqromb)
                
                ! Check for convergence
                if (abs(dqromb) <= EPS * abs(res)) return
            end if
            
            ! Carry forward the latest sum if we haven't converged yet
        end do

        ! Fallback: Non-convergence
        res = get_quiet_nan()

    end function integrate_romberg_simple

    !> @brief
    !> Integrates a context-aware function f(ctx, x) using Romberg integration.
    !> See `integrate_romberg_simple` for algorithm details.
    pure function integrate_romberg_context(ctx, func, a, b) result(res)
        type(fsps_context_t), intent(in) :: ctx
        procedure(func_ctx_interface) :: func
        real(sp), intent(in) :: a, b
        real(sp) :: res

        integer, parameter :: MAX_STEPS = 20
        integer, parameter :: K_ORDER = 5
        real(sp), parameter :: EPS = 1.0e-7_sp
        
        real(sp), dimension(MAX_STEPS + 1) :: h, s
        real(sp) :: dqromb, zero_h
        integer :: j

        zero_h = 0.0_sp
        h(1) = 1.0_sp

        call refine_trapezoid_context(ctx, func, a, b, s(1), 1)
        res = s(1)

        do j = 1, MAX_STEPS
            ! Copy previous estimate to current slot before refining
            s(j+1) = s(j)
            
            call refine_trapezoid_context(ctx, func, a, b, s(j+1), j+1)
            res = s(j+1)
            h(j+1) = 0.25_sp * h(j)

            if (j >= K_ORDER) then
                call polynomial_extrapolation(h(j-K_ORDER+1:j+1), s(j-K_ORDER+1:j+1), &
                                     zero_h, res, dqromb)
                if (abs(dqromb) <= EPS * abs(res)) return
            end if
        end do

        ! Fallback: Non-convergence
        res = get_quiet_nan()

    end function integrate_romberg_context

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> @brief
    !> Refines the trapezoidal approximation s_old to s_new by adding points 
    !> at the midpoints of the previous intervals.
    !> (Simple Function Version)
    pure subroutine refine_trapezoid_simple(func, a, b, s, n)
        procedure(func_interface) :: func
        real(sp), intent(in) :: a, b
        real(sp), intent(inout) :: s
        integer, intent(in) :: n

        real(sp) :: del, fsum
        integer :: it, i
        real(sp), allocatable :: x_points(:)

        if (n == 1) then
            s = 0.5_sp * (b - a) * sum(func([a, b]))
        else
            it = 2**(n - 2)
            del = (b - a) / real(it, sp)
            
            ! Modern array constructor replaces MYARTH
            ! Generates midpoints: a + 0.5*del, a + 1.5*del, ...
            x_points = [ (a + 0.5_sp * del + real(i - 1, sp) * del, i = 1, it) ]
            
            fsum = sum(func(x_points))
            s = 0.5_sp * (s + del * fsum)
        end if
    end subroutine refine_trapezoid_simple

    !> @brief
    !> Refines the trapezoidal approximation s_old to s_new.
    !> (Context Function Version)
    pure subroutine refine_trapezoid_context(ctx, func, a, b, s, n)
        type(fsps_context_t), intent(in) :: ctx
        procedure(func_ctx_interface) :: func
        real(sp), intent(in) :: a, b
        real(sp), intent(inout) :: s
        integer, intent(in) :: n

        real(sp) :: del, fsum
        integer :: it, i
        real(sp), allocatable :: x_points(:)

        if (n == 1) then
            s = 0.5_sp * (b - a) * sum(func(ctx, [a, b]))
        else
            it = 2**(n - 2)
            del = (b - a) / real(it, sp)
            x_points = [ (a + 0.5_sp * del + real(i - 1, sp) * del, i = 1, it) ]
            
            fsum = sum(func(ctx, x_points))
            s = 0.5_sp * (s + del * fsum)
        end if
    end subroutine refine_trapezoid_context

    !> @brief
    !> Polynomial interpolation/extrapolation (Neville's Algorithm).
    !> Given arrays xa and ya, returns value y at point x, and error estimate dy.
    pure subroutine polynomial_extrapolation(xa, ya, x, y, dy)
        real(sp), dimension(:), intent(in) :: xa, ya
        real(sp), intent(in) :: x
        real(sp), intent(out) :: y, dy

        integer :: m, n, ns, i
        real(sp), dimension(size(xa)) :: c, d, den, dist
        real(sp) :: w

        n = size(xa)
        c = ya
        d = ya
        dist = xa - x
        
        ! Find nearest neighbor index
        ns = minloc(abs(dist), 1)
        y = ya(ns)
        ns = ns - 1

        do m = 1, n - 1
            do i = 1, n - m
                den(i) = dist(i) - dist(i + m)
                w = c(i + 1) - d(i)
                
                ! Protect against division by zero (identical support points)
                if (den(i) == 0.0_sp) then
                    y = get_quiet_nan()
                    dy = get_quiet_nan()
                    return
                end if
                
                den(i) = w / den(i)
                d(i) = dist(i + m) * den(i)
                c(i) = dist(i) * den(i)
            end do
            
            if (2 * ns < n - m) then
                dy = c(ns + 1)
            else
                dy = d(ns)
                ns = ns - 1
            end if
            y = y + dy
        end do
    end subroutine polynomial_extrapolation

    !> @brief Helper to generate a Quiet NaN
    pure function get_quiet_nan() result(res)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real(sp) :: res
        res = ieee_value(0.0_sp, ieee_quiet_nan)
    end function get_quiet_nan

end module fsps_integration