module fsps_interpolation
    !> @brief
    !> Provides robust, type-safe utilities for 1D linear interpolation.
    !>
    !> @details
    !> This module contains routines for finding intervals in monotonic arrays
    !> (bisection search) and performing linear interpolation/extrapolation.
    !>
    !> It exposes a generic interface `interpolate_linear` that automatically
    !> handles both scalar and array query points.
    
    use fsps_precision, only: WP
    implicit none

    private
    
    ! Publicly expose the generic name and the utility
    public :: find_interval, interpolate_linear
    
    !> @brief
    !> Abstract interface for a generic linear interpolation function.
    interface interpolate_linear
        module procedure interpolate_linear_array
        module procedure interpolate_linear_scalar
    end interface

contains

    !> @brief
    !> Linearly interpolates a function y = f(x) at a single query point.
    !>
    !> @details
    !> Interpolates a value from the table (x_in, y_in) at point x_out.
    !> Uses linear interpolation. For query points outside the range of x_in,
    !> it uses the slope of the nearest boundary interval (linear extrapolation).
    !>
    !> @param[in] x_in   The x-coordinates of the data table (must be monotonic).
    !> @param[in] y_in   The y-coordinates of the data table.
    !> @param[in] x_out  The coordinate at which to interpolate.
    !>
    !> @return    y_out  The interpolated value.
    pure function interpolate_linear_scalar(x_in, y_in, x_out) result(y_out)
        real(WP), dimension(:), intent(in), contiguous :: x_in, y_in
        real(WP), intent(in) :: x_out
        real(WP) :: y_out

        integer :: idx
        real(WP) :: slope
        
        ! Return NaN if inputs are insufficient or mismatched
        if (size(x_in) < 2 .or. (size(x_in) /= size(y_in))) then
            y_out = get_quiet_nan()
            return
        end if

        ! Find the interval index and clamp to valid interpolation range.
        idx = find_interval(x_in, x_out)
        idx = max(1, min(idx, size(x_in) - 1))

        ! Calculate slope: (y2 - y1) / (x2 - x1)
        slope = (y_in(idx+1) - y_in(idx)) / (x_in(idx+1) - x_in(idx))

        ! Interpolate: y = y1 + slope * (x - x1)
        y_out = y_in(idx) + slope * (x_out - x_in(idx))

    end function interpolate_linear_scalar

    !> @brief
    !> Linearly interpolates a function y = f(x) at multiple query points.
    !>
    !> @details
    !> Interpolates values from the table (x_in, y_in) to the points defined in x_out.
    !> Uses linear interpolation. For query points outside the range of x_in,
    !> it uses the slope of the nearest boundary interval (linear extrapolation).
    !>
    !> @param[in] x_in   The x-coordinates of the data table (must be monotonic).
    !> @param[in] y_in   The y-coordinates of the data table.
    !> @param[in] x_out  The coordinates at which to interpolate.
    !>
    !> @return    y_out  The interpolated values at x_out.
    pure function interpolate_linear_array(x_in, y_in, x_out) result(y_out)
        real(WP), dimension(:), intent(in), contiguous :: x_in, y_in, x_out
        real(WP), dimension(size(x_out)) :: y_out

        integer :: i, idx, n
        real(WP) :: slope
        
        n = size(x_in)
        
        ! Check sizes once before looping
        if (n < 2 .or. (n /= size(y_in))) then
            ! Return array of NaNs on error
            y_out(:) = get_quiet_nan()
            return
        end if

        ! Loop over all query points.
        ! The contiguous attribute allows the compiler to vectorize this loop.
        do i = 1, size(x_out)
            ! Find the interval index and clamp to valid interpolation range.
            idx = find_interval(x_in, x_out(i))
            idx = max(1, min(idx, n - 1))

            ! Calculate slope for this interval: (y2 - y1) / (x2 - x1)
            slope = (y_in(idx+1) - y_in(idx)) / (x_in(idx+1) - x_in(idx))

            ! Interpolate: y = y1 + slope * (x - x1)
            y_out(i) = y_in(idx) + slope * (x_out(i) - x_in(idx))
        end do

    end function interpolate_linear_array

    !> @brief
    !> Finds the index `i` in an array such that `array(i) <= value < array(i+1)`.
    !>
    !> @details
    !> Uses a bisection (binary search) algorithm to find the interval containing the
    !> input value. It is designed to work efficiently with both monotonic increasing 
    !> and decreasing arrays. 
    !>
    !> This implementation is typically used for interpolation. For values exactly 
    !> matching the boundaries, it clamps the index to ensure valid interpolation intervals.
    !>
    !> @param[in] array  The 1D sorted array to search (monotonic).
    !> @param[in] value  The value to find within the array.
    !>
    !> @return    idx    The lower index of the interval. 
    !>                   Returns 1 if value < array(1).
    !>                   Returns n-1 if value > array(n).
    pure function find_interval(array, value) result(idx)
        real(WP), dimension(:), intent(in), contiguous :: array
        real(WP), intent(in) :: value
        integer :: idx

        integer :: n, lower, upper, mid
        logical :: is_ascending

        n = size(array)

        ! Safety check for empty or singleton arrays
        if (n < 2) then
            idx = 1
            return
        end if

        ! Determine sort order
        is_ascending = (array(n) >= array(1))

        ! Initialize bisection limits
        ! Note: 0 and n+1 are used to handle out-of-bound values gracefully
        lower = 0
        upper = n + 1

        ! Perform Binary Search
        do while (upper - lower > 1)
            mid = (upper + lower) / 2
            
            ! The .eqv. operator cleverly handles both ascending and descending logic:
            ! Ascending + Value >= Mid -> True  (Keep upper half -> lower = mid)
            ! Descending + Value >= Mid -> False (Keep lower half -> upper = mid)
            if (is_ascending .eqv. (value >= array(mid))) then
                lower = mid
            else
                upper = mid
            end if
        end do

        ! Match legacy locate behavior for out-of-bounds values.
        if (value == array(1)) then
            idx = 1
        else if (value == array(n)) then
            idx = n - 1
        else
            idx = lower
        end if

    end function find_interval

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------
    
    !> @brief
    !> Helper to generate a Quiet NaN (Not a Number)
    !>
    !> @details
    !> Wraps the IEEE_ARITHMETIC intrinsic to avoid module namespace pollution
    !> and ensures a safe return value for error conditions in PURE functions.
    pure function get_quiet_nan() result(res)
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        real(WP) :: res
        res = ieee_value(0.0_wp, ieee_quiet_nan)
    end function get_quiet_nan

end module fsps_interpolation