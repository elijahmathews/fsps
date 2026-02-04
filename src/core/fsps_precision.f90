module fsps_precision
    !> @brief
    !> Defines the precision parameter for FSPS.
    !>
    !> @details
    !> This module serves as the root dependency for precision definitions (WP)
    !> and is used by almost all other modules.

    ! Use the modern ISO standard for precision definitions
    use, intrinsic :: iso_fortran_env, only: real32, real64

    implicit none
    private

    public :: WP

! Default to Double (64-bit) unless the flag is set
#ifdef FORCE_SINGLE_PRECISION
    integer, parameter :: WP = real32
#else
    integer, parameter :: WP = real64
#endif

end module fsps_precision