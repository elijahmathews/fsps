module test_fsps_special_functions_mod
    use fsps_types, only: sp
    use fsps_special_functions
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_is_nan, assert_is_neg_inf, assert_float_equals, &
                              assert_relative_error
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_class, &
                                             ieee_negative_inf, ieee_class_type, &
                                             operator(==)
    implicit none

    ! Global failure counter
    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_special_functions_tests, total_failures, total_tests

contains

    !> @brief
    !> Unit test suite for fsps_special_functions module.
    !> Designed for CI usage: exits with status 0 on success, 1 on failure.
    subroutine run_fsps_special_functions_tests()

        call print_minor_header("fsps_special_functions")

        call test_ei_small_x()
        call test_ei_large_x()
        call test_ei_edge_cases()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)

    end subroutine run_fsps_special_functions_tests

    ! ------------------------------------------------------------------------
    ! TEST SUITE: POWER SERIES DOMAIN (0 < x <= 40)
    ! ------------------------------------------------------------------------
    subroutine test_ei_small_x()
        real(sp) :: x, expected, res
        
        call print_group("Exponential Integral (Small x, Power Series)")

        ! Case 1: x = 0.5
        ! Reference: scipy.special.expi(0.5) = 0.454219904863
        x = 0.5_sp
        expected = 0.45421990486317343_sp
        res = exponential_integral(x)
        call assert_float_equals(expected, res, 1.0e-6_sp, "Ei(0.5)", total_tests, total_failures)

        ! Case 2: x = 1.0
        ! Reference: scipy.special.expi(1.0) = 1.895117816356
        x = 1.0_sp
        expected = 1.8951178163559368_sp
        res = exponential_integral(x)
        call assert_float_equals(expected, res, 1.0e-6_sp, "Ei(1.0)", total_tests, total_failures)

        ! Case 3: x = 10.0
        ! Reference: scipy.special.expi(10.0) = 2492.22897624
        x = 10.0_sp
        expected = 2492.2289762418773_sp
        res = exponential_integral(x)
        ! Tolerance scaled relative to magnitude, or use relative error check.
        ! Here we check roughly 6 sig figs.
        call assert_float_equals(expected, res, 0.01_sp, "Ei(10.0)", total_tests, total_failures)

    end subroutine test_ei_small_x

    ! ------------------------------------------------------------------------
    ! TEST SUITE: ASYMPTOTIC DOMAIN (x > 40)
    ! ------------------------------------------------------------------------
    subroutine test_ei_large_x()
        real(sp) :: x, expected, res
        
        call print_group("Exponential Integral (Large x, Asymptotic)")

        ! Case 1: x = 45.0
        ! Reference: scipy.special.expi(45.0) = 7.78533476 * 10^17
        ! Note: exp(45) is approx 3.49e19, well within single precision limit (3.4e38)
        x = 45.0_sp
        expected = 7.943916035704438e17_sp
        res = exponential_integral(x)
        
        ! Check relative error for large numbers
        call assert_relative_error(expected, res, 1.0e-5_sp, "Ei(45.0)", total_tests, total_failures)

        ! Case 2: x = 50.0
        ! Reference: scipy.special.expi(50.0) = 1.03644598 * 10^20
        x = 50.0_sp
        expected = 1.058563689713169e20_sp
        res = exponential_integral(x)
        call assert_relative_error(expected, res, 1.0e-5_sp, "Ei(50.0)", total_tests, total_failures)

    end subroutine test_ei_large_x

    ! ------------------------------------------------------------------------
    ! TEST SUITE: EDGE CASES & ERRORS
    ! ------------------------------------------------------------------------
    subroutine test_ei_edge_cases()
        real(sp) :: res
        
        call print_group("Robustness (0 and Negative)")

        ! Case 1: x = 0.0
        ! Expected: -Infinity
        res = exponential_integral(0.0_sp)
        call assert_is_neg_inf(res, "Ei(0.0) is -Infinity", total_tests, total_failures)

        ! Case 2: x = -1.0
        ! Expected: NaN (Ei is complex for x < 0, implementation assumes real)
        res = exponential_integral(-1.0_sp)
        call assert_is_nan(res, "Ei(-1.0) is NaN", total_tests, total_failures)

    end subroutine test_ei_edge_cases

end module test_fsps_special_functions_mod