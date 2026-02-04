module test_utils_mod
    use fsps_constants, only: SP
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_class, ieee_negative_inf, ieee_class_type, operator(==)
    implicit none

    ! ANSI Color Codes
    character(len=4), parameter :: C_RESET    = char(27)//'[0m'
    character(len=4), parameter :: C_BOLD     = char(27)//'[1m'
    character(len=5), parameter :: C_RED      = char(27)//'[31m'
    character(len=5), parameter :: C_GREEN    = char(27)//'[32m'
    character(len=7), parameter :: C_RED_BOLD = char(27)//'[1;31m' ! Bold Red

    public :: assert_float_equals, assert_int_equals, assert_true
    public :: assert_is_nan, assert_is_neg_inf
    public :: print_major_header, print_minor_header, print_group, print_summary_line
    public :: C_RESET, C_BOLD, C_RED, C_GREEN, C_RED_BOLD

contains

    !> @brief Print a major section header in Bold
    subroutine print_major_header(text)
        character(len=*), intent(in) :: text

        print *
        print *, C_BOLD // "=========================================================" // C_RESET
        print *, C_BOLD // trim(text) // C_RESET
        print *, C_BOLD // "=========================================================" // C_RESET
    end subroutine

    !> @brief Print a minor section header in Bold
    subroutine print_minor_header(text)
        character(len=*), intent(in) :: text

        print *
        print *, C_BOLD // "---------------------------------------------------------" // C_RESET
        print *, C_BOLD // trim(text) // C_RESET
        print *, C_BOLD // "---------------------------------------------------------" // C_RESET
    end subroutine

    !> @brief Print a test group sub-header in Bold
    subroutine print_group(text)
        character(len=*), intent(in) :: text
        print *, C_BOLD // ">> Group: " // trim(text) // C_RESET
    end subroutine

    !> @brief Print a formatted summary line (Label: Pass / Total)
    !> Handles the coloring logic: Green numbers if all passed, Bold Red if failures exist.
    subroutine print_summary_line(label, n_pass, n_total)
        character(len=*), intent(in) :: label
        integer, intent(in) :: n_pass, n_total
        character(len=20) :: s_pass, s_total ! Buffers for integer conversion
        
        ! Write integers to strings to allow concatenation
        write(s_pass, '(I0)') n_pass
        write(s_total, '(I0)') n_total
        
        if (n_pass == n_total) then
            ! All passed: Green numbers
            print *, trim(label), ": ", C_GREEN // trim(s_pass) // C_RESET, "/", C_GREEN // trim(s_total) // C_RESET
        else
            ! Failures: Bold Red numbers
            print *, trim(label), ": ", C_RED_BOLD // trim(s_pass) // C_RESET, "/", C_RED_BOLD // trim(s_total) // C_RESET
        end if
    end subroutine

    subroutine assert_float_equals(expected, actual, tol, label, n_tests, n_fails)
        real(SP), intent(in) :: expected, actual, tol
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label
        
        n_tests = n_tests + 1
        if (abs(expected - actual) <= tol) then
            ! Pass: "PASS" is Green, label is normal
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            ! Fail: "FAIL" is Red Bold, label is normal
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
            print *, "         Expected: ", expected, " Got: ", actual, " Diff: ", abs(expected-actual)
        end if
    end subroutine assert_float_equals

    subroutine assert_int_equals(expected, actual, label, n_tests, n_fails)
        integer, intent(in) :: expected, actual
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label
        
        n_tests = n_tests + 1
        if (expected == actual) then
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
            print *, "         Expected: ", expected, " Got: ", actual
        end if
    end subroutine assert_int_equals

    subroutine assert_true(condition, label, n_tests, n_fails)
        logical, intent(in) :: condition
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label

        n_tests = n_tests + 1
        if (condition) then
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
        end if
    end subroutine assert_true

    subroutine assert_is_nan(val, label, n_tests, n_fails)
        real(SP), intent(in) :: val
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label

        n_tests = n_tests + 1
        if (ieee_is_nan(val)) then
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
            print *, "         Expected NaN, Got: ", val
        end if
    end subroutine assert_is_nan

    subroutine assert_is_neg_inf(val, label, n_tests, n_fails)
        real(SP), intent(in) :: val
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label
        type(ieee_class_type) :: cl

        cl = ieee_class(val)
        n_tests = n_tests + 1
        if (cl == ieee_negative_inf) then
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
            print *, "         Expected -Inf"
        end if
    end subroutine assert_is_neg_inf

    subroutine assert_relative_error(expected, actual, rel_tol, label, n_tests, n_fails)
        real(SP), intent(in) :: expected, actual, rel_tol
        integer, intent(inout) :: n_tests, n_fails
        character(len=*), intent(in) :: label
        real(SP) :: rel_diff
        
        rel_diff = abs((expected - actual) / expected)
        
        n_tests = n_tests + 1
        if (rel_diff <= rel_tol) then
            print *, "  [" // C_GREEN // "PASS" // C_RESET // "] " // label
        else
            n_fails = n_fails + 1
            print *, "  [" // C_RED_BOLD // "FAIL" // C_RESET // "] " // label
            print *, "         Expected:", expected, " Got:", actual, " RelDiff:", rel_diff
        end if
    end subroutine assert_relative_error

end module test_utils_mod