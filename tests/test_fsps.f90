program test_fsps
    !> @brief
    !> Master unit test runner for all FSPS modules.
    !> Aggregates results from individual test modules.

    ! IMF modules
    use test_fsps_imf_mod, only: run_fsps_imf_tests, &
                                 failures_imf => total_failures, &
                                 tests_imf => total_tests
    ! Math modules
    use test_fsps_integration_mod, only: run_fsps_integration_tests, &
                                         failures_integration => total_failures, &
                                         tests_integration => total_tests
    use test_fsps_interpolation_mod, only: run_fsps_interpolation_tests, &
                                           failures_interpolation => total_failures, &
                                           tests_interpolation => total_tests
    use test_fsps_special_functions_mod, only: run_fsps_special_functions_tests, &
                                           failures_special_functions => total_failures, &
                                           tests_special_functions => total_tests
    ! Test utilities
    use test_utils_mod, only: print_summary_line, print_major_header
    implicit none

    integer :: grand_total_failures = 0
    integer :: grand_total_tests = 0

    call print_major_header("FSPS UNIT TEST SUITE")

    ! Run all the tests
    call run_fsps_imf_tests()
    call run_fsps_integration_tests()
    call run_fsps_interpolation_tests()
    call run_fsps_special_functions_tests()

    ! --- Summary ---
    call print_major_header("FSPS UNIT TEST FINAL REPORT")
    call print_summary_line("fsps_imf", (tests_imf - failures_imf), tests_imf)
    call print_summary_line("fsps_integration", (tests_integration - failures_integration), tests_integration)
    call print_summary_line("fsps_interpolation", (tests_interpolation - failures_interpolation), tests_interpolation)
    call print_summary_line("fsps_special_functions", (tests_special_functions - failures_special_functions), tests_special_functions)
    
    grand_total_failures = failures_imf + &
                           failures_integration + &
                           failures_interpolation + &
                           failures_special_functions
    
    grand_total_tests = tests_imf + &
                        tests_integration + &
                        tests_interpolation + &
                        tests_special_functions

    print *
    call print_summary_line("Result", grand_total_tests, grand_total_tests)

    if (grand_total_failures == 0) then
        stop 0
    else
        stop 1
    end if

end program test_fsps