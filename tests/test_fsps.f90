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
    use test_fsps_dust_mod, only: run_fsps_dust_tests, &
                                  failures_dust => total_failures, &
                                  tests_dust => total_tests
    use test_fsps_gas_mod, only: run_fsps_gas_tests, &
                                 failures_gas => total_failures, &
                                 tests_gas => total_tests
    use test_fsps_stellar_modifications_mod, only: run_fsps_stellar_modifications_tests, &
                                 failures_stellar_modifications => total_failures, &
                                 tests_stellar_modifications => total_tests
    use test_fsps_cosmology_mod, only: run_fsps_cosmology_tests, &
                                  failures_cosmology => total_failures, &
                                  tests_cosmology => total_tests
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
    call run_fsps_dust_tests()
    call run_fsps_gas_tests()
    call run_fsps_stellar_modifications_tests()
    call run_fsps_cosmology_tests()

    ! --- Summary ---
    call print_major_header("FSPS UNIT TEST FINAL REPORT")
    call print_summary_line( &
        "fsps_imf", &
        (tests_imf - failures_imf), &
        tests_imf &
    )
    call print_summary_line( &
        "fsps_integration", &
        (tests_integration - failures_integration), &
        tests_integration &
    )
    call print_summary_line( &
        "fsps_interpolation", &
        (tests_interpolation - failures_interpolation), &
        tests_interpolation &
    )
    call print_summary_line( &
        "fsps_special_functions", &
        (tests_special_functions - failures_special_functions), &
        tests_special_functions &
    )
    call print_summary_line( &
        "fsps_dust", &
        (tests_dust - failures_dust), &
        tests_dust &
    )
    call print_summary_line( &
        "fsps_gas", &
        (tests_gas - failures_gas), &
        tests_gas &
    )
    call print_summary_line( &
        "fsps_stellar_modifications", &
        (tests_stellar_modifications - failures_stellar_modifications), &
        tests_stellar_modifications &
    )
    call print_summary_line( &
        "fsps_cosmology", &
        (tests_cosmology - failures_cosmology), &
        tests_cosmology &
    )
    
    grand_total_failures = failures_imf + &
                           failures_integration + &
                           failures_interpolation + &
                           failures_special_functions + &
                           failures_dust + &
                           failures_gas + &
                           failures_stellar_modifications + &
                           failures_cosmology
    
    grand_total_tests = tests_imf + &
                        tests_integration + &
                        tests_interpolation + &
                        tests_special_functions + &
                        tests_dust + &
                        tests_gas + &
                        tests_stellar_modifications + &
                        tests_cosmology

    print *
    call print_summary_line("Result", grand_total_tests - grand_total_failures, grand_total_tests)
    print *

    if (grand_total_failures == 0) then
        stop 0
    else
        stop 1
    end if

end program test_fsps