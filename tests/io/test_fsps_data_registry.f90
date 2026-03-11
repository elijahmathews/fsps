module test_fsps_data_registry_mod
    use fsps_data_backend, only: data_backend_t, backend_status_t, backend_status_ok
    use fsps_data_backend_hdf5, only: hdf5_backend_t
    use fsps_data_registry, only: create_data_backend
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, assert_true, assert_int_equals

    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_data_registry_tests, total_failures, total_tests

contains

    subroutine run_fsps_data_registry_tests()
        call print_minor_header("fsps_data_registry")

        call test_backend_status_helpers()
        call test_create_backend_hdf5_and_auto()
        call test_create_backend_replaces_existing_instance()
        call test_create_backend_legacy_error()
        call test_create_backend_unknown_error()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_data_registry_tests

    subroutine test_backend_status_helpers()
        type(backend_status_t) :: status

        call print_group("backend_status helper methods")

        call status%clear()
        call assert_int_equals(0, status%code, 'clear keeps code=0 when message absent', total_tests, total_failures)
        call assert_true(.not. allocated(status%message), &
                 'clear leaves message absent when already absent', total_tests, total_failures)

        call status%set_ok()
        call assert_int_equals(0, status%code, 'set_ok sets code=0', total_tests, total_failures)
        call assert_true(allocated(status%message), 'set_ok allocates message', total_tests, total_failures)
        call assert_true(backend_status_ok(status), 'backend_status_ok true for set_ok', total_tests, total_failures)

        call status%set_error(7777, 'sentinel backend error')
        call assert_int_equals(7777, status%code, 'set_error stores code', total_tests, total_failures)
        call assert_true(allocated(status%message), 'set_error stores message', total_tests, total_failures)
        call assert_true(.not. backend_status_ok(status), 'backend_status_ok false for error status', total_tests, total_failures)

        call status%set_error(8888, 'replacement backend error')
        call assert_int_equals(8888, status%code, 'set_error replaces prior error code', total_tests, total_failures)

        call status%set_ok()
        call assert_int_equals(0, status%code, 'set_ok clears prior error code', total_tests, total_failures)

        call status%set_ok()
        call assert_int_equals(0, status%code, 'set_ok remains stable on repeated calls', total_tests, total_failures)

        call status%clear()
        call assert_int_equals(0, status%code, 'clear resets code', total_tests, total_failures)
        call assert_true(.not. allocated(status%message), 'clear deallocates message', total_tests, total_failures)
    end subroutine test_backend_status_helpers

    subroutine test_create_backend_hdf5_and_auto()
        class(data_backend_t), allocatable :: backend
        type(backend_status_t) :: status

        call print_group("create_data_backend: valid modes")

        call create_data_backend('hdf5', backend, status)
        call assert_int_equals(0, status%code, 'hdf5 mode returns OK', total_tests, total_failures)
        call assert_true(allocated(backend), 'hdf5 mode allocates backend', total_tests, total_failures)
        if (allocated(backend)) deallocate(backend)

        call create_data_backend('fsds_hdf5', backend, status)
        call assert_int_equals(0, status%code, 'fsds_hdf5 mode returns OK', total_tests, total_failures)
        call assert_true(allocated(backend), 'fsds_hdf5 mode allocates backend', total_tests, total_failures)
        if (allocated(backend)) deallocate(backend)

        call create_data_backend('auto', backend, status)
        call assert_int_equals(0, status%code, 'auto mode returns OK', total_tests, total_failures)
        call assert_true(allocated(backend), 'auto mode allocates backend', total_tests, total_failures)
        if (allocated(backend)) deallocate(backend)

        call create_data_backend('fsds_auto', backend, status)
        call assert_int_equals(0, status%code, 'fsds_auto mode returns OK', total_tests, total_failures)
        call assert_true(allocated(backend), 'fsds_auto mode allocates backend', total_tests, total_failures)
        if (allocated(backend)) deallocate(backend)
    end subroutine test_create_backend_hdf5_and_auto

    subroutine test_create_backend_replaces_existing_instance()
        class(data_backend_t), allocatable :: backend
        type(backend_status_t) :: status

        call print_group("create_data_backend: replace preallocated backend")

        allocate(hdf5_backend_t :: backend)
        call assert_true(allocated(backend), 'manual pre-allocation succeeds', total_tests, total_failures)

        call create_data_backend('auto', backend, status)
        call assert_int_equals(0, status%code, 'auto mode succeeds when backend already allocated', total_tests, total_failures)
        call assert_true(allocated(backend), 'backend remains allocated after replacement', total_tests, total_failures)

        if (allocated(backend)) deallocate(backend)
    end subroutine test_create_backend_replaces_existing_instance

    subroutine test_create_backend_legacy_error()
        class(data_backend_t), allocatable :: backend
        type(backend_status_t) :: status

        call print_group("create_data_backend: legacy mode errors")

        call create_data_backend('legacy', backend, status)
        call assert_int_equals(4003, status%code, 'legacy mode returns removed-backend error', total_tests, total_failures)
        call assert_true(.not. allocated(backend), 'legacy mode does not allocate backend', total_tests, total_failures)

        call create_data_backend('fsds_legacy', backend, status)
        call assert_int_equals(4003, status%code, 'fsds_legacy mode returns removed-backend error', total_tests, total_failures)
        call assert_true(.not. allocated(backend), 'fsds_legacy mode does not allocate backend', total_tests, total_failures)
    end subroutine test_create_backend_legacy_error

    subroutine test_create_backend_unknown_error()
        class(data_backend_t), allocatable :: backend
        type(backend_status_t) :: status

        call print_group("create_data_backend: unknown mode errors")

        call create_data_backend('definitely_unknown_backend', backend, status)
        call assert_int_equals(4002, status%code, 'unknown mode returns unknown-mode error', total_tests, total_failures)
        call assert_true(.not. allocated(backend), 'unknown mode does not allocate backend', total_tests, total_failures)
    end subroutine test_create_backend_unknown_error

end module test_fsps_data_registry_mod
