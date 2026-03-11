module test_fsps_data_loader_mod
    use fsps_precision, only: WP
    use fsps_context_types, only: fsps_context_t
    use fsps_api, only: fsps_create, fsps_destroy
    use fsps_data_backend, only: backend_status_t
    use fsps_data_loader, only: resolve_data_backend_uri, fsps_data_open, fsps_data_close, fsps_data_query_axis
    use fsps_data_schema, only: axis_desc_t
    use test_fsps_data_backend_hdf5_mod, only: create_test_hdf5_fixture, remove_test_hdf5_fixture, make_backend_uri
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_true, assert_int_equals, assert_float_equals

    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_data_loader_tests, total_failures, total_tests

contains

    subroutine run_fsps_data_loader_tests()
        call print_minor_header("fsps_data_loader")

        call test_resolve_data_backend_uri_defaults()
        call test_resolve_data_backend_uri_dust_override()
        call test_open_query_close_hdf5_backend()
        call test_open_missing_file_error()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_data_loader_tests

    subroutine test_resolve_data_backend_uri_defaults()
        type(fsps_context_t) :: ctx
        character(len=:), allocatable :: uri

        call print_group("resolve_data_backend_uri defaults")

        call fsps_create(ctx)
        ctx%sps_home = '/tmp/fsps_data_loader_test_home'
        ctx%state%isoc_type = 'mist'
        ctx%state%spec_type = 'miles'
        ctx%state%str_dustem = 'DL07'

        call resolve_data_backend_uri(ctx, uri)

        call assert_true(index(uri, '|mist|miles|DL07') > 0, &
                 'default URI carries selected library tuple', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_resolve_data_backend_uri_defaults

    subroutine test_resolve_data_backend_uri_dust_override()
        type(fsps_context_t) :: ctx
        character(len=:), allocatable :: uri

        call print_group("resolve_data_backend_uri dust override")

        call fsps_create(ctx)
        ctx%sps_home = '/tmp/fsps_data_loader_test_home'
        ctx%state%isoc_type = 'mist'
        ctx%state%spec_type = 'miles'
        ctx%state%str_dustem = 'DL07'

        call resolve_data_backend_uri(ctx, uri, 'THEMIS')

        call assert_true(index(uri, '|mist|miles|THEMIS') > 0, &
                 'explicit dust override used in URI tuple', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_resolve_data_backend_uri_dust_override

    subroutine test_open_query_close_hdf5_backend()
        type(backend_status_t) :: status
        type(axis_desc_t) :: axis
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("data loader open/query/close")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'THEMIS')

        call fsps_data_open(uri, status)
        call assert_int_equals(0, status%code, 'fsps_data_open succeeds on fixture file', &
                       total_tests, total_failures)

        call fsps_data_query_axis('lambda', axis, status)
        call assert_int_equals(0, status%code, 'fsps_data_query_axis lambda succeeds', &
                       total_tests, total_failures)
        if (status%code == 0) then
            call assert_true(allocated(axis%values), 'lambda axis values allocated', total_tests, total_failures)
            if (allocated(axis%values)) then
                call assert_int_equals(3, size(axis%values), 'lambda axis size matches fixture', total_tests, total_failures)
                call assert_float_equals(200.0_wp, axis%values(2), 1.0e-12_wp, &
                                         'lambda axis midpoint matches fixture', total_tests, total_failures)
            end if
        end if

        call fsps_data_close(status)
        call assert_int_equals(0, status%code, 'fsps_data_close succeeds', total_tests, total_failures)

        call axis%clear()
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_open_query_close_hdf5_backend

    subroutine test_open_missing_file_error()
        type(backend_status_t) :: status

        call print_group("data loader open missing file")

        call fsps_data_open('/tmp/no_such_fsps_fixture.h5|mist|miles|DL07', status)
        call assert_true(status%code /= 0, 'open missing file returns non-zero status', total_tests, total_failures)

        call fsps_data_close(status)
    end subroutine test_open_missing_file_error

end module test_fsps_data_loader_mod
