module test_fsps_data_backend_hdf5_mod
    use fsps_precision, only: WP
    use fsps_data_backend, only: backend_status_t
    use fsps_data_backend_hdf5, only: hdf5_backend_t, resolve_axis_dataset_path
    use fsps_data_schema, only: axis_desc_t, library_manifest_t, dataset_desc_t, spectral_slice_t, spectral_grid_t
    use test_fsps_hdf5_fixture_support_mod, only: create_test_hdf5_fixture, remove_test_hdf5_fixture, make_backend_uri
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, assert_true
    use test_utils_mod, only: assert_int_equals, assert_float_equals
    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_data_backend_hdf5_tests, total_failures, total_tests

contains

    subroutine run_fsps_data_backend_hdf5_tests()
        call print_minor_header("fsps_data_backend_hdf5")

        call test_resolve_axis_dataset_path_spectral_axes()
        call test_resolve_axis_dataset_path_non_spectral_axis()
        call test_resolve_axis_dataset_path_empty_spec_type()
        call test_resolve_axis_dataset_path_extended()

        call test_open_uri_validation_and_open_state()
        call test_open_missing_file_error()
        call test_has_path_and_query_axis_behavior()
        call test_read_manifest_behavior_and_dust_mapping()

        call test_closed_backend_read_errors()
        call test_real_rank_reads_success()
        call test_int_rank_reads_success()
        call test_rank_mismatch_errors()
        call test_missing_dataset_errors()

        call test_spectral_closed_backend_and_input_validation()
        call test_read_spectral_slice_rank5_with_mask()
        call test_read_spectral_slice_rank4_and_rank_checks()
        call test_read_spectral_neighborhood_rank5_with_mask()
        call test_read_spectral_neighborhood_rank4_and_rank_checks()
        call test_valid_mask_path_missing_errors()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_data_backend_hdf5_tests

    subroutine test_resolve_axis_dataset_path_spectral_axes()
        character(len=:), allocatable :: path

        call print_group("Axis Path: Spectral Axis Resolution")

        path = resolve_axis_dataset_path('lambda', 'BaSeL')
        call assert_true(path == '/libraries/spectra/basel/axes/lambda', &
                         'lambda path uses library-local spectral axis', total_tests, total_failures)

        path = resolve_axis_dataset_path('z', 'c3k_afe+0.0')
        call assert_true(path == '/libraries/spectra/c3k_afe+0.0/axes/z', &
                 'z path uses library-local spectral axis', total_tests, total_failures)
    end subroutine test_resolve_axis_dataset_path_spectral_axes

    subroutine test_resolve_axis_dataset_path_non_spectral_axis()
        character(len=:), allocatable :: path

        call print_group("Axis Path: Global Fallback for Non-Spectral")

        path = resolve_axis_dataset_path('age', 'miles')
        call assert_true(path == '/axes/age', &
                         'non-spectral axis resolves to global axis path', total_tests, total_failures)
    end subroutine test_resolve_axis_dataset_path_non_spectral_axis

    subroutine test_resolve_axis_dataset_path_empty_spec_type()
        character(len=:), allocatable :: path

        call print_group("Axis Path: Empty Spec Fallback")

        path = resolve_axis_dataset_path('logt', '')
        call assert_true(path == '/axes/logt', &
                         'empty spec type resolves to global axis path', total_tests, total_failures)
    end subroutine test_resolve_axis_dataset_path_empty_spec_type

    subroutine test_resolve_axis_dataset_path_extended()
        character(len=:), allocatable :: path

        call print_group("Axis path resolver edge cases")

        path = resolve_axis_dataset_path('Alpha_Fe', 'MILES')
        call assert_true(path == '/libraries/spectra/miles/axes/alpha_fe', &
                         'alpha_fe resolves to spectral-local lower-case path', total_tests, total_failures)

        path = resolve_axis_dataset_path('  LOGT ', '  Miles  ')
        call assert_true(path == '/libraries/spectra/miles/axes/logt', &
                 'axis/spec names are trimmed and lowercased', total_tests, total_failures)

        path = resolve_axis_dataset_path('age', 'MILES')
        call assert_true(path == '/axes/age', 'non-spectral axis stays global', total_tests, total_failures)
    end subroutine test_resolve_axis_dataset_path_extended

    subroutine test_open_uri_validation_and_open_state()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Open/close URI validation and backend state")

        call backend%open('bad_uri', status)
        call assert_int_equals(2005, status%code, 'malformed URI returns 2005', total_tests, total_failures)
        call assert_true(.not. backend%is_open(), 'backend remains closed after malformed URI', total_tests, total_failures)

        call backend%open('/tmp/x.h5|mist|miles', status)
        call assert_int_equals(2005, status%code, 'URI with too few fields returns 2005', total_tests, total_failures)

        call backend%open('/tmp/x.h5||miles|themis', status)
        call assert_int_equals(2005, status%code, 'URI with empty middle field returns parser error', total_tests, total_failures)

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'MIST', 'MILES', 'themis')
        call backend%open(uri, status)
        call assert_int_equals(0, status%code, 'valid URI opens fixture file', total_tests, total_failures)
        call assert_true(backend%is_open(), 'backend open state is true after open', total_tests, total_failures)

        call backend%close(status)
        call assert_int_equals(0, status%code, 'close returns success', total_tests, total_failures)
        call assert_true(.not. backend%is_open(), 'backend open state false after close', total_tests, total_failures)

        call backend%close(status)
        call assert_int_equals(0, status%code, 'close is idempotent on already-closed backend', total_tests, total_failures)

        call remove_test_hdf5_fixture(file_path)
    end subroutine test_open_uri_validation_and_open_state

    subroutine test_open_missing_file_error()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: uri

        call print_group("Open missing file error path")

        uri = make_backend_uri('/tmp/this_file_does_not_exist_fsps_test.h5', 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        call assert_int_equals(2002, status%code, 'missing HDF5 file returns 2002', total_tests, total_failures)
        call assert_true(.not. backend%is_open(), 'backend remains closed when file open fails', total_tests, total_failures)
    end subroutine test_open_missing_file_error

    subroutine test_has_path_and_query_axis_behavior()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(axis_desc_t) :: axis
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("has_path and query_axis resolution behavior")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for query-axis tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call assert_true(backend%has_path('/libraries/spectra/miles/axes/lambda'), &
                         'has_path true for existing path', total_tests, total_failures)
        call assert_true(.not. backend%has_path('/libraries/spectra/miles/axes/not_real'), &
                         'has_path false for missing path', total_tests, total_failures)

        call backend%query_axis('lambda', axis, status)
        call assert_int_equals(0, status%code, 'query_axis lambda succeeds', total_tests, total_failures)
        call assert_int_equals(3, axis%n, 'lambda axis length is 3', total_tests, total_failures)
        call assert_true(axis%path == '/libraries/spectra/miles/axes/lambda', &
                 'lambda query uses spectral-local path', total_tests, total_failures)
        call assert_float_equals(200.0_wp, axis%values(2), 1.0e-12_wp, &
                     'lambda axis value(2) matches fixture', total_tests, total_failures)

        call backend%query_axis('age', axis, status)
        call assert_int_equals(0, status%code, 'query_axis age succeeds', total_tests, total_failures)
        call assert_true(axis%path == '/axes/age', 'age query uses global axis path', total_tests, total_failures)

        call backend%query_axis('nebular_lambda', axis, status)
        call assert_int_equals(0, status%code, 'query_axis nebular_lambda succeeds', total_tests, total_failures)
        call assert_true(axis%path == '/axes/nebular_lambda', &
                 'nebular_lambda query uses global axis path', total_tests, total_failures)
        call assert_float_equals(205.0_wp, axis%values(2), 1.0e-12_wp, &
                 'nebular_lambda axis sample value matches fixture', total_tests, total_failures)

        call backend%query_axis('nonexistent_axis', axis, status)
        call assert_int_equals(2202, status%code, 'missing axis returns 2202', total_tests, total_failures)

        call backend%close(status)
        call backend%query_axis('lambda', axis, status)
        call assert_int_equals(2201, status%code, 'query_axis on closed backend returns 2201', total_tests, total_failures)

        call remove_test_hdf5_fixture(file_path)
    end subroutine test_has_path_and_query_axis_behavior

    subroutine test_read_manifest_behavior_and_dust_mapping()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(library_manifest_t) :: manifest
        character(len=:), allocatable :: file_path, uri
        character(len=:), allocatable :: path12, path47
        logical :: ok

        call print_group("Manifest behavior and legacy-compatible dust mapping")

        call backend%read_manifest(manifest, status)
        call assert_int_equals(2101, status%code, 'read_manifest on closed backend returns 2101', total_tests, total_failures)

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'MIST', 'MILES', 'THEMIS')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for manifest tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call backend%read_manifest(manifest, status)
        call assert_int_equals(0, status%code, 'read_manifest succeeds', total_tests, total_failures)
        call assert_true(trim(manifest%isoc_name) == 'mist', 'manifest stores lower-cased isoc name', total_tests, total_failures)
        call assert_true(trim(manifest%spec_name) == 'miles', 'manifest stores lower-cased spec name', total_tests, total_failures)
        call assert_true(trim(manifest%dust_name) == 'THEMIS', &
                 'dust type THEMIS preserved in manifest', total_tests, total_failures)
        call assert_int_equals(61, size(manifest%datasets), 'manifest declares expected dataset count', total_tests, total_failures)

        path12 = dataset_path_for_role(manifest, 'nebular_wd_line_pos')
        path47 = dataset_path_for_role(manifest, 'dust_em_qpah')
        call assert_true(path12 == '/libraries/nebular/mist/WD/line_pos', &
                         'nebular WD path prefers isoc-specific location when present', total_tests, total_failures)
        call assert_true(index(path47, '/libraries/dust/emission/THEMIS/qpah') == 1, &
                         'dust emission path uses THEMIS branch', total_tests, total_failures)

        call backend%close(status)

        uri = make_backend_uri(file_path, 'padova', 'miles', 'not_themis')
        call backend%open(uri, status)
        call backend%read_manifest(manifest, status)
        call assert_int_equals(0, status%code, 'read_manifest succeeds for non-THEMIS dust input', total_tests, total_failures)
        call assert_true(trim(manifest%dust_name) == 'DL07', 'non-THEMIS dust input maps to DL07', total_tests, total_failures)

        path12 = dataset_path_for_role(manifest, 'nebular_wd_line_pos')
        path47 = dataset_path_for_role(manifest, 'dust_em_qpah')
        call assert_true(path12 == '/libraries/nebular/WD/line_pos', &
                         'nebular WD path falls back to global when isoc-specific path absent', total_tests, total_failures)
        call assert_true(index(path47, '/libraries/dust/emission/DL07/qpah') == 1, &
                         'dust emission path uses DL07 branch', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_read_manifest_behavior_and_dust_mapping

    subroutine test_closed_backend_read_errors()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        real(WP), allocatable :: r1(:), r2(:,:), r3(:,:,:), r4(:,:,:,:)
        integer, allocatable :: i1(:), i2(:,:), i3(:,:,:)

        call print_group("Closed-backend read error codes")

        call backend%read_real_1d('/tests/r1', r1, status)
        call assert_int_equals(2301, status%code, 'read_real_1d on closed backend -> 2301', total_tests, total_failures)

        call backend%read_real_2d('/tests/r2', r2, status)
        call assert_int_equals(2311, status%code, 'read_real_2d on closed backend -> 2311', total_tests, total_failures)

        call backend%read_real_3d('/tests/r3', r3, status)
        call assert_int_equals(2317, status%code, 'read_real_3d on closed backend -> 2317', total_tests, total_failures)

        call backend%read_real_4d('/tests/r4', r4, status)
        call assert_int_equals(2329, status%code, 'read_real_4d on closed backend -> 2329', total_tests, total_failures)

        call backend%read_int_1d('/tests/i1', i1, status)
        call assert_int_equals(2321, status%code, 'read_int_1d on closed backend -> 2321', total_tests, total_failures)

        call backend%read_int_2d('/tests/i2', i2, status)
        call assert_int_equals(2331, status%code, 'read_int_2d on closed backend -> 2331', total_tests, total_failures)

        call backend%read_int_3d('/tests/i3', i3, status)
        call assert_int_equals(2337, status%code, 'read_int_3d on closed backend -> 2337', total_tests, total_failures)
    end subroutine test_closed_backend_read_errors

    subroutine test_real_rank_reads_success()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok
        real(WP), allocatable :: r1(:), r2(:,:), r3(:,:,:), r4(:,:,:,:)

        call print_group("Real rank-1..4 reads from fixture")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for real-read tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call backend%read_real_1d('/tests/r1', r1, status)
        call assert_int_equals(0, status%code, 'read_real_1d succeeds', total_tests, total_failures)
        call assert_int_equals(3, size(r1), 'r1 size is 3', total_tests, total_failures)
        call assert_float_equals(3.0_wp, r1(3), 1.0e-12_wp, 'r1 sample value matches', total_tests, total_failures)

        call backend%read_real_2d('/tests/r2', r2, status)
        call assert_int_equals(0, status%code, 'read_real_2d succeeds', total_tests, total_failures)
        call assert_int_equals(2, size(r2,1), 'r2 dim1 is 2', total_tests, total_failures)
        call assert_int_equals(3, size(r2,2), 'r2 dim2 is 3', total_tests, total_failures)
        call assert_float_equals(6.0_wp, r2(2,3), 1.0e-12_wp, 'r2 sample value matches', total_tests, total_failures)

        call backend%read_real_3d('/tests/r3', r3, status)
        call assert_int_equals(0, status%code, 'read_real_3d succeeds', total_tests, total_failures)
        call assert_int_equals(2, size(r3,1), 'r3 dim1 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(r3,2), 'r3 dim2 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(r3,3), 'r3 dim3 is 2', total_tests, total_failures)
        call assert_float_equals(8.0_wp, r3(2,2,2), 1.0e-12_wp, 'r3 sample value matches', total_tests, total_failures)

        call backend%read_real_4d('/tests/r4', r4, status)
        call assert_int_equals(0, status%code, 'read_real_4d succeeds', total_tests, total_failures)
        call assert_int_equals(2, size(r4,1), 'r4 dim1 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(r4,2), 'r4 dim2 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(r4,3), 'r4 dim3 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(r4,4), 'r4 dim4 is 2', total_tests, total_failures)
        call assert_float_equals(16.0_wp, r4(2,2,2,2), 1.0e-12_wp, 'r4 sample value matches', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_real_rank_reads_success

    subroutine test_int_rank_reads_success()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok
        integer, allocatable :: i1(:), i2(:,:), i3(:,:,:)

        call print_group("Integer rank-1..3 reads from fixture")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for int-read tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call backend%read_int_1d('/tests/i1', i1, status)
        call assert_int_equals(0, status%code, 'read_int_1d succeeds', total_tests, total_failures)
        call assert_int_equals(3, size(i1), 'i1 size is 3', total_tests, total_failures)
        call assert_int_equals(6, i1(3), 'i1 sample value matches', total_tests, total_failures)

        call backend%read_int_2d('/tests/i2', i2, status)
        call assert_int_equals(0, status%code, 'read_int_2d succeeds', total_tests, total_failures)
        call assert_int_equals(2, size(i2,1), 'i2 dim1 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(i2,2), 'i2 dim2 is 2', total_tests, total_failures)
        call assert_int_equals(4, i2(2,2), 'i2 sample value matches', total_tests, total_failures)

        call backend%read_int_3d('/tests/i3', i3, status)
        call assert_int_equals(0, status%code, 'read_int_3d succeeds', total_tests, total_failures)
        call assert_int_equals(2, size(i3,1), 'i3 dim1 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(i3,2), 'i3 dim2 is 2', total_tests, total_failures)
        call assert_int_equals(2, size(i3,3), 'i3 dim3 is 2', total_tests, total_failures)
        call assert_int_equals(8, i3(2,2,2), 'i3 sample value matches', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_int_rank_reads_success

    subroutine test_rank_mismatch_errors()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok
        real(WP), allocatable :: r1(:), r2(:,:), r3(:,:,:), r4(:,:,:,:)
        integer, allocatable :: i1(:), i2(:,:), i3(:,:,:)

        call print_group("Rank mismatch error-path coverage")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for mismatch tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call backend%read_real_1d('/tests/r2', r1, status)
        call assert_int_equals(2304, status%code, 'rank-2 passed to read_real_1d -> 2304', total_tests, total_failures)

        call backend%read_real_2d('/tests/r3', r2, status)
        call assert_int_equals(2314, status%code, 'rank-3 passed to read_real_2d -> 2314', total_tests, total_failures)

        call backend%read_real_3d('/tests/r4', r3, status)
        call assert_int_equals(2320, status%code, 'rank-4 passed to read_real_3d -> 2320', total_tests, total_failures)

        call backend%read_real_4d('/tests/r3', r4, status)
        call assert_int_equals(2332, status%code, 'rank-3 passed to read_real_4d -> 2332', total_tests, total_failures)

        call backend%read_int_1d('/tests/i2', i1, status)
        call assert_int_equals(2324, status%code, 'rank-2 passed to read_int_1d -> 2324', total_tests, total_failures)

        call backend%read_int_2d('/tests/i3', i2, status)
        call assert_int_equals(2334, status%code, 'rank-3 passed to read_int_2d -> 2334', total_tests, total_failures)

        call backend%read_int_3d('/tests/i2', i3, status)
        call assert_int_equals(2340, status%code, 'rank-2 passed to read_int_3d -> 2340', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_rank_mismatch_errors

    subroutine test_missing_dataset_errors()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok
        real(WP), allocatable :: r1(:), r2(:,:), r3(:,:,:), r4(:,:,:,:)
        integer, allocatable :: i1(:), i2(:,:), i3(:,:,:)

        call print_group("Missing dataset open failures")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for missing-dataset tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call backend%read_real_1d('/tests/does_not_exist', r1, status)
        call assert_int_equals(2302, status%code, 'missing dataset in read_real_1d -> 2302', total_tests, total_failures)

        call backend%read_real_2d('/tests/does_not_exist', r2, status)
        call assert_int_equals(2312, status%code, 'missing dataset in read_real_2d -> 2312', total_tests, total_failures)

        call backend%read_real_3d('/tests/does_not_exist', r3, status)
        call assert_int_equals(2318, status%code, 'missing dataset in read_real_3d -> 2318', total_tests, total_failures)

        call backend%read_real_4d('/tests/does_not_exist', r4, status)
        call assert_int_equals(2330, status%code, 'missing dataset in read_real_4d -> 2330', total_tests, total_failures)

        call backend%read_int_1d('/tests/does_not_exist', i1, status)
        call assert_int_equals(2322, status%code, 'missing dataset in read_int_1d -> 2322', total_tests, total_failures)

        call backend%read_int_2d('/tests/does_not_exist', i2, status)
        call assert_int_equals(2332, status%code, 'missing dataset in read_int_2d -> 2332', total_tests, total_failures)

        call backend%read_int_3d('/tests/does_not_exist', i3, status)
        call assert_int_equals(2338, status%code, 'missing dataset in read_int_3d -> 2338', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_missing_dataset_errors

    subroutine test_spectral_closed_backend_and_input_validation()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_slice_t) :: slice
        type(spectral_grid_t) :: neighborhood
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Spectral API closed-backend and input guards")

        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'

        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2401, status%code, 'slice on closed backend -> 2401', total_tests, total_failures)

        call backend%read_spectral_neighborhood(dataset, 1, 1, 1, 1, neighborhood, status)
        call assert_int_equals(2501, status%code, 'neighborhood on closed backend -> 2501', total_tests, total_failures)

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for input-validation tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call dataset%clear()
        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2402, status%code, 'slice with unallocated path -> 2402', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec5'
        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2403, status%code, 'slice with unallocated dims_csv -> 2403', total_tests, total_failures)

        call dataset%clear()
        call backend%read_spectral_neighborhood(dataset, 1, 1, 1, 1, neighborhood, status)
        call assert_int_equals(2502, status%code, 'neighborhood with unallocated path -> 2502', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec5'
        call backend%read_spectral_neighborhood(dataset, 1, 1, 1, 1, neighborhood, status)
        call assert_int_equals(2503, status%code, 'neighborhood with unallocated dims_csv -> 2503', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_spectral_closed_backend_and_input_validation

    subroutine test_read_spectral_slice_rank5_with_mask()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_slice_t) :: slice
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Rank-5 spectral slice reads and bounds")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for rank-5 slice tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        dataset%has_valid_mask = .true.
        dataset%valid_mask_path = '/tests/spec5_valid'
        dataset%has_missing_value = .true.
        dataset%missing_value = -999.0_wp

        call backend%read_spectral_slice(dataset, 2, 1, slice, status)
        call assert_int_equals(0, status%code, 'rank-5 slice read succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(3, size(slice%flux,1), 'slice lambda size is 3', total_tests, total_failures)
            call assert_int_equals(2, size(slice%flux,2), 'slice logt size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(slice%flux,3), 'slice logg size is 2', total_tests, total_failures)
            call assert_float_equals(32111.0_wp, slice%flux(3,1,1), 1.0e-12_wp, &
                         'slice sample value matches formula', total_tests, total_failures)
            call assert_true(allocated(slice%valid), 'slice valid-mask allocated', total_tests, total_failures)
            if (allocated(slice%valid)) then
                call assert_true(.not. slice%valid(2,2), &
                     'slice valid-mask carries expected false cell', total_tests, total_failures)
            end if
            call assert_float_equals(-999.0_wp, slice%missing_value, 1.0e-12_wp, &
                         'slice missing_value propagated', total_tests, total_failures)
        end if

        call backend%read_spectral_slice(dataset, 0, 1, slice, status)
        call assert_int_equals(2410, status%code, 'slice z bound check -> 2410', total_tests, total_failures)

        call backend%read_spectral_slice(dataset, 1, 3, slice, status)
        call assert_int_equals(2411, status%code, 'slice afe bound check -> 2411', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_read_spectral_slice_rank5_with_mask

    subroutine test_read_spectral_slice_rank4_and_rank_checks()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_slice_t) :: slice
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Rank-4 spectral slice behavior and rank compatibility")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for rank-4 slice tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        dataset%path = '/tests/spec4'
        dataset%dims_csv = 'lambda,z,logt,logg'
        dataset%has_valid_mask = .true.
        dataset%valid_mask_path = '/tests/spec4_valid'

        call backend%read_spectral_slice(dataset, 2, 1, slice, status)
        call assert_int_equals(0, status%code, 'rank-4 slice read succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_float_equals(32021.0_wp, slice%flux(3,2,1), 1.0e-12_wp, &
                         'rank-4 slice sample value matches formula', total_tests, total_failures)
            call assert_true(allocated(slice%valid), 'rank-4 slice valid-mask allocated', total_tests, total_failures)
            if (allocated(slice%valid)) then
                call assert_true(.not. slice%valid(2,1), &
                     'rank-4 valid-mask carries expected false cell', total_tests, total_failures)
            end if
        end if

        call backend%read_spectral_slice(dataset, 2, 2, slice, status)
        call assert_int_equals(2415, status%code, 'rank-4 slice requires iafe=1 -> 2415', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec4'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2408, status%code, 'afe dims on rank-4 dataset -> 2408', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,logt,logg'
        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2409, status%code, 'no-afe dims on rank-5 dataset -> 2409', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_read_spectral_slice_rank4_and_rank_checks

    subroutine test_read_spectral_neighborhood_rank5_with_mask()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_grid_t) :: neighborhood
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Rank-5 spectral neighborhood reads and bounds")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for rank-5 neighborhood tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        dataset%has_valid_mask = .true.
        dataset%valid_mask_path = '/tests/spec5_valid'

        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 2, neighborhood, status)
        call assert_int_equals(0, status%code, 'rank-5 neighborhood read succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(3, size(neighborhood%flux,1), 'neighborhood lambda size is 3', total_tests, total_failures)
            call assert_int_equals(2, size(neighborhood%flux,2), 'neighborhood z size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(neighborhood%flux,3), 'neighborhood afe size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(neighborhood%flux,4), 'neighborhood logt size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(neighborhood%flux,5), 'neighborhood logg size is 2', total_tests, total_failures)
            call assert_float_equals(32122.0_wp, neighborhood%flux(3,2,1,2,2), 1.0e-12_wp, &
                         'neighborhood sample value matches formula', total_tests, total_failures)
            call assert_true(allocated(neighborhood%valid), 'neighborhood valid-mask allocated', total_tests, total_failures)
            if (allocated(neighborhood%valid)) then
                call assert_true(.not. neighborhood%valid(2,1,2,2), &
                     'neighborhood valid-mask carries expected false cell', total_tests, total_failures)
            end if
        end if

        call backend%read_spectral_neighborhood(dataset, 0, 1, 1, 1, neighborhood, status)
        call assert_int_equals(2508, status%code, 'invalid neighborhood z bounds -> 2508', total_tests, total_failures)

        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 3, neighborhood, status)
        call assert_int_equals(2510, status%code, 'invalid neighborhood afe bounds -> 2510', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_read_spectral_neighborhood_rank5_with_mask

    subroutine test_read_spectral_neighborhood_rank4_and_rank_checks()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_grid_t) :: neighborhood
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Rank-4 spectral neighborhood behavior and rank compatibility")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for rank-4 neighborhood tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        dataset%path = '/tests/spec4'
        dataset%dims_csv = 'lambda,z,logt,logg'
        dataset%has_valid_mask = .true.
        dataset%valid_mask_path = '/tests/spec4_valid'

        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 1, neighborhood, status)
        call assert_int_equals(0, status%code, 'rank-4 neighborhood read succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(1, size(neighborhood%flux,3), &
                           'rank-4 neighborhood expands afe as degenerate dim', total_tests, total_failures)
            call assert_float_equals(31021.0_wp, neighborhood%flux(3,1,1,2,1), 1.0e-12_wp, &
                         'rank-4 neighborhood sample value matches formula', total_tests, total_failures)
            call assert_true(allocated(neighborhood%valid), 'rank-4 neighborhood valid-mask allocated', total_tests, total_failures)
            if (allocated(neighborhood%valid)) then
                call assert_true(.not. neighborhood%valid(2,1,2,1), &
                     'rank-4 neighborhood mask carries expected false cell', total_tests, total_failures)
            end if
        end if

        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 2, neighborhood, status)
        call assert_int_equals(2515, status%code, &
                       'rank-4 neighborhood requires iafe bounds [1,1] -> 2515', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec4'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 1, neighborhood, status)
        call assert_int_equals(2509, status%code, 'afe dims on rank-4 neighborhood dataset -> 2509', total_tests, total_failures)

        call dataset%clear()
        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,logt,logg'
        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 1, neighborhood, status)
        call assert_int_equals(2514, status%code, 'no-afe dims on rank-5 neighborhood dataset -> 2514', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_read_spectral_neighborhood_rank4_and_rank_checks

    subroutine test_valid_mask_path_missing_errors()
        type(hdf5_backend_t) :: backend
        type(backend_status_t) :: status
        type(dataset_desc_t) :: dataset
        type(spectral_slice_t) :: slice
        type(spectral_grid_t) :: neighborhood
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("Validity-mask metadata guards")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for valid-mask guard tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        dataset%has_valid_mask = .true.

        call backend%read_spectral_slice(dataset, 1, 1, slice, status)
        call assert_int_equals(2601, status%code, 'slice with has_valid_mask but missing path -> 2601', total_tests, total_failures)

        call backend%read_spectral_neighborhood(dataset, 1, 2, 1, 2, neighborhood, status)
        call assert_int_equals(2701, status%code, &
                       'neighborhood with has_valid_mask but missing path -> 2701', total_tests, total_failures)

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_valid_mask_path_missing_errors

    function dataset_path_for_role(manifest, role) result(path)
        type(library_manifest_t), intent(in) :: manifest
        character(len=*), intent(in) :: role
        character(len=:), allocatable :: path
        integer :: i

        path = ''
        if (.not. allocated(manifest%datasets)) return
        do i = 1, size(manifest%datasets)
            if (allocated(manifest%datasets(i)%role)) then
                if (trim(manifest%datasets(i)%role) == trim(role)) then
                    if (allocated(manifest%datasets(i)%path)) then
                        path = manifest%datasets(i)%path
                        return
                    end if
                end if
            end if
        end do
    end function dataset_path_for_role

end module test_fsps_data_backend_hdf5_mod
