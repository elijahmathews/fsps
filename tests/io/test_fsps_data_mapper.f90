module test_fsps_data_mapper_mod
    use fsps_precision, only: WP
    use fsps_data_backend, only: backend_status_t
    use fsps_data_backend_hdf5, only: hdf5_backend_t
    use fsps_data_mapper, only: fsps_data_mapper_t
    use fsps_data_schema, only: dataset_desc_t, library_manifest_t, spectral_grid_t, spectral_slice_t, &
                                dust_emission_t, agn_dust_t, xrb_spectra_t
    use test_fsps_hdf5_fixture_support_mod, only: create_test_hdf5_fixture, remove_test_hdf5_fixture, make_backend_uri
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_true, assert_int_equals, assert_float_equals
    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_data_mapper_tests, total_failures, total_tests

contains

    subroutine run_fsps_data_mapper_tests()
        call print_minor_header("fsps_data_mapper")

        call test_find_bracketing_indices()
        call test_map_spectral_grid_success_and_validation()
        call test_map_spectral_slice_pick_logic()
        call test_map_auxiliary_roles_and_errors()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_data_mapper_tests

    subroutine test_find_bracketing_indices()
        type(fsps_data_mapper_t) :: mapper
        real(WP), dimension(4) :: axis_asc
        real(WP), dimension(4) :: axis_desc
        real(WP), dimension(1) :: axis_single
        integer :: ilo, ihi
        real(WP) :: whi

        call print_group("find_bracketing_indices edge cases")

        axis_asc = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
        axis_desc = [4.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
        axis_single = [5.0_wp]

        call mapper%find_bracketing_indices(axis_single, 5.0_wp, ilo, ihi, whi)
        call assert_int_equals(1, ilo, 'singleton axis -> i_lo=1', total_tests, total_failures)
        call assert_int_equals(1, ihi, 'singleton axis -> i_hi=1', total_tests, total_failures)
        call assert_float_equals(0.0_wp, whi, 1.0e-12_wp, 'singleton axis -> w_hi=0', total_tests, total_failures)

        call mapper%find_bracketing_indices(axis_asc, 2.5_wp, ilo, ihi, whi)
        call assert_int_equals(2, ilo, 'ascending interior lower index', total_tests, total_failures)
        call assert_int_equals(3, ihi, 'ascending interior upper index', total_tests, total_failures)
        call assert_float_equals(0.5_wp, whi, 1.0e-12_wp, 'ascending interior weight', total_tests, total_failures)

        call mapper%find_bracketing_indices(axis_asc, 0.1_wp, ilo, ihi, whi)
        call assert_int_equals(1, ilo, 'ascending low clamp lower index', total_tests, total_failures)
        call assert_int_equals(1, ihi, 'ascending low clamp upper index', total_tests, total_failures)

        call mapper%find_bracketing_indices(axis_desc, 2.5_wp, ilo, ihi, whi)
        call assert_int_equals(2, ilo, 'descending interior lower index', total_tests, total_failures)
        call assert_int_equals(3, ihi, 'descending interior upper index', total_tests, total_failures)
        call assert_float_equals(0.5_wp, whi, 1.0e-12_wp, 'descending interior weight', total_tests, total_failures)
    end subroutine test_find_bracketing_indices

    subroutine test_map_spectral_grid_success_and_validation()
        type(fsps_data_mapper_t) :: mapper
        type(hdf5_backend_t) :: backend
        type(dataset_desc_t) :: dataset
        type(spectral_grid_t) :: grid
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("map_spectral_grid success and validation")

        call dataset%clear()
        call mapper%map_spectral_grid(backend, dataset, grid, status)
        call assert_int_equals(1000, status%code, 'unallocated dims_csv -> 1000', total_tests, total_failures)

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for map_spectral_grid tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call dataset%clear()
        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,alpha_fe,logt,logg'
        dataset%has_valid_mask = .true.
        dataset%valid_mask_path = '/tests/spec5_valid'
        dataset%has_missing_value = .true.
        dataset%missing_value = -999.0_wp

        call mapper%map_spectral_grid(backend, dataset, grid, status)
        call assert_int_equals(0, status%code, 'map_spectral_grid succeeds for rank-5 dataset', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(3, size(grid%flux, 1), 'grid lambda size is 3', total_tests, total_failures)
            call assert_int_equals(2, size(grid%flux, 2), 'grid z size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(grid%flux, 3), 'grid afe size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(grid%flux, 4), 'grid logt size is 2', total_tests, total_failures)
            call assert_int_equals(2, size(grid%flux, 5), 'grid logg size is 2', total_tests, total_failures)
            call assert_float_equals(32222.0_wp, grid%flux(3, 2, 2, 2, 2), 1.0e-12_wp, &
                                     'grid sample value matches fixture formula', total_tests, total_failures)
            call assert_true(allocated(grid%valid), 'grid valid mask allocated', total_tests, total_failures)
            if (allocated(grid%valid)) then
                call assert_true(.not. grid%valid(2, 1, 2, 2), &
                                 'grid valid mask carries expected false cell', total_tests, total_failures)
            end if
            call assert_float_equals(-999.0_wp, grid%missing_value, 1.0e-12_wp, &
                                     'grid missing value propagated', total_tests, total_failures)
        end if

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_map_spectral_grid_success_and_validation

    subroutine test_map_spectral_slice_pick_logic()
        type(fsps_data_mapper_t) :: mapper
        type(hdf5_backend_t) :: backend
        type(dataset_desc_t) :: dataset
        type(spectral_slice_t) :: slice
        type(backend_status_t) :: status
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("map_spectral_slice nearest-index behavior")

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for map_spectral_slice tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        call dataset%clear()
        dataset%path = '/tests/spec4'
        dataset%dims_csv = 'lambda,z,logt,logg'
        call mapper%map_spectral_slice(backend, dataset, 0.019_wp, 0.99_wp, slice, status)
        call assert_int_equals(0, status%code, 'map_spectral_slice succeeds without afe axis', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(2, slice%iz, 'z pick chooses upper bracket when weight >= 0.5', total_tests, total_failures)
            call assert_int_equals(1, slice%iafe, 'no-afe dataset uses iafe=1', total_tests, total_failures)
            call assert_float_equals(32021.0_wp, slice%flux(3,2,1), 1.0e-12_wp, &
                                     'rank-4 slice sample value matches expected z=2,iafe=1', total_tests, total_failures)
        end if

        call dataset%clear()
        dataset%path = '/tests/spec5'
        dataset%dims_csv = 'lambda,z,afe,logt,logg'
        call mapper%map_spectral_slice(backend, dataset, 0.011_wp, 0.19_wp, slice, status)
        call assert_int_equals(0, status%code, 'map_spectral_slice succeeds with afe axis', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(1, slice%iz, 'z pick chooses lower bracket when weight < 0.5', total_tests, total_failures)
            call assert_int_equals(2, slice%iafe, 'afe pick chooses upper bracket when weight >= 0.5', total_tests, total_failures)
            call assert_float_equals(31211.0_wp, slice%flux(3,1,1), 1.0e-12_wp, &
                                     'rank-5 slice sample value matches expected z=1,iafe=2', total_tests, total_failures)
        end if

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_map_spectral_slice_pick_logic

    subroutine test_map_auxiliary_roles_and_errors()
        type(fsps_data_mapper_t) :: mapper
        type(hdf5_backend_t) :: backend
        type(library_manifest_t) :: manifest
        type(backend_status_t) :: status
        type(dust_emission_t) :: dust
        type(agn_dust_t) :: agn
        type(xrb_spectra_t) :: xrb
        character(len=:), allocatable :: file_path, uri
        logical :: ok

        call print_group("mapper role lookup errors and auxiliary table mapping")

        call manifest%clear()
        call mapper%map_dust_emission(backend, manifest, dust, status)
        call assert_int_equals(1301, status%code, 'empty manifest datasets -> 1301', total_tests, total_failures)

        call create_test_hdf5_fixture(file_path, ok)
        call assert_true(ok, 'fixture creation succeeds', total_tests, total_failures)
        if (.not. ok) return

        uri = make_backend_uri(file_path, 'mist', 'miles', 'themis')
        call backend%open(uri, status)
        if (status%code /= 0) then
            call assert_true(.false., 'backend open should succeed for auxiliary mapper tests', total_tests, total_failures)
            call remove_test_hdf5_fixture(file_path)
            return
        end if

        allocate(manifest%datasets(1))
        manifest%datasets(1)%role = 'dust_em_qpah'
        call mapper%map_dust_emission(backend, manifest, dust, status)
        call assert_int_equals(1302, status%code, 'role with missing path -> 1302', total_tests, total_failures)

        call manifest%clear()
        allocate(manifest%datasets(4))
        manifest%datasets(1)%role = 'dust_em_qpah'; manifest%datasets(1)%path = '/libraries/spectra/miles/axes/afe'
        manifest%datasets(2)%role = 'dust_em_umin'; manifest%datasets(2)%path = '/axes/z'
        manifest%datasets(3)%role = 'dust_em_lam';  manifest%datasets(3)%path = '/axes/lambda'
        manifest%datasets(4)%role = 'dust_em_spec'; manifest%datasets(4)%path = '/tests/r3'

        call mapper%map_dust_emission(backend, manifest, dust, status)
        call assert_int_equals(0, status%code, 'map_dust_emission succeeds with role-complete manifest', &
                       total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(2, size(dust%qpah), 'dust qpah size matches mapped axis', total_tests, total_failures)
            call assert_int_equals(2, size(dust%spec, 1), 'dust spec dim1 matches mapped 3D test dataset', &
                                   total_tests, total_failures)
        end if

        call manifest%clear()
        allocate(manifest%datasets(3))
        manifest%datasets(1)%role = 'agn_dust_tau';  manifest%datasets(1)%path = '/tests/r1'
        manifest%datasets(2)%role = 'agn_dust_lam';  manifest%datasets(2)%path = '/axes/lambda'
        manifest%datasets(3)%role = 'agn_dust_spec'; manifest%datasets(3)%path = '/tests/r2'

        call mapper%map_agn_dust(backend, manifest, agn, status)
        call assert_int_equals(0, status%code, 'map_agn_dust succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(3, size(agn%tau), 'agn tau size is 3', total_tests, total_failures)
            call assert_int_equals(2, size(agn%spec, 1), 'agn spec dim1 is 2', total_tests, total_failures)
            call assert_int_equals(3, size(agn%spec, 2), 'agn spec dim2 is 3', total_tests, total_failures)
        end if

        call manifest%clear()
        allocate(manifest%datasets(4))
        manifest%datasets(1)%role = 'xrb_lam';  manifest%datasets(1)%path = '/axes/lambda'
        manifest%datasets(2)%role = 'xrb_age';  manifest%datasets(2)%path = '/axes/age'
        manifest%datasets(3)%role = 'xrb_z';    manifest%datasets(3)%path = '/axes/z'
        manifest%datasets(4)%role = 'xrb_spec'; manifest%datasets(4)%path = '/tests/r3'

        call mapper%map_xrb(backend, manifest, xrb, status)
        call assert_int_equals(0, status%code, 'map_xrb succeeds', total_tests, total_failures)
        if (status%code == 0) then
            call assert_int_equals(3, size(xrb%lam), 'xrb lambda size is 3', total_tests, total_failures)
            call assert_int_equals(2, size(xrb%spec, 1), 'xrb spec dim1 is 2', total_tests, total_failures)
        end if

        call backend%close(status)
        call remove_test_hdf5_fixture(file_path)
    end subroutine test_map_auxiliary_roles_and_errors

end module test_fsps_data_mapper_mod
