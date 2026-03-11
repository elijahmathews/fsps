module test_fsps_data_schema_mod
    use fsps_precision, only: WP
    use fsps_data_schema, only: axis_desc_t, dataset_desc_t, library_manifest_t, spectral_grid_t, spectral_slice_t, &
                                isochrone_grid_t, nebular_grid_t, aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, &
                                dust_emission_t, agn_dust_t, dust_attenuation_t, xrb_spectra_t
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, assert_true, &
                              assert_int_equals, assert_float_equals

    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_data_schema_tests, total_failures, total_tests

contains

    subroutine run_fsps_data_schema_tests()
        call print_minor_header("fsps_data_schema")

        call test_axis_desc_clear()
        call test_dataset_desc_clear_and_has_axis()
        call test_library_manifest_clear()
        call test_spectral_grid_clear()
        call test_spectral_slice_clear()
        call test_isochrone_grid_clear()
        call test_nebular_grid_clear()
        call test_auxiliary_grid_clears()
        call test_dust_and_xrb_clears()
        call test_schema_copy_assignment_paths()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_data_schema_tests

    subroutine test_axis_desc_clear()
        type(axis_desc_t) :: axis

        call print_group("axis_desc clear")

        axis%name = 'lambda'
        axis%unit = 'Angstrom'
        axis%path = '/axes/lambda'
        axis%n = 3
        allocate(axis%values(3))
        axis%values = [100.0_wp, 200.0_wp, 300.0_wp]

        call axis%clear()

        call assert_true(.not. allocated(axis%name), 'axis name deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(axis%unit), 'axis unit deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(axis%path), 'axis path deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(axis%values), 'axis values deallocated', total_tests, total_failures)
        call assert_int_equals(0, axis%n, 'axis n reset', total_tests, total_failures)
    end subroutine test_axis_desc_clear

    subroutine test_dataset_desc_clear_and_has_axis()
        type(dataset_desc_t) :: ds

        call print_group("dataset_desc clear and has_axis")

        ds%role = 'spectral_base'
        ds%path = '/libraries/spectra/miles/base/spectral_grid_nd'
        ds%dims_csv = 'lambda,z,alpha_fe,logt,logg'
        ds%unit = 'Lsun/Hz/Msun'
        ds%transform = 'none'
        ds%dtype = 'float64'
        ds%representation = 'dense_nd'
        allocate(ds%shape(5))
        ds%shape = [3, 2, 2, 2, 2]
        ds%required = .false.
        ds%has_missing_value = .true.
        ds%missing_value = -99.0_wp
        ds%has_valid_mask = .true.
        ds%valid_mask_path = '/tests/spec5_valid'

        call assert_true(ds%has_axis('z'), 'has_axis z', total_tests, total_failures)
        call assert_true(ds%has_axis('Alpha_Fe'), 'has_axis alpha_fe case-insensitive', total_tests, total_failures)
        call assert_true(ds%has_axis('logg'), 'has_axis logg', total_tests, total_failures)
        call assert_true(.not. ds%has_axis('age'), 'has_axis age false', total_tests, total_failures)

        call ds%clear()

        call assert_true(.not. allocated(ds%role), 'dataset role deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(ds%path), 'dataset path deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(ds%dims_csv), 'dataset dims_csv deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(ds%shape), 'dataset shape deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(ds%valid_mask_path), 'dataset mask path deallocated', total_tests, total_failures)
        call assert_true(ds%required, 'dataset required reset true', total_tests, total_failures)
        call assert_true(.not. ds%has_missing_value, 'dataset has_missing_value reset false', total_tests, total_failures)
        call assert_true(.not. ds%has_valid_mask, 'dataset has_valid_mask reset false', total_tests, total_failures)
        call assert_float_equals(-1.0e30_wp, ds%missing_value, 1.0e-12_wp, &
                     'dataset missing_value reset', total_tests, total_failures)
        call assert_true(.not. ds%has_axis('z'), 'has_axis false when dims_csv not allocated', total_tests, total_failures)
    end subroutine test_dataset_desc_clear_and_has_axis

    subroutine test_library_manifest_clear()
        type(library_manifest_t) :: manifest

        call print_group("library_manifest clear")

        manifest%fsds_version = '1.0'
        manifest%isoc_name = 'mist'
        manifest%spec_name = 'miles'
        manifest%dust_name = 'THEMIS'
        manifest%producer = 'test'

        allocate(manifest%datasets(2))
        manifest%datasets(1)%role = 'spectral_base'
        manifest%datasets(2)%role = 'isochrone_track'

        allocate(manifest%axes(2))
        manifest%axes(1)%name = 'lambda'
        manifest%axes(2)%name = 'z'

        call manifest%clear()

        call assert_true(.not. allocated(manifest%datasets), 'manifest datasets deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%axes), 'manifest axes deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%fsds_version), 'manifest fsds_version deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%isoc_name), 'manifest isoc_name deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%spec_name), 'manifest spec_name deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%dust_name), 'manifest dust_name deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(manifest%producer), 'manifest producer deallocated', total_tests, total_failures)
    end subroutine test_library_manifest_clear

    subroutine test_spectral_grid_clear()
        type(spectral_grid_t) :: grid

        call print_group("spectral_grid clear")

        allocate(grid%flux(2,2,2,2,2))
        allocate(grid%valid(2,2,2,2))
        allocate(grid%axis_lambda(2), grid%axis_z(2), grid%axis_afe(2), grid%axis_logt(2), grid%axis_logg(2))
        grid%missing_value = -99.0_wp

        call grid%clear()

        call assert_true(.not. allocated(grid%flux), 'spectral flux deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%valid), 'spectral valid deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%axis_lambda), 'axis_lambda deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%axis_z), 'axis_z deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%axis_afe), 'axis_afe deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%axis_logt), 'axis_logt deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%axis_logg), 'axis_logg deallocated', total_tests, total_failures)
        call assert_float_equals(-1.0e30_wp, grid%missing_value, 1.0e-12_wp, &
                     'spectral missing_value reset', total_tests, total_failures)
    end subroutine test_spectral_grid_clear

    subroutine test_spectral_slice_clear()
        type(spectral_slice_t) :: slice

        call print_group("spectral_slice clear")

        allocate(slice%flux(2,2,2))
        allocate(slice%valid(2,2))
        slice%iz = 2
        slice%iafe = 3
        slice%missing_value = -123.0_wp

        call slice%clear()

        call assert_true(.not. allocated(slice%flux), 'slice flux deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(slice%valid), 'slice valid deallocated', total_tests, total_failures)
        call assert_int_equals(0, slice%iz, 'slice iz reset', total_tests, total_failures)
        call assert_int_equals(0, slice%iafe, 'slice iafe reset', total_tests, total_failures)
        call assert_float_equals(-1.0e30_wp, slice%missing_value, 1.0e-12_wp, &
                     'slice missing_value reset', total_tests, total_failures)
    end subroutine test_spectral_slice_clear

    subroutine test_isochrone_grid_clear()
        type(isochrone_grid_t) :: grid

        call print_group("isochrone_grid clear")

        allocate(grid%nmass(2,2), grid%timestep_logyr(2,2))
        allocate(grid%mini(2,2,2), grid%mact(2,2,2), grid%logl(2,2,2), grid%logt(2,2,2), grid%logg(2,2,2))
        allocate(grid%phase(2,2,2), grid%ffco(2,2,2), grid%lmdot(2,2,2))

        call grid%clear()

        call assert_true(.not. allocated(grid%nmass), 'nmass deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%timestep_logyr), 'timestep deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%mini), 'mini deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%mact), 'mact deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%logl), 'logl deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%logt), 'logt deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%logg), 'logg deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%phase), 'phase deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%ffco), 'ffco deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%lmdot), 'lmdot deallocated', total_tests, total_failures)
    end subroutine test_isochrone_grid_clear

    subroutine test_nebular_grid_clear()
        type(nebular_grid_t) :: grid

        call print_group("nebular_grid clear")

        allocate(grid%line_pos(2), grid%logz(2), grid%age(2), grid%logu(2))
        allocate(grid%cont(2,2,2,2), grid%line(2,2,2,2))
        grid%missing_value = -99.0_wp

        call grid%clear()

        call assert_true(.not. allocated(grid%line_pos), 'line_pos deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%logz), 'logz deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%age), 'age deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%logu), 'logu deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%cont), 'cont deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(grid%line), 'line deallocated', total_tests, total_failures)
        call assert_float_equals(-1.0e30_wp, grid%missing_value, 1.0e-12_wp, &
                     'nebular missing_value reset', total_tests, total_failures)
    end subroutine test_nebular_grid_clear

    subroutine test_auxiliary_grid_clears()
        type(aux_wmbasic_t) :: wmb
        type(aux_pagb_t) :: pagb
        type(aux_wr_t) :: wr
        type(aux_agb_t) :: agb

        call print_group("auxiliary grid clears")

        allocate(wmb%logt(2), wmb%z(2), wmb%lam(2), wmb%spec(2,2,2,2))
        call wmb%clear()
        call assert_true(.not. allocated(wmb%logt), 'wmb logt deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wmb%z), 'wmb z deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wmb%lam), 'wmb lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wmb%spec), 'wmb spec deallocated', total_tests, total_failures)

        allocate(pagb%logt(2), pagb%lam(2), pagb%spec(2,2,2))
        call pagb%clear()
        call assert_true(.not. allocated(pagb%logt), 'pagb logt deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(pagb%lam), 'pagb lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(pagb%spec), 'pagb spec deallocated', total_tests, total_failures)

        allocate(wr%logt_wn(2), wr%logt_wc(2), wr%z(2), wr%lam(2), wr%spec_wn(2,2,2), wr%spec_wc(2,2,2))
        call wr%clear()
        call assert_true(.not. allocated(wr%logt_wn), 'wr logt_wn deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wr%logt_wc), 'wr logt_wc deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wr%z), 'wr z deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wr%lam), 'wr lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wr%spec_wn), 'wr spec_wn deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(wr%spec_wc), 'wr spec_wc deallocated', total_tests, total_failures)

        allocate(agb%lam_o(2), agb%lam_c(2), agb%lam_car(2), agb%z_o(2))
        allocate(agb%logt_o(2,2), agb%logt_c(2), agb%logt_car(2))
        allocate(agb%spec_o(2,2), agb%spec_c(2,2), agb%spec_car(2,2))
        call agb%clear()
        call assert_true(.not. allocated(agb%lam_o), 'agb lam_o deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%lam_c), 'agb lam_c deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%lam_car), 'agb lam_car deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%z_o), 'agb z_o deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%logt_o), 'agb logt_o deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%logt_c), 'agb logt_c deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%logt_car), 'agb logt_car deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%spec_o), 'agb spec_o deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%spec_c), 'agb spec_c deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agb%spec_car), 'agb spec_car deallocated', total_tests, total_failures)
    end subroutine test_auxiliary_grid_clears

    subroutine test_dust_and_xrb_clears()
        type(dust_emission_t) :: dust
        type(agn_dust_t) :: agn
        type(dust_attenuation_t) :: att
        type(xrb_spectra_t) :: xrb

        call print_group("dust and xrb clears")

        allocate(dust%qpah(2), dust%umin(2), dust%lam(2), dust%spec(2,2,2))
        call dust%clear()
        call assert_true(.not. allocated(dust%qpah), 'dust qpah deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(dust%umin), 'dust umin deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(dust%lam), 'dust lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(dust%spec), 'dust spec deallocated', total_tests, total_failures)

        allocate(agn%tau(2), agn%lam(2), agn%spec(2,2))
        call agn%clear()
        call assert_true(.not. allocated(agn%tau), 'agn tau deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agn%lam), 'agn lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(agn%spec), 'agn spec deallocated', total_tests, total_failures)

        allocate(att%wg_lam(2), att%wg_spec(2,2,2,2), att%smc_lam(2), att%smc_ext(2))
        call att%clear()
        call assert_true(.not. allocated(att%wg_lam), 'att wg_lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(att%wg_spec), 'att wg_spec deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(att%smc_lam), 'att smc_lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(att%smc_ext), 'att smc_ext deallocated', total_tests, total_failures)

        allocate(xrb%lam(2), xrb%age(2), xrb%z(2), xrb%spec(2,2,2))
        call xrb%clear()
        call assert_true(.not. allocated(xrb%lam), 'xrb lam deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(xrb%age), 'xrb age deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(xrb%z), 'xrb z deallocated', total_tests, total_failures)
        call assert_true(.not. allocated(xrb%spec), 'xrb spec deallocated', total_tests, total_failures)
    end subroutine test_dust_and_xrb_clears

    subroutine test_schema_copy_assignment_paths()
        type(axis_desc_t) :: axis_a, axis_b
        type(dataset_desc_t) :: ds_a, ds_b
        type(library_manifest_t) :: mf_a, mf_b
        type(spectral_grid_t) :: sg_a, sg_b
        type(spectral_slice_t) :: ss_a, ss_b
        type(isochrone_grid_t) :: ig_a, ig_b
        type(nebular_grid_t) :: ng_a, ng_b
        type(aux_wmbasic_t) :: wmb_a, wmb_b
        type(aux_pagb_t) :: pagb_a, pagb_b
        type(aux_wr_t) :: wr_a, wr_b
        type(aux_agb_t) :: agb_a, agb_b
        type(dust_emission_t) :: de_a, de_b
        type(agn_dust_t) :: agn_a, agn_b
        type(dust_attenuation_t) :: att_a, att_b
        type(xrb_spectra_t) :: xrb_a, xrb_b

        call print_group("schema copy/assignment paths")

        axis_a%name = 'z'
        axis_a%n = 2
        allocate(axis_a%values(2)); axis_a%values = [0.01_wp, 0.02_wp]
        axis_b = axis_a
        call assert_int_equals(2, axis_b%n, 'axis assignment preserves size', total_tests, total_failures)

        ds_a%dims_csv = 'lambda,z,logt,logg'
        allocate(ds_a%shape(4)); ds_a%shape = [3,2,2,2]
        ds_b = ds_a
        call assert_true(ds_b%has_axis('z'), 'dataset assignment preserves dims_csv', total_tests, total_failures)

        allocate(mf_a%datasets(1)); mf_a%datasets(1)%role = 'spectral_base'
        allocate(mf_a%axes(1)); mf_a%axes(1)%name = 'lambda'
        mf_b = mf_a
        call assert_int_equals(1, size(mf_b%datasets), 'manifest assignment copies datasets', total_tests, total_failures)

        allocate(sg_a%flux(2,2,1,2,2)); sg_a%flux = 1.0_wp
        sg_b = sg_a
        call assert_int_equals(2, size(sg_b%flux, 1), 'spectral grid assignment copies flux', total_tests, total_failures)

        allocate(ss_a%flux(2,2,2)); ss_a%flux = 2.0_wp
        ss_b = ss_a
        call assert_int_equals(2, size(ss_b%flux, 1), 'spectral slice assignment copies flux', total_tests, total_failures)

        allocate(ig_a%nmass(2,2)); ig_a%nmass = 1
        ig_b = ig_a
        call assert_int_equals(2, size(ig_b%nmass, 1), 'isochrone assignment copies arrays', total_tests, total_failures)

        allocate(ng_a%line_pos(2)); ng_a%line_pos = [1.0_wp, 2.0_wp]
        ng_b = ng_a
        call assert_int_equals(2, size(ng_b%line_pos), 'nebular assignment copies arrays', total_tests, total_failures)

        allocate(wmb_a%logt(2)); wmb_a%logt = [3.5_wp, 3.6_wp]
        wmb_b = wmb_a
        call assert_int_equals(2, size(wmb_b%logt), 'wmbasic assignment copies arrays', total_tests, total_failures)

        allocate(pagb_a%lam(2)); pagb_a%lam = [100.0_wp, 200.0_wp]
        pagb_b = pagb_a
        call assert_int_equals(2, size(pagb_b%lam), 'pagb assignment copies arrays', total_tests, total_failures)

        allocate(wr_a%lam(2)); wr_a%lam = [100.0_wp, 200.0_wp]
        wr_b = wr_a
        call assert_int_equals(2, size(wr_b%lam), 'wr assignment copies arrays', total_tests, total_failures)

        allocate(agb_a%lam_o(2)); agb_a%lam_o = [1.0_wp, 2.0_wp]
        agb_b = agb_a
        call assert_int_equals(2, size(agb_b%lam_o), 'agb assignment copies arrays', total_tests, total_failures)

        allocate(de_a%lam(2)); de_a%lam = [1.0_wp, 2.0_wp]
        de_b = de_a
        call assert_int_equals(2, size(de_b%lam), 'dust emission assignment copies arrays', total_tests, total_failures)

        allocate(agn_a%tau(2)); agn_a%tau = [1.0_wp, 2.0_wp]
        agn_b = agn_a
        call assert_int_equals(2, size(agn_b%tau), 'agn dust assignment copies arrays', total_tests, total_failures)

        allocate(att_a%smc_ext(2)); att_a%smc_ext = [0.1_wp, 0.2_wp]
        att_b = att_a
        call assert_int_equals(2, size(att_b%smc_ext), 'dust attenuation assignment copies arrays', total_tests, total_failures)

        allocate(xrb_a%lam(2)); xrb_a%lam = [1.0_wp, 2.0_wp]
        xrb_b = xrb_a
        call assert_int_equals(2, size(xrb_b%lam), 'xrb assignment copies arrays', total_tests, total_failures)
    end subroutine test_schema_copy_assignment_paths

end module test_fsps_data_schema_mod
