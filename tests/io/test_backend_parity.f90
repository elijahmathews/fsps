program test_backend_parity
    use fsps_precision, only: WP
    use fsps_data_backend, only: data_backend_t, backend_status_t, backend_status_ok
    use fsps_data_registry, only: create_data_backend
    use fsps_data_schema, only: dataset_desc_t, library_manifest_t, spectral_grid_t, isochrone_grid_t, nebular_grid_t, &
                                aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, dust_emission_t, agn_dust_t, &
                                dust_attenuation_t, xrb_spectra_t
    use fsps_data_mapper, only: fsps_data_mapper_t
    implicit none

    class(data_backend_t), allocatable :: backend_legacy, backend_hdf5
    type(backend_status_t) :: status
    type(library_manifest_t) :: manifest_legacy, manifest_hdf5
    type(fsps_data_mapper_t) :: mapper
    type(dataset_desc_t) :: spectral_dataset_legacy, spectral_dataset_hdf5
    type(spectral_grid_t) :: spec_legacy, spec_hdf5
    type(isochrone_grid_t) :: iso_legacy, iso_hdf5
    type(nebular_grid_t) :: neb_wd_legacy, neb_wd_hdf5, neb_nd_legacy, neb_nd_hdf5
    type(aux_wmbasic_t) :: wmb_legacy, wmb_hdf5
    type(aux_pagb_t) :: pagb_legacy, pagb_hdf5
    type(aux_wr_t) :: wr_legacy, wr_hdf5
    type(aux_agb_t) :: agb_legacy, agb_hdf5
    type(dust_emission_t) :: dem_legacy, dem_hdf5
    type(agn_dust_t) :: agn_legacy, agn_hdf5
    type(dust_attenuation_t) :: datt_legacy, datt_hdf5
    type(xrb_spectra_t) :: xrb_legacy, xrb_hdf5

    character(len=1024) :: sps_home, hdf5_path
    character(len=64) :: isoc_type, spec_type, dust_type
    character(len=:), allocatable :: legacy_uri, hdf5_uri
    integer :: argc
    integer :: n_checks, n_failures

    argc = command_argument_count()
    if (argc < 4) then
        write(*, '(A)') 'Usage: test_backend_parity <sps_home> <hdf5_path> <isoc_type> <spec_type> [dust_type]'
        error stop 2
    end if

    call get_command_argument(1, sps_home)
    call get_command_argument(2, hdf5_path)
    call get_command_argument(3, isoc_type)
    call get_command_argument(4, spec_type)
    if (argc >= 5) then
        call get_command_argument(5, dust_type)
    else
        dust_type = 'DL07'
    end if

    legacy_uri = trim(sps_home)//'|'//trim(isoc_type)//'|'//trim(spec_type)//'|'//trim(dust_type)
    hdf5_uri = trim(hdf5_path)//'|'//trim(isoc_type)//'|'//trim(spec_type)//'|'//trim(dust_type)

    n_checks = 0
    n_failures = 0

    call create_data_backend('legacy', backend_legacy, status)
    call require_ok('create legacy backend', status)

    call create_data_backend('hdf5', backend_hdf5, status)
    call require_ok('create hdf5 backend', status)

    call backend_legacy%open(legacy_uri, status)
    call require_ok('open legacy backend', status)

    call backend_hdf5%open(hdf5_uri, status)
    call require_ok('open hdf5 backend', status)

    call backend_legacy%read_manifest(manifest_legacy, status)
    call require_ok('read legacy manifest', status)

    call backend_hdf5%read_manifest(manifest_hdf5, status)
    call require_ok('read hdf5 manifest', status)

    call find_dataset_by_role(manifest_legacy, 'spectral_base', spectral_dataset_legacy, status)
    call require_ok('find legacy spectral_base', status)

    call find_dataset_by_role(manifest_hdf5, 'spectral_base', spectral_dataset_hdf5, status)
    call require_ok('find hdf5 spectral_base', status)

    call mapper%map_spectral_grid(backend_legacy, spectral_dataset_legacy, spec_legacy, status)
    call require_ok('map legacy spectral grid', status)
    call mapper%map_spectral_grid(backend_hdf5, spectral_dataset_hdf5, spec_hdf5, status)
    call require_ok('map hdf5 spectral grid', status)
    call compare_real_1d('spectral.axis_lambda', spec_legacy%axis_lambda, spec_hdf5%axis_lambda, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_1d('spectral.axis_z', spec_legacy%axis_z, spec_hdf5%axis_z, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_1d('spectral.axis_afe', spec_legacy%axis_afe, spec_hdf5%axis_afe, 1.0e-8_wp, n_checks, n_failures)
    call compare_axis_size('spectral.axis_logt.size', spec_legacy%axis_logt, spec_hdf5%axis_logt, n_checks, n_failures)
    call compare_axis_size('spectral.axis_logg.size', spec_legacy%axis_logg, spec_hdf5%axis_logg, n_checks, n_failures)
    call compare_real_5d('spectral.flux', spec_legacy%flux, spec_hdf5%flux, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_isochrone_grid(backend_legacy, manifest_legacy, iso_legacy, status)
    call require_ok('map legacy isochrone', status)
    call mapper%map_isochrone_grid(backend_hdf5, manifest_hdf5, iso_hdf5, status)
    call require_ok('map hdf5 isochrone', status)
    call compare_int_2d('iso.nmass', iso_legacy%nmass, iso_hdf5%nmass, n_checks, n_failures)
    call compare_real_2d('iso.timestep', iso_legacy%timestep_logyr, iso_hdf5%timestep_logyr, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_3d('iso.mini', iso_legacy%mini, iso_hdf5%mini, 2.0e-5_wp, n_checks, n_failures)
    call compare_real_3d('iso.mact', iso_legacy%mact, iso_hdf5%mact, 2.0e-5_wp, n_checks, n_failures)
    call compare_real_3d('iso.logl', iso_legacy%logl, iso_hdf5%logl, 1.0e-6_wp, n_checks, n_failures)
    call compare_real_3d('iso.logt', iso_legacy%logt, iso_hdf5%logt, 1.0e-6_wp, n_checks, n_failures)
    call compare_real_3d('iso.logg', iso_legacy%logg, iso_hdf5%logg, 1.0e-6_wp, n_checks, n_failures)
    call compare_real_3d('iso.phase', iso_legacy%phase, iso_hdf5%phase, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_3d('iso.ffco', iso_legacy%ffco, iso_hdf5%ffco, 1.0e-6_wp, n_checks, n_failures)
    call compare_real_3d('iso.lmdot', iso_legacy%lmdot, iso_hdf5%lmdot, 1.0e-6_wp, n_checks, n_failures)

    call mapper%map_nebular_grid(backend_legacy, manifest_legacy, 'WD', neb_wd_legacy, status)
    call require_ok('map legacy nebular WD', status)
    call mapper%map_nebular_grid(backend_hdf5, manifest_hdf5, 'WD', neb_wd_hdf5, status)
    call require_ok('map hdf5 nebular WD', status)
    call compare_real_1d('nebular_wd.line_pos', neb_wd_legacy%line_pos, neb_wd_hdf5%line_pos, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_4d('nebular_wd.cont', neb_wd_legacy%cont, neb_wd_hdf5%cont, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_4d('nebular_wd.line', neb_wd_legacy%line, neb_wd_hdf5%line, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_nebular_grid(backend_legacy, manifest_legacy, 'ND', neb_nd_legacy, status)
    call require_ok('map legacy nebular ND', status)
    call mapper%map_nebular_grid(backend_hdf5, manifest_hdf5, 'ND', neb_nd_hdf5, status)
    call require_ok('map hdf5 nebular ND', status)
    call compare_real_4d('nebular_nd.cont', neb_nd_legacy%cont, neb_nd_hdf5%cont, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_4d('nebular_nd.line', neb_nd_legacy%line, neb_nd_hdf5%line, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_wmbasic(backend_legacy, manifest_legacy, wmb_legacy, status)
    call require_ok('map legacy wmbasic', status)
    call mapper%map_wmbasic(backend_hdf5, manifest_hdf5, wmb_hdf5, status)
    call require_ok('map hdf5 wmbasic', status)
    call compare_real_1d('wmb.logt', wmb_legacy%logt, wmb_hdf5%logt, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_1d('wmb.z', wmb_legacy%z, wmb_hdf5%z, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_1d('wmb.lam', wmb_legacy%lam, wmb_hdf5%lam, 1.0e-8_wp, n_checks, n_failures)
    call compare_real_4d('wmb.spec', wmb_legacy%spec, wmb_hdf5%spec, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_pagb(backend_legacy, manifest_legacy, pagb_legacy, status)
    call require_ok('map legacy pagb', status)
    call mapper%map_pagb(backend_hdf5, manifest_hdf5, pagb_hdf5, status)
    call require_ok('map hdf5 pagb', status)
    call compare_real_3d('pagb.spec', pagb_legacy%spec, pagb_hdf5%spec, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_wr(backend_legacy, manifest_legacy, wr_legacy, status)
    call require_ok('map legacy wr', status)
    call mapper%map_wr(backend_hdf5, manifest_hdf5, wr_hdf5, status)
    call require_ok('map hdf5 wr', status)
    call compare_real_3d('wr.spec_wn', wr_legacy%spec_wn, wr_hdf5%spec_wn, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_3d('wr.spec_wc', wr_legacy%spec_wc, wr_hdf5%spec_wc, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_agb(backend_legacy, manifest_legacy, agb_legacy, status)
    call require_ok('map legacy agb', status)
    call mapper%map_agb(backend_hdf5, manifest_hdf5, agb_hdf5, status)
    call require_ok('map hdf5 agb', status)
    call compare_real_2d('agb.spec_o', agb_legacy%spec_o, agb_hdf5%spec_o, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_2d('agb.spec_c', agb_legacy%spec_c, agb_hdf5%spec_c, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_2d('agb.spec_car', agb_legacy%spec_car, agb_hdf5%spec_car, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_dust_emission(backend_legacy, manifest_legacy, dem_legacy, status)
    call require_ok('map legacy dust emission', status)
    call mapper%map_dust_emission(backend_hdf5, manifest_hdf5, dem_hdf5, status)
    call require_ok('map hdf5 dust emission', status)
    call compare_real_3d('dust_em.spec', dem_legacy%spec, dem_hdf5%spec, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_agn_dust(backend_legacy, manifest_legacy, agn_legacy, status)
    call require_ok('map legacy agn dust', status)
    call mapper%map_agn_dust(backend_hdf5, manifest_hdf5, agn_hdf5, status)
    call require_ok('map hdf5 agn dust', status)
    call compare_real_2d('agn_dust.spec', agn_legacy%spec, agn_hdf5%spec, 1.0e-5_wp, n_checks, n_failures)

    call mapper%map_dust_attenuation(backend_legacy, manifest_legacy, datt_legacy, status)
    call require_ok('map legacy dust attenuation', status)
    call mapper%map_dust_attenuation(backend_hdf5, manifest_hdf5, datt_hdf5, status)
    call require_ok('map hdf5 dust attenuation', status)
    call compare_real_4d('dust_att.wg_spec', datt_legacy%wg_spec, datt_hdf5%wg_spec, 1.0e-5_wp, n_checks, n_failures)
    call compare_real_1d('dust_att.smc_ext', datt_legacy%smc_ext, datt_hdf5%smc_ext, 1.0e-8_wp, n_checks, n_failures)

    call mapper%map_xrb(backend_legacy, manifest_legacy, xrb_legacy, status)
    call require_ok('map legacy xrb', status)
    call mapper%map_xrb(backend_hdf5, manifest_hdf5, xrb_hdf5, status)
    call require_ok('map hdf5 xrb', status)
    call compare_real_3d('xrb.spec', xrb_legacy%spec, xrb_hdf5%spec, 1.0e-5_wp, n_checks, n_failures)

    call backend_legacy%close(status)
    call require_ok('close legacy backend', status)
    call backend_hdf5%close(status)
    call require_ok('close hdf5 backend', status)

    write(*, '(A,1x,I0,1x,A,1x,I0)') 'PARITY SUMMARY:', n_checks - n_failures, 'passed /', n_checks
    if (n_failures > 0) then
        error stop 1
    end if

contains

    subroutine require_ok(label, st)
        character(len=*), intent(in) :: label
        type(backend_status_t), intent(in) :: st

        if (.not. backend_status_ok(st)) then
            if (allocated(st%message)) then
                write(*, '(A,1x,A,1x,I0,2A)') 'ERROR:', trim(label), st%code, ' ', trim(st%message)
            else
                write(*, '(A,1x,A,1x,I0)') 'ERROR:', trim(label), st%code
            end if
            error stop 1
        end if
    end subroutine require_ok

    subroutine find_dataset_by_role(manifest, role, dataset, st)
        type(library_manifest_t), intent(in) :: manifest
        character(len=*), intent(in) :: role
        type(dataset_desc_t), intent(inout) :: dataset
        type(backend_status_t), intent(out) :: st
        integer :: i

        call st%set_ok()
        call dataset%clear()

        if (.not. allocated(manifest%datasets)) then
            call st%set_error(9801, 'Manifest has no datasets.')
            return
        end if

        do i = 1, size(manifest%datasets)
            if (.not. allocated(manifest%datasets(i)%role)) cycle
            if (trim(manifest%datasets(i)%role) /= trim(role)) cycle
            dataset = manifest%datasets(i)
            return
        end do

        call st%set_error(9802, 'Role not found: '//trim(role))
    end subroutine find_dataset_by_role

    subroutine compare_int_2d(label, a, b, checks, fails)
        character(len=*), intent(in) :: label
        integer, intent(in) :: a(:,:), b(:,:)
        integer, intent(inout) :: checks, fails
        integer :: max_abs

        checks = checks + 1
        if (size(a,1) /= size(b,1) .or. size(a,2) /= size(b,2)) then
            fails = fails + 1
            write(*, '(A,1x,A,1x,A,2(I0,1x),A,2(I0,1x))') 'FAIL', trim(label), 'shape', &
                size(a,1), size(a,2), 'vs', size(b,1), size(b,2)
            return
        end if

        max_abs = maxval(abs(a - b))
        if (max_abs /= 0) then
            fails = fails + 1
            write(*, '(A,1x,A,1x,A,1x,I0)') 'FAIL', trim(label), 'max_abs=', max_abs
        else
            write(*, '(A,1x,A)') 'PASS', trim(label)
        end if
    end subroutine compare_int_2d

    subroutine compare_real_1d(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:), b(:)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        call compare_real_flat(label, a, b, tol, checks, fails)
    end subroutine compare_real_1d

    subroutine compare_axis_size(label, a, b, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:), b(:)
        integer, intent(inout) :: checks, fails

        checks = checks + 1
        if (size(a) /= size(b)) then
            fails = fails + 1
            write(*, '(A,1x,A,1x,A,1x,I0,1x,A,1x,I0)') 'FAIL', trim(label), 'size', size(a), 'vs', size(b)
        else
            write(*, '(A,1x,A)') 'PASS', trim(label)
        end if
    end subroutine compare_axis_size

    subroutine compare_real_2d(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:,:), b(:,:)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        call compare_real_flat(label, a, b, tol, checks, fails)
    end subroutine compare_real_2d

    subroutine compare_real_3d(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:,:,:), b(:,:,:)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        call compare_real_flat(label, a, b, tol, checks, fails)
    end subroutine compare_real_3d

    subroutine compare_real_4d(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        call compare_real_flat(label, a, b, tol, checks, fails)
    end subroutine compare_real_4d

    subroutine compare_real_5d(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(:,:,:,:,:), b(:,:,:,:,:)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        call compare_real_flat(label, a, b, tol, checks, fails)
    end subroutine compare_real_5d

    subroutine compare_real_flat(label, a, b, tol, checks, fails)
        character(len=*), intent(in) :: label
        real(WP), intent(in) :: a(..), b(..)
        real(WP), intent(in) :: tol
        integer, intent(inout) :: checks, fails
        real(WP), allocatable :: aa(:), bb(:), denom(:)
        real(WP) :: max_abs, max_rel

        checks = checks + 1

        if (size(a) /= size(b)) then
            fails = fails + 1
            write(*, '(A,1x,A,1x,A,1x,I0,1x,A,1x,I0)') 'FAIL', trim(label), 'size', size(a), 'vs', size(b)
            return
        end if

        select rank (a)
        rank (1)
            select rank (b)
            rank (1)
                aa = reshape(a, [size(a)])
                bb = reshape(b, [size(b)])
            end select
        rank (2)
            select rank (b)
            rank (2)
                aa = reshape(a, [size(a)])
                bb = reshape(b, [size(b)])
            end select
        rank (3)
            select rank (b)
            rank (3)
                aa = reshape(a, [size(a)])
                bb = reshape(b, [size(b)])
            end select
        rank (4)
            select rank (b)
            rank (4)
                aa = reshape(a, [size(a)])
                bb = reshape(b, [size(b)])
            end select
        rank (5)
            select rank (b)
            rank (5)
                aa = reshape(a, [size(a)])
                bb = reshape(b, [size(b)])
            end select
        rank default
            fails = fails + 1
            write(*, '(A,1x,A,1x,A)') 'FAIL', trim(label), 'unsupported rank'
            return
        end select

        max_abs = maxval(abs(aa - bb))
        allocate(denom(size(aa)))
        denom = max(abs(aa), 1.0e-30_wp)
        max_rel = maxval(abs(aa - bb)/denom)

        if (max_abs > tol) then
            fails = fails + 1
            write(*, '(A,1x,A,1x,A,1x,ES12.5,1x,A,1x,ES12.5)') 'FAIL', trim(label), 'max_abs=', max_abs, 'max_rel=', max_rel
        else
            write(*, '(A,1x,A,1x,A,1x,ES12.5,1x,A,1x,ES12.5)') 'PASS', trim(label), 'max_abs=', max_abs, 'max_rel=', max_rel
        end if
    end subroutine compare_real_flat

end program test_backend_parity
