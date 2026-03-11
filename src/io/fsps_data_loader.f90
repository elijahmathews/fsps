module fsps_data_loader
    use fsps_context_types, only: fsps_context_t
    use fsps_data_backend, only: backend_status_t, data_backend_t, backend_status_ok
    use fsps_data_registry, only: create_data_backend
    use fsps_data_mapper, only: fsps_data_mapper_t
    use fsps_data_schema, only: spectral_grid_t, isochrone_grid_t, nebular_grid_t, &
                                aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, &
                                dust_emission_t, agn_dust_t, dust_attenuation_t, &
                                xrb_spectra_t, library_manifest_t, dataset_desc_t, axis_desc_t

    implicit none
    private

    public :: resolve_data_backend_uri
    public :: fsps_data_open
    public :: fsps_data_close
    public :: fsps_data_load_spectral_library
    public :: fsps_data_load_isochrones
    public :: fsps_data_load_nebular
    public :: fsps_data_query_axis
    public :: fsps_data_load_wmbasic
    public :: fsps_data_load_pagb
    public :: fsps_data_load_wr
    public :: fsps_data_load_agb
    public :: fsps_data_load_dust_emission
    public :: fsps_data_load_agn_dust
    public :: fsps_data_load_dust_attenuation
    public :: fsps_data_load_xrb

    class(data_backend_t), allocatable, save :: fsps_backend
    type(fsps_data_mapper_t), save :: fsps_mapper
    type(library_manifest_t), save :: fsps_manifest
    logical, save :: fsps_backend_open = .false.

contains

    subroutine resolve_data_backend_uri(ctx, uri, dust_name)
        type(fsps_context_t), intent(in) :: ctx
        character(len=:), allocatable, intent(out) :: uri
        character(len=*), intent(in), optional :: dust_name

        character(len=1024) :: hdf5_file_path
        character(len=1024) :: hdf5_file_path_env
        character(len=64) :: dust_part
        integer :: env_stat

        hdf5_file_path = trim(ctx%sps_home)//'/data/fsps_data_v1.h5'
        hdf5_file_path_env = ''
        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=hdf5_file_path_env, status=env_stat)
        if (env_stat == 0 .and. len_trim(hdf5_file_path_env) > 0) then
            hdf5_file_path = trim(hdf5_file_path_env)
        end if

        if (present(dust_name)) then
            dust_part = trim(dust_name)
        else
            dust_part = trim(ctx%state%str_dustem)
        end if

        uri = trim(hdf5_file_path)//'|'//trim(ctx%state%isoc_type)//'|'// &
              trim(ctx%state%spec_type)//'|'//trim(dust_part)
    end subroutine resolve_data_backend_uri

    subroutine fsps_data_open(uri, status)
        character(len=*), intent(in) :: uri
        type(backend_status_t), intent(out) :: status
        type(backend_status_t) :: close_status

        call status%set_ok()

        if (fsps_backend_open) then
            call fsps_data_close(close_status)
            if (close_status%code /= 0) then
                call status%set_error(close_status%code, trim(close_status%message))
                return
            end if
        end if

        call create_data_backend(fsps_backend, status)
        if (.not. backend_status_ok(status)) return

        call fsps_backend%open(trim(uri), status)
        if (.not. backend_status_ok(status)) then
            if (allocated(fsps_backend)) deallocate(fsps_backend)
            return
        end if

        call fsps_manifest%clear()
        call fsps_backend%read_manifest(fsps_manifest, status)
        if (.not. backend_status_ok(status)) then
            call fsps_backend%close(close_status)
            if (allocated(fsps_backend)) deallocate(fsps_backend)
            call fsps_manifest%clear()
            return
        end if

        fsps_backend_open = .true.
    end subroutine fsps_data_open

    subroutine fsps_data_close(status)
        type(backend_status_t), intent(out) :: status
        type(backend_status_t) :: close_status

        call status%set_ok()

        if (allocated(fsps_backend)) then
            call fsps_backend%close(close_status)
            if (.not. backend_status_ok(close_status)) then
                call status%set_error(close_status%code, trim(close_status%message))
            end if
            deallocate(fsps_backend)
        end if

        call fsps_manifest%clear()
        fsps_backend_open = .false.
    end subroutine fsps_data_close

    subroutine fsps_data_load_spectral_library(role, grid, status)
        character(len=*), intent(in) :: role
        type(spectral_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status
        type(dataset_desc_t) :: dataset

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call find_dataset_by_role(role, dataset, status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_spectral_grid(fsps_backend, dataset, grid, status)
    end subroutine fsps_data_load_spectral_library

    subroutine fsps_data_load_isochrones(grid, status)
        type(isochrone_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_isochrone_grid(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_isochrones

    subroutine fsps_data_load_nebular(component, grid, status)
        character(len=*), intent(in) :: component
        type(nebular_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_nebular_grid(fsps_backend, fsps_manifest, component, grid, status)
    end subroutine fsps_data_load_nebular

    subroutine fsps_data_query_axis(axis_name, axis_desc, status)
        character(len=*), intent(in) :: axis_name
        type(axis_desc_t), intent(inout) :: axis_desc
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call axis_desc%clear()
        call fsps_backend%query_axis(axis_name, axis_desc, status)
    end subroutine fsps_data_query_axis

    subroutine fsps_data_load_wmbasic(grid, status)
        type(aux_wmbasic_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_wmbasic(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_wmbasic

    subroutine fsps_data_load_pagb(grid, status)
        type(aux_pagb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_pagb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_pagb

    subroutine fsps_data_load_wr(grid, status)
        type(aux_wr_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_wr(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_wr

    subroutine fsps_data_load_agb(grid, status)
        type(aux_agb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_agb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_agb

    subroutine fsps_data_load_dust_emission(grid, status)
        type(dust_emission_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_dust_emission(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_dust_emission

    subroutine fsps_data_load_agn_dust(grid, status)
        type(agn_dust_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_agn_dust(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_agn_dust

    subroutine fsps_data_load_dust_attenuation(grid, status)
        type(dust_attenuation_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_dust_attenuation(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_dust_attenuation

    subroutine fsps_data_load_xrb(grid, status)
        type(xrb_spectra_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_xrb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_xrb

    subroutine ensure_backend_ready(status)
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        if (.not. fsps_backend_open .or. .not. allocated(fsps_backend)) then
            call status%set_error(9101, 'FSPS data backend is not open.')
            return
        end if
    end subroutine ensure_backend_ready

    subroutine find_dataset_by_role(role, dataset, status)
        character(len=*), intent(in) :: role
        type(dataset_desc_t), intent(out) :: dataset
        type(backend_status_t), intent(out) :: status
        integer :: i

        call status%set_ok()
        call dataset%clear()

        if (.not. allocated(fsps_manifest%datasets)) then
            call status%set_error(9102, 'FSPS manifest has no datasets.')
            return
        end if

        do i = 1, size(fsps_manifest%datasets)
            if (.not. allocated(fsps_manifest%datasets(i)%role)) cycle
            if (trim(fsps_manifest%datasets(i)%role) /= trim(role)) cycle
            dataset = fsps_manifest%datasets(i)
            return
        end do

        call status%set_error(9103, 'Dataset role not found in manifest: '//trim(role))
    end subroutine find_dataset_by_role

end module fsps_data_loader
