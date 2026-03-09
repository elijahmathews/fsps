module fsps_data_mapper
    !> @brief Mapping layer between backend I/O and FSPS dense spectral containers.
    !>
    !> @details
    !> This module performs schema-aware translation from backend reads into
    !> canonical Fortran types. It owns non-I/O logic such as:
    !> - Parsing `dims_csv`
    !> - Handling degenerate alpha-axis insertion (`afe=1`) for 4D inputs
    !> - Bracketing index computations on numeric axes

    use fsps_precision, only: WP
    use fsps_strings, only: to_lower
    use fsps_data_schema, only: axis_desc_t, dataset_desc_t, library_manifest_t, &
                                spectral_grid_t, spectral_slice_t, isochrone_grid_t, nebular_grid_t, &
                                aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, &
                                dust_emission_t, agn_dust_t, dust_attenuation_t, xrb_spectra_t
    use fsps_data_backend, only: data_backend_t, backend_status_t, backend_status_ok

    implicit none
    private

    public :: fsps_data_mapper_t

    !> @brief Stateless mapper object.
    type :: fsps_data_mapper_t
    contains
        !> @brief Build a canonical dense 5D spectral grid from backend slices.
        procedure, public :: map_spectral_grid
        !> @brief Read one spectral slice at requested physical coordinates.
        procedure, public :: map_spectral_slice
        !> @brief Compute lower/upper bracketing indices and interpolation weight.
        procedure, public :: find_bracketing_indices
        !> @brief Build dense isochrone tables from manifest roles.
        procedure, public :: map_isochrone_grid
        !> @brief Build dense nebular tables from manifest roles.
        procedure, public :: map_nebular_grid
        !> @brief Build dense WMBasic auxiliary tables from manifest roles.
        procedure, public :: map_wmbasic
        !> @brief Build dense Post-AGB auxiliary tables from manifest roles.
        procedure, public :: map_pagb
        !> @brief Build dense Wolf-Rayet auxiliary tables from manifest roles.
        procedure, public :: map_wr
        !> @brief Build dense AGB auxiliary tables from manifest roles.
        procedure, public :: map_agb
        !> @brief Build dense dust emission tables from manifest roles.
        procedure, public :: map_dust_emission
        !> @brief Build dense AGN dust tables from manifest roles.
        procedure, public :: map_agn_dust
        !> @brief Build dense dust attenuation tables from manifest roles.
        procedure, public :: map_dust_attenuation
        !> @brief Build dense XRB spectral tables from manifest roles.
        procedure, public :: map_xrb
    end type fsps_data_mapper_t

contains

    !> @brief Populate a canonical dense spectral grid.
    !>
    !> @details
    !> If `dataset%dims_csv` omits `afe` (4D input), this routine allocates the
    !> output as `(lambda, z, 1, logt, logg)` and writes all data into `afe=1`.
    !>
    !> @param[inout] self    Mapper instance.
    !> @param[inout] backend Active backend implementation.
    !> @param[in]    dataset Spectral dataset descriptor.
    !> @param[inout] grid    Output canonical dense spectral grid.
    !> @param[out]   status  Operation status.
    subroutine map_spectral_grid(self, backend, dataset, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(dataset_desc_t), intent(in) :: dataset
        type(spectral_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        type(axis_desc_t) :: ax_lambda, ax_z, ax_afe, ax_logt, ax_logg
        type(spectral_slice_t) :: slice
        logical :: has_afe
        integer :: n_lambda, nz, nafe, n_logt, n_logg
        integer :: iz, iafe, iafe_backend

        call status%set_ok()
        call grid%clear()

        if (.not. allocated(dataset%dims_csv)) then
            call status%set_error(1000, 'Dataset descriptor has unallocated dims_csv in map_spectral_grid.')
            return
        end if

        call require_axis(backend, 'lambda', ax_lambda, status)
        if (.not. backend_status_ok(status)) return
        call require_axis(backend, 'z', ax_z, status)
        if (.not. backend_status_ok(status)) return
        call require_axis(backend, 'logt', ax_logt, status)
        if (.not. backend_status_ok(status)) return
        call require_axis(backend, 'logg', ax_logg, status)
        if (.not. backend_status_ok(status)) return

        has_afe = dataset%has_axis('afe') .or. dataset%has_axis('alpha_fe')

        if (has_afe) then
            call require_axis(backend, 'afe', ax_afe, status)
            if (.not. backend_status_ok(status)) then
                ! Allow alternate axis spelling in backend registry.
                call require_axis(backend, 'alpha_fe', ax_afe, status)
                if (.not. backend_status_ok(status)) return
            end if
            nafe = ax_afe%n
        else
            nafe = 1
        end if

        n_lambda = ax_lambda%n
        nz = ax_z%n
        n_logt = ax_logt%n
        n_logg = ax_logg%n

        if (n_lambda <= 0 .or. nz <= 0 .or. nafe <= 0 .or. n_logt <= 0 .or. n_logg <= 0) then
            call status%set_error(1001, 'Invalid axis sizes while mapping spectral grid.')
            return
        end if

        allocate(grid%axis_lambda(n_lambda))
        allocate(grid%axis_z(nz))
        allocate(grid%axis_afe(nafe))
        allocate(grid%axis_logt(n_logt))
        allocate(grid%axis_logg(n_logg))

        grid%axis_lambda = ax_lambda%values
        grid%axis_z = ax_z%values
        if (has_afe) then
            grid%axis_afe = ax_afe%values
        else
            grid%axis_afe(1) = 0.0_wp
        end if
        grid%axis_logt = ax_logt%values
        grid%axis_logg = ax_logg%values

        if (dataset%has_missing_value) then
            grid%missing_value = dataset%missing_value
        else
            grid%missing_value = -1.0e99_wp
        end if

        allocate(grid%flux(n_lambda, nz, nafe, n_logt, n_logg))
        grid%flux = grid%missing_value

        if (dataset%has_valid_mask) then
            allocate(grid%valid(nz, nafe, n_logt, n_logg))
            grid%valid = .false.
        end if

        do iz = 1, nz
            do iafe = 1, nafe
                call slice%clear()

                if (has_afe) then
                    iafe_backend = iafe
                else
                    iafe_backend = 1
                end if

                call backend%read_spectral_slice(dataset, iz, iafe_backend, slice, status)
                if (.not. backend_status_ok(status)) return

                if (.not. allocated(slice%flux)) then
                    call status%set_error(1002, 'Backend returned unallocated spectral slice.')
                    return
                end if

                if (size(slice%flux, 1) /= n_lambda .or. size(slice%flux, 2) /= n_logt .or. size(slice%flux, 3) /= n_logg) then
                    call status%set_error(1003, 'Spectral slice shape mismatch in mapper.')
                    return
                end if

                grid%flux(:, iz, iafe, :, :) = slice%flux

                if (allocated(grid%valid)) then
                    if (allocated(slice%valid)) then
                        if (size(slice%valid, 1) /= n_logt .or. size(slice%valid, 2) /= n_logg) then
                            call status%set_error(1004, 'Slice validity mask shape mismatch in mapper.')
                            return
                        end if
                        grid%valid(iz, iafe, :, :) = slice%valid
                    else
                        ! Conservative default when dataset declares a mask but backend slice does not provide one.
                        grid%valid(iz, iafe, :, :) = .false.
                    end if
                end if
            end do
        end do
    end subroutine map_spectral_grid

    !> @brief Read one spectral slice at requested physical coordinates.
    !>
    !> @details
    !> The mapper computes bracketing indices on `z` and optional `afe` axes,
    !> then chooses the nearest index for a single-slice read.
    !>
    !> If `afe` is absent in `dataset%dims_csv`, the backend is always called
    !> with `iafe=1`, and output metadata is set accordingly.
    !>
    !> @param[inout] self      Mapper instance.
    !> @param[inout] backend   Active backend implementation.
    !> @param[in]    dataset   Spectral dataset descriptor.
    !> @param[in]    z_value   Requested metallicity coordinate.
    !> @param[in]    afe_value Requested alpha coordinate.
    !> @param[inout] slice     Output spectral slice.
    !> @param[out]   status    Operation status.
    subroutine map_spectral_slice(self, backend, dataset, z_value, afe_value, slice, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(dataset_desc_t), intent(in) :: dataset
        real(WP), intent(in) :: z_value
        real(WP), intent(in) :: afe_value
        type(spectral_slice_t), intent(inout) :: slice
        type(backend_status_t), intent(out) :: status

        type(axis_desc_t) :: ax_z, ax_afe
        logical :: has_afe
        integer :: iz_lo, iz_hi, ia_lo, ia_hi
        real(WP) :: wz_hi, wa_hi
        integer :: iz_pick, ia_pick

        call status%set_ok()
        call slice%clear()

        if (.not. allocated(dataset%dims_csv)) then
            call status%set_error(1200, 'Dataset descriptor has unallocated dims_csv in map_spectral_slice.')
            return
        end if

        call require_axis(backend, 'z', ax_z, status)
        if (.not. backend_status_ok(status)) return

        call self%find_bracketing_indices(ax_z%values, z_value, iz_lo, iz_hi, wz_hi)
        iz_pick = merge(iz_hi, iz_lo, wz_hi >= 0.5_wp)

        has_afe = dataset%has_axis('afe') .or. dataset%has_axis('alpha_fe')
        if (has_afe) then
            call require_axis(backend, 'afe', ax_afe, status)
            if (.not. backend_status_ok(status)) then
                call require_axis(backend, 'alpha_fe', ax_afe, status)
                if (.not. backend_status_ok(status)) return
            end if
            call self%find_bracketing_indices(ax_afe%values, afe_value, ia_lo, ia_hi, wa_hi)
            ia_pick = merge(ia_hi, ia_lo, wa_hi >= 0.5_wp)
        else
            ia_pick = 1
        end if

        call backend%read_spectral_slice(dataset, iz_pick, ia_pick, slice, status)
        if (.not. backend_status_ok(status)) return

        slice%iz = iz_pick
        if (has_afe) then
            slice%iafe = ia_pick
        else
            slice%iafe = 1
        end if

        if (dataset%has_missing_value) then
            slice%missing_value = dataset%missing_value
        else
            slice%missing_value = -1.0e99_wp
        end if
    end subroutine map_spectral_slice

    !> @brief Compute bracketing indices and upper weight on a monotonic axis.
    !>
    !> @param[in]  axis_values Monotonic coordinate array.
    !> @param[in]  x           Query value.
    !> @param[out] i_lo        Lower bracket index (1-based).
    !> @param[out] i_hi        Upper bracket index (1-based).
    !> @param[out] w_hi        Fraction toward upper bracket in `[0,1]`.
    subroutine find_bracketing_indices(self, axis_values, x, i_lo, i_hi, w_hi)
        class(fsps_data_mapper_t), intent(inout) :: self
        real(WP), intent(in) :: axis_values(:)
        real(WP), intent(in) :: x
        integer, intent(out) :: i_lo
        integer, intent(out) :: i_hi
        real(WP), intent(out) :: w_hi

        integer :: n, lo, hi, mid
        logical :: ascending

        n = size(axis_values)
        if (n <= 1) then
            i_lo = 1
            i_hi = 1
            w_hi = 0.0_wp
            return
        end if

        ascending = (axis_values(n) >= axis_values(1))

        if (ascending) then
            if (x <= axis_values(1)) then
                i_lo = 1
                i_hi = 1
                w_hi = 0.0_wp
                return
            end if
            if (x >= axis_values(n)) then
                i_lo = n
                i_hi = n
                w_hi = 0.0_wp
                return
            end if
        else
            if (x >= axis_values(1)) then
                i_lo = 1
                i_hi = 1
                w_hi = 0.0_wp
                return
            end if
            if (x <= axis_values(n)) then
                i_lo = n
                i_hi = n
                w_hi = 0.0_wp
                return
            end if
        end if

        lo = 1
        hi = n

        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (ascending) then
                if (x >= axis_values(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            else
                if (x <= axis_values(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            end if
        end do

        i_lo = lo
        i_hi = hi

        if (abs(axis_values(i_hi) - axis_values(i_lo)) <= tiny(1.0_wp)) then
            w_hi = 0.0_wp
        else
            w_hi = (x - axis_values(i_lo)) / (axis_values(i_hi) - axis_values(i_lo))
            w_hi = max(0.0_wp, min(1.0_wp, w_hi))
        end if

    end subroutine find_bracketing_indices

    !> @brief Require one named axis from backend.
    subroutine require_axis(backend, axis_name, axis_desc, status)
        class(data_backend_t), intent(inout) :: backend
        character(len=*), intent(in) :: axis_name
        type(axis_desc_t), intent(inout) :: axis_desc
        type(backend_status_t), intent(out) :: status

        call axis_desc%clear()
        call backend%query_axis(axis_name, axis_desc, status)
        if (.not. backend_status_ok(status)) return

        if (axis_desc%n <= 0 .or. .not. allocated(axis_desc%values)) then
            call status%set_error(1101, 'Axis query returned empty axis: '//trim(axis_name))
            return
        end if

        if (size(axis_desc%values) /= axis_desc%n) then
            call status%set_error(1102, 'Axis size mismatch for axis: '//trim(axis_name))
            return
        end if
    end subroutine require_axis

    !> @brief Extract dataset path by role from manifest.
    subroutine get_path_by_role(manifest, role, path, status)
        type(library_manifest_t), intent(in) :: manifest
        character(len=*), intent(in) :: role
        character(len=:), allocatable, intent(out) :: path
        type(backend_status_t), intent(out) :: status

        integer :: i

        call status%set_ok()
        if (allocated(path)) deallocate(path)

        if (.not. allocated(manifest%datasets)) then
            call status%set_error(1301, 'Manifest has no datasets in get_path_by_role.')
            return
        end if

        do i = 1, size(manifest%datasets)
            if (.not. allocated(manifest%datasets(i)%role)) cycle
            if (trim(manifest%datasets(i)%role) /= trim(role)) cycle
            if (.not. allocated(manifest%datasets(i)%path)) then
                call status%set_error(1302, 'Dataset path missing for role: '//trim(role))
                return
            end if
            path = trim(manifest%datasets(i)%path)
            return
        end do

        call status%set_error(1303, 'Role not found in manifest: '//trim(role))
    end subroutine get_path_by_role

    !> @brief Populate dense isochrone tables from manifest roles.
    subroutine map_isochrone_grid(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(isochrone_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'isoc_nmass', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_int_2d(path, grid%nmass, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_timestep', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%timestep_logyr, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_mini', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%mini, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_mact', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%mact, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_logl', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%logl, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_logt', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%logt, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_logg', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%logg, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_phase', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%phase, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_ffco', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%ffco, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'isoc_lmdot', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%lmdot, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_isochrone_grid

    !> @brief Populate dense nebular tables from manifest roles.
    subroutine map_nebular_grid(self, backend, manifest, component, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        character(len=*), intent(in) :: component
        type(nebular_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: role_prefix
        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        role_prefix = 'nebular_'//trim(to_lower(trim(component)))//'_'

        call get_path_by_role(manifest, trim(role_prefix)//'line_pos', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%line_pos, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, trim(role_prefix)//'logz', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logz, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, trim(role_prefix)//'age', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%age, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, trim(role_prefix)//'logu', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logu, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, trim(role_prefix)//'cont', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_4d(path, grid%cont, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, trim(role_prefix)//'line', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_4d(path, grid%line, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_nebular_grid

    !> @brief Populate dense WMBasic auxiliary tables from manifest roles.
    subroutine map_wmbasic(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(aux_wmbasic_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'wmb_logt', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wmb_z', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%z, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wmb_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wmb_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_4d(path, grid%spec, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_wmbasic

    !> @brief Populate dense Post-AGB auxiliary tables from manifest roles.
    subroutine map_pagb(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(aux_pagb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'pagb_logt', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'pagb_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'pagb_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%spec, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_pagb

    !> @brief Populate dense Wolf-Rayet auxiliary tables from manifest roles.
    subroutine map_wr(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(aux_wr_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'wr_logt_wn', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt_wn, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wr_logt_wc', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt_wc, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wr_z', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%z, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wr_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wr_spec_wn', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%spec_wn, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'wr_spec_wc', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%spec_wc, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_wr

    !> @brief Populate dense AGB auxiliary tables from manifest roles.
    subroutine map_agb(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(aux_agb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'agb_z_o', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%z_o, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_logt_c', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt_c, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_logt_car', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%logt_car, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_lam_o', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam_o, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_lam_c', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam_c, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_lam_car', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam_car, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_logt_o', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%logt_o, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_spec_o', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%spec_o, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_spec_c', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%spec_c, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agb_spec_car', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%spec_car, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_agb

    !> @brief Populate dense dust emission tables from manifest roles.
    subroutine map_dust_emission(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(dust_emission_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'dust_em_qpah', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%qpah, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_em_umin', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%umin, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_em_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_em_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%spec, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_dust_emission

    !> @brief Populate dense AGN dust tables from manifest roles.
    subroutine map_agn_dust(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(agn_dust_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'agn_dust_tau', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%tau, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agn_dust_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'agn_dust_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_2d(path, grid%spec, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_agn_dust

    !> @brief Populate dense dust attenuation tables from manifest roles.
    subroutine map_dust_attenuation(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(dust_attenuation_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'dust_att_wg_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%wg_lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_att_wg_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_4d(path, grid%wg_spec, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_att_smc_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%smc_lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'dust_att_smc_ext', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%smc_ext, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_dust_attenuation

    !> @brief Populate dense XRB spectral tables from manifest roles.
    subroutine map_xrb(self, backend, manifest, grid, status)
        class(fsps_data_mapper_t), intent(inout) :: self
        class(data_backend_t), intent(inout) :: backend
        type(library_manifest_t), intent(in) :: manifest
        type(xrb_spectra_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()
        call grid%clear()

        call get_path_by_role(manifest, 'xrb_lam', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%lam, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'xrb_age', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%age, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'xrb_z', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_1d(path, grid%z, status)
        if (.not. backend_status_ok(status)) return

        call get_path_by_role(manifest, 'xrb_spec', path, status)
        if (.not. backend_status_ok(status)) return
        call backend%read_real_3d(path, grid%spec, status)
        if (.not. backend_status_ok(status)) return
    end subroutine map_xrb

end module fsps_data_mapper
