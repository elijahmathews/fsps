module fsps_data_backend_hdf5
    !> @brief Concrete HDF5 backend implementing `data_backend_t`.
    !>
    !> @details
    !> This backend reads dense FSDS datasets from an HDF5 file using the
    !> standard Fortran HDF5 interface (`use hdf5`).
    !>
    !> Important behaviors:
    !> - Uses hyperslab selections for spectral slice and neighborhood reads.
    !> - Does not materialize full spectral 5D cubes unless explicitly requested.
    !> - Returns errors through `backend_status_t` and avoids program termination.

    use fsps_precision, only: WP
    use fsps_data_schema, only: axis_desc_t, dataset_desc_t, library_manifest_t, &
                                spectral_grid_t, spectral_slice_t
    use fsps_data_backend, only: data_backend_t, backend_status_t
    use fsps_strings, only: to_lower
    use hdf5

    implicit none
    private

    public :: hdf5_backend_t
    public :: resolve_axis_dataset_path

    !> @brief HDF5-backed implementation of the abstract FSDS backend.
    type, extends(data_backend_t) :: hdf5_backend_t
        private
        integer(HID_T) :: file_id = -1_HID_T
        logical :: is_opened = .false.
        logical :: hdf5_initialized = .false.
        character(len=:), allocatable :: source_uri
        character(len=:), allocatable :: data_uri
        character(len=:), allocatable :: isoc_type
        character(len=:), allocatable :: spec_type
        character(len=:), allocatable :: dust_type
    contains
        procedure, public :: open => hdf5_backend_open
        procedure, public :: close => hdf5_backend_close
        procedure, public :: is_open => hdf5_backend_is_open

        procedure, public :: read_manifest => hdf5_backend_read_manifest
        procedure, public :: has_path => hdf5_backend_has_path
        procedure, public :: query_axis => hdf5_backend_query_axis

        procedure, public :: read_real_1d => hdf5_backend_read_real_1d
        procedure, public :: read_real_2d => hdf5_backend_read_real_2d
        procedure, public :: read_real_3d => hdf5_backend_read_real_3d
        procedure, public :: read_real_4d => hdf5_backend_read_real_4d
        procedure, public :: read_int_1d => hdf5_backend_read_int_1d
        procedure, public :: read_int_2d => hdf5_backend_read_int_2d
        procedure, public :: read_int_3d => hdf5_backend_read_int_3d

        procedure, public :: read_spectral_slice => hdf5_backend_read_spectral_slice
        procedure, public :: read_spectral_neighborhood => hdf5_backend_read_spectral_neighborhood
    end type hdf5_backend_t

contains

    pure function resolve_axis_dataset_path(axis_name, spec_type) result(path)
        character(len=*), intent(in) :: axis_name
        character(len=*), intent(in) :: spec_type
        character(len=:), allocatable :: path
        character(len=:), allocatable :: spec_name

        character(len=:), allocatable :: axis_l, spec_l

        axis_l = trim(to_lower(adjustl(trim(axis_name))))
        spec_l = trim(to_lower(adjustl(trim(spec_type))))

        select case (axis_l)
        case ('lambda', 'z', 'logt', 'logg', 'afe', 'alpha_fe')
            if (len_trim(spec_l) > 0) then
                path = '/libraries/spectra/'//trim(spec_l)//'/axes/'//trim(axis_l)
            else
                path = '/axes/'//trim(axis_l)
            end if
        case default
            path = '/axes/'//trim(axis_l)
        end select
    end function resolve_axis_dataset_path

    !> @brief Open an HDF5 source file in read-only mode.
    subroutine hdf5_backend_open(self, source_uri, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: source_uri
        type(backend_status_t), intent(out) :: status

        integer :: hdferr
        integer :: sep_pos_1, sep_pos_2, sep_pos_3

        call status%set_ok()

        if (self%is_opened) then
            call self%close(status)
            if (status%code /= 0) return
        end if

        sep_pos_1 = index(trim(source_uri), '|')
        if (sep_pos_1 <= 1) then
            call status%set_error(2005, 'HDF5 backend URI must be DATA_URI|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_2 = index(trim(source_uri(sep_pos_1+1:len_trim(source_uri))), '|')
        if (sep_pos_2 <= 1) then
            call status%set_error(2005, 'HDF5 backend URI must be DATA_URI|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_2 = sep_pos_2 + sep_pos_1
        sep_pos_3 = index(trim(source_uri(sep_pos_2+1:len_trim(source_uri))), '|')
        if (sep_pos_3 <= 1) then
            call status%set_error(2005, 'HDF5 backend URI must be DATA_URI|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_3 = sep_pos_3 + sep_pos_2
        if (sep_pos_3 >= len_trim(source_uri)) then
            call status%set_error(2005, 'HDF5 backend URI must be DATA_URI|isoc_type|spec_type|dust_type.')
            return
        end if

        if (allocated(self%source_uri)) deallocate(self%source_uri)
        if (allocated(self%data_uri)) deallocate(self%data_uri)
        if (allocated(self%isoc_type)) deallocate(self%isoc_type)
        if (allocated(self%spec_type)) deallocate(self%spec_type)
        if (allocated(self%dust_type)) deallocate(self%dust_type)
        self%source_uri = trim(source_uri)
        self%data_uri = trim(source_uri(1:sep_pos_1-1))
        self%isoc_type = trim(to_lower(trim(source_uri(sep_pos_1+1:sep_pos_2-1))))
        self%spec_type = trim(to_lower(trim(source_uri(sep_pos_2+1:sep_pos_3-1))))
        self%dust_type = trim(source_uri(sep_pos_3+1:len_trim(source_uri)))
        if (trim(to_lower(self%dust_type)) == 'themis') then
            self%dust_type = 'THEMIS'
        else
            self%dust_type = 'DL07'
        end if

        if (len_trim(self%data_uri) == 0 .or. len_trim(self%isoc_type) == 0 .or. &
            len_trim(self%spec_type) == 0 .or. len_trim(self%dust_type) == 0) then
            call status%set_error(2006, 'HDF5 backend URI contains empty DATA_URI, isoc_type, spec_type, or dust_type.')
            return
        end if

        call h5open_f(hdferr)
        if (hdferr < 0) then
            call status%set_error(2001, 'HDF5 initialization failed in open().')
            return
        end if
        self%hdf5_initialized = .true.

        call h5fopen_f(trim(self%data_uri), H5F_ACC_RDONLY_F, self%file_id, hdferr)
        if (hdferr < 0) then
            call status%set_error(2002, 'Failed to open HDF5 file: '//trim(self%data_uri))
            self%file_id = -1_HID_T
            call h5close_f(hdferr)
            self%hdf5_initialized = .false.
            return
        end if

        self%is_opened = .true.
    end subroutine hdf5_backend_open

    !> @brief Close the HDF5 file and finalize HDF5 interface.
    subroutine hdf5_backend_close(self, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: hdferr

        call status%set_ok()

        if (self%is_opened .and. self%file_id >= 0_HID_T) then
            call h5fclose_f(self%file_id, hdferr)
            if (hdferr /= 0) then
                call status%set_error(2003, 'Failed to close HDF5 file identifier.')
                return
            end if
            self%file_id = -1_HID_T
        end if

        self%is_opened = .false.

        if (self%hdf5_initialized) then
            call h5close_f(hdferr)
            if (hdferr /= 0) then
                call status%set_error(2004, 'Failed to finalize HDF5 interface.')
                return
            end if
            self%hdf5_initialized = .false.
        end if

        if (allocated(self%source_uri)) deallocate(self%source_uri)
        if (allocated(self%data_uri)) deallocate(self%data_uri)
        if (allocated(self%isoc_type)) deallocate(self%isoc_type)
        if (allocated(self%spec_type)) deallocate(self%spec_type)
        if (allocated(self%dust_type)) deallocate(self%dust_type)
    end subroutine hdf5_backend_close

    !> @brief Return `.true.` if backend has an open file.
    logical function hdf5_backend_is_open(self)
        class(hdf5_backend_t), intent(in) :: self
        hdf5_backend_is_open = self%is_opened
    end function hdf5_backend_is_open

    !> @brief Read top-level manifest metadata.
    !>
    !> @details
    !> This interface is intentionally conservative for the initial contract:
    !> it validates file openness and populates core string fields with defaults.
    !> Rich manifest population can be expanded as schema implementation grows.
    subroutine hdf5_backend_read_manifest(self, manifest, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(library_manifest_t), intent(inout) :: manifest
        type(backend_status_t), intent(out) :: status
        character(len=:), allocatable :: neb_root, neb_candidate

        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2101, 'read_manifest called on closed backend.')
            return
        end if

        call manifest%clear()

        manifest%fsds_version = 'unknown'
        manifest%producer = 'unknown'
        if (allocated(self%isoc_type)) then
            manifest%isoc_name = trim(self%isoc_type)
        else
            manifest%isoc_name = ''
        end if
        if (allocated(self%spec_type)) then
            manifest%spec_name = trim(self%spec_type)
        else
            manifest%spec_name = ''
        end if
        if (allocated(self%dust_type)) then
            manifest%dust_name = trim(self%dust_type)
        else
            manifest%dust_name = ''
        end if

        neb_candidate = '/libraries/nebular/'//trim(self%isoc_type)
        if (self%has_path(trim(neb_candidate)//'/WD/line_pos')) then
            neb_root = trim(neb_candidate)
        else
            neb_root = '/libraries/nebular'
        end if

        allocate(manifest%datasets(61))
        call manifest%datasets(1)%clear()
        manifest%datasets(1)%role = 'spectral_base'
        manifest%datasets(1)%dims_csv = 'lambda,z,logt,logg'
        manifest%datasets(1)%path = '/libraries/spectra/'//trim(self%spec_type)//'/base/spectral_grid_nd'
        manifest%datasets(1)%unit = 'Lsun/Hz/Msun'
        manifest%datasets(1)%dtype = 'float32'
        manifest%datasets(1)%representation = 'dense_nd'

        call manifest%datasets(2)%clear()
        manifest%datasets(2)%role = 'isoc_nmass'
        manifest%datasets(2)%dims_csv = 'nt,nz'
        manifest%datasets(2)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/nmass'
        manifest%datasets(2)%dtype = 'int32'
        manifest%datasets(2)%representation = 'dense_nd'

        call manifest%datasets(3)%clear()
        manifest%datasets(3)%role = 'isoc_timestep'
        manifest%datasets(3)%dims_csv = 'nt,nz'
        manifest%datasets(3)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/timestep'
        manifest%datasets(3)%dtype = 'float64'
        manifest%datasets(3)%representation = 'dense_nd'

        call manifest%datasets(4)%clear()
        manifest%datasets(4)%role = 'isoc_mini'
        manifest%datasets(4)%dims_csv = 'nm,nt,nz'
        manifest%datasets(4)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/mini'
        manifest%datasets(4)%dtype = 'float64'
        manifest%datasets(4)%representation = 'dense_nd'

        call manifest%datasets(5)%clear()
        manifest%datasets(5)%role = 'isoc_mact'
        manifest%datasets(5)%dims_csv = 'nm,nt,nz'
        manifest%datasets(5)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/mact'
        manifest%datasets(5)%dtype = 'float64'
        manifest%datasets(5)%representation = 'dense_nd'

        call manifest%datasets(6)%clear()
        manifest%datasets(6)%role = 'isoc_logl'
        manifest%datasets(6)%dims_csv = 'nm,nt,nz'
        manifest%datasets(6)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/logl'
        manifest%datasets(6)%dtype = 'float64'
        manifest%datasets(6)%representation = 'dense_nd'

        call manifest%datasets(7)%clear()
        manifest%datasets(7)%role = 'isoc_logt'
        manifest%datasets(7)%dims_csv = 'nm,nt,nz'
        manifest%datasets(7)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/logt'
        manifest%datasets(7)%dtype = 'float64'
        manifest%datasets(7)%representation = 'dense_nd'

        call manifest%datasets(8)%clear()
        manifest%datasets(8)%role = 'isoc_logg'
        manifest%datasets(8)%dims_csv = 'nm,nt,nz'
        manifest%datasets(8)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/logg'
        manifest%datasets(8)%dtype = 'float64'
        manifest%datasets(8)%representation = 'dense_nd'

        call manifest%datasets(9)%clear()
        manifest%datasets(9)%role = 'isoc_phase'
        manifest%datasets(9)%dims_csv = 'nm,nt,nz'
        manifest%datasets(9)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/phase'
        manifest%datasets(9)%dtype = 'float64'
        manifest%datasets(9)%representation = 'dense_nd'

        call manifest%datasets(10)%clear()
        manifest%datasets(10)%role = 'isoc_ffco'
        manifest%datasets(10)%dims_csv = 'nm,nt,nz'
        manifest%datasets(10)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/ffco'
        manifest%datasets(10)%dtype = 'float64'
        manifest%datasets(10)%representation = 'dense_nd'

        call manifest%datasets(11)%clear()
        manifest%datasets(11)%role = 'isoc_lmdot'
        manifest%datasets(11)%dims_csv = 'nm,nt,nz'
        manifest%datasets(11)%path = '/libraries/isochrones/'//trim(self%isoc_type)//'/tracks/lmdot'
        manifest%datasets(11)%dtype = 'float64'
        manifest%datasets(11)%representation = 'dense_nd'

        call manifest%datasets(12)%clear()
        manifest%datasets(12)%role = 'nebular_wd_line_pos'
        manifest%datasets(12)%dims_csv = 'line'
        manifest%datasets(12)%path = trim(neb_root)//'/WD/line_pos'
        manifest%datasets(12)%dtype = 'float64'
        manifest%datasets(12)%representation = 'dense_nd'

        call manifest%datasets(13)%clear()
        manifest%datasets(13)%role = 'nebular_wd_logz'
        manifest%datasets(13)%dims_csv = 'z'
        manifest%datasets(13)%path = trim(neb_root)//'/WD/logz'
        manifest%datasets(13)%dtype = 'float64'
        manifest%datasets(13)%representation = 'dense_nd'

        call manifest%datasets(14)%clear()
        manifest%datasets(14)%role = 'nebular_wd_age'
        manifest%datasets(14)%dims_csv = 'age'
        manifest%datasets(14)%path = trim(neb_root)//'/WD/age'
        manifest%datasets(14)%dtype = 'float64'
        manifest%datasets(14)%representation = 'dense_nd'

        call manifest%datasets(15)%clear()
        manifest%datasets(15)%role = 'nebular_wd_logu'
        manifest%datasets(15)%dims_csv = 'u'
        manifest%datasets(15)%path = trim(neb_root)//'/WD/logu'
        manifest%datasets(15)%dtype = 'float64'
        manifest%datasets(15)%representation = 'dense_nd'

        call manifest%datasets(16)%clear()
        manifest%datasets(16)%role = 'nebular_wd_cont'
        manifest%datasets(16)%dims_csv = 'lam,z,age,u'
        manifest%datasets(16)%path = trim(neb_root)//'/WD/cont'
        manifest%datasets(16)%dtype = 'float32'
        manifest%datasets(16)%representation = 'dense_nd'

        call manifest%datasets(17)%clear()
        manifest%datasets(17)%role = 'nebular_wd_line'
        manifest%datasets(17)%dims_csv = 'line,z,age,u'
        manifest%datasets(17)%path = trim(neb_root)//'/WD/lines'
        manifest%datasets(17)%dtype = 'float32'
        manifest%datasets(17)%representation = 'dense_nd'

        call manifest%datasets(18)%clear()
        manifest%datasets(18)%role = 'nebular_nd_line_pos'
        manifest%datasets(18)%dims_csv = 'line'
        manifest%datasets(18)%path = trim(neb_root)//'/ND/line_pos'
        manifest%datasets(18)%dtype = 'float64'
        manifest%datasets(18)%representation = 'dense_nd'

        call manifest%datasets(19)%clear()
        manifest%datasets(19)%role = 'nebular_nd_logz'
        manifest%datasets(19)%dims_csv = 'z'
        manifest%datasets(19)%path = trim(neb_root)//'/ND/logz'
        manifest%datasets(19)%dtype = 'float64'
        manifest%datasets(19)%representation = 'dense_nd'

        call manifest%datasets(20)%clear()
        manifest%datasets(20)%role = 'nebular_nd_age'
        manifest%datasets(20)%dims_csv = 'age'
        manifest%datasets(20)%path = trim(neb_root)//'/ND/age'
        manifest%datasets(20)%dtype = 'float64'
        manifest%datasets(20)%representation = 'dense_nd'

        call manifest%datasets(21)%clear()
        manifest%datasets(21)%role = 'nebular_nd_logu'
        manifest%datasets(21)%dims_csv = 'u'
        manifest%datasets(21)%path = trim(neb_root)//'/ND/logu'
        manifest%datasets(21)%dtype = 'float64'
        manifest%datasets(21)%representation = 'dense_nd'

        call manifest%datasets(22)%clear()
        manifest%datasets(22)%role = 'nebular_nd_cont'
        manifest%datasets(22)%dims_csv = 'lam,z,age,u'
        manifest%datasets(22)%path = trim(neb_root)//'/ND/cont'
        manifest%datasets(22)%dtype = 'float32'
        manifest%datasets(22)%representation = 'dense_nd'

        call manifest%datasets(23)%clear()
        manifest%datasets(23)%role = 'nebular_nd_line'
        manifest%datasets(23)%dims_csv = 'line,z,age,u'
        manifest%datasets(23)%path = trim(neb_root)//'/ND/lines'
        manifest%datasets(23)%dtype = 'float32'
        manifest%datasets(23)%representation = 'dense_nd'

        call manifest%datasets(24)%clear()
        manifest%datasets(24)%role = 'wmb_logt'
        manifest%datasets(24)%dims_csv = 'logt'
        manifest%datasets(24)%path = '/libraries/auxiliary/wmbasic/logt'
        manifest%datasets(24)%dtype = 'float64'
        manifest%datasets(24)%representation = 'dense_nd'

        call manifest%datasets(25)%clear()
        manifest%datasets(25)%role = 'wmb_z'
        manifest%datasets(25)%dims_csv = 'z'
        manifest%datasets(25)%path = '/libraries/auxiliary/wmbasic/z'
        manifest%datasets(25)%dtype = 'float64'
        manifest%datasets(25)%representation = 'dense_nd'

        call manifest%datasets(26)%clear()
        manifest%datasets(26)%role = 'wmb_lam'
        manifest%datasets(26)%dims_csv = 'lam'
        manifest%datasets(26)%path = '/libraries/auxiliary/wmbasic/lam'
        manifest%datasets(26)%dtype = 'float64'
        manifest%datasets(26)%representation = 'dense_nd'

        call manifest%datasets(27)%clear()
        manifest%datasets(27)%role = 'wmb_spec'
        manifest%datasets(27)%dims_csv = 'lam,logt,logg,z'
        manifest%datasets(27)%path = '/libraries/auxiliary/wmbasic/spec'
        manifest%datasets(27)%dtype = 'float32'
        manifest%datasets(27)%representation = 'dense_nd'

        call manifest%datasets(28)%clear()
        manifest%datasets(28)%role = 'pagb_logt'
        manifest%datasets(28)%dims_csv = 'logt'
        manifest%datasets(28)%path = '/libraries/auxiliary/pagb/logt'
        manifest%datasets(28)%dtype = 'float64'
        manifest%datasets(28)%representation = 'dense_nd'

        call manifest%datasets(29)%clear()
        manifest%datasets(29)%role = 'pagb_lam'
        manifest%datasets(29)%dims_csv = 'lam'
        manifest%datasets(29)%path = '/libraries/auxiliary/pagb/lam'
        manifest%datasets(29)%dtype = 'float64'
        manifest%datasets(29)%representation = 'dense_nd'

        call manifest%datasets(30)%clear()
        manifest%datasets(30)%role = 'pagb_spec'
        manifest%datasets(30)%dims_csv = 'lam,logt,z'
        manifest%datasets(30)%path = '/libraries/auxiliary/pagb/spec'
        manifest%datasets(30)%dtype = 'float32'
        manifest%datasets(30)%representation = 'dense_nd'

        call manifest%datasets(31)%clear()
        manifest%datasets(31)%role = 'wr_logt_wn'
        manifest%datasets(31)%dims_csv = 'logt'
        manifest%datasets(31)%path = '/libraries/auxiliary/wr/logt_wn'
        manifest%datasets(31)%dtype = 'float64'
        manifest%datasets(31)%representation = 'dense_nd'

        call manifest%datasets(32)%clear()
        manifest%datasets(32)%role = 'wr_logt_wc'
        manifest%datasets(32)%dims_csv = 'logt'
        manifest%datasets(32)%path = '/libraries/auxiliary/wr/logt_wc'
        manifest%datasets(32)%dtype = 'float64'
        manifest%datasets(32)%representation = 'dense_nd'

        call manifest%datasets(33)%clear()
        manifest%datasets(33)%role = 'wr_z'
        manifest%datasets(33)%dims_csv = 'z'
        manifest%datasets(33)%path = '/libraries/auxiliary/wr/z'
        manifest%datasets(33)%dtype = 'float64'
        manifest%datasets(33)%representation = 'dense_nd'

        call manifest%datasets(34)%clear()
        manifest%datasets(34)%role = 'wr_lam'
        manifest%datasets(34)%dims_csv = 'lam'
        manifest%datasets(34)%path = '/libraries/auxiliary/wr/lam'
        manifest%datasets(34)%dtype = 'float64'
        manifest%datasets(34)%representation = 'dense_nd'

        call manifest%datasets(35)%clear()
        manifest%datasets(35)%role = 'wr_spec_wn'
        manifest%datasets(35)%dims_csv = 'lam,logt,z'
        manifest%datasets(35)%path = '/libraries/auxiliary/wr/spec_wn'
        manifest%datasets(35)%dtype = 'float32'
        manifest%datasets(35)%representation = 'dense_nd'

        call manifest%datasets(36)%clear()
        manifest%datasets(36)%role = 'wr_spec_wc'
        manifest%datasets(36)%dims_csv = 'lam,logt,z'
        manifest%datasets(36)%path = '/libraries/auxiliary/wr/spec_wc'
        manifest%datasets(36)%dtype = 'float32'
        manifest%datasets(36)%representation = 'dense_nd'

        call manifest%datasets(37)%clear()
        manifest%datasets(37)%role = 'agb_z_o'
        manifest%datasets(37)%dims_csv = 'z'
        manifest%datasets(37)%path = '/libraries/auxiliary/agb/z_o'
        manifest%datasets(37)%dtype = 'float64'
        manifest%datasets(37)%representation = 'dense_nd'

        call manifest%datasets(38)%clear()
        manifest%datasets(38)%role = 'agb_logt_c'
        manifest%datasets(38)%dims_csv = 'logt'
        manifest%datasets(38)%path = '/libraries/auxiliary/agb/logt_c'
        manifest%datasets(38)%dtype = 'float64'
        manifest%datasets(38)%representation = 'dense_nd'

        call manifest%datasets(39)%clear()
        manifest%datasets(39)%role = 'agb_logt_car'
        manifest%datasets(39)%dims_csv = 'logt'
        manifest%datasets(39)%path = '/libraries/auxiliary/agb/logt_car'
        manifest%datasets(39)%dtype = 'float64'
        manifest%datasets(39)%representation = 'dense_nd'

        call manifest%datasets(40)%clear()
        manifest%datasets(40)%role = 'agb_lam_o'
        manifest%datasets(40)%dims_csv = 'lam'
        manifest%datasets(40)%path = '/libraries/auxiliary/agb/lam_o'
        manifest%datasets(40)%dtype = 'float64'
        manifest%datasets(40)%representation = 'dense_nd'

        call manifest%datasets(41)%clear()
        manifest%datasets(41)%role = 'agb_lam_c'
        manifest%datasets(41)%dims_csv = 'lam'
        manifest%datasets(41)%path = '/libraries/auxiliary/agb/lam_c'
        manifest%datasets(41)%dtype = 'float64'
        manifest%datasets(41)%representation = 'dense_nd'

        call manifest%datasets(42)%clear()
        manifest%datasets(42)%role = 'agb_lam_car'
        manifest%datasets(42)%dims_csv = 'lam'
        manifest%datasets(42)%path = '/libraries/auxiliary/agb/lam_car'
        manifest%datasets(42)%dtype = 'float64'
        manifest%datasets(42)%representation = 'dense_nd'

        call manifest%datasets(43)%clear()
        manifest%datasets(43)%role = 'agb_logt_o'
        manifest%datasets(43)%dims_csv = 'z,logt'
        manifest%datasets(43)%path = '/libraries/auxiliary/agb/logt_o'
        manifest%datasets(43)%dtype = 'float64'
        manifest%datasets(43)%representation = 'dense_nd'

        call manifest%datasets(44)%clear()
        manifest%datasets(44)%role = 'agb_spec_o'
        manifest%datasets(44)%dims_csv = 'lam,logt'
        manifest%datasets(44)%path = '/libraries/auxiliary/agb/spec_o'
        manifest%datasets(44)%dtype = 'float32'
        manifest%datasets(44)%representation = 'dense_nd'

        call manifest%datasets(45)%clear()
        manifest%datasets(45)%role = 'agb_spec_c'
        manifest%datasets(45)%dims_csv = 'lam,logt'
        manifest%datasets(45)%path = '/libraries/auxiliary/agb/spec_c'
        manifest%datasets(45)%dtype = 'float32'
        manifest%datasets(45)%representation = 'dense_nd'

        call manifest%datasets(46)%clear()
        manifest%datasets(46)%role = 'agb_spec_car'
        manifest%datasets(46)%dims_csv = 'lam,logt'
        manifest%datasets(46)%path = '/libraries/auxiliary/agb/spec_car'
        manifest%datasets(46)%dtype = 'float32'
        manifest%datasets(46)%representation = 'dense_nd'

        call manifest%datasets(47)%clear()
        manifest%datasets(47)%role = 'dust_em_qpah'
        manifest%datasets(47)%dims_csv = 'qpah'
        manifest%datasets(47)%path = '/libraries/dust/emission/'//trim(self%dust_type)//'/qpah'
        manifest%datasets(47)%dtype = 'float64'
        manifest%datasets(47)%representation = 'dense_nd'

        call manifest%datasets(48)%clear()
        manifest%datasets(48)%role = 'dust_em_umin'
        manifest%datasets(48)%dims_csv = 'umin'
        manifest%datasets(48)%path = '/libraries/dust/emission/'//trim(self%dust_type)//'/umin'
        manifest%datasets(48)%dtype = 'float64'
        manifest%datasets(48)%representation = 'dense_nd'

        call manifest%datasets(49)%clear()
        manifest%datasets(49)%role = 'dust_em_lam'
        manifest%datasets(49)%dims_csv = 'lam'
        manifest%datasets(49)%path = '/libraries/dust/emission/'//trim(self%dust_type)//'/lam'
        manifest%datasets(49)%dtype = 'float64'
        manifest%datasets(49)%representation = 'dense_nd'

        call manifest%datasets(50)%clear()
        manifest%datasets(50)%role = 'dust_em_spec'
        manifest%datasets(50)%dims_csv = 'lam,qpah,umin'
        manifest%datasets(50)%path = '/libraries/dust/emission/'//trim(self%dust_type)//'/spec'
        manifest%datasets(50)%dtype = 'float32'
        manifest%datasets(50)%representation = 'dense_nd'

        call manifest%datasets(51)%clear()
        manifest%datasets(51)%role = 'agn_dust_tau'
        manifest%datasets(51)%dims_csv = 'tau'
        manifest%datasets(51)%path = '/libraries/dust/agn/tau'
        manifest%datasets(51)%dtype = 'float64'
        manifest%datasets(51)%representation = 'dense_nd'

        call manifest%datasets(52)%clear()
        manifest%datasets(52)%role = 'agn_dust_lam'
        manifest%datasets(52)%dims_csv = 'lam'
        manifest%datasets(52)%path = '/libraries/dust/agn/lam'
        manifest%datasets(52)%dtype = 'float64'
        manifest%datasets(52)%representation = 'dense_nd'

        call manifest%datasets(53)%clear()
        manifest%datasets(53)%role = 'agn_dust_spec'
        manifest%datasets(53)%dims_csv = 'lam,tau'
        manifest%datasets(53)%path = '/libraries/dust/agn/spec'
        manifest%datasets(53)%dtype = 'float32'
        manifest%datasets(53)%representation = 'dense_nd'

        call manifest%datasets(54)%clear()
        manifest%datasets(54)%role = 'dust_att_wg_lam'
        manifest%datasets(54)%dims_csv = 'lam'
        manifest%datasets(54)%path = '/libraries/dust/attenuation/wg_lam'
        manifest%datasets(54)%dtype = 'float64'
        manifest%datasets(54)%representation = 'dense_nd'

        call manifest%datasets(55)%clear()
        manifest%datasets(55)%role = 'dust_att_wg_spec'
        manifest%datasets(55)%dims_csv = 'lam,geom,tau,type'
        manifest%datasets(55)%path = '/libraries/dust/attenuation/wg_spec'
        manifest%datasets(55)%dtype = 'float32'
        manifest%datasets(55)%representation = 'dense_nd'

        call manifest%datasets(56)%clear()
        manifest%datasets(56)%role = 'dust_att_smc_lam'
        manifest%datasets(56)%dims_csv = 'lam'
        manifest%datasets(56)%path = '/libraries/dust/attenuation/smc_lam'
        manifest%datasets(56)%dtype = 'float64'
        manifest%datasets(56)%representation = 'dense_nd'

        call manifest%datasets(57)%clear()
        manifest%datasets(57)%role = 'dust_att_smc_ext'
        manifest%datasets(57)%dims_csv = 'lam'
        manifest%datasets(57)%path = '/libraries/dust/attenuation/smc_ext'
        manifest%datasets(57)%dtype = 'float64'
        manifest%datasets(57)%representation = 'dense_nd'

        call manifest%datasets(58)%clear()
        manifest%datasets(58)%role = 'xrb_lam'
        manifest%datasets(58)%dims_csv = 'lam'
        manifest%datasets(58)%path = '/libraries/xrb/lam'
        manifest%datasets(58)%dtype = 'float64'
        manifest%datasets(58)%representation = 'dense_nd'

        call manifest%datasets(59)%clear()
        manifest%datasets(59)%role = 'xrb_age'
        manifest%datasets(59)%dims_csv = 'age'
        manifest%datasets(59)%path = '/libraries/xrb/age'
        manifest%datasets(59)%dtype = 'float64'
        manifest%datasets(59)%representation = 'dense_nd'

        call manifest%datasets(60)%clear()
        manifest%datasets(60)%role = 'xrb_z'
        manifest%datasets(60)%dims_csv = 'z'
        manifest%datasets(60)%path = '/libraries/xrb/z'
        manifest%datasets(60)%dtype = 'float64'
        manifest%datasets(60)%representation = 'dense_nd'

        call manifest%datasets(61)%clear()
        manifest%datasets(61)%role = 'xrb_spec'
        manifest%datasets(61)%dims_csv = 'lam,age,z'
        manifest%datasets(61)%path = '/libraries/xrb/spec'
        manifest%datasets(61)%dtype = 'float32'
        manifest%datasets(61)%representation = 'dense_nd'
    end subroutine hdf5_backend_read_manifest

    !> @brief Check whether a path exists in the HDF5 file.
    logical function hdf5_backend_has_path(self, path)
        class(hdf5_backend_t), intent(in) :: self
        character(len=*), intent(in) :: path

        integer :: hdferr
        logical :: exists

        if (.not. self%is_opened) then
            hdf5_backend_has_path = .false.
            return
        end if

        exists = .false.
        call h5lexists_f(self%file_id, trim(path), exists, hdferr)
        if (hdferr < 0) then
            hdf5_backend_has_path = .false.
        else
            hdf5_backend_has_path = exists
        end if
    end function hdf5_backend_has_path

    !> @brief Read one axis descriptor.
    subroutine hdf5_backend_query_axis(self, axis_name, axis_desc, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: axis_name
        type(axis_desc_t), intent(inout) :: axis_desc
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path
        real(WP), allocatable :: vals(:)

        call status%set_ok()
        call axis_desc%clear()

        if (.not. self%is_opened) then
            call status%set_error(2201, 'query_axis called on closed backend.')
            return
        end if

        path = resolve_axis_dataset_path(axis_name, self%spec_type)
        if (.not. self%has_path(path)) then
            path = '/axes/'//trim(axis_name)
            if (.not. self%has_path(path)) then
                call status%set_error(2202, 'Axis path not found: '//trim(path))
                return
            end if
        end if

        call self%read_real_1d(path, vals, status)
        if (status%code /= 0) return

        axis_desc%name = trim(axis_name)
        axis_desc%path = path
        axis_desc%unit = ''
        axis_desc%n = size(vals)
        allocate(axis_desc%values(axis_desc%n))
        axis_desc%values = vals
    end subroutine hdf5_backend_query_axis

    !> @brief Read a real 1-D dataset.
    subroutine hdf5_backend_read_real_1d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T) :: npoints
        integer(HSIZE_T), dimension(1) :: read_dims
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2301, 'read_real_1d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr < 0) then
            call status%set_error(2302, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr < 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2303, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 1) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2304, 'Dataset is not rank-1: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_npoints_f(space_id, npoints, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2305, 'Failed to read rank-1 size: '//trim(dataset_path))
            return
        end if

        if (npoints <= 0_HSIZE_T) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2305, 'Rank-1 dataset has zero size: '//trim(dataset_path))
            return
        end if

        h5_real_type = get_h5_real_type()

        read_dims(1) = npoints
        allocate(values(int(npoints)))
        call h5dread_f(dset_id, h5_real_type, values, read_dims, hdferr)
        if (hdferr < 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2306, 'Failed to read rank-1 real dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_real_1d

    !> @brief Read a real 2-D dataset.
    subroutine hdf5_backend_read_real_2d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2311, 'read_real_2d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr < 0) then
            call status%set_error(2312, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr < 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2313, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 2) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2314, 'Dataset is not rank-2: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2315, 'Failed to read rank-2 dimensions: '//trim(dataset_path))
            return
        end if

        h5_real_type = get_h5_real_type()

        allocate(values(dims(1), dims(2)))
        call h5dread_f(dset_id, h5_real_type, values, dims, hdferr)
        if (hdferr < 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2316, 'Failed to read rank-2 real dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_real_2d

    !> @brief Read a real 3-D dataset.
    subroutine hdf5_backend_read_real_3d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2317, 'read_real_3d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr < 0) then
            call status%set_error(2318, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr < 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2319, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 3) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2320, 'Dataset is not rank-3: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2327, 'Failed to read rank-3 dimensions: '//trim(dataset_path))
            return
        end if

        h5_real_type = get_h5_real_type()

        allocate(values(dims(1), dims(2), dims(3)))
        call h5dread_f(dset_id, h5_real_type, values, dims, hdferr)
        if (hdferr < 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2328, 'Failed to read rank-3 real dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_real_3d

    !> @brief Read a real 4-D dataset.
    subroutine hdf5_backend_read_real_4d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:,:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2329, 'read_real_4d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr < 0) then
            call status%set_error(2330, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr < 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2331, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 4) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2332, 'Dataset is not rank-4: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2333, 'Failed to read rank-4 dimensions: '//trim(dataset_path))
            return
        end if

        h5_real_type = get_h5_real_type()

        allocate(values(dims(1), dims(2), dims(3), dims(4)))
        call h5dread_f(dset_id, h5_real_type, values, dims, hdferr)
        if (hdferr /= 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2334, 'Failed to read rank-4 real dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_real_4d

    !> @brief Read an integer 1-D dataset.
    subroutine hdf5_backend_read_int_1d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2321, 'read_int_1d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2322, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2323, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 1) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2324, 'Dataset is not rank-1: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2325, 'Failed to read rank-1 dimensions: '//trim(dataset_path))
            return
        end if

        allocate(values(dims(1)))
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2326, 'Failed to read rank-1 integer dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_int_1d

    !> @brief Read an integer 2-D dataset.
    subroutine hdf5_backend_read_int_2d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2331, 'read_int_2d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2332, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2333, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 2) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2334, 'Dataset is not rank-2: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2335, 'Failed to read rank-2 dimensions: '//trim(dataset_path))
            return
        end if

        allocate(values(dims(1), dims(2)))
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2336, 'Failed to read rank-2 integer dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_int_2d

    !> @brief Read an integer 3-D dataset.
    subroutine hdf5_backend_read_int_3d(self, dataset_path, values, status)
        class(hdf5_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:,:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
        integer :: rank, hdferr

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(2337, 'read_int_3d called on closed backend.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset_path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2338, 'Failed to open dataset: '//trim(dataset_path))
            return
        end if

        call h5dget_space_f(dset_id, space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2339, 'Failed to get dataspace: '//trim(dataset_path))
            return
        end if

        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (hdferr /= 0 .or. rank /= 3) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2340, 'Dataset is not rank-3: '//trim(dataset_path))
            return
        end if

        allocate(dims(rank), maxdims(rank))
        call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2341, 'Failed to read rank-3 dimensions: '//trim(dataset_path))
            return
        end if

        allocate(values(dims(1), dims(2), dims(3)))
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) then
            if (allocated(values)) deallocate(values)
            call h5sclose_f(space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2342, 'Failed to read rank-3 integer dataset: '//trim(dataset_path))
            return
        end if

        call h5sclose_f(space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_int_3d

    !> @brief Read one spectral hyperslab at fixed `(iz, iafe)`.
    !>
    !> @details
    !> Uses HDF5 dataspace hyperslab selection, never full-cube read.
    subroutine hdf5_backend_read_spectral_slice(self, dataset, iz, iafe, slice, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        integer, intent(in) :: iz
        integer, intent(in) :: iafe
        type(spectral_slice_t), intent(inout) :: slice
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, file_space_id, mem_space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T), allocatable :: dims(:), offset(:), count(:)
        type(axis_desc_t) :: ax_lambda, ax_z, ax_logt, ax_logg, ax_afe
        real(WP), allocatable :: tmp4(:,:,:,:), tmp5(:,:,:,:,:)
        character(len=:), allocatable :: axis_path
        logical :: has_afe
        integer :: rank, hdferr

        call status%set_ok()
        call slice%clear()

        if (.not. self%is_opened) then
            call status%set_error(2401, 'read_spectral_slice called on closed backend.')
            return
        end if

        if (.not. allocated(dataset%path)) then
            call status%set_error(2402, 'Dataset path is unallocated in read_spectral_slice.')
            return
        end if
        if (.not. allocated(dataset%dims_csv)) then
            call status%set_error(2403, 'Dataset dims_csv is unallocated in read_spectral_slice.')
            return
        end if

        has_afe = dataset%has_axis('afe') .or. dataset%has_axis('alpha_fe')
        h5_real_type = get_h5_real_type()

        call h5dopen_f(self%file_id, trim(dataset%path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2404, 'Failed to open spectral dataset: '//trim(dataset%path))
            return
        end if

        call h5dget_space_f(dset_id, file_space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2405, 'Failed to get spectral dataspace.')
            return
        end if

        call h5sget_simple_extent_ndims_f(file_space_id, rank, hdferr)
        if (hdferr /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2406, 'Failed to get spectral rank.')
            return
        end if

        allocate(dims(rank), offset(rank), count(rank))

        axis_path = resolve_axis_dataset_path('lambda', self%spec_type)
        if (.not. self%has_path(axis_path)) axis_path = '/axes/lambda'
        call self%read_real_1d(axis_path, ax_lambda%values, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2407, 'Failed to query lambda axis for spectral dimensions.')
            return
        end if
        ax_lambda%n = size(ax_lambda%values)

        axis_path = resolve_axis_dataset_path('z', self%spec_type)
        if (.not. self%has_path(axis_path)) axis_path = '/axes/z'
        call self%read_real_1d(axis_path, ax_z%values, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2407, 'Failed to query z axis for spectral dimensions.')
            return
        end if
        ax_z%n = size(ax_z%values)

        axis_path = resolve_axis_dataset_path('logt', self%spec_type)
        if (.not. self%has_path(axis_path)) axis_path = '/axes/logt'
        call self%read_real_1d(axis_path, ax_logt%values, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2407, 'Failed to query logt axis for spectral dimensions.')
            return
        end if
        ax_logt%n = size(ax_logt%values)

        axis_path = resolve_axis_dataset_path('logg', self%spec_type)
        if (.not. self%has_path(axis_path)) axis_path = '/axes/logg'
        call self%read_real_1d(axis_path, ax_logg%values, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2407, 'Failed to query logg axis for spectral dimensions.')
            return
        end if
        ax_logg%n = size(ax_logg%values)

        if (has_afe) then
            axis_path = resolve_axis_dataset_path('afe', self%spec_type)
            if (.not. self%has_path(axis_path)) axis_path = '/axes/afe'
            call self%read_real_1d(axis_path, ax_afe%values, status)
            if (status%code /= 0) then
                axis_path = resolve_axis_dataset_path('alpha_fe', self%spec_type)
                if (.not. self%has_path(axis_path)) axis_path = '/axes/alpha_fe'
                call self%read_real_1d(axis_path, ax_afe%values, status)
                if (status%code /= 0) then
                    call h5sclose_f(file_space_id, hdferr)
                    call h5dclose_f(dset_id, hdferr)
                    call status%set_error(2407, 'Failed to query afe axis for spectral dimensions.')
                    return
                end if
            end if
            ax_afe%n = size(ax_afe%values)
            dims = [int(ax_lambda%n,HSIZE_T), int(ax_z%n,HSIZE_T), int(ax_afe%n,HSIZE_T), &
                    int(ax_logt%n,HSIZE_T), int(ax_logg%n,HSIZE_T)]
        else
            dims = [int(ax_lambda%n,HSIZE_T), int(ax_z%n,HSIZE_T), int(ax_logt%n,HSIZE_T), int(ax_logg%n,HSIZE_T)]
        end if

        if (has_afe .and. rank /= 5) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2408, 'Expected rank-5 spectral dataset for afe-enabled data.')
            return
        end if
        if ((.not. has_afe) .and. rank /= 4) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2409, 'Expected rank-4 spectral dataset for afe-degenerate data.')
            return
        end if

        if (iz < 1 .or. iz > int(dims(2))) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2410, 'z index out of bounds in read_spectral_slice.')
            return
        end if

        if (has_afe) then
            if (iafe < 1 .or. iafe > int(dims(3))) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2411, 'afe index out of bounds in read_spectral_slice.')
                return
            end if

            offset = [0_HSIZE_T, int(iz-1,HSIZE_T), int(iafe-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [dims(1), 1_HSIZE_T, 1_HSIZE_T, dims(4), dims(5)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2412, 'Failed to select spectral hyperslab (rank-5).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2413, 'Failed to create memory dataspace (rank-5).')
                return
            end if

            allocate(tmp5(count(1), count(2), count(3), count(4), count(5)))
            call h5dread_f(dset_id, h5_real_type, tmp5, count, hdferr, file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2414, 'Failed to read spectral hyperslab (rank-5).')
                return
            end if

            allocate(slice%flux(count(1), count(4), count(5)))
            slice%flux = tmp5(:, 1, 1, :, :)

            call h5sclose_f(mem_space_id, hdferr)

        else
            if (iafe /= 1) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2415, 'afe index must be 1 for rank-4 spectral dataset.')
                return
            end if

            offset = [0_HSIZE_T, int(iz-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [dims(1), 1_HSIZE_T, dims(3), dims(4)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2416, 'Failed to select spectral hyperslab (rank-4).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2417, 'Failed to create memory dataspace (rank-4).')
                return
            end if

            allocate(tmp4(count(1), count(2), count(3), count(4)))
            call h5dread_f(dset_id, h5_real_type, tmp4, count, hdferr, file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2418, 'Failed to read spectral hyperslab (rank-4).')
                return
            end if

            allocate(slice%flux(count(1), count(3), count(4)))
            slice%flux = tmp4(:, 1, :, :)

            call h5sclose_f(mem_space_id, hdferr)
        end if

        slice%iz = iz
        slice%iafe = iafe
        if (dataset%has_missing_value) then
            slice%missing_value = dataset%missing_value
        else
            slice%missing_value = -1.0e99_wp
        end if

        call maybe_read_slice_valid_mask(self, dataset, has_afe, iz, iafe, size(slice%flux,2), &
                                         size(slice%flux,3), slice%valid, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            return
        end if

        call h5sclose_f(file_space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_spectral_slice

    !> @brief Read a spectral neighborhood hyperslab over `(z, afe)`.
    subroutine hdf5_backend_read_spectral_neighborhood(self, dataset, iz_lo, iz_hi, iafe_lo, iafe_hi, neighborhood, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        integer, intent(in) :: iz_lo
        integer, intent(in) :: iz_hi
        integer, intent(in) :: iafe_lo
        integer, intent(in) :: iafe_hi
        type(spectral_grid_t), intent(inout) :: neighborhood
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, file_space_id, mem_space_id
        integer(HID_T) :: h5_real_type
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:), offset(:), count(:)
        real(WP), allocatable :: tmp4(:,:,:,:)
        logical :: has_afe
        integer :: rank, zcount, acount, hdferr

        call status%set_ok()
        call neighborhood%clear()

        if (.not. self%is_opened) then
            call status%set_error(2501, 'read_spectral_neighborhood called on closed backend.')
            return
        end if

        if (.not. allocated(dataset%path)) then
            call status%set_error(2502, 'Dataset path is unallocated in read_spectral_neighborhood.')
            return
        end if
        if (.not. allocated(dataset%dims_csv)) then
            call status%set_error(2503, 'Dataset dims_csv is unallocated in read_spectral_neighborhood.')
            return
        end if

        has_afe = dataset%has_axis('afe') .or. dataset%has_axis('alpha_fe')
        h5_real_type = get_h5_real_type()

        call h5dopen_f(self%file_id, trim(dataset%path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2504, 'Failed to open spectral dataset: '//trim(dataset%path))
            return
        end if

        call h5dget_space_f(dset_id, file_space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2505, 'Failed to get spectral dataspace.')
            return
        end if

        call h5sget_simple_extent_ndims_f(file_space_id, rank, hdferr)
        if (hdferr /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2506, 'Failed to get spectral rank.')
            return
        end if

        allocate(dims(rank), maxdims(rank), offset(rank), count(rank))
        call h5sget_simple_extent_dims_f(file_space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2507, 'Failed to read spectral dimensions.')
            return
        end if

        if (iz_lo < 1 .or. iz_hi < iz_lo .or. iz_hi > int(dims(2))) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2508, 'Invalid z neighborhood bounds.')
            return
        end if
        zcount = iz_hi - iz_lo + 1

        if (has_afe) then
            if (rank /= 5) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2509, 'Expected rank-5 spectral dataset for afe-enabled neighborhood.')
                return
            end if
            if (iafe_lo < 1 .or. iafe_hi < iafe_lo .or. iafe_hi > int(dims(3))) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2510, 'Invalid afe neighborhood bounds.')
                return
            end if
            acount = iafe_hi - iafe_lo + 1

            offset = [0_HSIZE_T, int(iz_lo-1,HSIZE_T), int(iafe_lo-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [dims(1), int(zcount,HSIZE_T), int(acount,HSIZE_T), dims(4), dims(5)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2511, 'Failed to select neighborhood hyperslab (rank-5).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2512, 'Failed to create neighborhood memspace (rank-5).')
                return
            end if

            allocate(neighborhood%flux(count(1), count(2), count(3), count(4), count(5)))
            call h5dread_f(dset_id, h5_real_type, neighborhood%flux, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2513, 'Failed to read neighborhood hyperslab (rank-5).')
                return
            end if
            call h5sclose_f(mem_space_id, hdferr)

        else
            if (rank /= 4) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2514, 'Expected rank-4 spectral dataset for afe-degenerate neighborhood.')
                return
            end if
            if (iafe_lo /= 1 .or. iafe_hi /= 1) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2515, 'afe neighborhood must be [1,1] for rank-4 spectral data.')
                return
            end if
            acount = 1

            offset = [0_HSIZE_T, int(iz_lo-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [dims(1), int(zcount,HSIZE_T), dims(3), dims(4)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2516, 'Failed to select neighborhood hyperslab (rank-4).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2517, 'Failed to create neighborhood memspace (rank-4).')
                return
            end if

            allocate(tmp4(count(1), count(2), count(3), count(4)))
            call h5dread_f(dset_id, h5_real_type, tmp4, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2518, 'Failed to read neighborhood hyperslab (rank-4).')
                return
            end if

            allocate(neighborhood%flux(count(1), count(2), 1, count(3), count(4)))
            neighborhood%flux(:, :, 1, :, :) = tmp4
            call h5sclose_f(mem_space_id, hdferr)
        end if

        if (dataset%has_missing_value) then
            neighborhood%missing_value = dataset%missing_value
        else
            neighborhood%missing_value = -1.0e99_wp
        end if

        call maybe_read_neighborhood_valid_mask(self, dataset, has_afe, iz_lo, iz_hi, iafe_lo, iafe_hi, &
                                                size(neighborhood%flux, 4), size(neighborhood%flux, 5), &
                                                neighborhood%valid, status)
        if (status%code /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            return
        end if

        call h5sclose_f(file_space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine hdf5_backend_read_spectral_neighborhood

    !> @brief Read optional 2-D validity mask for one spectral slice.
    subroutine maybe_read_slice_valid_mask(self, dataset, has_afe, iz, iafe, n_logt, n_logg, valid, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        logical, intent(in) :: has_afe
        integer, intent(in) :: iz, iafe, n_logt, n_logg
        logical, allocatable, intent(inout) :: valid(:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, file_space_id, mem_space_id
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:), offset(:), count(:)
        integer, allocatable :: ibuf3(:,:,:), ibuf2(:,:)
        integer :: rank, hdferr

        call status%set_ok()
        if (allocated(valid)) deallocate(valid)

        if (.not. dataset%has_valid_mask) return
        if (.not. allocated(dataset%valid_mask_path)) then
            call status%set_error(2601, 'Dataset declares has_valid_mask but valid_mask_path is missing.')
            return
        end if

        call h5dopen_f(self%file_id, trim(dataset%valid_mask_path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2602, 'Failed to open validity mask dataset: '//trim(dataset%valid_mask_path))
            return
        end if

        call h5dget_space_f(dset_id, file_space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2603, 'Failed to get validity mask dataspace.')
            return
        end if

        call h5sget_simple_extent_ndims_f(file_space_id, rank, hdferr)
        if (hdferr /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2604, 'Failed to read validity mask rank.')
            return
        end if

        allocate(dims(rank), maxdims(rank), offset(rank), count(rank))
        call h5sget_simple_extent_dims_f(file_space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2605, 'Failed to read validity mask dimensions.')
            return
        end if

        if (has_afe) then
            if (rank /= 4) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2606, 'Expected rank-4 validity mask for afe-enabled data.')
                return
            end if
            offset = [int(iz-1,HSIZE_T), int(iafe-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [1_HSIZE_T, 1_HSIZE_T, int(n_logt,HSIZE_T), int(n_logg,HSIZE_T)]
            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2607, 'Failed selecting validity mask hyperslab (rank-4).')
                return
            end if
            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2608, 'Failed creating mask memspace (rank-4).')
                return
            end if
            allocate(ibuf3(1, n_logt, n_logg))
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ibuf3, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2609, 'Failed reading validity mask hyperslab (rank-4).')
                return
            end if
            allocate(valid(n_logt, n_logg))
            valid = (ibuf3(1, :, :) /= 0)
            call h5sclose_f(mem_space_id, hdferr)
        else
            if (rank /= 3) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2610, 'Expected rank-3 validity mask for afe-degenerate data.')
                return
            end if
            offset = [int(iz-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [1_HSIZE_T, int(n_logt,HSIZE_T), int(n_logg,HSIZE_T)]
            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2611, 'Failed selecting validity mask hyperslab (rank-3).')
                return
            end if
            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2612, 'Failed creating mask memspace (rank-3).')
                return
            end if
            allocate(ibuf2(n_logt, n_logg))
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ibuf2, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2613, 'Failed reading validity mask hyperslab (rank-3).')
                return
            end if
            allocate(valid(n_logt, n_logg))
            valid = (ibuf2 /= 0)
            call h5sclose_f(mem_space_id, hdferr)
        end if

        call h5sclose_f(file_space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine maybe_read_slice_valid_mask

    !> @brief Read optional 4-D validity mask for a spectral neighborhood.
    subroutine maybe_read_neighborhood_valid_mask(self, dataset, has_afe, iz_lo, iz_hi, iafe_lo, iafe_hi, &
                                                  n_logt, n_logg, valid, status)
        class(hdf5_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        logical, intent(in) :: has_afe
        integer, intent(in) :: iz_lo, iz_hi, iafe_lo, iafe_hi, n_logt, n_logg
        logical, allocatable, intent(inout) :: valid(:,:,:,:)
        type(backend_status_t), intent(out) :: status

        integer(HID_T) :: dset_id, file_space_id, mem_space_id
        integer(HSIZE_T), allocatable :: dims(:), maxdims(:), offset(:), count(:)
        integer, allocatable :: ibuf4(:,:,:,:), ibuf3(:,:,:)
        integer :: rank, zcount, acount, hdferr

        call status%set_ok()
        if (allocated(valid)) deallocate(valid)

        if (.not. dataset%has_valid_mask) return
        if (.not. allocated(dataset%valid_mask_path)) then
            call status%set_error(2701, 'Dataset declares has_valid_mask but valid_mask_path is missing.')
            return
        end if

        zcount = iz_hi - iz_lo + 1
        acount = iafe_hi - iafe_lo + 1

        call h5dopen_f(self%file_id, trim(dataset%valid_mask_path), dset_id, hdferr)
        if (hdferr /= 0) then
            call status%set_error(2702, 'Failed to open neighborhood validity mask dataset.')
            return
        end if

        call h5dget_space_f(dset_id, file_space_id, hdferr)
        if (hdferr /= 0) then
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2703, 'Failed to get neighborhood mask dataspace.')
            return
        end if

        call h5sget_simple_extent_ndims_f(file_space_id, rank, hdferr)
        if (hdferr /= 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2704, 'Failed to read neighborhood mask rank.')
            return
        end if

        allocate(dims(rank), maxdims(rank), offset(rank), count(rank))
        call h5sget_simple_extent_dims_f(file_space_id, dims, maxdims, hdferr)
        if (hdferr < 0) then
            call h5sclose_f(file_space_id, hdferr)
            call h5dclose_f(dset_id, hdferr)
            call status%set_error(2705, 'Failed to read neighborhood mask dimensions.')
            return
        end if

        if (has_afe) then
            if (rank /= 4) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2706, 'Expected rank-4 mask for afe-enabled neighborhood.')
                return
            end if

            offset = [int(iz_lo-1,HSIZE_T), int(iafe_lo-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [int(zcount,HSIZE_T), int(acount,HSIZE_T), int(n_logt,HSIZE_T), int(n_logg,HSIZE_T)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2707, 'Failed selecting neighborhood mask hyperslab (rank-4).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2708, 'Failed creating neighborhood mask memspace (rank-4).')
                return
            end if

            allocate(ibuf4(zcount, acount, n_logt, n_logg))
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ibuf4, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2709, 'Failed reading neighborhood mask hyperslab (rank-4).')
                return
            end if

            allocate(valid(zcount, acount, n_logt, n_logg))
            valid = (ibuf4 /= 0)
            call h5sclose_f(mem_space_id, hdferr)

        else
            if (rank /= 3) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2710, 'Expected rank-3 mask for afe-degenerate neighborhood.')
                return
            end if

            offset = [int(iz_lo-1,HSIZE_T), 0_HSIZE_T, 0_HSIZE_T]
            count  = [int(zcount,HSIZE_T), int(n_logt,HSIZE_T), int(n_logg,HSIZE_T)]

            call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, offset, count, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2711, 'Failed selecting neighborhood mask hyperslab (rank-3).')
                return
            end if

            call h5screate_simple_f(rank, count, mem_space_id, hdferr)
            if (hdferr /= 0) then
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2712, 'Failed creating neighborhood mask memspace (rank-3).')
                return
            end if

            allocate(ibuf3(zcount, n_logt, n_logg))
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, ibuf3, count, hdferr, &
                           file_space_id=file_space_id, mem_space_id=mem_space_id)
            if (hdferr /= 0) then
                call h5sclose_f(mem_space_id, hdferr)
                call h5sclose_f(file_space_id, hdferr)
                call h5dclose_f(dset_id, hdferr)
                call status%set_error(2713, 'Failed reading neighborhood mask hyperslab (rank-3).')
                return
            end if

            allocate(valid(zcount, 1, n_logt, n_logg))
            valid(:, 1, :, :) = (ibuf3 /= 0)
            call h5sclose_f(mem_space_id, hdferr)
        end if

        call h5sclose_f(file_space_id, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine maybe_read_neighborhood_valid_mask

    !> @brief Map FSPS working precision to a matching native HDF5 real type.
    !>
    !> @details
    !> Returns `H5T_NATIVE_DOUBLE` when `WP` is double precision, otherwise
    !> returns `H5T_NATIVE_REAL` for single precision builds.
    pure function get_h5_real_type() result(h5_type)
        integer(HID_T) :: h5_type

        if (kind(1.0_wp) == kind(1.0d0)) then
            h5_type = H5T_NATIVE_DOUBLE
        else
            h5_type = H5T_NATIVE_REAL
        end if
    end function get_h5_real_type

end module fsps_data_backend_hdf5
