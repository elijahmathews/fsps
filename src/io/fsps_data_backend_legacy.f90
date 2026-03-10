module fsps_data_backend_legacy
    !> @brief Legacy FSPS backend implementing `data_backend_t` for spectral reads.

    use fsps_precision, only: WP
    use fsps_constants, only: NDIM_LOGT, NDIM_LOGG, NDIM_WMB_LOGT, NDIM_WMB_LOGG, NDIM_PAGB, NDIM_WR, &
                              N_AGB_O, N_AGB_C, N_AGB_CAR, BASEL_STR, NM, NLINES, NEMLINE, &
                              NLAM_NEBCONT, NEBNZ, NEBNAGE, NEBNIP, NAGNDUST, NAGNDUST_SPEC, &
                                                            NSPEC_XRB, NT_XRB, NZ_XRB, &
                              NT_MIST, NZ_MIST, NT_PADOVA, NZ_PADOVA, NT_PARSEC, NZ_PARSEC, &
                              NT_BASTI, NZ_BASTI, NT_GENEVA, NZ_GENEVA
    use fsps_data_schema, only: axis_desc_t, dataset_desc_t, library_manifest_t, &
                                spectral_grid_t, spectral_slice_t, isochrone_grid_t, nebular_grid_t, &
                                aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, dust_emission_t, &
                                agn_dust_t, dust_attenuation_t, xrb_spectra_t
    use fsps_data_backend, only: data_backend_t, backend_status_t
    use fsps_interpolation, only: interpolate_linear
    use fsps_strings, only: to_lower
    use, intrinsic :: iso_fortran_env, only: real32, file_storage_size, error_unit, iostat_end

    implicit none
    private

    integer, parameter :: NLAMWR = 1963
    integer, parameter :: NSPEC_AGB = 6146
    integer, parameter :: NSPEC_ARINGER = 9032
    integer, parameter :: NT_BPASS = 43
    integer, parameter :: NZ_BPASS = 12
    real(WP), parameter :: QPAH_ARR_DL07(7) = [ &
        0.47_wp, 1.12_wp, 1.77_wp, 2.50_wp, 3.19_wp, 3.90_wp, 4.58_wp &
    ]
    real(WP), parameter :: UMIN_ARR_DL07(22) = [ &
        0.1_wp,  0.15_wp, 0.2_wp,  0.3_wp,  0.4_wp, 0.5_wp, &
        0.7_wp,  0.8_wp,  1.0_wp,  1.2_wp,  1.5_wp, 2.0_wp, &
        2.5_wp,  3.0_wp,  4.0_wp,  5.0_wp,  7.0_wp, 8.0_wp, &
        12.0_wp, 15.0_wp, 20.0_wp, 25.0_wp &
    ]
    real(WP), parameter :: QPAH_ARR_THEMIS(11) = (100.0_wp / 2.2_wp) * [ &
        0.02_wp, 0.06_wp, 0.10_wp, 0.14_wp, 0.17_wp, 0.20_wp, &
        0.24_wp, 0.28_wp, 0.32_wp, 0.36_wp, 0.40_wp &
    ]
    real(WP), parameter :: UMIN_ARR_THEMIS(37) = [ &
        0.1_wp,  0.12_wp, 0.15_wp, 0.17_wp, 0.2_wp,  0.25_wp, &
        0.3_wp,  0.35_wp, 0.4_wp,  0.5_wp,  0.6_wp,  0.7_wp, &
        0.8_wp,  1.0_wp,  1.2_wp,  1.5_wp,  1.7_wp,  2.0_wp, &
        2.5_wp,  3.0_wp,  3.5_wp,  4.0_wp,  5.0_wp,  6.0_wp, &
        7.0_wp,  8.0_wp,  10.0_wp, 12.0_wp, 15.0_wp, 17.0_wp, &
        20.0_wp, 25.0_wp, 30.0_wp, 35.0_wp, 40.0_wp, 50.0_wp, &
        80.0_wp &
    ]

    public :: legacy_backend_t

    type, extends(data_backend_t) :: legacy_backend_t
        private
        logical :: is_opened = .false.
        character(len=:), allocatable :: source_uri
        character(len=:), allocatable :: sps_home
        character(len=:), allocatable :: isoc_type
        character(len=:), allocatable :: spec_type
        character(len=:), allocatable :: dust_type
        integer :: nspec = 0
        real(WP), allocatable :: zlegend(:)
        real(WP), allocatable :: spec_lambda(:)
        type(isochrone_grid_t) :: iso_cache
        type(nebular_grid_t) :: neb_cache_wd
        type(nebular_grid_t) :: neb_cache_nd
        type(aux_wmbasic_t) :: wmb_cache
        type(aux_pagb_t) :: pagb_cache
        type(aux_wr_t) :: wr_cache
        type(aux_agb_t) :: agb_cache
        type(dust_emission_t) :: dust_em_cache
        type(agn_dust_t) :: agn_dust_cache
        type(dust_attenuation_t) :: dust_att_cache
        type(xrb_spectra_t) :: xrb_cache
    contains
        procedure, public :: open => legacy_backend_open
        procedure, public :: close => legacy_backend_close
        procedure, public :: is_open => legacy_backend_is_open

        procedure, public :: read_manifest => legacy_backend_read_manifest
        procedure, public :: has_path => legacy_backend_has_path
        procedure, public :: query_axis => legacy_backend_query_axis

        procedure, public :: read_real_1d => legacy_backend_read_real_1d
        procedure, public :: read_real_2d => legacy_backend_read_real_2d
        procedure, public :: read_real_3d => legacy_backend_read_real_3d
        procedure, public :: read_real_4d => legacy_backend_read_real_4d
        procedure, public :: read_int_1d => legacy_backend_read_int_1d
        procedure, public :: read_int_2d => legacy_backend_read_int_2d
        procedure, public :: read_int_3d => legacy_backend_read_int_3d

        procedure, public :: read_spectral_slice => legacy_backend_read_spectral_slice
        procedure, public :: read_spectral_neighborhood => legacy_backend_read_spectral_neighborhood
    end type legacy_backend_t

contains

    subroutine legacy_backend_open(self, source_uri, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: source_uri
        type(backend_status_t), intent(out) :: status

        integer :: sep_pos_1, sep_pos_2, sep_pos_3

        call status%set_ok()
        call clear_open_state(self)

        sep_pos_1 = index(trim(source_uri), '|')
        if (sep_pos_1 <= 1) then
            call status%set_error(3001, 'Legacy backend URI must be SPS_HOME|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_2 = index(trim(source_uri(sep_pos_1+1:len_trim(source_uri))), '|')
        if (sep_pos_2 <= 1) then
            call status%set_error(3001, 'Legacy backend URI must be SPS_HOME|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_2 = sep_pos_2 + sep_pos_1
        sep_pos_3 = index(trim(source_uri(sep_pos_2+1:len_trim(source_uri))), '|')
        if (sep_pos_3 <= 1) then
            call status%set_error(3001, 'Legacy backend URI must be SPS_HOME|isoc_type|spec_type|dust_type.')
            return
        end if
        sep_pos_3 = sep_pos_3 + sep_pos_2
        if (sep_pos_3 >= len_trim(source_uri)) then
            call status%set_error(3001, 'Legacy backend URI must be SPS_HOME|isoc_type|spec_type|dust_type.')
            return
        end if

        self%source_uri = trim(source_uri)
        self%sps_home = trim(source_uri(1:sep_pos_1-1))
        self%isoc_type = trim(to_lower(trim(source_uri(sep_pos_1+1:sep_pos_2-1))))
        self%spec_type = trim(to_lower(trim(source_uri(sep_pos_2+1:sep_pos_3-1))))
        self%dust_type = trim(source_uri(sep_pos_3+1:len_trim(source_uri)))
        if (trim(to_lower(self%dust_type)) == 'themis') then
            self%dust_type = 'THEMIS'
        else
            self%dust_type = 'DL07'
        end if

        if (len_trim(self%sps_home) == 0 .or. len_trim(self%isoc_type) == 0 .or. &
            len_trim(self%spec_type) == 0 .or. len_trim(self%dust_type) == 0) then
            call status%set_error(3002, 'Legacy backend URI contains empty SPS_HOME, isoc_type, spec_type, or dust_type.')
            return
        end if

        call load_legacy_iso_zlegend(self, status)
        if (status%code /= 0) return

        call load_legacy_isochrones(self, status)
        if (status%code /= 0) return

        call load_legacy_spec_lambda(self, status)
        if (status%code /= 0) return

        call load_legacy_wmbasic(self, status)
        if (status%code /= 0) return

        call load_legacy_pagb(self, status)
        if (status%code /= 0) return

        call load_legacy_wr(self, status)
        if (status%code /= 0) return

        call load_legacy_agb(self, status)
        if (status%code /= 0) return

        call load_legacy_dust_emission(self, status)
        if (status%code /= 0) return

        call load_legacy_agn_dust(self, status)
        if (status%code /= 0) return

        call load_legacy_dust_attenuation(self, status)
        if (status%code /= 0) return

        call load_legacy_xrb(self, status)
        if (status%code /= 0) return

        call load_legacy_nebular_grid(self, .true., self%neb_cache_wd, status)
        if (status%code /= 0) return

        call load_legacy_nebular_grid(self, .false., self%neb_cache_nd, status)
        if (status%code /= 0) return

        self%is_opened = .true.
    end subroutine legacy_backend_open

    subroutine legacy_backend_close(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        call status%set_ok()
        call clear_open_state(self)
    end subroutine legacy_backend_close

    logical function legacy_backend_is_open(self)
        class(legacy_backend_t), intent(in) :: self
        legacy_backend_is_open = self%is_opened
    end function legacy_backend_is_open

    subroutine legacy_backend_read_manifest(self, manifest, status)
        class(legacy_backend_t), intent(inout) :: self
        type(library_manifest_t), intent(inout) :: manifest
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(3010, 'read_manifest called on closed legacy backend.')
            return
        end if

        call manifest%clear()

        manifest%fsds_version = 'legacy-adapter'
        manifest%producer = 'fsps-legacy'
        manifest%isoc_name = self%isoc_type
        manifest%spec_name = self%spec_type
        manifest%dust_name = self%dust_type

        allocate(manifest%datasets(61))
        call manifest%datasets(1)%clear()
        manifest%datasets(1)%role = 'spectral_base'
        manifest%datasets(1)%path = 'dummy'
        manifest%datasets(1)%dims_csv = 'lambda,z,logt,logg'
        manifest%datasets(1)%unit = 'Lsun/Hz/Msun'
        manifest%datasets(1)%dtype = 'float32'
        manifest%datasets(1)%representation = 'dense_nd'

        call manifest%datasets(2)%clear()
        manifest%datasets(2)%role = 'isoc_nmass'
        manifest%datasets(2)%path = 'legacy:isoc_nmass'
        manifest%datasets(2)%dims_csv = 'nt,nz'
        manifest%datasets(2)%dtype = 'int32'
        manifest%datasets(2)%representation = 'dense_nd'

        call manifest%datasets(3)%clear()
        manifest%datasets(3)%role = 'isoc_timestep'
        manifest%datasets(3)%path = 'legacy:isoc_timestep'
        manifest%datasets(3)%dims_csv = 'nt,nz'
        manifest%datasets(3)%dtype = 'float64'
        manifest%datasets(3)%representation = 'dense_nd'

        call manifest%datasets(4)%clear()
        manifest%datasets(4)%role = 'isoc_mini'
        manifest%datasets(4)%path = 'legacy:isoc_mini'
        manifest%datasets(4)%dims_csv = 'nm,nt,nz'
        manifest%datasets(4)%dtype = 'float64'
        manifest%datasets(4)%representation = 'dense_nd'

        call manifest%datasets(5)%clear()
        manifest%datasets(5)%role = 'isoc_mact'
        manifest%datasets(5)%path = 'legacy:isoc_mact'
        manifest%datasets(5)%dims_csv = 'nm,nt,nz'
        manifest%datasets(5)%dtype = 'float64'
        manifest%datasets(5)%representation = 'dense_nd'

        call manifest%datasets(6)%clear()
        manifest%datasets(6)%role = 'isoc_logl'
        manifest%datasets(6)%path = 'legacy:isoc_logl'
        manifest%datasets(6)%dims_csv = 'nm,nt,nz'
        manifest%datasets(6)%dtype = 'float64'
        manifest%datasets(6)%representation = 'dense_nd'

        call manifest%datasets(7)%clear()
        manifest%datasets(7)%role = 'isoc_logt'
        manifest%datasets(7)%path = 'legacy:isoc_logt'
        manifest%datasets(7)%dims_csv = 'nm,nt,nz'
        manifest%datasets(7)%dtype = 'float64'
        manifest%datasets(7)%representation = 'dense_nd'

        call manifest%datasets(8)%clear()
        manifest%datasets(8)%role = 'isoc_logg'
        manifest%datasets(8)%path = 'legacy:isoc_logg'
        manifest%datasets(8)%dims_csv = 'nm,nt,nz'
        manifest%datasets(8)%dtype = 'float64'
        manifest%datasets(8)%representation = 'dense_nd'

        call manifest%datasets(9)%clear()
        manifest%datasets(9)%role = 'isoc_phase'
        manifest%datasets(9)%path = 'legacy:isoc_phase'
        manifest%datasets(9)%dims_csv = 'nm,nt,nz'
        manifest%datasets(9)%dtype = 'float64'
        manifest%datasets(9)%representation = 'dense_nd'

        call manifest%datasets(10)%clear()
        manifest%datasets(10)%role = 'isoc_ffco'
        manifest%datasets(10)%path = 'legacy:isoc_ffco'
        manifest%datasets(10)%dims_csv = 'nm,nt,nz'
        manifest%datasets(10)%dtype = 'float64'
        manifest%datasets(10)%representation = 'dense_nd'

        call manifest%datasets(11)%clear()
        manifest%datasets(11)%role = 'isoc_lmdot'
        manifest%datasets(11)%path = 'legacy:isoc_lmdot'
        manifest%datasets(11)%dims_csv = 'nm,nt,nz'
        manifest%datasets(11)%dtype = 'float64'
        manifest%datasets(11)%representation = 'dense_nd'

        call manifest%datasets(12)%clear()
        manifest%datasets(12)%role = 'nebular_wd_line_pos'
        manifest%datasets(12)%path = 'legacy:nebular:wd:line_pos'
        manifest%datasets(12)%dims_csv = 'line'
        manifest%datasets(12)%dtype = 'float64'
        manifest%datasets(12)%representation = 'dense_nd'

        call manifest%datasets(13)%clear()
        manifest%datasets(13)%role = 'nebular_wd_logz'
        manifest%datasets(13)%path = 'legacy:nebular:wd:logz'
        manifest%datasets(13)%dims_csv = 'z'
        manifest%datasets(13)%dtype = 'float64'
        manifest%datasets(13)%representation = 'dense_nd'

        call manifest%datasets(14)%clear()
        manifest%datasets(14)%role = 'nebular_wd_age'
        manifest%datasets(14)%path = 'legacy:nebular:wd:age'
        manifest%datasets(14)%dims_csv = 'age'
        manifest%datasets(14)%dtype = 'float64'
        manifest%datasets(14)%representation = 'dense_nd'

        call manifest%datasets(15)%clear()
        manifest%datasets(15)%role = 'nebular_wd_logu'
        manifest%datasets(15)%path = 'legacy:nebular:wd:logu'
        manifest%datasets(15)%dims_csv = 'u'
        manifest%datasets(15)%dtype = 'float64'
        manifest%datasets(15)%representation = 'dense_nd'

        call manifest%datasets(16)%clear()
        manifest%datasets(16)%role = 'nebular_wd_cont'
        manifest%datasets(16)%path = 'legacy:nebular:wd:cont'
        manifest%datasets(16)%dims_csv = 'lam,z,age,u'
        manifest%datasets(16)%dtype = 'float32'
        manifest%datasets(16)%representation = 'dense_nd'

        call manifest%datasets(17)%clear()
        manifest%datasets(17)%role = 'nebular_wd_line'
        manifest%datasets(17)%path = 'legacy:nebular:wd:line'
        manifest%datasets(17)%dims_csv = 'line,z,age,u'
        manifest%datasets(17)%dtype = 'float32'
        manifest%datasets(17)%representation = 'dense_nd'

        call manifest%datasets(18)%clear()
        manifest%datasets(18)%role = 'nebular_nd_line_pos'
        manifest%datasets(18)%path = 'legacy:nebular:nd:line_pos'
        manifest%datasets(18)%dims_csv = 'line'
        manifest%datasets(18)%dtype = 'float64'
        manifest%datasets(18)%representation = 'dense_nd'

        call manifest%datasets(19)%clear()
        manifest%datasets(19)%role = 'nebular_nd_logz'
        manifest%datasets(19)%path = 'legacy:nebular:nd:logz'
        manifest%datasets(19)%dims_csv = 'z'
        manifest%datasets(19)%dtype = 'float64'
        manifest%datasets(19)%representation = 'dense_nd'

        call manifest%datasets(20)%clear()
        manifest%datasets(20)%role = 'nebular_nd_age'
        manifest%datasets(20)%path = 'legacy:nebular:nd:age'
        manifest%datasets(20)%dims_csv = 'age'
        manifest%datasets(20)%dtype = 'float64'
        manifest%datasets(20)%representation = 'dense_nd'

        call manifest%datasets(21)%clear()
        manifest%datasets(21)%role = 'nebular_nd_logu'
        manifest%datasets(21)%path = 'legacy:nebular:nd:logu'
        manifest%datasets(21)%dims_csv = 'u'
        manifest%datasets(21)%dtype = 'float64'
        manifest%datasets(21)%representation = 'dense_nd'

        call manifest%datasets(22)%clear()
        manifest%datasets(22)%role = 'nebular_nd_cont'
        manifest%datasets(22)%path = 'legacy:nebular:nd:cont'
        manifest%datasets(22)%dims_csv = 'lam,z,age,u'
        manifest%datasets(22)%dtype = 'float32'
        manifest%datasets(22)%representation = 'dense_nd'

        call manifest%datasets(23)%clear()
        manifest%datasets(23)%role = 'nebular_nd_line'
        manifest%datasets(23)%path = 'legacy:nebular:nd:line'
        manifest%datasets(23)%dims_csv = 'line,z,age,u'
        manifest%datasets(23)%dtype = 'float32'
        manifest%datasets(23)%representation = 'dense_nd'

        call manifest%datasets(24)%clear()
        manifest%datasets(24)%role = 'wmb_logt'
        manifest%datasets(24)%path = 'legacy:wmbasic:logt'
        manifest%datasets(24)%dims_csv = 'logt'
        manifest%datasets(24)%dtype = 'float64'
        manifest%datasets(24)%representation = 'dense_nd'

        call manifest%datasets(25)%clear()
        manifest%datasets(25)%role = 'wmb_z'
        manifest%datasets(25)%path = 'legacy:wmbasic:z'
        manifest%datasets(25)%dims_csv = 'z'
        manifest%datasets(25)%dtype = 'float64'
        manifest%datasets(25)%representation = 'dense_nd'

        call manifest%datasets(26)%clear()
        manifest%datasets(26)%role = 'wmb_lam'
        manifest%datasets(26)%path = 'legacy:wmbasic:lam'
        manifest%datasets(26)%dims_csv = 'lam'
        manifest%datasets(26)%dtype = 'float64'
        manifest%datasets(26)%representation = 'dense_nd'

        call manifest%datasets(27)%clear()
        manifest%datasets(27)%role = 'wmb_spec'
        manifest%datasets(27)%path = 'legacy:wmbasic:spec'
        manifest%datasets(27)%dims_csv = 'lam,logt,logg,z'
        manifest%datasets(27)%dtype = 'float32'
        manifest%datasets(27)%representation = 'dense_nd'

        call manifest%datasets(28)%clear()
        manifest%datasets(28)%role = 'pagb_logt'
        manifest%datasets(28)%path = 'legacy:pagb:logt'
        manifest%datasets(28)%dims_csv = 'logt'
        manifest%datasets(28)%dtype = 'float64'
        manifest%datasets(28)%representation = 'dense_nd'

        call manifest%datasets(29)%clear()
        manifest%datasets(29)%role = 'pagb_lam'
        manifest%datasets(29)%path = 'legacy:pagb:lam'
        manifest%datasets(29)%dims_csv = 'lam'
        manifest%datasets(29)%dtype = 'float64'
        manifest%datasets(29)%representation = 'dense_nd'

        call manifest%datasets(30)%clear()
        manifest%datasets(30)%role = 'pagb_spec'
        manifest%datasets(30)%path = 'legacy:pagb:spec'
        manifest%datasets(30)%dims_csv = 'lam,logt,z'
        manifest%datasets(30)%dtype = 'float32'
        manifest%datasets(30)%representation = 'dense_nd'

        call manifest%datasets(31)%clear()
        manifest%datasets(31)%role = 'wr_logt_wn'
        manifest%datasets(31)%path = 'legacy:wr:logt_wn'
        manifest%datasets(31)%dims_csv = 'logt'
        manifest%datasets(31)%dtype = 'float64'
        manifest%datasets(31)%representation = 'dense_nd'

        call manifest%datasets(32)%clear()
        manifest%datasets(32)%role = 'wr_logt_wc'
        manifest%datasets(32)%path = 'legacy:wr:logt_wc'
        manifest%datasets(32)%dims_csv = 'logt'
        manifest%datasets(32)%dtype = 'float64'
        manifest%datasets(32)%representation = 'dense_nd'

        call manifest%datasets(33)%clear()
        manifest%datasets(33)%role = 'wr_z'
        manifest%datasets(33)%path = 'legacy:wr:z'
        manifest%datasets(33)%dims_csv = 'z'
        manifest%datasets(33)%dtype = 'float64'
        manifest%datasets(33)%representation = 'dense_nd'

        call manifest%datasets(34)%clear()
        manifest%datasets(34)%role = 'wr_lam'
        manifest%datasets(34)%path = 'legacy:wr:lam'
        manifest%datasets(34)%dims_csv = 'lam'
        manifest%datasets(34)%dtype = 'float64'
        manifest%datasets(34)%representation = 'dense_nd'

        call manifest%datasets(35)%clear()
        manifest%datasets(35)%role = 'wr_spec_wn'
        manifest%datasets(35)%path = 'legacy:wr:spec_wn'
        manifest%datasets(35)%dims_csv = 'lam,logt,z'
        manifest%datasets(35)%dtype = 'float32'
        manifest%datasets(35)%representation = 'dense_nd'

        call manifest%datasets(36)%clear()
        manifest%datasets(36)%role = 'wr_spec_wc'
        manifest%datasets(36)%path = 'legacy:wr:spec_wc'
        manifest%datasets(36)%dims_csv = 'lam,logt,z'
        manifest%datasets(36)%dtype = 'float32'
        manifest%datasets(36)%representation = 'dense_nd'

        call manifest%datasets(37)%clear()
        manifest%datasets(37)%role = 'agb_z_o'
        manifest%datasets(37)%path = 'legacy:agb:z_o'
        manifest%datasets(37)%dims_csv = 'z'
        manifest%datasets(37)%dtype = 'float64'
        manifest%datasets(37)%representation = 'dense_nd'

        call manifest%datasets(38)%clear()
        manifest%datasets(38)%role = 'agb_logt_c'
        manifest%datasets(38)%path = 'legacy:agb:logt_c'
        manifest%datasets(38)%dims_csv = 'logt'
        manifest%datasets(38)%dtype = 'float64'
        manifest%datasets(38)%representation = 'dense_nd'

        call manifest%datasets(39)%clear()
        manifest%datasets(39)%role = 'agb_logt_car'
        manifest%datasets(39)%path = 'legacy:agb:logt_car'
        manifest%datasets(39)%dims_csv = 'logt'
        manifest%datasets(39)%dtype = 'float64'
        manifest%datasets(39)%representation = 'dense_nd'

        call manifest%datasets(40)%clear()
        manifest%datasets(40)%role = 'agb_lam_o'
        manifest%datasets(40)%path = 'legacy:agb:lam_o'
        manifest%datasets(40)%dims_csv = 'lam'
        manifest%datasets(40)%dtype = 'float64'
        manifest%datasets(40)%representation = 'dense_nd'

        call manifest%datasets(41)%clear()
        manifest%datasets(41)%role = 'agb_lam_c'
        manifest%datasets(41)%path = 'legacy:agb:lam_c'
        manifest%datasets(41)%dims_csv = 'lam'
        manifest%datasets(41)%dtype = 'float64'
        manifest%datasets(41)%representation = 'dense_nd'

        call manifest%datasets(42)%clear()
        manifest%datasets(42)%role = 'agb_lam_car'
        manifest%datasets(42)%path = 'legacy:agb:lam_car'
        manifest%datasets(42)%dims_csv = 'lam'
        manifest%datasets(42)%dtype = 'float64'
        manifest%datasets(42)%representation = 'dense_nd'

        call manifest%datasets(43)%clear()
        manifest%datasets(43)%role = 'agb_logt_o'
        manifest%datasets(43)%path = 'legacy:agb:logt_o'
        manifest%datasets(43)%dims_csv = 'z,logt'
        manifest%datasets(43)%dtype = 'float64'
        manifest%datasets(43)%representation = 'dense_nd'

        call manifest%datasets(44)%clear()
        manifest%datasets(44)%role = 'agb_spec_o'
        manifest%datasets(44)%path = 'legacy:agb:spec_o'
        manifest%datasets(44)%dims_csv = 'lam,logt'
        manifest%datasets(44)%dtype = 'float32'
        manifest%datasets(44)%representation = 'dense_nd'

        call manifest%datasets(45)%clear()
        manifest%datasets(45)%role = 'agb_spec_c'
        manifest%datasets(45)%path = 'legacy:agb:spec_c'
        manifest%datasets(45)%dims_csv = 'lam,logt'
        manifest%datasets(45)%dtype = 'float32'
        manifest%datasets(45)%representation = 'dense_nd'

        call manifest%datasets(46)%clear()
        manifest%datasets(46)%role = 'agb_spec_car'
        manifest%datasets(46)%path = 'legacy:agb:spec_car'
        manifest%datasets(46)%dims_csv = 'lam,logt'
        manifest%datasets(46)%dtype = 'float32'
        manifest%datasets(46)%representation = 'dense_nd'

        call manifest%datasets(47)%clear()
        manifest%datasets(47)%role = 'dust_em_qpah'
        manifest%datasets(47)%path = 'legacy:dust:em:qpah'
        manifest%datasets(47)%dims_csv = 'qpah'
        manifest%datasets(47)%dtype = 'float64'
        manifest%datasets(47)%representation = 'dense_nd'

        call manifest%datasets(48)%clear()
        manifest%datasets(48)%role = 'dust_em_umin'
        manifest%datasets(48)%path = 'legacy:dust:em:umin'
        manifest%datasets(48)%dims_csv = 'umin'
        manifest%datasets(48)%dtype = 'float64'
        manifest%datasets(48)%representation = 'dense_nd'

        call manifest%datasets(49)%clear()
        manifest%datasets(49)%role = 'dust_em_lam'
        manifest%datasets(49)%path = 'legacy:dust:em:lam'
        manifest%datasets(49)%dims_csv = 'lam'
        manifest%datasets(49)%dtype = 'float64'
        manifest%datasets(49)%representation = 'dense_nd'

        call manifest%datasets(50)%clear()
        manifest%datasets(50)%role = 'dust_em_spec'
        manifest%datasets(50)%path = 'legacy:dust:em:spec'
        manifest%datasets(50)%dims_csv = 'lam,qpah,umin'
        manifest%datasets(50)%dtype = 'float32'
        manifest%datasets(50)%representation = 'dense_nd'

        call manifest%datasets(51)%clear()
        manifest%datasets(51)%role = 'agn_dust_tau'
        manifest%datasets(51)%path = 'legacy:dust:agn:tau'
        manifest%datasets(51)%dims_csv = 'tau'
        manifest%datasets(51)%dtype = 'float64'
        manifest%datasets(51)%representation = 'dense_nd'

        call manifest%datasets(52)%clear()
        manifest%datasets(52)%role = 'agn_dust_lam'
        manifest%datasets(52)%path = 'legacy:dust:agn:lam'
        manifest%datasets(52)%dims_csv = 'lam'
        manifest%datasets(52)%dtype = 'float64'
        manifest%datasets(52)%representation = 'dense_nd'

        call manifest%datasets(53)%clear()
        manifest%datasets(53)%role = 'agn_dust_spec'
        manifest%datasets(53)%path = 'legacy:dust:agn:spec'
        manifest%datasets(53)%dims_csv = 'lam,tau'
        manifest%datasets(53)%dtype = 'float32'
        manifest%datasets(53)%representation = 'dense_nd'

        call manifest%datasets(54)%clear()
        manifest%datasets(54)%role = 'dust_att_wg_lam'
        manifest%datasets(54)%path = 'legacy:dust:att:wg_lam'
        manifest%datasets(54)%dims_csv = 'lam'
        manifest%datasets(54)%dtype = 'float64'
        manifest%datasets(54)%representation = 'dense_nd'

        call manifest%datasets(55)%clear()
        manifest%datasets(55)%role = 'dust_att_wg_spec'
        manifest%datasets(55)%path = 'legacy:dust:att:wg_spec'
        manifest%datasets(55)%dims_csv = 'lam,geom,tau,type'
        manifest%datasets(55)%dtype = 'float32'
        manifest%datasets(55)%representation = 'dense_nd'

        call manifest%datasets(56)%clear()
        manifest%datasets(56)%role = 'dust_att_smc_lam'
        manifest%datasets(56)%path = 'legacy:dust:att:smc_lam'
        manifest%datasets(56)%dims_csv = 'lam'
        manifest%datasets(56)%dtype = 'float64'
        manifest%datasets(56)%representation = 'dense_nd'

        call manifest%datasets(57)%clear()
        manifest%datasets(57)%role = 'dust_att_smc_ext'
        manifest%datasets(57)%path = 'legacy:dust:att:smc_ext'
        manifest%datasets(57)%dims_csv = 'lam'
        manifest%datasets(57)%dtype = 'float64'
        manifest%datasets(57)%representation = 'dense_nd'

        call manifest%datasets(58)%clear()
        manifest%datasets(58)%role = 'xrb_lam'
        manifest%datasets(58)%path = 'legacy:xrb:lam'
        manifest%datasets(58)%dims_csv = 'lam'
        manifest%datasets(58)%dtype = 'float64'
        manifest%datasets(58)%representation = 'dense_nd'

        call manifest%datasets(59)%clear()
        manifest%datasets(59)%role = 'xrb_age'
        manifest%datasets(59)%path = 'legacy:xrb:age'
        manifest%datasets(59)%dims_csv = 'age'
        manifest%datasets(59)%dtype = 'float64'
        manifest%datasets(59)%representation = 'dense_nd'

        call manifest%datasets(60)%clear()
        manifest%datasets(60)%role = 'xrb_z'
        manifest%datasets(60)%path = 'legacy:xrb:z'
        manifest%datasets(60)%dims_csv = 'z'
        manifest%datasets(60)%dtype = 'float64'
        manifest%datasets(60)%representation = 'dense_nd'

        call manifest%datasets(61)%clear()
        manifest%datasets(61)%role = 'xrb_spec'
        manifest%datasets(61)%path = 'legacy:xrb:spec'
        manifest%datasets(61)%dims_csv = 'lam,age,z'
        manifest%datasets(61)%dtype = 'float32'
        manifest%datasets(61)%representation = 'dense_nd'
    end subroutine legacy_backend_read_manifest

    logical function legacy_backend_has_path(self, path)
        class(legacy_backend_t), intent(in) :: self
        character(len=*), intent(in) :: path

        if (.not. self%is_opened) then
            legacy_backend_has_path = .false.
            return
        end if

        select case (trim(path))
          case ('dummy', 'legacy:isoc_nmass', 'legacy:isoc_timestep', 'legacy:isoc_mini', 'legacy:isoc_mact', &
              'legacy:isoc_logl', 'legacy:isoc_logt', 'legacy:isoc_logg', 'legacy:isoc_phase', 'legacy:isoc_ffco', &
              'legacy:isoc_lmdot', 'legacy:nebular:wd:line_pos', 'legacy:nebular:wd:logz', 'legacy:nebular:wd:age', &
              'legacy:nebular:wd:logu', 'legacy:nebular:wd:cont', 'legacy:nebular:wd:line', &
              'legacy:nebular:nd:line_pos', 'legacy:nebular:nd:logz', 'legacy:nebular:nd:age', &
                            'legacy:nebular:nd:logu', 'legacy:nebular:nd:cont', 'legacy:nebular:nd:line', &
                            'legacy:wmbasic:logt', 'legacy:wmbasic:z', 'legacy:wmbasic:lam', 'legacy:wmbasic:spec', &
                            'legacy:pagb:logt', 'legacy:pagb:lam', 'legacy:pagb:spec', &
                            'legacy:wr:logt_wn', 'legacy:wr:logt_wc', 'legacy:wr:z', 'legacy:wr:lam', &
                            'legacy:wr:spec_wn', 'legacy:wr:spec_wc', &
                            'legacy:agb:z_o', 'legacy:agb:logt_c', 'legacy:agb:logt_car', &
                            'legacy:agb:lam_o', 'legacy:agb:lam_c', 'legacy:agb:lam_car', &
                            'legacy:agb:logt_o', 'legacy:agb:spec_o', 'legacy:agb:spec_c', 'legacy:agb:spec_car', &
                            'legacy:dust:em:qpah', 'legacy:dust:em:umin', 'legacy:dust:em:lam', 'legacy:dust:em:spec', &
                            'legacy:dust:agn:tau', 'legacy:dust:agn:lam', 'legacy:dust:agn:spec', &
                            'legacy:dust:att:wg_lam', 'legacy:dust:att:wg_spec', &
                            'legacy:dust:att:smc_lam', 'legacy:dust:att:smc_ext', &
                            'legacy:xrb:lam', 'legacy:xrb:age', 'legacy:xrb:z', 'legacy:xrb:spec')
            legacy_backend_has_path = .true.
        case default
            legacy_backend_has_path = .false.
        end select
    end function legacy_backend_has_path

    subroutine legacy_backend_query_axis(self, axis_name, axis_desc, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: axis_name
        type(axis_desc_t), intent(inout) :: axis_desc
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: axis
        character(len=:), allocatable :: path
        integer :: i
        real(WP), allocatable :: bpass_time(:), bpass_mass(:,:)

        call axis_desc%clear()
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(3020, 'query_axis called on closed legacy backend.')
            return
        end if

        axis = trim(to_lower(trim(axis_name)))
        axis_desc%name = axis

        select case (axis)
        case ('lambda')
            call build_lambda_file_path(self, path, status)
            if (status%code /= 0) return
            call read_real_column_file(path, axis_desc%values, status)
            if (status%code /= 0) return
            axis_desc%path = path
            axis_desc%unit = 'Angstrom'
            axis_desc%n = size(axis_desc%values)
            self%nspec = axis_desc%n

        case ('z')
            call build_zlegend_file_path(self, path, status)
            if (status%code /= 0) return
            call read_real_column_file(path, axis_desc%values, status)
            if (status%code /= 0) return
            axis_desc%path = path
            axis_desc%unit = 'Z'
            axis_desc%n = size(axis_desc%values)
            if (allocated(self%zlegend)) deallocate(self%zlegend)
            allocate(self%zlegend(axis_desc%n))
            self%zlegend = axis_desc%values

        case ('logt')
            axis_desc%unit = 'dex'
            if (trim(self%spec_type) == 'bpass') then
                call read_bpass_mass_table(self, bpass_time, bpass_mass, status)
                if (status%code /= 0) return
                axis_desc%n = size(bpass_time)
                axis_desc%path = trim(self%sps_home)//'/data/isochrones/BPASS/bpass.mass'
                allocate(axis_desc%values(axis_desc%n))
                axis_desc%values = bpass_time
                if (allocated(bpass_time)) deallocate(bpass_time)
                if (allocated(bpass_mass)) deallocate(bpass_mass)
            else
                axis_desc%n = NDIM_LOGT
                axis_desc%path = 'legacy:dummy:logt'
                allocate(axis_desc%values(axis_desc%n))
                do i = 1, axis_desc%n
                    axis_desc%values(i) = real(i, WP)
                end do
            end if

        case ('logg')
            axis_desc%unit = 'dex'
            if (trim(self%spec_type) == 'bpass') then
                axis_desc%n = 1
                axis_desc%path = 'legacy:dummy:bpass_logg'
                allocate(axis_desc%values(1))
                axis_desc%values(1) = 0.0_wp
            else
                axis_desc%n = NDIM_LOGG
                axis_desc%path = 'legacy:dummy:logg'
                allocate(axis_desc%values(axis_desc%n))
                do i = 1, axis_desc%n
                    axis_desc%values(i) = real(i, WP)
                end do
            end if

        case ('afe', 'alpha_fe')
            call status%set_error(3021, 'Legacy backend does not support an afe axis.')

        case default
            call status%set_error(3022, 'Unknown axis requested from legacy backend: '//trim(axis_name))
        end select
    end subroutine legacy_backend_query_axis

    subroutine legacy_backend_read_real_1d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        select case (trim(dataset_path))
        case ('legacy:nebular:wd:line_pos')
            if (.not. allocated(self%neb_cache_wd%line_pos)) then
                call status%set_error(3030, 'WD nebular line_pos cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%line_pos)))
            values = self%neb_cache_wd%line_pos
        case ('legacy:nebular:wd:logz')
            if (.not. allocated(self%neb_cache_wd%logz)) then
                call status%set_error(3030, 'WD nebular logz cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%logz)))
            values = self%neb_cache_wd%logz
        case ('legacy:nebular:wd:age')
            if (.not. allocated(self%neb_cache_wd%age)) then
                call status%set_error(3030, 'WD nebular age cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%age)))
            values = self%neb_cache_wd%age
        case ('legacy:nebular:wd:logu')
            if (.not. allocated(self%neb_cache_wd%logu)) then
                call status%set_error(3030, 'WD nebular logu cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%logu)))
            values = self%neb_cache_wd%logu
        case ('legacy:nebular:nd:line_pos')
            if (.not. allocated(self%neb_cache_nd%line_pos)) then
                call status%set_error(3030, 'ND nebular line_pos cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%line_pos)))
            values = self%neb_cache_nd%line_pos
        case ('legacy:nebular:nd:logz')
            if (.not. allocated(self%neb_cache_nd%logz)) then
                call status%set_error(3030, 'ND nebular logz cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%logz)))
            values = self%neb_cache_nd%logz
        case ('legacy:nebular:nd:age')
            if (.not. allocated(self%neb_cache_nd%age)) then
                call status%set_error(3030, 'ND nebular age cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%age)))
            values = self%neb_cache_nd%age
        case ('legacy:nebular:nd:logu')
            if (.not. allocated(self%neb_cache_nd%logu)) then
                call status%set_error(3030, 'ND nebular logu cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%logu)))
            values = self%neb_cache_nd%logu
        case ('legacy:wmbasic:logt')
            if (.not. allocated(self%wmb_cache%logt)) then
                call status%set_error(3030, 'WMBasic logt cache is not allocated.')
                return
            end if
            allocate(values(size(self%wmb_cache%logt)))
            values = self%wmb_cache%logt
        case ('legacy:wmbasic:z')
            if (.not. allocated(self%wmb_cache%z)) then
                call status%set_error(3030, 'WMBasic z cache is not allocated.')
                return
            end if
            allocate(values(size(self%wmb_cache%z)))
            values = self%wmb_cache%z
        case ('legacy:wmbasic:lam')
            if (.not. allocated(self%wmb_cache%lam)) then
                call status%set_error(3030, 'WMBasic lam cache is not allocated.')
                return
            end if
            allocate(values(size(self%wmb_cache%lam)))
            values = self%wmb_cache%lam
        case ('legacy:pagb:logt')
            if (.not. allocated(self%pagb_cache%logt)) then
                call status%set_error(3030, 'Post-AGB logt cache is not allocated.')
                return
            end if
            allocate(values(size(self%pagb_cache%logt)))
            values = self%pagb_cache%logt
        case ('legacy:pagb:lam')
            if (.not. allocated(self%pagb_cache%lam)) then
                call status%set_error(3030, 'Post-AGB lam cache is not allocated.')
                return
            end if
            allocate(values(size(self%pagb_cache%lam)))
            values = self%pagb_cache%lam
        case ('legacy:wr:logt_wn')
            if (.not. allocated(self%wr_cache%logt_wn)) then
                call status%set_error(3030, 'WR WN logt cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%logt_wn)))
            values = self%wr_cache%logt_wn
        case ('legacy:wr:logt_wc')
            if (.not. allocated(self%wr_cache%logt_wc)) then
                call status%set_error(3030, 'WR WC logt cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%logt_wc)))
            values = self%wr_cache%logt_wc
        case ('legacy:wr:z')
            if (.not. allocated(self%wr_cache%z)) then
                call status%set_error(3030, 'WR z cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%z)))
            values = self%wr_cache%z
        case ('legacy:wr:lam')
            if (.not. allocated(self%wr_cache%lam)) then
                call status%set_error(3030, 'WR lam cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%lam)))
            values = self%wr_cache%lam
        case ('legacy:agb:z_o')
            if (.not. allocated(self%agb_cache%z_o)) then
                call status%set_error(3030, 'AGB z_o cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%z_o)))
            values = self%agb_cache%z_o
        case ('legacy:agb:logt_c')
            if (.not. allocated(self%agb_cache%logt_c)) then
                call status%set_error(3030, 'AGB logt_c cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%logt_c)))
            values = self%agb_cache%logt_c
        case ('legacy:agb:logt_car')
            if (.not. allocated(self%agb_cache%logt_car)) then
                call status%set_error(3030, 'AGB logt_car cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%logt_car)))
            values = self%agb_cache%logt_car
        case ('legacy:agb:lam_o')
            if (.not. allocated(self%agb_cache%lam_o)) then
                call status%set_error(3030, 'AGB lam_o cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%lam_o)))
            values = self%agb_cache%lam_o
        case ('legacy:agb:lam_c')
            if (.not. allocated(self%agb_cache%lam_c)) then
                call status%set_error(3030, 'AGB lam_c cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%lam_c)))
            values = self%agb_cache%lam_c
        case ('legacy:agb:lam_car')
            if (.not. allocated(self%agb_cache%lam_car)) then
                call status%set_error(3030, 'AGB lam_car cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%lam_car)))
            values = self%agb_cache%lam_car
        case ('legacy:dust:em:qpah')
            if (.not. allocated(self%dust_em_cache%qpah)) then
                call status%set_error(3030, 'Dust emission qpah cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_em_cache%qpah)))
            values = self%dust_em_cache%qpah
        case ('legacy:dust:em:umin')
            if (.not. allocated(self%dust_em_cache%umin)) then
                call status%set_error(3030, 'Dust emission umin cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_em_cache%umin)))
            values = self%dust_em_cache%umin
        case ('legacy:dust:em:lam')
            if (.not. allocated(self%dust_em_cache%lam)) then
                call status%set_error(3030, 'Dust emission lam cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_em_cache%lam)))
            values = self%dust_em_cache%lam
        case ('legacy:dust:agn:tau')
            if (.not. allocated(self%agn_dust_cache%tau)) then
                call status%set_error(3030, 'AGN dust tau cache is not allocated.')
                return
            end if
            allocate(values(size(self%agn_dust_cache%tau)))
            values = self%agn_dust_cache%tau
        case ('legacy:dust:agn:lam')
            if (.not. allocated(self%agn_dust_cache%lam)) then
                call status%set_error(3030, 'AGN dust lam cache is not allocated.')
                return
            end if
            allocate(values(size(self%agn_dust_cache%lam)))
            values = self%agn_dust_cache%lam
        case ('legacy:dust:att:wg_lam')
            if (.not. allocated(self%dust_att_cache%wg_lam)) then
                call status%set_error(3030, 'Dust attenuation WG lambda cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_att_cache%wg_lam)))
            values = self%dust_att_cache%wg_lam
        case ('legacy:dust:att:smc_lam')
            if (.not. allocated(self%dust_att_cache%smc_lam)) then
                call status%set_error(3030, 'Dust attenuation SMC lambda cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_att_cache%smc_lam)))
            values = self%dust_att_cache%smc_lam
        case ('legacy:dust:att:smc_ext')
            if (.not. allocated(self%dust_att_cache%smc_ext)) then
                call status%set_error(3030, 'Dust attenuation SMC extinction cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_att_cache%smc_ext)))
            values = self%dust_att_cache%smc_ext
        case ('legacy:xrb:lam')
            if (.not. allocated(self%xrb_cache%lam)) then
                call status%set_error(3030, 'XRB lambda cache is not allocated.')
                return
            end if
            allocate(values(size(self%xrb_cache%lam)))
            values = self%xrb_cache%lam
        case ('legacy:xrb:age')
            if (.not. allocated(self%xrb_cache%age)) then
                call status%set_error(3030, 'XRB age cache is not allocated.')
                return
            end if
            allocate(values(size(self%xrb_cache%age)))
            values = self%xrb_cache%age
        case ('legacy:xrb:z')
            if (.not. allocated(self%xrb_cache%z)) then
                call status%set_error(3030, 'XRB metallicity cache is not allocated.')
                return
            end if
            allocate(values(size(self%xrb_cache%z)))
            values = self%xrb_cache%z
        case default
            call status%set_error(3030, 'Legacy backend read_real_1d path not implemented: '//trim(dataset_path))
        end select
    end subroutine legacy_backend_read_real_1d

    subroutine legacy_backend_read_real_2d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        select case (trim(dataset_path))
        case ('legacy:isoc_timestep')
            if (.not. allocated(self%iso_cache%timestep_logyr)) then
                call status%set_error(3031, 'Isochrone timestep cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%timestep_logyr, 1), size(self%iso_cache%timestep_logyr, 2)))
            values = self%iso_cache%timestep_logyr
        case ('legacy:agb:logt_o')
            if (.not. allocated(self%agb_cache%logt_o)) then
                call status%set_error(3031, 'AGB logt_o cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%logt_o, 1), size(self%agb_cache%logt_o, 2)))
            values = self%agb_cache%logt_o
        case ('legacy:agb:spec_o')
            if (.not. allocated(self%agb_cache%spec_o)) then
                call status%set_error(3031, 'AGB spec_o cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%spec_o, 1), size(self%agb_cache%spec_o, 2)))
            values = self%agb_cache%spec_o
        case ('legacy:agb:spec_c')
            if (.not. allocated(self%agb_cache%spec_c)) then
                call status%set_error(3031, 'AGB spec_c cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%spec_c, 1), size(self%agb_cache%spec_c, 2)))
            values = self%agb_cache%spec_c
        case ('legacy:agb:spec_car')
            if (.not. allocated(self%agb_cache%spec_car)) then
                call status%set_error(3031, 'AGB spec_car cache is not allocated.')
                return
            end if
            allocate(values(size(self%agb_cache%spec_car, 1), size(self%agb_cache%spec_car, 2)))
            values = self%agb_cache%spec_car
        case ('legacy:dust:agn:spec')
            if (.not. allocated(self%agn_dust_cache%spec)) then
                call status%set_error(3031, 'AGN dust spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%agn_dust_cache%spec, 1), size(self%agn_dust_cache%spec, 2)))
            values = self%agn_dust_cache%spec
        case default
            call status%set_error(3031, 'Legacy backend read_real_2d path not implemented: '//trim(dataset_path))
        end select
    end subroutine legacy_backend_read_real_2d

    subroutine legacy_backend_read_real_3d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:,:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        select case (trim(dataset_path))
        case ('legacy:isoc_mini')
            if (.not. allocated(self%iso_cache%mini)) then
                call status%set_error(3034, 'Isochrone mini cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%mini,1), size(self%iso_cache%mini,2), size(self%iso_cache%mini,3)))
            values = self%iso_cache%mini
        case ('legacy:isoc_mact')
            if (.not. allocated(self%iso_cache%mact)) then
                call status%set_error(3034, 'Isochrone mact cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%mact,1), size(self%iso_cache%mact,2), size(self%iso_cache%mact,3)))
            values = self%iso_cache%mact
        case ('legacy:isoc_logl')
            if (.not. allocated(self%iso_cache%logl)) then
                call status%set_error(3034, 'Isochrone logl cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%logl,1), size(self%iso_cache%logl,2), size(self%iso_cache%logl,3)))
            values = self%iso_cache%logl
        case ('legacy:isoc_logt')
            if (.not. allocated(self%iso_cache%logt)) then
                call status%set_error(3034, 'Isochrone logt cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%logt,1), size(self%iso_cache%logt,2), size(self%iso_cache%logt,3)))
            values = self%iso_cache%logt
        case ('legacy:isoc_logg')
            if (.not. allocated(self%iso_cache%logg)) then
                call status%set_error(3034, 'Isochrone logg cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%logg,1), size(self%iso_cache%logg,2), size(self%iso_cache%logg,3)))
            values = self%iso_cache%logg
        case ('legacy:isoc_phase')
            if (.not. allocated(self%iso_cache%phase)) then
                call status%set_error(3034, 'Isochrone phase cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%phase,1), size(self%iso_cache%phase,2), size(self%iso_cache%phase,3)))
            values = self%iso_cache%phase
        case ('legacy:isoc_ffco')
            if (.not. allocated(self%iso_cache%ffco)) then
                call status%set_error(3034, 'Isochrone ffco cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%ffco,1), size(self%iso_cache%ffco,2), size(self%iso_cache%ffco,3)))
            values = self%iso_cache%ffco
        case ('legacy:isoc_lmdot')
            if (.not. allocated(self%iso_cache%lmdot)) then
                call status%set_error(3034, 'Isochrone lmdot cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%lmdot,1), size(self%iso_cache%lmdot,2), size(self%iso_cache%lmdot,3)))
            values = self%iso_cache%lmdot
        case ('legacy:pagb:spec')
            if (.not. allocated(self%pagb_cache%spec)) then
                call status%set_error(3034, 'Post-AGB spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%pagb_cache%spec,1), size(self%pagb_cache%spec,2), size(self%pagb_cache%spec,3)))
            values = self%pagb_cache%spec
        case ('legacy:wr:spec_wn')
            if (.not. allocated(self%wr_cache%spec_wn)) then
                call status%set_error(3034, 'WR WN spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%spec_wn,1), size(self%wr_cache%spec_wn,2), size(self%wr_cache%spec_wn,3)))
            values = self%wr_cache%spec_wn
        case ('legacy:wr:spec_wc')
            if (.not. allocated(self%wr_cache%spec_wc)) then
                call status%set_error(3034, 'WR WC spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%wr_cache%spec_wc,1), size(self%wr_cache%spec_wc,2), size(self%wr_cache%spec_wc,3)))
            values = self%wr_cache%spec_wc
        case ('legacy:dust:em:spec')
            if (.not. allocated(self%dust_em_cache%spec)) then
                call status%set_error(3034, 'Dust emission spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_em_cache%spec,1), size(self%dust_em_cache%spec,2), size(self%dust_em_cache%spec,3)))
            values = self%dust_em_cache%spec
        case ('legacy:xrb:spec')
            if (.not. allocated(self%xrb_cache%spec)) then
                call status%set_error(3034, 'XRB spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%xrb_cache%spec,1), size(self%xrb_cache%spec,2), size(self%xrb_cache%spec,3)))
            values = self%xrb_cache%spec
        case default
            call status%set_error(3034, 'Legacy backend read_real_3d path not implemented: '//trim(dataset_path))
        end select
    end subroutine legacy_backend_read_real_3d

    subroutine legacy_backend_read_real_4d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        real(WP), allocatable, intent(out) :: values(:,:,:,:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        select case (trim(dataset_path))
        case ('legacy:nebular:wd:cont')
            if (.not. allocated(self%neb_cache_wd%cont)) then
                call status%set_error(3036, 'WD nebular continuum cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%cont,1), size(self%neb_cache_wd%cont,2), &
                           size(self%neb_cache_wd%cont,3), size(self%neb_cache_wd%cont,4)))
            values = self%neb_cache_wd%cont
        case ('legacy:nebular:wd:line')
            if (.not. allocated(self%neb_cache_wd%line)) then
                call status%set_error(3036, 'WD nebular line cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_wd%line,1), size(self%neb_cache_wd%line,2), &
                           size(self%neb_cache_wd%line,3), size(self%neb_cache_wd%line,4)))
            values = self%neb_cache_wd%line
        case ('legacy:nebular:nd:cont')
            if (.not. allocated(self%neb_cache_nd%cont)) then
                call status%set_error(3036, 'ND nebular continuum cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%cont,1), size(self%neb_cache_nd%cont,2), &
                           size(self%neb_cache_nd%cont,3), size(self%neb_cache_nd%cont,4)))
            values = self%neb_cache_nd%cont
        case ('legacy:nebular:nd:line')
            if (.not. allocated(self%neb_cache_nd%line)) then
                call status%set_error(3036, 'ND nebular line cache is not allocated.')
                return
            end if
            allocate(values(size(self%neb_cache_nd%line,1), size(self%neb_cache_nd%line,2), &
                           size(self%neb_cache_nd%line,3), size(self%neb_cache_nd%line,4)))
            values = self%neb_cache_nd%line
        case ('legacy:wmbasic:spec')
            if (.not. allocated(self%wmb_cache%spec)) then
                call status%set_error(3036, 'WMBasic spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%wmb_cache%spec,1), size(self%wmb_cache%spec,2), &
                           size(self%wmb_cache%spec,3), size(self%wmb_cache%spec,4)))
            values = self%wmb_cache%spec
        case ('legacy:dust:att:wg_spec')
            if (.not. allocated(self%dust_att_cache%wg_spec)) then
                call status%set_error(3036, 'Dust attenuation WG spec cache is not allocated.')
                return
            end if
            allocate(values(size(self%dust_att_cache%wg_spec,1), size(self%dust_att_cache%wg_spec,2), &
                           size(self%dust_att_cache%wg_spec,3), size(self%dust_att_cache%wg_spec,4)))
            values = self%dust_att_cache%wg_spec
        case default
            call status%set_error(3036, 'Legacy backend read_real_4d path not implemented: '//trim(dataset_path))
        end select
    end subroutine legacy_backend_read_real_4d

    subroutine legacy_backend_read_int_1d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_error(3032, 'Legacy backend read_int_1d is not implemented.')
    end subroutine legacy_backend_read_int_1d

    subroutine legacy_backend_read_int_2d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:,:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        select case (trim(dataset_path))
        case ('legacy:isoc_nmass')
            if (.not. allocated(self%iso_cache%nmass)) then
                call status%set_error(3033, 'Isochrone nmass cache is not allocated.')
                return
            end if
            allocate(values(size(self%iso_cache%nmass, 1), size(self%iso_cache%nmass, 2)))
            values = self%iso_cache%nmass
        case default
            call status%set_error(3033, 'Legacy backend read_int_2d path not implemented: '//trim(dataset_path))
        end select
    end subroutine legacy_backend_read_int_2d

    subroutine legacy_backend_read_int_3d(self, dataset_path, values, status)
        class(legacy_backend_t), intent(inout) :: self
        character(len=*), intent(in) :: dataset_path
        integer, allocatable, intent(out) :: values(:,:,:)
        type(backend_status_t), intent(out) :: status

        if (allocated(values)) deallocate(values)
        call status%set_error(3035, 'Legacy backend read_int_3d is not yet implemented.')
    end subroutine legacy_backend_read_int_3d

    subroutine legacy_backend_read_spectral_slice(self, dataset, iz, iafe, slice, status)
        class(legacy_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        integer, intent(in) :: iz
        integer, intent(in) :: iafe
        type(spectral_slice_t), intent(inout) :: slice
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: file_path
        real(real32), allocatable :: buffer_32(:,:,:)
        real(WP), allocatable :: bpass_cube(:,:,:)
        integer :: u_file, io_stat, rec_len, file_unit_size

        call slice%clear()
        call status%set_ok()

        if (.not. self%is_opened) then
            call status%set_error(3040, 'read_spectral_slice called on closed legacy backend.')
            return
        end if

        if (iafe /= 1) then
            call status%set_error(3041, 'Legacy backend supports only iafe=1.')
            return
        end if

        if (self%nspec <= 0) then
            call status%set_error(3042, 'Legacy backend nspec unknown; query lambda axis before slice reads.')
            return
        end if

        if (.not. allocated(self%zlegend)) then
            call status%set_error(3043, 'Legacy backend zlegend unknown; query z axis before slice reads.')
            return
        end if

        if (iz < 1 .or. iz > size(self%zlegend)) then
            call status%set_error(3044, 'Requested iz index out of range in legacy backend.')
            return
        end if

        call build_binary_file_path(self, iz, file_path, status)
        if (status%code /= 0) return

        if (trim(self%spec_type) == 'bpass') then
            allocate(bpass_cube(self%nspec, NT_BPASS, size(self%zlegend)))

            open(newunit=u_file, file=trim(file_path), status='old', access='stream', form='unformatted', &
                 action='read', iostat=io_stat)
            if (io_stat /= 0) then
                deallocate(bpass_cube)
                call status%set_error(3045, 'Failed to open legacy BPASS spectral binary: '//trim(file_path))
                return
            end if

            read(u_file, iostat=io_stat) bpass_cube
            close(u_file)
            if (io_stat /= 0) then
                deallocate(bpass_cube)
                call status%set_error(3046, 'Failed to read legacy BPASS spectral record: '//trim(file_path))
                return
            end if

            allocate(slice%flux(self%nspec, NT_BPASS, 1))
            slice%flux(:, :, 1) = bpass_cube(:, :, iz)
            deallocate(bpass_cube)
        else
            file_unit_size = file_storage_size / 8
            rec_len = (self%nspec * NDIM_LOGT * NDIM_LOGG * 4) / file_unit_size

            allocate(buffer_32(self%nspec, NDIM_LOGT, NDIM_LOGG))

            open(newunit=u_file, file=trim(file_path), status='old', access='direct', recl=rec_len, &
                 form='unformatted', action='read', iostat=io_stat)
            if (io_stat /= 0) then
                deallocate(buffer_32)
                call status%set_error(3045, 'Failed to open legacy spectral binary: '//trim(file_path))
                return
            end if

            read(u_file, rec=1, iostat=io_stat) buffer_32
            close(u_file)
            if (io_stat /= 0) then
                deallocate(buffer_32)
                call status%set_error(3046, 'Failed to read legacy spectral record: '//trim(file_path))
                return
            end if

            allocate(slice%flux(self%nspec, NDIM_LOGT, NDIM_LOGG))
            slice%flux = real(buffer_32, WP)
            deallocate(buffer_32)
        end if

        slice%iz = iz
        slice%iafe = 1
    end subroutine legacy_backend_read_spectral_slice

    subroutine legacy_backend_read_spectral_neighborhood(self, dataset, iz_lo, iz_hi, iafe_lo, iafe_hi, neighborhood, status)
        class(legacy_backend_t), intent(inout) :: self
        type(dataset_desc_t), intent(in) :: dataset
        integer, intent(in) :: iz_lo
        integer, intent(in) :: iz_hi
        integer, intent(in) :: iafe_lo
        integer, intent(in) :: iafe_hi
        type(spectral_grid_t), intent(inout) :: neighborhood
        type(backend_status_t), intent(out) :: status

        call neighborhood%clear()
        call status%set_error(3050, 'Legacy backend read_spectral_neighborhood is not implemented.')
    end subroutine legacy_backend_read_spectral_neighborhood

    subroutine clear_open_state(self)
        class(legacy_backend_t), intent(inout) :: self

        self%is_opened = .false.
        self%nspec = 0
        call self%iso_cache%clear()
        call self%neb_cache_wd%clear()
        call self%neb_cache_nd%clear()
        call self%wmb_cache%clear()
        call self%pagb_cache%clear()
        call self%wr_cache%clear()
        call self%agb_cache%clear()
        call self%dust_em_cache%clear()
        call self%agn_dust_cache%clear()
        call self%dust_att_cache%clear()
        call self%xrb_cache%clear()
        if (allocated(self%source_uri)) deallocate(self%source_uri)
        if (allocated(self%sps_home)) deallocate(self%sps_home)
        if (allocated(self%isoc_type)) deallocate(self%isoc_type)
        if (allocated(self%spec_type)) deallocate(self%spec_type)
        if (allocated(self%dust_type)) deallocate(self%dust_type)
        if (allocated(self%spec_lambda)) deallocate(self%spec_lambda)
        if (allocated(self%zlegend)) deallocate(self%zlegend)
    end subroutine clear_open_state

    subroutine load_legacy_spec_lambda(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: path

        call status%set_ok()

        call build_lambda_file_path(self, path, status)
        if (status%code /= 0) return

        call read_real_column_file(path, self%spec_lambda, status)
        if (status%code /= 0) return

        self%nspec = size(self%spec_lambda)
    end subroutine load_legacy_spec_lambda

    subroutine load_legacy_wmbasic(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer, parameter :: NZWMB = 12
        integer, parameter :: NSPEC_WMB = 5508

        integer :: u_file, u_spec, io_stat
        integer :: iz, i
        character(len=1024) :: file_path
        character(len=6) :: z_str
        real(WP) :: lam_i
        real(WP) :: g1(NDIM_WMB_LOGT), g2(NDIM_WMB_LOGT), g3(NDIM_WMB_LOGT)

        call status%set_ok()
        call self%wmb_cache%clear()

        allocate(self%wmb_cache%logt(NDIM_WMB_LOGT))
        allocate(self%wmb_cache%z(NZWMB))
        allocate(self%wmb_cache%lam(NSPEC_WMB))
        allocate(self%wmb_cache%spec(NSPEC_WMB, NDIM_WMB_LOGT, NDIM_WMB_LOGG, NZWMB))

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/WMBASIC.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3101, 'Failed to open WMBasic Teff file: '//trim(file_path))
            return
        end if
        do i = 1, NDIM_WMB_LOGT
            read(u_file, *, iostat=io_stat) self%wmb_cache%logt(i)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3102, 'Failed to read WMBasic Teff file: '//trim(file_path))
                return
            end if
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/WMBASIC_zlegend.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3103, 'Failed to open WMBasic zlegend file: '//trim(file_path))
            return
        end if

        do iz = 1, NZWMB
            read(u_file, *, iostat=io_stat) self%wmb_cache%z(iz)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3104, 'Failed to read WMBasic zlegend file: '//trim(file_path))
                return
            end if

            write(z_str, '(F6.4)') self%wmb_cache%z(iz)
            file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/WMBASIC_z'//trim(adjustl(z_str))//'.spec'
            open(newunit=u_spec, file=trim(file_path), status='old', action='read', iostat=io_stat)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3105, 'Failed to open WMBasic spectrum file: '//trim(file_path))
                return
            end if

            do i = 1, NSPEC_WMB
                read(u_spec, *, iostat=io_stat) lam_i, g1, g2, g3
                if (io_stat /= 0) then
                    close(u_spec)
                    close(u_file)
                    call status%set_error(3106, 'Failed to read WMBasic spectrum file: '//trim(file_path))
                    return
                end if

                if (iz == 1) self%wmb_cache%lam(i) = lam_i
                self%wmb_cache%spec(i, :, 1, iz) = g1
                self%wmb_cache%spec(i, :, 2, iz) = g2
                self%wmb_cache%spec(i, :, 3, iz) = g3
            end do

            close(u_spec)
        end do

        close(u_file)
    end subroutine load_legacy_wmbasic

    subroutine load_legacy_pagb(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer, parameter :: NSPEC_PAGB = 9281

        integer :: u_file, io_stat
        integer :: i
        character(len=1024) :: file_path
        real(WP) :: tval

        call status%set_ok()
        call self%pagb_cache%clear()

        allocate(self%pagb_cache%logt(NDIM_PAGB))
        allocate(self%pagb_cache%lam(NSPEC_PAGB))
        allocate(self%pagb_cache%spec(NSPEC_PAGB, NDIM_PAGB, 2))

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/ipagb.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3111, 'Failed to open Post-AGB Teff file: '//trim(file_path))
            return
        end if
        do i = 1, NDIM_PAGB
            read(u_file, *, iostat=io_stat) tval
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3112, 'Failed to read Post-AGB Teff file: '//trim(file_path))
                return
            end if
            self%pagb_cache%logt(i) = log10(tval)
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/ipagb_halo.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3113, 'Failed to open Post-AGB halo file: '//trim(file_path))
            return
        end if
        do i = 1, NSPEC_PAGB
            read(u_file, *, iostat=io_stat) self%pagb_cache%lam(i), self%pagb_cache%spec(i, :, 1)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3114, 'Failed to read Post-AGB halo file: '//trim(file_path))
                return
            end if
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/ipagb_solar.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3115, 'Failed to open Post-AGB solar file: '//trim(file_path))
            return
        end if
        do i = 1, NSPEC_PAGB
            read(u_file, *, iostat=io_stat) tval, self%pagb_cache%spec(i, :, 2)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3116, 'Failed to read Post-AGB solar file: '//trim(file_path))
                return
            end if
        end do
        close(u_file)
    end subroutine load_legacy_pagb

    subroutine load_legacy_wr(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer, parameter :: NZ_WR = 5

        integer :: u_file, io_stat
        integer :: i, j
        character(len=1024) :: file_path
        real(WP) :: dummy_teff, zval

        call status%set_ok()
        call self%wr_cache%clear()

        allocate(self%wr_cache%logt_wn(NDIM_WR))
        allocate(self%wr_cache%logt_wc(NDIM_WR))
        allocate(self%wr_cache%z(NZ_WR))
        allocate(self%wr_cache%lam(NLAMWR))
        allocate(self%wr_cache%spec_wn(NLAMWR, NDIM_WR, NZ_WR))
        allocate(self%wr_cache%spec_wc(NLAMWR, NDIM_WR, NZ_WR))

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/CMFGEN_WN.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3117, 'Failed to open WR WN Teff file: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) self%wr_cache%logt_wn
        close(u_file)
        if (io_stat /= 0) then
            call status%set_error(3118, 'Failed to read WR WN Teff file: '//trim(file_path))
            return
        end if

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/CMFGEN_WC.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3119, 'Failed to open WR WC Teff file: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) self%wr_cache%logt_wc
        close(u_file)
        if (io_stat /= 0) then
            call status%set_error(3120, 'Failed to read WR WC Teff file: '//trim(file_path))
            return
        end if

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/CMFGEN_WN_Zall.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3121, 'Failed to open WR WN spectra file: '//trim(file_path))
            return
        end if

        read(u_file, *, iostat=io_stat) self%wr_cache%lam
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3122, 'Failed to read WR WN wavelength grid: '//trim(file_path))
            return
        end if

        do j = 1, NZ_WR
            do i = 1, NDIM_WR
                read(u_file, *, iostat=io_stat) dummy_teff, zval
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3123, 'Failed to read WR WN metadata row: '//trim(file_path))
                    return
                end if
                if (i == 1) self%wr_cache%z(j) = zval

                read(u_file, *, iostat=io_stat) self%wr_cache%spec_wn(:, i, j)
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3124, 'Failed to read WR WN spectrum row: '//trim(file_path))
                    return
                end if
            end do
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/Hot_spectra/CMFGEN_WC_Zall.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3125, 'Failed to open WR WC spectra file: '//trim(file_path))
            return
        end if

        read(u_file, *, iostat=io_stat) self%wr_cache%lam
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3126, 'Failed to read WR WC wavelength grid: '//trim(file_path))
            return
        end if

        do j = 1, NZ_WR
            do i = 1, NDIM_WR
                read(u_file, *, iostat=io_stat) dummy_teff, zval
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3127, 'Failed to read WR WC metadata row: '//trim(file_path))
                    return
                end if
                if (i == 1) self%wr_cache%z(j) = zval

                read(u_file, *, iostat=io_stat) self%wr_cache%spec_wc(:, i, j)
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3128, 'Failed to read WR WC spectrum row: '//trim(file_path))
                    return
                end if
            end do
        end do
        close(u_file)
    end subroutine load_legacy_wr

    subroutine load_legacy_agb(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat
        integer :: i
        character(len=1024) :: file_path
        character(len=1) :: char_dummy
        real(WP) :: dummy_val
        real(WP) :: tagb_logt_o_raw(22, N_AGB_O)
        real(WP) :: tagb_logz_o(22)
        real(WP) :: temp_row_o(N_AGB_O), temp_row_c(N_AGB_C), temp_row_car(N_AGB_CAR)

        call status%set_ok()
        call self%agb_cache%clear()

        allocate(self%agb_cache%lam_o(NSPEC_AGB))
        allocate(self%agb_cache%lam_c(NSPEC_AGB))
        allocate(self%agb_cache%lam_car(NSPEC_ARINGER))
        allocate(self%agb_cache%z_o(22))
        allocate(self%agb_cache%logt_o(22, N_AGB_O))
        allocate(self%agb_cache%logt_c(N_AGB_C))
        allocate(self%agb_cache%logt_car(N_AGB_CAR))
        allocate(self%agb_cache%spec_o(NSPEC_AGB, N_AGB_O))
        allocate(self%agb_cache%spec_c(NSPEC_AGB, N_AGB_C))
        allocate(self%agb_cache%spec_car(NSPEC_ARINGER, N_AGB_CAR))

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Orich.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3129, 'Failed to open AGB Orich Teff file: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) char_dummy
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3130, 'Failed to read AGB Orich Teff header: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) dummy_val, tagb_logz_o
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3131, 'Failed to read AGB Orich metallicity grid: '//trim(file_path))
            return
        end if
        do i = 1, N_AGB_O
            read(u_file, *, iostat=io_stat) dummy_val, tagb_logt_o_raw(:, i)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3132, 'Failed to read AGB Orich Teff rows: '//trim(file_path))
                return
            end if
        end do
        close(u_file)

        self%agb_cache%z_o = tagb_logz_o
        self%agb_cache%logt_o = log10(tagb_logt_o_raw)

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Crich.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3133, 'Failed to open AGB Crich Teff file: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) char_dummy
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3134, 'Failed to read AGB Crich Teff header: '//trim(file_path))
            return
        end if
        do i = 1, N_AGB_C
            read(u_file, *, iostat=io_stat) dummy_val, self%agb_cache%logt_c(i)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3135, 'Failed to read AGB Crich Teff rows: '//trim(file_path))
                return
            end if
        end do
        close(u_file)
        self%agb_cache%logt_c = log10(self%agb_cache%logt_c)

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Orich.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3136, 'Failed to open AGB Orich spectra file: '//trim(file_path))
            return
        end if
        do i = 1, NSPEC_AGB
            read(u_file, *, iostat=io_stat) self%agb_cache%lam_o(i), temp_row_o
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3137, 'Failed to read AGB Orich spectra file: '//trim(file_path))
                return
            end if
            self%agb_cache%spec_o(i, :) = temp_row_o
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Crich.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3138, 'Failed to open AGB Crich spectra file: '//trim(file_path))
            return
        end if
        do i = 1, NSPEC_AGB
            read(u_file, *, iostat=io_stat) self%agb_cache%lam_c(i), temp_row_c
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3139, 'Failed to read AGB Crich spectra file: '//trim(file_path))
                return
            end if
            self%agb_cache%spec_c(i, :) = temp_row_c
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Crich_Aringer.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3140, 'Failed to open AGB Aringer Teff file: '//trim(file_path))
            return
        end if
        read(u_file, *, iostat=io_stat) self%agb_cache%logt_car
        close(u_file)
        if (io_stat /= 0) then
            call status%set_error(3141, 'Failed to read AGB Aringer Teff file: '//trim(file_path))
            return
        end if
        self%agb_cache%logt_car = log10(self%agb_cache%logt_car)

        file_path = trim(self%sps_home)//'/data/spectra/AGB_spectra/Crich_Aringer.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3142, 'Failed to open AGB Aringer spectra file: '//trim(file_path))
            return
        end if
        do i = 1, NSPEC_ARINGER
            read(u_file, *, iostat=io_stat) self%agb_cache%lam_car(i), temp_row_car
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3143, 'Failed to read AGB Aringer spectra file: '//trim(file_path))
                return
            end if
            self%agb_cache%spec_car(i, :) = temp_row_car
        end do
        close(u_file)
    end subroutine load_legacy_agb

    subroutine load_legacy_dust_emission(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, k, i_spec
        integer :: ndim_dust, numin_dust, nqpah
        character(len=1024) :: file_path
        character(len=1) :: qpah_char
        character(len=1024) :: header_line
        real(WP), allocatable :: raw_lam(:)
        real(WP), allocatable :: raw_spec(:,:)

        call status%set_ok()
        call self%dust_em_cache%clear()

        if (trim(self%dust_type) == 'THEMIS') then
            ndim_dust = 576
            numin_dust = 37
            nqpah = 11
            allocate(self%dust_em_cache%qpah(nqpah))
            allocate(self%dust_em_cache%umin(numin_dust))
            self%dust_em_cache%qpah = QPAH_ARR_THEMIS
            self%dust_em_cache%umin = UMIN_ARR_THEMIS
        else
            ndim_dust = 1001
            numin_dust = 22
            nqpah = 7
            allocate(self%dust_em_cache%qpah(nqpah))
            allocate(self%dust_em_cache%umin(numin_dust))
            self%dust_em_cache%qpah = QPAH_ARR_DL07
            self%dust_em_cache%umin = UMIN_ARR_DL07
        end if

        allocate(self%dust_em_cache%lam(ndim_dust))
        allocate(self%dust_em_cache%spec(ndim_dust, nqpah, numin_dust*2))
        self%dust_em_cache%lam = 0.0_wp
        self%dust_em_cache%spec = 0.0_wp

        allocate(raw_lam(ndim_dust))
        allocate(raw_spec(ndim_dust, numin_dust*2))

        do k = 1, nqpah
            if (k - 1 == 10) then
                file_path = trim(self%sps_home)//'/data/dust/dustem/'//trim(self%dust_type)//'_MW3.1_100.dat'
            else
                write(qpah_char, '(I1)') k - 1
                file_path = trim(self%sps_home)//'/data/dust/dustem/'//trim(self%dust_type)//'_MW3.1_'//qpah_char//'0.dat'
            end if

            open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
            if (io_stat /= 0) then
                call status%set_error(3144, 'Failed to open dust emission file: '//trim(file_path))
                return
            end if

            read(u_file, '(A)', iostat=io_stat) header_line
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3145, 'Failed to read dust emission header: '//trim(file_path))
                return
            end if
            read(u_file, '(A)', iostat=io_stat) header_line
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3145, 'Failed to read dust emission header: '//trim(file_path))
                return
            end if

            do i_spec = 1, ndim_dust
                read(u_file, *, iostat=io_stat) raw_lam(i_spec), raw_spec(i_spec, :)
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3146, 'Failed to read dust emission table: '//trim(file_path))
                    return
                end if
            end do
            close(u_file)

            raw_lam = raw_lam * 1.0e4_wp
            if (k == 1) then
                self%dust_em_cache%lam = raw_lam
            else
                if (any(abs(raw_lam - self%dust_em_cache%lam) > 0.0_wp)) then
                    call status%set_error(3147, 'Dust emission wavelength grids are inconsistent across qpah files.')
                    return
                end if
            end if

            self%dust_em_cache%spec(:, k, :) = raw_spec
        end do
    end subroutine load_legacy_dust_emission

    subroutine load_legacy_agn_dust(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, i
        character(len=1024) :: file_path
        character(len=1024) :: header_line

        call status%set_ok()
        call self%agn_dust_cache%clear()

        allocate(self%agn_dust_cache%tau(NAGNDUST))
        allocate(self%agn_dust_cache%lam(NAGNDUST_SPEC))
        allocate(self%agn_dust_cache%spec(NAGNDUST_SPEC, NAGNDUST))

        file_path = trim(self%sps_home)//'/data/dust/Nenkova08_y010_torusg_n10_q2.0.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3148, 'Failed to open AGN dust file: '//trim(file_path))
            return
        end if

        do i = 1, 3
            read(u_file, '(A)', iostat=io_stat) header_line
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3149, 'Failed to read AGN dust header: '//trim(file_path))
                return
            end if
        end do

        read(u_file, *, iostat=io_stat) self%agn_dust_cache%tau
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3150, 'Failed to read AGN dust tau grid: '//trim(file_path))
            return
        end if

        do i = 1, NAGNDUST_SPEC
            read(u_file, *, iostat=io_stat) self%agn_dust_cache%lam(i), self%agn_dust_cache%spec(i, :)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3151, 'Failed to read AGN dust spectra table: '//trim(file_path))
                return
            end if
        end do
        close(u_file)
    end subroutine load_legacy_agn_dust

    subroutine load_legacy_dust_attenuation(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer, parameter :: NWG_LAM = 25
        integer, parameter :: NWG_GEOM = 18
        integer, parameter :: NWG_TAU = 6
        integer, parameter :: NSMC_LAM = 30

        integer :: u_file, io_stat, i, j
        character(len=1024) :: file_path
        character(len=1024) :: header_line
        real(WP) :: d1, lam_tmp

        call status%set_ok()
        call self%dust_att_cache%clear()

        allocate(self%dust_att_cache%wg_lam(NWG_LAM))
        allocate(self%dust_att_cache%wg_spec(NWG_LAM, NWG_GEOM, NWG_TAU, 2))
        allocate(self%dust_att_cache%smc_lam(NSMC_LAM))
        allocate(self%dust_att_cache%smc_ext(NSMC_LAM))

        file_path = trim(self%sps_home)//'/data/dust/alldirty_h.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3152, 'Failed to open dust attenuation file: '//trim(file_path))
            return
        end if

        do i = 1, 3
            read(u_file, '(A)', iostat=io_stat) header_line
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3153, 'Failed to read dust attenuation header: '//trim(file_path))
                return
            end if
        end do

        do i = 1, NWG_GEOM
            do j = 1, NWG_LAM
                read(u_file, *, iostat=io_stat) self%dust_att_cache%wg_lam(j), d1, self%dust_att_cache%wg_spec(j, i, :, 1)
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3154, 'Failed to read WG shell attenuation table: '//trim(file_path))
                    return
                end if
            end do
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/dust/alldirty_c.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3155, 'Failed to open dust attenuation file: '//trim(file_path))
            return
        end if

        do i = 1, 3
            read(u_file, '(A)', iostat=io_stat) header_line
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3156, 'Failed to read dust attenuation header: '//trim(file_path))
                return
            end if
        end do

        do i = 1, NWG_GEOM
            do j = 1, NWG_LAM
                read(u_file, *, iostat=io_stat) lam_tmp, d1, self%dust_att_cache%wg_spec(j, i, :, 2)
                if (io_stat /= 0) then
                    close(u_file)
                    call status%set_error(3157, 'Failed to read WG cloudy attenuation table: '//trim(file_path))
                    return
                end if
                if (abs(lam_tmp - self%dust_att_cache%wg_lam(j)) > 1.0e-12_wp) then
                    close(u_file)
                    call status%set_error(3157, 'WG attenuation wavelength grids are inconsistent between shell and cloudy tables.')
                    return
                end if
            end do
        end do
        close(u_file)

        file_path = trim(self%sps_home)//'/data/dust/Gordon03_table4.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3158, 'Failed to open SMC attenuation file: '//trim(file_path))
            return
        end if

        do i = 1, NSMC_LAM
            read(u_file, *, iostat=io_stat) self%dust_att_cache%smc_lam(NSMC_LAM-i+1), d1, self%dust_att_cache%smc_ext(NSMC_LAM-i+1)
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3159, 'Failed to read SMC attenuation table: '//trim(file_path))
                return
            end if
        end do
        close(u_file)

        self%dust_att_cache%smc_lam = self%dust_att_cache%smc_lam * 1.0e4_wp
    end subroutine load_legacy_dust_attenuation

    subroutine load_legacy_xrb(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, i, j
        character(len=1024) :: file_path
        character(len=5), dimension(NZ_XRB) :: zz_str
        real(WP), allocatable :: tspec(:)

        call status%set_ok()
        call self%xrb_cache%clear()

        allocate(self%xrb_cache%lam(NSPEC_XRB))
        allocate(self%xrb_cache%age(NT_XRB))
        allocate(self%xrb_cache%z(NZ_XRB))
        allocate(self%xrb_cache%spec(NSPEC_XRB, NT_XRB, NZ_XRB))
        allocate(tspec(NSPEC_XRB))

        file_path = trim(self%sps_home)//'/data/spectra/xrb/xsp.lambda'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            deallocate(tspec)
            call status%set_error(3160, 'Failed to open XRB lambda file: '//trim(file_path))
            return
        end if

        do i = 1, NSPEC_XRB
            read(u_file, *, iostat=io_stat) self%xrb_cache%lam(i)
            if (io_stat /= 0) then
                close(u_file)
                deallocate(tspec)
                call status%set_error(3161, 'Failed to read XRB lambda file: '//trim(file_path))
                return
            end if
        end do
        close(u_file)

        self%xrb_cache%age = log10((/1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 8.0_wp, 10.0_wp, 12.6_wp, 16.0_wp, 20.0_wp/)) + 6.0_wp
        self%xrb_cache%z = (/-1.3_wp, -1.0_wp, -0.8_wp, -0.7_wp, -0.5_wp, -0.4_wp, -0.3_wp, -0.2_wp, 0.0_wp, 0.2_wp, 0.3_wp/)
        zz_str = (/'-1.30', '-1.00', '-0.80', '-0.70', '-0.50', '-0.40', '-0.30', '-0.20', '+0.00', '+0.20', '+0.30'/)

        do j = 1, NZ_XRB
            file_path = trim(self%sps_home)//'/data/spectra/xrb/xsp_feh'//zz_str(j)//'.spec'
            open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
            if (io_stat /= 0) then
                deallocate(tspec)
                call status%set_error(3162, 'Failed to open XRB spec file: '//trim(file_path))
                return
            end if

            do i = 1, NT_XRB
                read(u_file, *, iostat=io_stat) tspec
                if (io_stat /= 0) then
                    close(u_file)
                    deallocate(tspec)
                    call status%set_error(3163, 'Failed to read XRB spec file: '//trim(file_path))
                    return
                end if
                self%xrb_cache%spec(:, i, j) = tspec
            end do
            close(u_file)
        end do

        deallocate(tspec)
    end subroutine load_legacy_xrb

    subroutine load_legacy_nebular_grid(self, use_cloudy, grid, status)
        class(legacy_backend_t), intent(inout) :: self
        logical, intent(in) :: use_cloudy
        type(nebular_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, i, j, k
        character(len=1024) :: file_path_cont, file_path_lines
        character(len=32) :: suffix
        real(WP), allocatable :: raw_lam(:), raw_spec(:)
        real(WP) :: val_logz, val_age, val_logu
        real(WP), parameter :: NEBULAR_FLOOR = 10.0_wp**(-95.0_wp)

        call status%set_ok()
        call grid%clear()

        if (.not. allocated(self%spec_lambda)) then
            call status%set_error(3090, 'Legacy spec_lambda must be loaded before nebular grids.')
            return
        end if

        if (use_cloudy) then
            suffix = '_WD_'
        else
            suffix = '_ND_'
        end if

        file_path_cont = trim(self%sps_home)//'/data/nebular/ZAU'//trim(suffix)//trim(self%isoc_type)//'.cont'
        file_path_lines = trim(self%sps_home)//'/data/nebular/ZAU'//trim(suffix)//trim(self%isoc_type)//'.lines'

        allocate(grid%logz(NEBNZ), grid%age(NEBNAGE), grid%logu(NEBNIP), grid%line_pos(NEMLINE))
        allocate(grid%cont(self%nspec, NEBNZ, NEBNAGE, NEBNIP))
        allocate(grid%line(NEMLINE, NEBNZ, NEBNAGE, NEBNIP))

        open(newunit=u_file, file=trim(file_path_cont), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3091, 'Failed to open nebular continuum file: '//trim(file_path_cont))
            return
        end if

        read(u_file, *, iostat=io_stat)
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3092, 'Failed reading nebular continuum header: '//trim(file_path_cont))
            return
        end if

        allocate(raw_lam(NLAM_NEBCONT), raw_spec(NLAM_NEBCONT))
        read(u_file, *, iostat=io_stat) raw_lam
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3093, 'Failed reading nebular wavelength grid: '//trim(file_path_cont))
            return
        end if

        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    read(u_file, *, iostat=io_stat) val_logz, val_age, val_logu
                    if (io_stat /= 0) then
                        close(u_file)
                        call status%set_error(3094, 'Failed reading nebular continuum metadata: '//trim(file_path_cont))
                        return
                    end if

                    read(u_file, *, iostat=io_stat) raw_spec
                    if (io_stat /= 0) then
                        close(u_file)
                        call status%set_error(3095, 'Failed reading nebular continuum payload: '//trim(file_path_cont))
                        return
                    end if

                    grid%logz(i) = val_logz
                    grid%age(j) = val_age
                    grid%logu(k) = val_logu
                    grid%cont(:, i, j, k) = interpolate_linear(raw_lam, log10(raw_spec + NEBULAR_FLOOR), self%spec_lambda)
                end do
            end do
        end do
        close(u_file)

        open(newunit=u_file, file=trim(file_path_lines), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3096, 'Failed to open nebular lines file: '//trim(file_path_lines))
            return
        end if

        read(u_file, *, iostat=io_stat)
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3097, 'Failed reading nebular line header: '//trim(file_path_lines))
            return
        end if

        read(u_file, *, iostat=io_stat) grid%line_pos
        if (io_stat /= 0) then
            close(u_file)
            call status%set_error(3098, 'Failed reading nebular line positions: '//trim(file_path_lines))
            return
        end if

        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    read(u_file, *, iostat=io_stat) val_logz, val_age, val_logu
                    if (io_stat /= 0) then
                        close(u_file)
                        call status%set_error(3099, 'Failed reading nebular line metadata: '//trim(file_path_lines))
                        return
                    end if

                    read(u_file, *, iostat=io_stat) grid%line(:, i, j, k)
                    if (io_stat /= 0) then
                        close(u_file)
                        call status%set_error(3100, 'Failed reading nebular line payload: '//trim(file_path_lines))
                        return
                    end if

                    grid%age(j) = val_age
                end do
            end do
        end do
        close(u_file)

        grid%age = log10(grid%age)
        grid%line = log10(grid%line + NEBULAR_FLOOR)

        if (allocated(raw_lam)) deallocate(raw_lam)
        if (allocated(raw_spec)) deallocate(raw_spec)
    end subroutine load_legacy_nebular_grid

    subroutine load_legacy_iso_zlegend(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        character(len=1024) :: file_path
        character(len=5) :: mist_str
        real(WP), allocatable :: zvals(:)
        real(WP) :: val_read, zsol
        integer :: nz, u_file, io_stat, i

        call status%set_ok()

        select case (trim(self%isoc_type))
        case ('mist')
            nz = NZ_MIST
            zsol = 0.0142_wp
            file_path = trim(self%sps_home)//'/data/isochrones/MIST/zlegend.dat'
        case ('pdva')
            nz = NZ_PADOVA
            zsol = 0.019_wp
            file_path = trim(self%sps_home)//'/data/isochrones/Padova/Padova2007/zlegend.dat'
        case ('prsc')
            nz = NZ_PARSEC
            zsol = 0.01524_wp
            file_path = trim(self%sps_home)//'/data/isochrones/PARSEC/zlegend.dat'
        case ('bsti')
            nz = NZ_BASTI
            zsol = 0.020_wp
            file_path = trim(self%sps_home)//'/data/isochrones/BaSTI/zlegend.dat'
        case ('gnva')
            nz = NZ_GENEVA
            zsol = 0.020_wp
            file_path = trim(self%sps_home)//'/data/isochrones/Geneva/zlegend.dat'
        case ('bpss')
            nz = NZ_BPASS
            zsol = 0.020_wp
            file_path = trim(self%sps_home)//'/data/isochrones/BPASS/zlegend.dat'
        case default
            call status%set_error(3080, 'Unsupported isoc_type in legacy backend: '//trim(self%isoc_type))
            return
        end select

        allocate(zvals(nz))
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            deallocate(zvals)
            call status%set_error(3081, 'Failed to open isochrone zlegend file: '//trim(file_path))
            return
        end if

        do i = 1, nz
            if (trim(self%isoc_type) == 'mist') then
                read(u_file, '(A5)', iostat=io_stat) mist_str
                if (io_stat /= 0) exit
                read(mist_str(2:5), '(F4.2)', iostat=io_stat) val_read
                if (io_stat /= 0) exit
                if (mist_str(1:1) == 'm') then
                    zvals(i) = 10.0_wp**(-val_read) * zsol
                else
                    zvals(i) = 10.0_wp**(val_read) * zsol
                end if
            else
                read(u_file, *, iostat=io_stat) zvals(i)
                if (io_stat /= 0) exit
            end if
        end do
        close(u_file)

        if (io_stat /= 0) then
            deallocate(zvals)
            call status%set_error(3082, 'Failed reading isochrone zlegend values: '//trim(file_path))
            return
        end if

        if (allocated(self%zlegend)) deallocate(self%zlegend)
        allocate(self%zlegend(nz))
        self%zlegend = zvals
        deallocate(zvals)
    end subroutine load_legacy_iso_zlegend

    subroutine load_legacy_isochrones(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        integer :: nt, nz, z_idx
        real(WP) :: zsol

        call status%set_ok()
        call self%iso_cache%clear()

        select case (trim(self%isoc_type))
        case ('mist')
            nt = NT_MIST
            nz = NZ_MIST
            zsol = 0.0142_wp
        case ('pdva')
            nt = NT_PADOVA
            nz = NZ_PADOVA
            zsol = 0.019_wp
        case ('prsc')
            nt = NT_PARSEC
            nz = NZ_PARSEC
            zsol = 0.01524_wp
        case ('bsti')
            nt = NT_BASTI
            nz = NZ_BASTI
            zsol = 0.020_wp
        case ('gnva')
            nt = NT_GENEVA
            nz = NZ_GENEVA
            zsol = 0.020_wp
        case ('bpss')
            nt = NT_BPASS
            nz = NZ_BPASS
            zsol = 0.020_wp
        case default
            call status%set_error(3083, 'Unsupported isoc_type for legacy isochrone loading: '//trim(self%isoc_type))
            return
        end select

        if (.not. allocated(self%zlegend)) then
            call status%set_error(3084, 'Isochrone zlegend must be loaded before loading tracks.')
            return
        end if

        allocate(self%iso_cache%nmass(nt, nz))
        allocate(self%iso_cache%timestep_logyr(nt, nz))
        allocate(self%iso_cache%mini(NM, nt, nz))
        allocate(self%iso_cache%mact(NM, nt, nz))
        allocate(self%iso_cache%logl(NM, nt, nz))
        allocate(self%iso_cache%logt(NM, nt, nz))
        allocate(self%iso_cache%logg(NM, nt, nz))
        allocate(self%iso_cache%phase(NM, nt, nz))
        allocate(self%iso_cache%ffco(NM, nt, nz))
        allocate(self%iso_cache%lmdot(NM, nt, nz))

        self%iso_cache%nmass = 0
        self%iso_cache%timestep_logyr = 0.0_wp
        self%iso_cache%mini = 0.0_wp
        self%iso_cache%mact = 0.0_wp
        self%iso_cache%logl = 0.0_wp
        self%iso_cache%logt = 0.0_wp
        self%iso_cache%logg = 0.0_wp
        self%iso_cache%phase = 0.0_wp
        self%iso_cache%ffco = 0.0_wp
        self%iso_cache%lmdot = -99.0_wp

        if (trim(self%isoc_type) == 'bpss') then
            call load_legacy_bpss_isochrones(self, status)
            return
        end if

        do z_idx = 1, nz
            call read_one_legacy_iso_file(self, z_idx, nt, zsol, status)
            if (status%code /= 0) return
        end do
    end subroutine load_legacy_isochrones

    subroutine load_legacy_bpss_isochrones(self, status)
        class(legacy_backend_t), intent(inout) :: self
        type(backend_status_t), intent(out) :: status

        real(WP), allocatable :: time_full(:), mass_ssp(:,:)
        integer :: it, iz

        call status%set_ok()
        call read_bpass_mass_table(self, time_full, mass_ssp, status)
        if (status%code /= 0) return

        if (size(time_full) /= NT_BPASS .or. size(mass_ssp, 1) /= NT_BPASS .or. size(mass_ssp, 2) /= NZ_BPASS) then
            if (allocated(time_full)) deallocate(time_full)
            if (allocated(mass_ssp)) deallocate(mass_ssp)
            call status%set_error(3085, 'Unexpected BPASS mass table dimensions while loading BPSS isochrones.')
            return
        end if

        do it = 1, NT_BPASS
            do iz = 1, NZ_BPASS
                self%iso_cache%nmass(it, iz) = 1
                self%iso_cache%timestep_logyr(it, iz) = time_full(it)
                self%iso_cache%mini(1, it, iz) = mass_ssp(it, iz)
                self%iso_cache%mact(1, it, iz) = mass_ssp(it, iz)
                self%iso_cache%logl(1, it, iz) = 0.0_wp
                self%iso_cache%logt(1, it, iz) = 0.0_wp
                self%iso_cache%logg(1, it, iz) = 0.0_wp
                self%iso_cache%phase(1, it, iz) = 0.0_wp
                self%iso_cache%ffco(1, it, iz) = 0.0_wp
                self%iso_cache%lmdot(1, it, iz) = -99.0_wp
            end do
        end do

        if (allocated(time_full)) deallocate(time_full)
        if (allocated(mass_ssp)) deallocate(mass_ssp)
    end subroutine load_legacy_bpss_isochrones

    subroutine read_one_legacy_iso_file(self, z_idx, nt, zsol, status)
        class(legacy_backend_t), intent(inout) :: self
        integer, intent(in) :: z_idx
        integer, intent(in) :: nt
        real(WP), intent(in) :: zsol
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, i_line
        integer :: current_track, current_step
        character(len=1024) :: file_path
        character(len=4096) :: line_buf
        logical :: is_mist
        real(WP) :: r_logage, r_mini, r_mact, r_logl, r_logt, r_logg
        real(WP) :: r_ffco, r_phase, r_lmdot

        call status%set_ok()
        call get_isochrone_filename(self, z_idx, zsol, file_path, status)
        if (status%code /= 0) return

        is_mist = (trim(self%isoc_type) == 'mist')

        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3085, 'Cannot open isochrone file: '//trim(file_path))
            return
        end if

        current_track = 0
        current_step = 0

        line_loop: do i_line = 1, NLINES
            read(u_file, '(A)', iostat=io_stat) line_buf
            if (io_stat == iostat_end) exit line_loop
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3086, 'Error reading isochrone file line in: '//trim(file_path))
                return
            end if

            line_buf = adjustl(line_buf)
            if (len_trim(line_buf) == 0) cycle

            if (line_buf(1:1) == '#') then
                current_step = 0
                cycle
            end if

            if (is_mist) then
                read(line_buf, *, iostat=io_stat) r_logage, r_mini, r_mact, r_logl, r_logt, r_logg, r_ffco, r_phase, r_lmdot
            else
                read(line_buf, *, iostat=io_stat) r_logage, r_mini, r_mact, r_logl, r_logt, r_logg, r_ffco, r_phase
                r_lmdot = -99.0_wp
            end if

            if (io_stat /= 0) then
                current_step = 0
                cycle line_loop
            end if

            if (current_step == 0) then
                current_track = current_track + 1
                if (current_track > nt) then
                    close(u_file)
                    call status%set_error(3087, 'Tracks exceed expected NT in file: '//trim(file_path))
                    return
                end if
            end if

            current_step = current_step + 1
            if (current_step > NM) then
                close(u_file)
                call status%set_error(3088, 'Mass steps exceed NM in file: '//trim(file_path))
                return
            end if

            self%iso_cache%mini(current_step, current_track, z_idx) = r_mini
            self%iso_cache%mact(current_step, current_track, z_idx) = r_mact
            self%iso_cache%logl(current_step, current_track, z_idx) = r_logl
            self%iso_cache%logt(current_step, current_track, z_idx) = r_logt
            self%iso_cache%logg(current_step, current_track, z_idx) = r_logg
            self%iso_cache%ffco(current_step, current_track, z_idx) = r_ffco
            self%iso_cache%phase(current_step, current_track, z_idx) = r_phase
            self%iso_cache%lmdot(current_step, current_track, z_idx) = r_lmdot

            if (current_step == 1) then
                self%iso_cache%timestep_logyr(current_track, z_idx) = r_logage
            end if
            self%iso_cache%nmass(current_track, z_idx) = current_step
        end do line_loop

        close(u_file)
    end subroutine read_one_legacy_iso_file

    subroutine get_isochrone_filename(self, z_idx, zsol, file_path, status)
        class(legacy_backend_t), intent(in) :: self
        integer, intent(in) :: z_idx
        real(WP), intent(in) :: zsol
        character(len=*), intent(out) :: file_path
        type(backend_status_t), intent(out) :: status

        character(len=6) :: z_str
        real(WP) :: val_log, z_val

        call status%set_ok()

        z_val = self%zlegend(z_idx)

        if (trim(self%isoc_type) == 'mist') then
            if (z_val < tiny(0.0_wp)) then
                val_log = -99.0_wp
            else
                val_log = log10(z_val / zsol)
            end if

            if (val_log < -0.001_wp) then
                write(z_str, '("m", F4.2)') abs(val_log)
            else
                write(z_str, '("p", F4.2)') abs(val_log)
            end if
            file_path = trim(self%sps_home)//'/data/isochrones/MIST/isoc_z'//trim(z_str)//'.dat'
            return
        end if

        write(z_str, '(F6.4)') z_val
        select case (trim(self%isoc_type))
        case ('pdva')
            file_path = trim(self%sps_home)//'/data/isochrones/Padova/Padova2007/isoc_z'//z_str//'.dat'
        case ('prsc')
            file_path = trim(self%sps_home)//'/data/isochrones/PARSEC/isoc_z'//z_str//'.dat'
        case ('bsti')
            file_path = trim(self%sps_home)//'/data/isochrones/BaSTI/isoc_z'//z_str//'.dat'
        case ('gnva')
            file_path = trim(self%sps_home)//'/data/isochrones/Geneva/isoc_z'//z_str//'.dat'
        case default
            call status%set_error(3089, 'Unknown isoc_type for filename construction: '//trim(self%isoc_type))
        end select
    end subroutine get_isochrone_filename

    subroutine read_real_column_file(file_path, values, status)
        character(len=*), intent(in) :: file_path
        real(WP), allocatable, intent(out) :: values(:)
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, n, i
        real(WP) :: x

        if (allocated(values)) deallocate(values)
        call status%set_ok()

        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3060, 'Failed to open axis file: '//trim(file_path))
            return
        end if

        n = 0
        do
            read(u_file, *, iostat=io_stat) x
            if (io_stat == iostat_end) exit
            if (io_stat /= 0) then
                close(u_file)
                call status%set_error(3061, 'Failed to parse axis file: '//trim(file_path))
                return
            end if
            n = n + 1
        end do

        if (n <= 0) then
            close(u_file)
            call status%set_error(3062, 'Axis file is empty: '//trim(file_path))
            return
        end if

        rewind(u_file)
        allocate(values(n))
        do i = 1, n
            read(u_file, *, iostat=io_stat) values(i)
            if (io_stat /= 0) then
                close(u_file)
                if (allocated(values)) deallocate(values)
                call status%set_error(3063, 'Failed to read axis values from: '//trim(file_path))
                return
            end if
        end do

        close(u_file)
    end subroutine read_real_column_file

    subroutine build_lambda_file_path(self, file_path, status)
        class(legacy_backend_t), intent(in) :: self
        character(len=:), allocatable, intent(out) :: file_path
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        select case (trim(self%spec_type))
        case ('basel')
            file_path = trim(self%sps_home)//'/data/spectra/BaSeL3.1/basel.lambda'
        case ('miles')
            file_path = trim(self%sps_home)//'/data/spectra/MILES/miles.lambda'
        case ('bpass')
            file_path = trim(self%sps_home)//'/data/isochrones/BPASS/bpass.lambda'
        case ('ckc14')
            file_path = trim(self%sps_home)//'/data/spectra/CKC14/ckc14.lambda'
        case default
            if (len_trim(self%spec_type) >= 3 .and. self%spec_type(1:3) == 'c3k') then
                file_path = trim(self%sps_home)//'/data/spectra/C3K/'//trim(self%spec_type)//'.lambda'
            else
                call status%set_error(3070, 'Unknown legacy spec_type for lambda path: '//trim(self%spec_type))
            end if
        end select
    end subroutine build_lambda_file_path

    subroutine build_zlegend_file_path(self, file_path, status)
        class(legacy_backend_t), intent(in) :: self
        character(len=:), allocatable, intent(out) :: file_path
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        select case (trim(self%spec_type))
        case ('basel')
            file_path = trim(self%sps_home)//'/data/spectra/BaSeL3.1/zlegend.dat'
        case ('miles')
            file_path = trim(self%sps_home)//'/data/spectra/MILES/zlegend.dat'
        case ('bpass')
            file_path = trim(self%sps_home)//'/data/isochrones/BPASS/zlegend.dat'
        case ('ckc14')
            file_path = trim(self%sps_home)//'/data/spectra/CKC14/zlegend.dat'
        case default
            if (len_trim(self%spec_type) >= 3 .and. self%spec_type(1:3) == 'c3k') then
                file_path = trim(self%sps_home)//'/data/spectra/C3K/zlegend.dat'
            else
                call status%set_error(3071, 'Unknown legacy spec_type for zlegend path: '//trim(self%spec_type))
            end if
        end select
    end subroutine build_zlegend_file_path

    subroutine build_binary_file_path(self, iz, file_path, status)
        class(legacy_backend_t), intent(in) :: self
        integer, intent(in) :: iz
        character(len=:), allocatable, intent(out) :: file_path
        type(backend_status_t), intent(out) :: status

        character(len=6) :: z_str

        call status%set_ok()

        if (.not. allocated(self%zlegend)) then
            call status%set_error(3072, 'zlegend is not loaded in legacy backend.')
            return
        end if

        write(z_str, '(F6.4)') self%zlegend(iz)

        select case (trim(self%spec_type))
        case ('basel')
            file_path = trim(self%sps_home)//'/data/spectra/BaSeL3.1/basel_'//trim(BASEL_STR)//'_z'//trim(z_str)//'.spectra.bin'
        case ('miles')
            file_path = trim(self%sps_home)//'/data/spectra/MILES/imiles_z'//z_str//'.spectra.bin'
        case ('bpass')
            file_path = trim(self%sps_home)//'/data/isochrones/BPASS/bpass_v2.2_salpeter100.ssp.bin'
        case ('ckc14')
            file_path = trim(self%sps_home)//'/data/spectra/CKC14/'//trim(self%spec_type)//'_z'//z_str//'.spectra.bin'
        case default
            if (len_trim(self%spec_type) >= 3 .and. self%spec_type(1:3) == 'c3k') then
                file_path = trim(self%sps_home)//'/data/spectra/C3K/'//trim(self%spec_type)//'_z'//z_str//'.spectra.bin'
            else
                call status%set_error(3073, 'Unknown legacy spec_type for binary path: '//trim(self%spec_type))
            end if
        end select
    end subroutine build_binary_file_path

    subroutine read_bpass_mass_table(self, time_full, mass_ssp, status)
        class(legacy_backend_t), intent(in) :: self
        real(WP), allocatable, intent(out) :: time_full(:)
        real(WP), allocatable, intent(out) :: mass_ssp(:,:)
        type(backend_status_t), intent(out) :: status

        integer :: u_file, io_stat, it, nz
        character(len=1024) :: file_path

        if (allocated(time_full)) deallocate(time_full)
        if (allocated(mass_ssp)) deallocate(mass_ssp)
        call status%set_ok()

        nz = NZ_BPASS
        file_path = trim(self%sps_home)//'/data/isochrones/BPASS/bpass.mass'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call status%set_error(3164, 'Failed to open BPASS mass file: '//trim(file_path))
            return
        end if

        allocate(time_full(NT_BPASS))
        allocate(mass_ssp(NT_BPASS, nz))

        do it = 1, NT_BPASS
            read(u_file, *, iostat=io_stat) time_full(it), mass_ssp(it, :)
            if (io_stat /= 0) then
                close(u_file)
                if (allocated(time_full)) deallocate(time_full)
                if (allocated(mass_ssp)) deallocate(mass_ssp)
                call status%set_error(3165, 'Failed to parse BPASS mass row in: '//trim(file_path))
                return
            end if
        end do

        close(u_file)
    end subroutine read_bpass_mass_table

end module fsps_data_backend_legacy
