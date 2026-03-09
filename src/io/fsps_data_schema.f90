module fsps_data_schema
    !> @brief Core schema types for FSPS standardized data ingestion.
    !>
    !> @details
    !> This module defines pure data containers (no backend logic) used by
    !> metadata-driven ingestion. The design assumes dense N-D spectral grids
    !> on disk and in memory. Sparse physical libraries must be padded by
    !> external authoring/conversion tools before Fortran ingestion.
    !>
    !> Key memory-layout rule:
    !> - Spectral arrays are always ordered with wavelength first:
    !>   `flux(lambda, z, afe, logt, logg)`.
    !>   This keeps the lambda axis contiguous in memory, matching FSPS usage.

    use fsps_precision, only: WP
    use fsps_strings, only: to_lower

    implicit none
    private

    public :: axis_desc_t
    public :: dataset_desc_t
    public :: library_manifest_t
    public :: spectral_grid_t
    public :: spectral_slice_t
    public :: isochrone_grid_t
    public :: nebular_grid_t
    public :: aux_wmbasic_t
    public :: aux_pagb_t
    public :: aux_wr_t
    public :: aux_agb_t
    public :: dust_emission_t
    public :: agn_dust_t
    public :: dust_attenuation_t
    public :: xrb_spectra_t

    !> @brief Axis descriptor.
    !>
    !> @details
    !> Describes one independent coordinate axis in the schema.
    !> Examples: `lambda`, `z`, `afe`, `logt`, `logg`.
    type :: axis_desc_t
        !> Canonical axis name (e.g., "lambda", "z", "afe").
        character(len=:), allocatable :: name
        !> Unit string (e.g., "Angstrom", "dex").
        character(len=:), allocatable :: unit
        !> Fully-qualified backend path for this axis.
        character(len=:), allocatable :: path
        !> Number of axis points.
        integer :: n = 0
        !> Axis values in working precision.
        real(WP), allocatable :: values(:)
    contains
        procedure, public :: clear => axis_desc_clear
    end type axis_desc_t

    !> @brief Dataset descriptor.
    !>
    !> @details
    !> Metadata for one logical dataset in the FSDS container.
    type :: dataset_desc_t
        !> Logical role (e.g., "spectral_cube", "isochrone_track").
        character(len=:), allocatable :: role
        !> Fully-qualified backend path.
        character(len=:), allocatable :: path
        !> Ordered dimension labels (CSV), e.g. "lambda,z,afe,logt,logg".
        character(len=:), allocatable :: dims_csv
        !> Physical unit string.
        character(len=:), allocatable :: unit
        !> Transform tag (e.g., "none", "log10_plus_floor:1e-95").
        character(len=:), allocatable :: transform
        !> Storage dtype string (e.g., "float32", "float64").
        character(len=:), allocatable :: dtype
        !> Representation type. Enforced default is "dense_nd".
        character(len=:), allocatable :: representation
        !> Raw shape in storage order.
        integer, allocatable :: shape(:)
        !> Whether this dataset is required for the selected model tuple.
        logical :: required = .true.
        !> Whether a finite missing-value sentinel is defined.
        logical :: has_missing_value = .false.
        !> Finite sentinel for invalid nodes (e.g., -1.0e30).
        real(WP) :: missing_value = -1.0e30_wp
        !> Whether a validity mask is available.
        logical :: has_valid_mask = .false.
        !> Path to validity mask (if present).
        character(len=:), allocatable :: valid_mask_path
    contains
        procedure, public :: clear => dataset_desc_clear
        procedure, public :: has_axis => dataset_desc_has_axis
    end type dataset_desc_t

    !> @brief Manifest for one resolved library configuration.
    !>
    !> @details
    !> Aggregates axes and dataset descriptors selected for a concrete
    !> `(isochrone, spectral, dust)` tuple.
    type :: library_manifest_t
        !> FSDS schema version string.
        character(len=:), allocatable :: fsds_version
        !> Selected isochrone library identifier.
        character(len=:), allocatable :: isoc_name
        !> Selected spectral library identifier.
        character(len=:), allocatable :: spec_name
        !> Selected dust library identifier.
        character(len=:), allocatable :: dust_name
        !> Declared producer tag.
        character(len=:), allocatable :: producer
        !> Dataset descriptors.
        type(dataset_desc_t), allocatable :: datasets(:)
        !> Axis descriptors.
        type(axis_desc_t), allocatable :: axes(:)
    contains
        procedure, public :: clear => library_manifest_clear
    end type library_manifest_t

    !> @brief Dense spectral grid container.
    !>
    !> @details
    !> Canonical in-memory spectral representation.
    !>
    !> Array order is fixed to:
    !> `flux(lambda, z, afe, logt, logg)`.
    !>
    !> The `afe` dimension may be degenerate (`size=1`) for legacy libraries.
    type :: spectral_grid_t
        !> Spectral flux cube, dense and padded as needed.
        real(WP), allocatable :: flux(:,:,:,:,:)
        !> Optional validity mask for parameter-space nodes:
        !> `valid(z, afe, logt, logg)`.
        logical, allocatable :: valid(:,:,:,:)
        !> Axis vectors.
        real(WP), allocatable :: axis_lambda(:)
        real(WP), allocatable :: axis_z(:)
        real(WP), allocatable :: axis_afe(:)
        real(WP), allocatable :: axis_logt(:)
        real(WP), allocatable :: axis_logg(:)
        !> Finite sentinel used for invalid nodes in `flux`.
        real(WP) :: missing_value = -1.0e30_wp
    contains
        procedure, public :: clear => spectral_grid_clear
    end type spectral_grid_t

    !> @brief Hyperslab slice container for one `(z, afe)` location.
    !>
    !> @details
    !> Dense runtime working set with array order:
    !> `flux(lambda, logt, logg)`.
    type :: spectral_slice_t
        !> Spectral slab for one metallicity/alpha pair.
        real(WP), allocatable :: flux(:,:,:)
        !> Optional validity mask: `valid(logt, logg)`.
        logical, allocatable :: valid(:,:)
        !> Source index in z axis (1-based).
        integer :: iz = 0
        !> Source index in afe axis (1-based).
        integer :: iafe = 0
        !> Finite sentinel used for invalid nodes in `flux`.
        real(WP) :: missing_value = -1.0e30_wp
    contains
        procedure, public :: clear => spectral_slice_clear
    end type spectral_slice_t

    !> @brief Dense isochrone grid container.
    !>
    !> @details
    !> Canonical in-memory container for isochrone tables loaded from FSDS.
    !> Arrays are expected to be loaded in storage order `(nm, nt, nz)` and can
    !> be transposed later by callers that need legacy `(nz, nt, nm)` layout.
    type :: isochrone_grid_t
        integer, allocatable :: nmass(:,:)
        real(WP), allocatable :: timestep_logyr(:,:)
        real(WP), allocatable :: mini(:,:,:)
        real(WP), allocatable :: mact(:,:,:)
        real(WP), allocatable :: logl(:,:,:)
        real(WP), allocatable :: logt(:,:,:)
        real(WP), allocatable :: logg(:,:,:)
        real(WP), allocatable :: phase(:,:,:)
        real(WP), allocatable :: ffco(:,:,:)
        real(WP), allocatable :: lmdot(:,:,:)
    contains
        procedure, public :: clear => isochrone_grid_clear
    end type isochrone_grid_t

    !> @brief Dense nebular emission grid container.
    !>
    !> @details
    !> Canonical in-memory nebular representation with array orders:
    !> - `cont(lambda, z, age, u)`
    !> - `line(line, z, age, u)`
    type :: nebular_grid_t
        real(WP), allocatable :: line_pos(:)
        real(WP), allocatable :: logz(:)
        real(WP), allocatable :: age(:)
        real(WP), allocatable :: logu(:)
        real(WP), allocatable :: cont(:,:,:,:)
        real(WP), allocatable :: line(:,:,:,:)
        real(WP) :: missing_value = -1.0e30_wp
    contains
        procedure, public :: clear => nebular_grid_clear
    end type nebular_grid_t

    !> @brief Dense auxiliary WMBasic grid container.
    !>
    !> @details
    !> Canonical in-memory WMBasic representation with array order:
    !> - `spec(lam, logt, logg, z)`
    type :: aux_wmbasic_t
        real(WP), allocatable :: logt(:), z(:), lam(:)
        real(WP), allocatable :: spec(:,:,:,:)
    contains
        procedure, public :: clear => aux_wmbasic_clear
    end type aux_wmbasic_t

    !> @brief Dense auxiliary Post-AGB grid container.
    !>
    !> @details
    !> Canonical in-memory Post-AGB representation with array order:
    !> - `spec(lam, logt, z)`
    type :: aux_pagb_t
        real(WP), allocatable :: logt(:), lam(:)
        real(WP), allocatable :: spec(:,:,:)
    contains
        procedure, public :: clear => aux_pagb_clear
    end type aux_pagb_t

    !> @brief Dense auxiliary Wolf-Rayet grid container.
    !> @details Arrays ordered as: spec(lam, logt, z)
    type :: aux_wr_t
        real(WP), allocatable :: logt_wn(:), logt_wc(:), z(:), lam(:)
        real(WP), allocatable :: spec_wn(:,:,:), spec_wc(:,:,:)
    contains
        procedure, public :: clear => aux_wr_clear
    end type aux_wr_t

    !> @brief Dense auxiliary AGB grid container.
    !> @details Arrays ordered as: logt_o(z, logt), spec(lam, logt)
    type :: aux_agb_t
        real(WP), allocatable :: lam_o(:), lam_c(:), lam_car(:)
        real(WP), allocatable :: z_o(:), logt_o(:,:), logt_c(:), logt_car(:)
        real(WP), allocatable :: spec_o(:,:), spec_c(:,:), spec_car(:,:)
    contains
        procedure, public :: clear => aux_agb_clear
    end type aux_agb_t

    !> @brief Dense dust emission template grid container.
    !> @details Arrays ordered as: spec(lam, qpah, umin_cols)
    type :: dust_emission_t
        real(WP), allocatable :: qpah(:), umin(:), lam(:)
        real(WP), allocatable :: spec(:,:,:)
    contains
        procedure, public :: clear => dust_emission_clear
    end type dust_emission_t

    !> @brief Dense AGN dust template grid container.
    !> @details Arrays ordered as: spec(lam, tau)
    type :: agn_dust_t
        real(WP), allocatable :: tau(:), lam(:)
        real(WP), allocatable :: spec(:,:)
    contains
        procedure, public :: clear => agn_dust_clear
    end type agn_dust_t

    !> @brief Dense dust attenuation template grid container.
    !> @details Arrays ordered as: wg_spec(lam, geom, tau, type)
    type :: dust_attenuation_t
        real(WP), allocatable :: wg_lam(:)
        real(WP), allocatable :: wg_spec(:,:,:,:)
        real(WP), allocatable :: smc_lam(:)
        real(WP), allocatable :: smc_ext(:)
    contains
        procedure, public :: clear => dust_attenuation_clear
    end type dust_attenuation_t

    !> @brief Dense X-ray binary spectral grid container.
    !> @details Arrays ordered as: spec(lam, age, z)
    type :: xrb_spectra_t
        real(WP), allocatable :: lam(:)
        real(WP), allocatable :: age(:)
        real(WP), allocatable :: z(:)
        real(WP), allocatable :: spec(:,:,:)
    contains
        procedure, public :: clear => xrb_spectra_clear
    end type xrb_spectra_t

contains

    !> @brief Reset an axis descriptor to an uninitialized state.
    subroutine axis_desc_clear(self)
        class(axis_desc_t), intent(inout) :: self

        if (allocated(self%name)) deallocate(self%name)
        if (allocated(self%unit)) deallocate(self%unit)
        if (allocated(self%path)) deallocate(self%path)
        if (allocated(self%values)) deallocate(self%values)
        self%n = 0
    end subroutine axis_desc_clear

    !> @brief Reset a dataset descriptor to an uninitialized state.
    subroutine dataset_desc_clear(self)
        class(dataset_desc_t), intent(inout) :: self

        if (allocated(self%role)) deallocate(self%role)
        if (allocated(self%path)) deallocate(self%path)
        if (allocated(self%dims_csv)) deallocate(self%dims_csv)
        if (allocated(self%unit)) deallocate(self%unit)
        if (allocated(self%transform)) deallocate(self%transform)
        if (allocated(self%dtype)) deallocate(self%dtype)
        if (allocated(self%representation)) deallocate(self%representation)
        if (allocated(self%shape)) deallocate(self%shape)
        if (allocated(self%valid_mask_path)) deallocate(self%valid_mask_path)

        self%required = .true.
        self%has_missing_value = .false.
        self%missing_value = -1.0e30_wp
        self%has_valid_mask = .false.
    end subroutine dataset_desc_clear

    !> @brief Return `.true.` if `dims_csv` contains a named axis.
    !>
    !> @details
    !> Matching is case-insensitive and whitespace-insensitive for each CSV token.
    !> Returns `.false.` if `dims_csv` is not allocated.
    pure logical function dataset_desc_has_axis(self, axis_name)
        class(dataset_desc_t), intent(in) :: self
        character(len=*), intent(in) :: axis_name

        integer :: p0, p1, n
        character(len=:), allocatable :: token
        character(len=:), allocatable :: target

        dataset_desc_has_axis = .false.
        if (.not. allocated(self%dims_csv)) return

        n = len_trim(self%dims_csv)
        if (n <= 0) return

        target = trim(to_lower(trim(axis_name)))

        p0 = 1
        do
            p1 = index(self%dims_csv(p0:n), ',')
            if (p1 == 0) then
                token = trim(to_lower(trim(self%dims_csv(p0:n))))
                if (token == target) dataset_desc_has_axis = .true.
                exit
            else
                token = trim(to_lower(trim(self%dims_csv(p0:p0+p1-2))))
                if (token == target) then
                    dataset_desc_has_axis = .true.
                    exit
                end if
                p0 = p0 + p1
            end if
            if (p0 > n) exit
        end do
    end function dataset_desc_has_axis

    !> @brief Reset a library manifest and all nested descriptors.
    subroutine library_manifest_clear(self)
        class(library_manifest_t), intent(inout) :: self
        integer :: i

        if (allocated(self%datasets)) then
            do i = 1, size(self%datasets)
                call self%datasets(i)%clear()
            end do
            deallocate(self%datasets)
        end if

        if (allocated(self%axes)) then
            do i = 1, size(self%axes)
                call self%axes(i)%clear()
            end do
            deallocate(self%axes)
        end if

        if (allocated(self%fsds_version)) deallocate(self%fsds_version)
        if (allocated(self%isoc_name)) deallocate(self%isoc_name)
        if (allocated(self%spec_name)) deallocate(self%spec_name)
        if (allocated(self%dust_name)) deallocate(self%dust_name)
        if (allocated(self%producer)) deallocate(self%producer)
    end subroutine library_manifest_clear

    !> @brief Reset a dense spectral grid container.
    subroutine spectral_grid_clear(self)
        class(spectral_grid_t), intent(inout) :: self

        if (allocated(self%flux)) deallocate(self%flux)
        if (allocated(self%valid)) deallocate(self%valid)
        if (allocated(self%axis_lambda)) deallocate(self%axis_lambda)
        if (allocated(self%axis_z)) deallocate(self%axis_z)
        if (allocated(self%axis_afe)) deallocate(self%axis_afe)
        if (allocated(self%axis_logt)) deallocate(self%axis_logt)
        if (allocated(self%axis_logg)) deallocate(self%axis_logg)

        self%missing_value = -1.0e30_wp
    end subroutine spectral_grid_clear

    !> @brief Reset a spectral hyperslab container.
    subroutine spectral_slice_clear(self)
        class(spectral_slice_t), intent(inout) :: self

        if (allocated(self%flux)) deallocate(self%flux)
        if (allocated(self%valid)) deallocate(self%valid)

        self%iz = 0
        self%iafe = 0
        self%missing_value = -1.0e30_wp
    end subroutine spectral_slice_clear

    !> @brief Reset a dense isochrone grid container.
    subroutine isochrone_grid_clear(self)
        class(isochrone_grid_t), intent(inout) :: self

        if (allocated(self%nmass)) deallocate(self%nmass)
        if (allocated(self%timestep_logyr)) deallocate(self%timestep_logyr)
        if (allocated(self%mini)) deallocate(self%mini)
        if (allocated(self%mact)) deallocate(self%mact)
        if (allocated(self%logl)) deallocate(self%logl)
        if (allocated(self%logt)) deallocate(self%logt)
        if (allocated(self%logg)) deallocate(self%logg)
        if (allocated(self%phase)) deallocate(self%phase)
        if (allocated(self%ffco)) deallocate(self%ffco)
        if (allocated(self%lmdot)) deallocate(self%lmdot)
    end subroutine isochrone_grid_clear

    !> @brief Reset a dense nebular grid container.
    subroutine nebular_grid_clear(self)
        class(nebular_grid_t), intent(inout) :: self

        if (allocated(self%line_pos)) deallocate(self%line_pos)
        if (allocated(self%logz)) deallocate(self%logz)
        if (allocated(self%age)) deallocate(self%age)
        if (allocated(self%logu)) deallocate(self%logu)
        if (allocated(self%cont)) deallocate(self%cont)
        if (allocated(self%line)) deallocate(self%line)

        self%missing_value = -1.0e30_wp
    end subroutine nebular_grid_clear

    !> @brief Reset a dense WMBasic grid container.
    subroutine aux_wmbasic_clear(self)
        class(aux_wmbasic_t), intent(inout) :: self

        if (allocated(self%logt)) deallocate(self%logt)
        if (allocated(self%z)) deallocate(self%z)
        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%spec)) deallocate(self%spec)
    end subroutine aux_wmbasic_clear

    !> @brief Reset a dense Post-AGB grid container.
    subroutine aux_pagb_clear(self)
        class(aux_pagb_t), intent(inout) :: self

        if (allocated(self%logt)) deallocate(self%logt)
        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%spec)) deallocate(self%spec)
    end subroutine aux_pagb_clear

    !> @brief Reset a dense Wolf-Rayet grid container.
    subroutine aux_wr_clear(self)
        class(aux_wr_t), intent(inout) :: self

        if (allocated(self%logt_wn)) deallocate(self%logt_wn)
        if (allocated(self%logt_wc)) deallocate(self%logt_wc)
        if (allocated(self%z)) deallocate(self%z)
        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%spec_wn)) deallocate(self%spec_wn)
        if (allocated(self%spec_wc)) deallocate(self%spec_wc)
    end subroutine aux_wr_clear

    !> @brief Reset a dense AGB grid container.
    subroutine aux_agb_clear(self)
        class(aux_agb_t), intent(inout) :: self

        if (allocated(self%lam_o)) deallocate(self%lam_o)
        if (allocated(self%lam_c)) deallocate(self%lam_c)
        if (allocated(self%lam_car)) deallocate(self%lam_car)
        if (allocated(self%z_o)) deallocate(self%z_o)
        if (allocated(self%logt_o)) deallocate(self%logt_o)
        if (allocated(self%logt_c)) deallocate(self%logt_c)
        if (allocated(self%logt_car)) deallocate(self%logt_car)
        if (allocated(self%spec_o)) deallocate(self%spec_o)
        if (allocated(self%spec_c)) deallocate(self%spec_c)
        if (allocated(self%spec_car)) deallocate(self%spec_car)
    end subroutine aux_agb_clear

    !> @brief Reset a dense dust emission template grid container.
    subroutine dust_emission_clear(self)
        class(dust_emission_t), intent(inout) :: self

        if (allocated(self%qpah)) deallocate(self%qpah)
        if (allocated(self%umin)) deallocate(self%umin)
        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%spec)) deallocate(self%spec)
    end subroutine dust_emission_clear

    !> @brief Reset a dense AGN dust template grid container.
    subroutine agn_dust_clear(self)
        class(agn_dust_t), intent(inout) :: self

        if (allocated(self%tau)) deallocate(self%tau)
        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%spec)) deallocate(self%spec)
    end subroutine agn_dust_clear

    !> @brief Reset a dense dust attenuation template grid container.
    subroutine dust_attenuation_clear(self)
        class(dust_attenuation_t), intent(inout) :: self

        if (allocated(self%wg_lam)) deallocate(self%wg_lam)
        if (allocated(self%wg_spec)) deallocate(self%wg_spec)
        if (allocated(self%smc_lam)) deallocate(self%smc_lam)
        if (allocated(self%smc_ext)) deallocate(self%smc_ext)
    end subroutine dust_attenuation_clear

    !> @brief Reset a dense X-ray binary spectral grid container.
    subroutine xrb_spectra_clear(self)
        class(xrb_spectra_t), intent(inout) :: self

        if (allocated(self%lam)) deallocate(self%lam)
        if (allocated(self%age)) deallocate(self%age)
        if (allocated(self%z)) deallocate(self%z)
        if (allocated(self%spec)) deallocate(self%spec)
    end subroutine xrb_spectra_clear

end module fsps_data_schema
