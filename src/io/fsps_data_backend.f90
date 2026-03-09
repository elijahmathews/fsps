module fsps_data_backend
    !> @brief Abstract backend API for FSDS data access.
    !>
    !> @details
    !> This module defines the backend contract consumed by FSPS ingestion.
    !> Concrete implementations (HDF5, legacy adapters, etc.) must extend
    !> `data_backend_t` and implement all deferred bindings.
    !>
    !> Design goals:
    !> - Backend-agnostic orchestration in `fsps_initialization`.
    !> - Metadata-first workflow.
    !> - Explicit support for hyperslab slicing to avoid full 5D loads.

    use fsps_precision, only: WP
    use fsps_data_schema, only: axis_desc_t, dataset_desc_t, library_manifest_t, &
                                spectral_grid_t, spectral_slice_t

    implicit none
    private

    public :: backend_status_t
    public :: data_backend_t
    public :: backend_status_ok

    !> @brief Status container for backend operations.
    type :: backend_status_t
        !> Zero indicates success; nonzero indicates an error.
        integer :: code = 0
        !> Human-readable status/error message.
        character(len=:), allocatable :: message
    contains
        procedure, public :: clear => backend_status_clear
        procedure, public :: set_error => backend_status_set_error
        procedure, public :: set_ok => backend_status_set_ok
    end type backend_status_t

    !> @brief Abstract data backend.
    !>
    !> @details
    !> Implementations are responsible for opening data sources,
    !> exposing schema metadata, and reading dense arrays/slices.
    type, abstract :: data_backend_t
    contains
        ! Lifecycle
        procedure(open_backend_if), deferred :: open
        procedure(close_backend_if), deferred :: close
        procedure(is_open_if), deferred :: is_open

        ! Metadata
        procedure(read_manifest_if), deferred :: read_manifest
        procedure(has_path_if), deferred :: has_path
        procedure(query_axis_if), deferred :: query_axis

        ! Generic scalar/array reads
        procedure(read_real_1d_if), deferred :: read_real_1d
        procedure(read_real_2d_if), deferred :: read_real_2d
        procedure(read_real_3d_if), deferred :: read_real_3d
        procedure(read_real_4d_if), deferred :: read_real_4d
        procedure(read_int_1d_if), deferred :: read_int_1d
        procedure(read_int_2d_if), deferred :: read_int_2d
        procedure(read_int_3d_if), deferred :: read_int_3d

        ! Spectral hyperslab reads
        procedure(read_spectral_slice_if), deferred :: read_spectral_slice
        procedure(read_spectral_neighborhood_if), deferred :: read_spectral_neighborhood
    end type data_backend_t

    ! ---------------------------------------------------------------------
    ! Abstract interfaces
    ! ---------------------------------------------------------------------

    abstract interface

        !> @brief Open a backend data source.
        subroutine open_backend_if(self, source_uri, status)
            import :: data_backend_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: source_uri
            type(backend_status_t), intent(out) :: status
        end subroutine open_backend_if

        !> @brief Close a backend data source and release backend-local resources.
        subroutine close_backend_if(self, status)
            import :: data_backend_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            type(backend_status_t), intent(out) :: status
        end subroutine close_backend_if

        !> @brief Query whether the backend currently has an open source.
        logical function is_open_if(self)
            import :: data_backend_t
            class(data_backend_t), intent(in) :: self
        end function is_open_if

        !> @brief Read the resolved library manifest from the source.
        subroutine read_manifest_if(self, manifest, status)
            import :: data_backend_t, library_manifest_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            type(library_manifest_t), intent(inout) :: manifest
            type(backend_status_t), intent(out) :: status
        end subroutine read_manifest_if

        !> @brief Test whether a source path exists.
        logical function has_path_if(self, path)
            import :: data_backend_t
            class(data_backend_t), intent(in) :: self
            character(len=*), intent(in) :: path
        end function has_path_if

        !> @brief Query one axis descriptor by name.
        subroutine query_axis_if(self, axis_name, axis_desc, status)
            import :: data_backend_t, axis_desc_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: axis_name
            type(axis_desc_t), intent(inout) :: axis_desc
            type(backend_status_t), intent(out) :: status
        end subroutine query_axis_if

        !> @brief Read a real 1-D dataset.
        subroutine read_real_1d_if(self, dataset_path, values, status)
            import :: data_backend_t, WP, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            real(WP), allocatable, intent(out) :: values(:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_real_1d_if

        !> @brief Read a real 2-D dataset.
        subroutine read_real_2d_if(self, dataset_path, values, status)
            import :: data_backend_t, WP, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            real(WP), allocatable, intent(out) :: values(:,:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_real_2d_if

        !> @brief Read a real 3-D dataset.
        subroutine read_real_3d_if(self, dataset_path, values, status)
            import :: data_backend_t, WP, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            real(WP), allocatable, intent(out) :: values(:,:,:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_real_3d_if

        !> @brief Read a real 4-D dataset.
        subroutine read_real_4d_if(self, dataset_path, values, status)
            import :: data_backend_t, WP, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            real(WP), allocatable, intent(out) :: values(:,:,:,:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_real_4d_if

        !> @brief Read an integer 1-D dataset.
        subroutine read_int_1d_if(self, dataset_path, values, status)
            import :: data_backend_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            integer, allocatable, intent(out) :: values(:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_int_1d_if

        !> @brief Read an integer 2-D dataset.
        subroutine read_int_2d_if(self, dataset_path, values, status)
            import :: data_backend_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            integer, allocatable, intent(out) :: values(:,:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_int_2d_if

        !> @brief Read an integer 3-D dataset.
        subroutine read_int_3d_if(self, dataset_path, values, status)
            import :: data_backend_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            character(len=*), intent(in) :: dataset_path
            integer, allocatable, intent(out) :: values(:,:,:)
            type(backend_status_t), intent(out) :: status
        end subroutine read_int_3d_if

        !> @brief Read one dense spectral slice at fixed `(iz, iafe)`.
        !>
        !> @details
        !> Reads the hyperslab corresponding to
        !> `flux(:, iz, iafe, :, :)` into `slice%flux(lambda,logt,logg)`.
        subroutine read_spectral_slice_if(self, dataset, iz, iafe, slice, status)
            import :: data_backend_t, dataset_desc_t, spectral_slice_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            type(dataset_desc_t), intent(in) :: dataset
            integer, intent(in) :: iz
            integer, intent(in) :: iafe
            type(spectral_slice_t), intent(inout) :: slice
            type(backend_status_t), intent(out) :: status
        end subroutine read_spectral_slice_if

        !> @brief Read a dense spectral neighborhood in `(z, afe)`.
        !>
        !> @details
        !> Reads a bounded hyperslab block for interpolation neighborhoods,
        !> typically 2x2 in `(z, afe)`:
        !> `flux(:, iz_lo:iz_hi, iafe_lo:iafe_hi, :, :)`.
        !> Output is placed in `neighborhood%flux(lambda,z,afe,logt,logg)`.
        subroutine read_spectral_neighborhood_if(self, dataset, iz_lo, iz_hi, iafe_lo, iafe_hi, neighborhood, status)
            import :: data_backend_t, dataset_desc_t, spectral_grid_t, backend_status_t
            class(data_backend_t), intent(inout) :: self
            type(dataset_desc_t), intent(in) :: dataset
            integer, intent(in) :: iz_lo
            integer, intent(in) :: iz_hi
            integer, intent(in) :: iafe_lo
            integer, intent(in) :: iafe_hi
            type(spectral_grid_t), intent(inout) :: neighborhood
            type(backend_status_t), intent(out) :: status
        end subroutine read_spectral_neighborhood_if

    end interface

contains

    !> @brief True if backend operation succeeded.
    pure logical function backend_status_ok(status)
        type(backend_status_t), intent(in) :: status
        backend_status_ok = (status%code == 0)
    end function backend_status_ok

    !> @brief Reset a status object to "OK".
    subroutine backend_status_clear(self)
        class(backend_status_t), intent(inout) :: self
        self%code = 0
        if (allocated(self%message)) deallocate(self%message)
    end subroutine backend_status_clear

    !> @brief Set an error status.
    subroutine backend_status_set_error(self, code, message)
        class(backend_status_t), intent(inout) :: self
        integer, intent(in) :: code
        character(len=*), intent(in) :: message

        self%code = code
        if (allocated(self%message)) deallocate(self%message)
        self%message = trim(message)
    end subroutine backend_status_set_error

    !> @brief Set status to success.
    subroutine backend_status_set_ok(self)
        class(backend_status_t), intent(inout) :: self
        self%code = 0
        if (allocated(self%message)) deallocate(self%message)
        self%message = 'OK'
    end subroutine backend_status_set_ok

end module fsps_data_backend
