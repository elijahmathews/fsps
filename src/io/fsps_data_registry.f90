module fsps_data_registry
    !> @brief Backend factory/registry for FSDS data access.
    !>
    !> @details
    !> This module maps a user/backend mode request (`legacy`, `hdf5`, `auto`)
    !> to a concrete backend instance. Legacy backend support has been removed;
    !> only the HDF5 backend is available.

    use fsps_strings, only: to_lower
    use fsps_data_backend, only: data_backend_t, backend_status_t
    use fsps_data_backend_hdf5, only: hdf5_backend_t

    implicit none
    private

    public :: create_data_backend

contains

    !> @brief Create a concrete backend instance based on mode.
    !>
    !> @param[in]  mode    Requested mode: `legacy`, `hdf5`, or `auto`.
    !> @param[out] backend Allocatable polymorphic backend instance.
    !> @param[out] status  Factory status.
    subroutine create_data_backend(mode, backend, status)
        character(len=*), intent(in) :: mode
        class(data_backend_t), allocatable, intent(out) :: backend
        type(backend_status_t), intent(out) :: status

        character(len=:), allocatable :: mode_l

        call status%set_ok()

        if (allocated(backend)) deallocate(backend)

        mode_l = trim(to_lower(trim(mode)))

        select case (mode_l)
        case ('legacy', 'fsds_legacy')
            call status%set_error(4003, 'Legacy backend has been removed. Please use FSDS HDF5 format.')

        case ('hdf5', 'fsds_hdf5')
            allocate(hdf5_backend_t :: backend)

        case ('auto', 'fsds_auto')
            allocate(hdf5_backend_t :: backend)

        case default
            call status%set_error(4002, 'Unknown data backend mode: '//trim(mode))
        end select
    end subroutine create_data_backend

end module fsps_data_registry
