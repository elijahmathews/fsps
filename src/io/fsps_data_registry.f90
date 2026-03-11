module fsps_data_registry
    !> @brief Backend factory/registry for FSDS data access.
    !>
    !> @details
    !> FSPS is HDF5-only. This factory always returns the HDF5 backend.

    use fsps_data_backend, only: data_backend_t, backend_status_t
    use fsps_data_backend_hdf5, only: hdf5_backend_t

    implicit none
    private

    public :: create_data_backend

contains

    !> @brief Create the concrete HDF5 backend instance.
    !>
    !> @param[out] backend Allocatable polymorphic backend instance.
    !> @param[out] status  Factory status.
    subroutine create_data_backend(backend, status)
        class(data_backend_t), allocatable, intent(out) :: backend
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        if (allocated(backend)) deallocate(backend)

        allocate(hdf5_backend_t :: backend)
    end subroutine create_data_backend

end module fsps_data_registry
