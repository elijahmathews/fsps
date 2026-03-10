#include "fsps_build_config.h"

module fsps_data_registry
    !> @brief Backend factory/registry for FSDS data access.
    !>
    !> @details
    !> This module maps a user/backend mode request (`legacy`, `hdf5`, `auto`)
    !> to a concrete backend instance. Compile-time HDF5 support is controlled
    !> by Meson-generated preprocessor definitions in `fsps_build_config.h`.

    use fsps_strings, only: to_lower
    use fsps_data_backend, only: data_backend_t, backend_status_t
    use fsps_data_backend_legacy, only: legacy_backend_t
#if FSPS_HAS_HDF5 == 1
    use fsps_data_backend_hdf5, only: hdf5_backend_t
#endif

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
            allocate(legacy_backend_t :: backend)

        case ('hdf5', 'fsds_hdf5')
#if FSPS_HAS_HDF5 == 1
            allocate(hdf5_backend_t :: backend)
#else
            call status%set_error(4001, 'HDF5 backend requested, but FSPS was compiled without HDF5 support.')
#endif

        case ('auto', 'fsds_auto')
#if FSPS_HAS_HDF5 == 1
            allocate(hdf5_backend_t :: backend)
#else
            allocate(legacy_backend_t :: backend)
#endif

        case default
            call status%set_error(4002, 'Unknown data backend mode: '//trim(mode))
        end select
    end subroutine create_data_backend

end module fsps_data_registry
