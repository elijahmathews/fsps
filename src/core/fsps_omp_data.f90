module fsps_omp_data
    !> @brief
    !> OpenMP target data of a host array given by its address and size in bytes.
    !>
    !> @details
    !> The data directives here name only a plain byte array (an explicit-shape dummy), never a
    !> component of a derived type. With LLVM flang, `target exit data map(delete: s%x)` also
    !> deletes the parent structure s: for s inside the context mapped as a whole, the whole
    !> context mapping disappeared (seen on the Radeon 780M, flang 22), and map(release:) instead
    !> left stale mappings. Mapping the bytes leaves parent structures alone with every compiler.
    !> Device code then reaches such data only through kernel arguments (a present lookup by
    !> address), not through the device copy of a structure. Callers pass loc(x) and
    !> OMP_BYTES(x), a macro defined where used: int(storage_size(x) / 8, 8) * size(x, kind=8).
    !> No-ops without OMP_TARGET.

    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_f_pointer
    implicit none
    private

    public :: omp_data_to, omp_data_alloc, omp_data_delete, omp_data_update_to, omp_data_update_from

contains

    !> target enter data map(to:) of nbytes bytes at host address addr.
    subroutine omp_data_to(addr, nbytes)
        integer(8), intent(in) :: addr, nbytes
#ifdef OMP_TARGET
        integer(1), pointer, contiguous :: b(:)

        if (addr == 0 .or. nbytes <= 0) return
        call c_f_pointer(transfer(addr, c_null_ptr), b, [nbytes])
        call bytes_to(b, nbytes)
#endif
    end subroutine omp_data_to

    !> target enter data map(alloc:) of nbytes bytes at host address addr.
    subroutine omp_data_alloc(addr, nbytes)
        integer(8), intent(in) :: addr, nbytes
#ifdef OMP_TARGET
        integer(1), pointer, contiguous :: b(:)

        if (addr == 0 .or. nbytes <= 0) return
        call c_f_pointer(transfer(addr, c_null_ptr), b, [nbytes])
        call bytes_alloc(b, nbytes)
#endif
    end subroutine omp_data_alloc

    !> target exit data map(delete:) of nbytes bytes at host address addr.
    subroutine omp_data_delete(addr, nbytes)
        integer(8), intent(in) :: addr, nbytes
#ifdef OMP_TARGET
        integer(1), pointer, contiguous :: b(:)

        if (addr == 0 .or. nbytes <= 0) return
        call c_f_pointer(transfer(addr, c_null_ptr), b, [nbytes])
        call bytes_delete(b, nbytes)
#endif
    end subroutine omp_data_delete

    !> target update to() of nbytes bytes at host address addr.
    subroutine omp_data_update_to(addr, nbytes)
        integer(8), intent(in) :: addr, nbytes
#ifdef OMP_TARGET
        integer(1), pointer, contiguous :: b(:)

        if (addr == 0 .or. nbytes <= 0) return
        call c_f_pointer(transfer(addr, c_null_ptr), b, [nbytes])
        call bytes_update_to(b, nbytes)
#endif
    end subroutine omp_data_update_to

    !> target update from() of nbytes bytes at host address addr.
    subroutine omp_data_update_from(addr, nbytes)
        integer(8), intent(in) :: addr, nbytes
#ifdef OMP_TARGET
        integer(1), pointer, contiguous :: b(:)

        if (addr == 0 .or. nbytes <= 0) return
        call c_f_pointer(transfer(addr, c_null_ptr), b, [nbytes])
        call bytes_update_from(b, nbytes)
#endif
    end subroutine omp_data_update_from

#ifdef OMP_TARGET
    ! The directives on an explicit-shape dummy: only the data is mapped (no descriptor of a
    ! pointer or assumed-shape dummy, which would leave an entry for a stack address)

    subroutine bytes_to(x, n)
        integer(8), intent(in) :: n
        integer(1), intent(in) :: x(n)
        !$omp target enter data map(to: x)
    end subroutine bytes_to

    subroutine bytes_alloc(x, n)
        integer(8), intent(in) :: n
        integer(1), intent(in) :: x(n)
        !$omp target enter data map(alloc: x)
    end subroutine bytes_alloc

    subroutine bytes_delete(x, n)
        integer(8), intent(in) :: n
        integer(1), intent(in) :: x(n)
        !$omp target exit data map(delete: x)
    end subroutine bytes_delete

    subroutine bytes_update_to(x, n)
        integer(8), intent(in) :: n
        integer(1), intent(in) :: x(n)
        !$omp target update to(x)
    end subroutine bytes_update_to

    subroutine bytes_update_from(x, n)
        integer(8), intent(in) :: n
        integer(1), intent(inout) :: x(n)
        !$omp target update from(x)
    end subroutine bytes_update_from
#endif

end module fsps_omp_data
