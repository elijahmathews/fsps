module fsps_omp_check
    !> @brief
    !> Debug checks of the OpenMP target data before the batch kernels (FSPS_OMP_CHECK_PRE=1).
    !>
    !> @details
    !> A kernel that touches memory never mapped to the device faults on GPUs with separate
    !> memory (an SQ memory violation on AMD) or, on some systems, silently reads host memory.
    !> With FSPS_OMP_CHECK_PRE set, the batch checks before its first kernel that (1) every array
    !> it passes to kernels is present on the device (omp_target_is_present), and (2) every
    !> context component that device code reads through the device copy of ctx (inside declare
    !> target routines, where no per-kernel map refreshes it) points there to the component's
    !> mapped device copy, i.e. the descriptor was attached; a failure stops with the name of the
    !> array instead of a memory fault. The component list is the set of ctx%state components
    !> referenced in device code (update it with new device code). No-ops without OMP_TARGET.

    use fsps_context_types, only: fsps_context_t
    use, intrinsic :: iso_c_binding, only: c_ptr, c_loc, c_associated, c_intptr_t
#ifdef OMP_TARGET
    use omp_lib, only: omp_target_is_present, omp_get_mapped_ptr, omp_get_default_device
#endif
    implicit none
    private

    public :: omp_check_enabled, omp_check_present, omp_check_context

    integer, parameter :: N_COMP = 33
    integer, save :: enabled = -1

contains

    !> True if FSPS_OMP_CHECK_PRE is set to a value other than 0 (read once).
    logical function omp_check_enabled()
        character(len=8) :: v
        integer :: st

        if (enabled < 0) then
            call get_environment_variable('FSPS_OMP_CHECK_PRE', v, status=st)
            enabled = 0
            if (st == 0 .and. len_trim(v) > 0 .and. v /= '0') enabled = 1
        end if
        omp_check_enabled = enabled == 1
    end function omp_check_enabled

    !> Stops with a message if the host address addr (loc() of array name, used in routine where)
    !> is not mapped to the default device.
    subroutine omp_check_present(addr, name, where, n)
        integer(8), intent(in)        :: addr
        character(len=*), intent(in)  :: name, where
        integer, intent(in), optional :: n   !< elements; zero-size arrays are not mapped
#ifdef OMP_TARGET
        type(c_ptr) :: p

        if (addr == 0) return
        if (present(n)) then
            if (n == 0) return
        end if
        p = transfer(int(addr, c_intptr_t), p)
        if (omp_target_is_present(p, omp_get_default_device()) == 0) then
            write(0, '(5A)') 'FSPS_OMP_CHECK_PRE: ', name, ' is not mapped to the device before the kernels of ', &
                where, ' (it would fault there)'
            error stop 'FSPS_OMP_CHECK_PRE: unmapped array'
        end if
#endif
    end subroutine omp_check_present

    !> The address that device code finds for component id of ctx%state, read through the device
    !> copy of ctx (as a declare target routine reads it).
    integer(8) function device_component_address(ctx, id) result(a)
#ifdef OMP_TARGET
        !$omp declare target
#endif
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in) :: id

        a = 0
        select case (id)
        case (1); a = loc(ctx%state%agb_logt_c)
        case (2); a = loc(ctx%state%agb_logt_o)
        case (3); a = loc(ctx%state%agb_spec_c)
        case (4); a = loc(ctx%state%agb_spec_c_smooth)
        case (5); a = loc(ctx%state%agb_spec_car)
        case (6); a = loc(ctx%state%agb_spec_car_smooth)
        case (7); a = loc(ctx%state%agb_spec_o)
        case (8); a = loc(ctx%state%agb_spec_o_smooth)
        case (9); a = loc(ctx%state%agndust_spec)
        case (10); a = loc(ctx%state%dustem2_dustem)
        case (11); a = loc(ctx%state%flux_dagb)
        case (12); a = loc(ctx%state%g03smcextn)
        case (13); a = loc(ctx%state%gaussnebarr)
        case (14); a = loc(ctx%state%indexdefined)
        case (15); a = loc(ctx%state%magvega)
        case (16); a = loc(ctx%state%nebem_cont)
        case (17); a = loc(ctx%state%pagb_spec)
        case (18); a = loc(ctx%state%sfh_w_tmp1)
        case (19); a = loc(ctx%state%sfh_w_tmp2)
        case (20); a = loc(ctx%state%spec_lambda)
        case (21); a = loc(ctx%state%spec_nu)
        case (22); a = loc(ctx%state%speclib)
        case (23); a = loc(ctx%state%speclib_agb_smooth)
        case (24); a = loc(ctx%state%ssp_basis_lbol)
        case (25); a = loc(ctx%state%ssp_basis_mass)
        case (26); a = loc(ctx%state%ssp_basis_spec)
        case (27); a = loc(ctx%state%time_full)
        case (28); a = loc(ctx%state%wgdust)
        case (29); a = loc(ctx%state%wmb_spec)
        case (30); a = loc(ctx%state%wrc_spec)
        case (31); a = loc(ctx%state%wrn_spec)
        case (32); a = loc(ctx%state%xnebem_cont)
        case (33); a = loc(ctx%state%zlegend)
        end select
    end function device_component_address

    !> Checks the device copy of ctx before the kernels of where: each component used by device
    !> code is mapped and attached (see the module description).
    subroutine omp_check_context(ctx, where)
        type(fsps_context_t), intent(in), target :: ctx
        character(len=*), intent(in)             :: where
#ifdef OMP_TARGET
        integer(8) :: dbase(N_COMP)
        integer :: i

        call omp_check_present(loc(ctx), 'ctx', where)
        dbase = 0
        !$omp target map(from: dbase)
        do i = 1, N_COMP
            dbase(i) = device_component_address(ctx, i)
        end do
        !$omp end target
        if (associated(ctx%state%agb_logt_c)) call check_attached(1, 'agb_logt_c', c_loc(ctx%state%agb_logt_c), where, dbase)
        if (associated(ctx%state%agb_logt_o)) call check_attached(2, 'agb_logt_o', c_loc(ctx%state%agb_logt_o), where, dbase)
        if (associated(ctx%state%agb_spec_c)) call check_attached(3, 'agb_spec_c', c_loc(ctx%state%agb_spec_c), where, dbase)
        if (associated(ctx%state%agb_spec_c_smooth)) call check_attached(4, 'agb_spec_c_smooth', c_loc(ctx%state%agb_spec_c_smooth), where, dbase)
        if (associated(ctx%state%agb_spec_car)) call check_attached(5, 'agb_spec_car', c_loc(ctx%state%agb_spec_car), where, dbase)
        if (associated(ctx%state%agb_spec_car_smooth)) call check_attached(6, 'agb_spec_car_smooth', c_loc(ctx%state%agb_spec_car_smooth), where, dbase)
        if (associated(ctx%state%agb_spec_o)) call check_attached(7, 'agb_spec_o', c_loc(ctx%state%agb_spec_o), where, dbase)
        if (associated(ctx%state%agb_spec_o_smooth)) call check_attached(8, 'agb_spec_o_smooth', c_loc(ctx%state%agb_spec_o_smooth), where, dbase)
        if (associated(ctx%state%agndust_spec)) call check_attached(9, 'agndust_spec', c_loc(ctx%state%agndust_spec), where, dbase)
        if (associated(ctx%state%dustem2_dustem)) call check_attached(10, 'dustem2_dustem', c_loc(ctx%state%dustem2_dustem), where, dbase)
        if (associated(ctx%state%flux_dagb)) call check_attached(11, 'flux_dagb', c_loc(ctx%state%flux_dagb), where, dbase)
        if (associated(ctx%state%g03smcextn)) call check_attached(12, 'g03smcextn', c_loc(ctx%state%g03smcextn), where, dbase)
        if (associated(ctx%state%gaussnebarr)) call check_attached(13, 'gaussnebarr', c_loc(ctx%state%gaussnebarr), where, dbase)
        if (associated(ctx%state%indexdefined)) call check_attached(14, 'indexdefined', c_loc(ctx%state%indexdefined), where, dbase)
        if (associated(ctx%state%magvega)) call check_attached(15, 'magvega', c_loc(ctx%state%magvega), where, dbase)
        if (associated(ctx%state%nebem_cont)) call check_attached(16, 'nebem_cont', c_loc(ctx%state%nebem_cont), where, dbase)
        if (associated(ctx%state%pagb_spec)) call check_attached(17, 'pagb_spec', c_loc(ctx%state%pagb_spec), where, dbase)
        if (allocated(ctx%state%sfh_w_tmp1)) call check_attached(18, 'sfh_w_tmp1', c_loc(ctx%state%sfh_w_tmp1), where, dbase)
        if (allocated(ctx%state%sfh_w_tmp2)) call check_attached(19, 'sfh_w_tmp2', c_loc(ctx%state%sfh_w_tmp2), where, dbase)
        if (associated(ctx%state%spec_lambda)) call check_attached(20, 'spec_lambda', c_loc(ctx%state%spec_lambda), where, dbase)
        if (associated(ctx%state%spec_nu)) call check_attached(21, 'spec_nu', c_loc(ctx%state%spec_nu), where, dbase)
        if (associated(ctx%state%speclib)) call check_attached(22, 'speclib', c_loc(ctx%state%speclib), where, dbase)
        if (associated(ctx%state%speclib_agb_smooth)) call check_attached(23, 'speclib_agb_smooth', c_loc(ctx%state%speclib_agb_smooth), where, dbase)
        if (allocated(ctx%state%ssp_basis_lbol)) call check_attached(24, 'ssp_basis_lbol', c_loc(ctx%state%ssp_basis_lbol), where, dbase)
        if (allocated(ctx%state%ssp_basis_mass)) call check_attached(25, 'ssp_basis_mass', c_loc(ctx%state%ssp_basis_mass), where, dbase)
        if (allocated(ctx%state%ssp_basis_spec)) call check_attached(26, 'ssp_basis_spec', c_loc(ctx%state%ssp_basis_spec), where, dbase)
        if (associated(ctx%state%time_full)) call check_attached(27, 'time_full', c_loc(ctx%state%time_full), where, dbase)
        if (associated(ctx%state%wgdust)) call check_attached(28, 'wgdust', c_loc(ctx%state%wgdust), where, dbase)
        if (associated(ctx%state%wmb_spec)) call check_attached(29, 'wmb_spec', c_loc(ctx%state%wmb_spec), where, dbase)
        if (associated(ctx%state%wrc_spec)) call check_attached(30, 'wrc_spec', c_loc(ctx%state%wrc_spec), where, dbase)
        if (associated(ctx%state%wrn_spec)) call check_attached(31, 'wrn_spec', c_loc(ctx%state%wrn_spec), where, dbase)
        if (associated(ctx%state%xnebem_cont)) call check_attached(32, 'xnebem_cont', c_loc(ctx%state%xnebem_cont), where, dbase)
        if (associated(ctx%state%zlegend)) call check_attached(33, 'zlegend', c_loc(ctx%state%zlegend), where, dbase)
#endif
    end subroutine omp_check_context

#ifdef OMP_TARGET
    subroutine check_attached(id, name, host, where, dbase)
        integer, intent(in)              :: id
        character(len=*), intent(in)     :: name, where
        type(c_ptr), intent(in)          :: host
        integer(8), intent(in)           :: dbase(:)
        type(c_ptr) :: dev

        call omp_check_present(int(transfer(host, 0_c_intptr_t), 8), 'ctx%state%' // name, where)
        dev = omp_get_mapped_ptr(host, omp_get_default_device())
        if (int(transfer(dev, 0_c_intptr_t), 8) /= dbase(id)) then
            write(0, '(5A,Z0,A,Z0,A,Z0,A)') 'FSPS_OMP_CHECK_PRE: device copy of ctx%state%', name, &
                ' is not attached before the kernels of ', where, ': device code reads 0x', dbase(id), &
                ', mapped data is at 0x', transfer(dev, 0_c_intptr_t), ' (host 0x', transfer(host, 0_c_intptr_t), ')'
            error stop 'FSPS_OMP_CHECK_PRE: unattached component'
        end if
    end subroutine check_attached
#endif

end module fsps_omp_check
