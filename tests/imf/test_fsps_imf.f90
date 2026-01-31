module test_fsps_imf_mod
    use fsps_types, only: sp, &
                          chab_mc, chab_sigma2, chab_ind, &
                          vd_sigma2, vd_ah, vd_ind, vd_al, vd_nc
    use fsps_context_types, only: fsps_context_t
    use fsps_imf
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_float_equals, assert_true
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    implicit none

    ! Global failure counter
    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_imf_tests, total_failures, total_tests

contains

    !> @brief
    !> Unit test suite for fsps_imf module.
    !> Designed for CI usage: exits with status 0 on success, 1 on failure.
    subroutine run_fsps_imf_tests()

        call print_minor_header("fsps_imf")

        call test_salpeter()
        call test_chabrier()
        call test_kroupa()
        call test_imf_weighting()
        call test_compute_weights_normalization()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)

    end subroutine run_fsps_imf_tests

    ! ------------------------------------------------------------------------
    ! TEST SUITE: SALPETER IMF (Type 0)
    ! ------------------------------------------------------------------------
    subroutine test_salpeter()
        type(fsps_context_t), allocatable :: ctx
        real(sp), dimension(3) :: m = [0.5_sp, 1.0_sp, 2.0_sp]
        real(sp), dimension(3) :: res
        real(sp) :: expected
        
        call print_group("Salpeter IMF")

        allocate(ctx)
        ! Setup Context
        ctx%imf_type_val = 0
        ctx%state%salp_ind = 2.35_sp ! Standard Salpeter

        ! Test Number Density (dN/dM)
        res = get_imf_value(ctx, m, mass_weighted=.false.)
        
        ! Check m=1.0 (should be 1.0^-2.35 = 1.0)
        call assert_float_equals(1.0_sp, res(2), 1.0e-6_sp, "Salpeter(1.0)", total_tests, total_failures)
        
        ! Check m=2.0 (should be 2.0^-2.35)
        expected = 2.0_sp**(-2.35_sp)
        call assert_float_equals(expected, res(3), 1.0e-6_sp, "Salpeter(2.0)", total_tests, total_failures)

        ! Test Mass Weighting (M * dN/dM)
        res = get_imf_value(ctx, m, mass_weighted=.true.)
        expected = 2.0_sp * (2.0_sp**(-2.35_sp))
        call assert_float_equals(expected, res(3), 1.0e-6_sp, "Salpeter Mass-Weighted", total_tests, total_failures)

        deallocate(ctx)
    end subroutine test_salpeter

    ! ------------------------------------------------------------------------
    ! TEST SUITE: CHABRIER IMF (Type 1)
    ! ------------------------------------------------------------------------
    subroutine test_chabrier()
        type(fsps_context_t), allocatable :: ctx
        real(sp), dimension(2) :: m = [0.1_sp, 2.0_sp] ! Below and above 1.0 break
        real(sp), dimension(2) :: res
        real(sp) :: log_m, log_mc, term, expected
        
        call print_group("Chabrier IMF")

        allocate(ctx)
        ctx%imf_type_val = 1
        
        ! Note: chab_mc, chab_sigma2, chab_ind are from fsps_types (constants)

        res = get_imf_value(ctx, m, mass_weighted=.false.)

        ! Case 1: Low Mass (Log Normal) at m=0.1
        log_m = log10(0.1_sp)
        log_mc = log10(chab_mc)
        term = (log_m - log_mc)**2 / (2.0_sp * chab_sigma2)
        ! Formula: exp(-term) / m
        expected = exp(-term) / 0.1_sp
        call assert_float_equals(expected, res(1), 1.0e-5_sp, "Chabrier Low Mass (<1)", total_tests, total_failures)

        ! Case 2: High Mass (Power Law) at m=2.0
        ! Formula: exp(-term_at_mc) * m^(-ind) / m  (Continuity matching)
        term = log_mc**2 / (2.0_sp * chab_sigma2)
        expected = exp(-term) * 2.0_sp**(-chab_ind) / 2.0_sp
        call assert_float_equals(expected, res(2), 1.0e-5_sp, "Chabrier High Mass (>1)", total_tests, total_failures)

        deallocate(ctx)
    end subroutine test_chabrier

    ! ------------------------------------------------------------------------
    ! TEST SUITE: KROUPA IMF (Type 2)
    ! ------------------------------------------------------------------------
    subroutine test_kroupa()
        type(fsps_context_t), allocatable :: ctx
        real(sp), dimension(3) :: m = [0.2_sp, 0.7_sp, 2.0_sp] ! One in each regime
        real(sp), dimension(3) :: res
        real(sp) :: expected
        real(sp), dimension(3) :: alpha
        
        call print_group("Kroupa IMF")

        allocate(ctx)
        ctx%imf_type_val = 2
        ! Standard Kroupa exponents
        alpha = [1.3_sp, 2.3_sp, 2.3_sp]
        ctx%state%imf_alpha = alpha

        res = get_imf_value(ctx, m, mass_weighted=.false.)

        ! Regime 1: 0.08 <= m < 0.5 (m=0.2) -> m^-1.3
        expected = 0.2_sp**(-1.3_sp)
        call assert_float_equals(expected, res(1), 1.0e-5_sp, "Kroupa Low (0.2)", total_tests, total_failures)

        ! Regime 2: 0.5 <= m < 1.0 (m=0.7) -> C1 * m^-2.3
        ! C1 = 0.5^(-1.3 + 2.3) = 0.5^1.0 = 0.5
        expected = 0.5_sp * 0.7_sp**(-2.3_sp)
        call assert_float_equals(expected, res(2), 1.0e-5_sp, "Kroupa Mid (0.7)", total_tests, total_failures)

        ! Regime 3: m >= 1.0 (m=2.0) -> C2 * m^-2.3
        ! C2 = 0.5^(-1.3 + 2.3) = 0.5
        expected = 0.5_sp * 2.0_sp**(-2.3_sp)
        call assert_float_equals(expected, res(3), 1.0e-5_sp, "Kroupa High (2.0)", total_tests, total_failures)

        deallocate(ctx)
    end subroutine test_kroupa

    ! ------------------------------------------------------------------------
    ! TEST SUITE: COMPUTE_IMF_WEIGHTS (Normalization)
    ! ------------------------------------------------------------------------
    subroutine test_compute_weights_normalization()
        type(fsps_context_t), allocatable :: ctx
        real(sp), allocatable :: mini(:), weights(:)
        integer :: nmass
        
        call print_group("compute_imf_weights (Normalization)")

        allocate(ctx)
        
        ! Setup a simple Salpeter context
        ctx%imf_type_val = 0
        ctx%state%salp_ind = 2.35_sp
        
        ! Define limits covering the isochrone points
        ctx%state%imf_lower_limit = 0.1_sp
        ctx%state%imf_upper_limit = 100.0_sp
        ctx%state%imf_lower_bound = 0.1_sp

        ! Create a mock isochrone mass grid (log-spaced)
        nmass = 10
        allocate(mini(nmass))
        allocate(weights(nmass))
        
        ! Fill mini with values from 0.5 to 10.0
        ! (Well within the limits to avoid edge effect complexity in this test)
        mini = [0.5, 1.0, 1.5, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0, 80.0]

        ! Run computation
        call compute_imf_weights(ctx, mini, weights, nmass)
        
        call assert_true(all(weights > 0.0_sp), "All weights positive", total_tests, total_failures)
        call assert_true(weights(1) > weights(nmass), "Salpeter decreases with mass", total_tests, total_failures)

        deallocate(ctx, mini, weights)
    end subroutine test_compute_weights_normalization

    ! ------------------------------------------------------------------------
    ! TEST SUITE: MASS WEIGHTING FLAG
    ! ------------------------------------------------------------------------
    subroutine test_imf_weighting()
        type(fsps_context_t), allocatable :: ctx
        real(sp), dimension(1) :: m = [2.0_sp]
        real(sp), dimension(1) :: val_num, val_mass
        
        call print_group("Mass Weighting Flag")

        allocate(ctx)
        ctx%imf_type_val = 0 ! Salpeter
        ctx%state%salp_ind = 2.0_sp ! Simple square law

        ! Get dN/dM
        val_num = get_imf_value(ctx, m, mass_weighted=.false.)
        ! Get M * dN/dM
        val_mass = get_imf_value(ctx, m, mass_weighted=.true.)

        ! Check relationship
        call assert_float_equals(val_num(1) * m(1), val_mass(1), 1.0e-6_sp, &
                                 "Mass Weighted = Mass * Number Weighted", &
                                 total_tests, total_failures)

        deallocate(ctx)
    end subroutine test_imf_weighting

end module test_fsps_imf_mod