module test_fsps_dust_mod
    use fsps_types, only: sp, params
    use fsps_context_types, only: fsps_context_t
    use fsps_dust
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_float_equals, assert_true, assert_int_equals
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none

    ! Global failure counter
    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_dust_tests, total_failures, total_tests

    ! Local constants for expected-value calculations
    real(sp), parameter :: V_BAND_ANGSTROMS = 5500.0_sp
    real(sp), parameter :: UV_BUMP_CENTER = 2175.0_sp
    real(sp), parameter :: KC13_BUMP_WIDTH = 350.0_sp
    real(sp), parameter :: KC13_R_V_BASE = 4.05_sp

    real(sp), parameter :: CALZ_LAM_UV_MIN = 1200.0_sp
    real(sp), parameter :: CALZ_LAM_BREAK  = 6300.0_sp
    real(sp), parameter :: CALZ_LAM_IR_MAX = 22000.0_sp
    real(sp), parameter :: CALZ_R_V   = 4.05_sp
    real(sp), parameter :: CALZ_SCALE = 2.659_sp
    real(sp), parameter :: CALZ_UV_COEFFS(0:3) = [-2.156_sp, 1.509_sp, -0.198_sp, 0.011_sp]
    real(sp), parameter :: CALZ_OPT_COEFFS(0:1) = [-1.857_sp, 1.040_sp]

    real(sp), parameter :: EPS = 1.0e-6_sp
    real(sp), parameter :: SMALL_DELTA = 1.0e-4_sp
    real(sp), parameter :: LESS_SMALL_DELTA = 1.0e-2_sp

contains

    !> @brief
    !> Unit test suite for fsps_dust module.
    subroutine run_fsps_dust_tests()

        call print_minor_header("fsps_dust")

        call test_power_law()
        call test_ccm89()
        call test_calzetti()
        call test_wg_smc()
        call test_kriek_conroy()
        call test_reddy()
        call test_robustness()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)

    end subroutine run_fsps_dust_tests

    ! ------------------------------------------------------------------------
    ! TEST SUITE: POWER LAW (Type 0)
    ! ------------------------------------------------------------------------
    subroutine test_power_law()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(2) :: w2
        real(sp), dimension(3) :: w3
        real(sp), dimension(1) :: w1
        real(sp), dimension(2) :: res2
        real(sp), dimension(3) :: res3
        real(sp), dimension(1) :: res1

        call print_group("Power Law (Type 0)")

        allocate(ctx)

        ! Test 1.1: V-band normalization
        settings%dust_index = -0.7_sp
        w1 = [5500.0_sp]
        res1 = compute_attenuation_curve(w1, 0, settings, ctx)
        call assert_float_equals(1.0_sp, res1(1), 0.0_sp, "V-band normalization (5500A)", total_tests, total_failures)

        ! Test 1.2: Grey screen (index=0)
        settings%dust_index = 0.0_sp
        w3 = [1000.0_sp, 2000.0_sp, 5000.0_sp]
        res3 = compute_attenuation_curve(w3, 0, settings, ctx)
        call assert_true(all(abs(res3 - 1.0_sp) <= EPS), "Grey screen (index=0)", total_tests, total_failures)

        ! Test 1.3: Blue-tilted slope
        settings%dust_index = -1.0_sp
        w2 = [2750.0_sp, 5500.0_sp]
        res2 = compute_attenuation_curve(w2, 0, settings, ctx)
        call assert_float_equals(2.0_sp, res2(1), EPS, "Blue-tilted slope (2750A)", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_power_law

    ! ------------------------------------------------------------------------
    ! TEST SUITE: CCM89 Milky Way (Type 1)
    ! ------------------------------------------------------------------------
    subroutine test_ccm89()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(5) :: w_regions
        real(sp), dimension(5) :: res_regions
        real(sp), dimension(2) :: w_break
        real(sp), dimension(2) :: res_break
        real(sp), dimension(2) :: w_fuv
        real(sp), dimension(2) :: res_fuv
        real(sp), dimension(1) :: w_bump
        real(sp), dimension(1) :: res_uvb0, res_uvb1

        call print_group("CCM89 Milky Way (Type 1)")

        allocate(ctx)

        settings%mwr = 3.1_sp
        settings%uvb = 1.0_sp

        ! Test 2.1: Region coverage (IR, Optical, NUV, MUV, FUV)
        w_regions = [1.0e4_sp/0.7_sp, 1.0e4_sp/2.0_sp, 1.0e4_sp/4.0_sp, &
                     1.0e4_sp/6.5_sp, 1.0e4_sp/9.0_sp]
        res_regions = compute_attenuation_curve(w_regions, 1, settings, ctx)
        call assert_true(all(ieee_is_finite(res_regions)) .and. all(res_regions > 0.0_sp), &
                         "Region coverage (finite and > 0)", total_tests, total_failures)

        ! Test 2.2: Smoothing hack continuity at 3030A
        w_break = [3030.30_sp, 3030.31_sp]
        res_break = compute_attenuation_curve(w_break, 1, settings, ctx)
        call assert_true(abs(res_break(1) - res_break(2)) <= SMALL_DELTA, &
                         "Smoothing continuity (3030A)", total_tests, total_failures)

        ! Test 2.3: Far-UV clamping (lambda < 833A)
        w_fuv = [800.0_sp, 700.0_sp]
        res_fuv = compute_attenuation_curve(w_fuv, 1, settings, ctx)
        call assert_float_equals(res_fuv(1), res_fuv(2), 0.0_sp, "FUV clamp equality", total_tests, total_failures)

        ! Test 2.4: UV bump scaling
        w_bump = [2175.0_sp]
        settings%uvb = 0.0_sp
        res_uvb0 = compute_attenuation_curve(w_bump, 1, settings, ctx)
        settings%uvb = 1.0_sp
        res_uvb1 = compute_attenuation_curve(w_bump, 1, settings, ctx)
        call assert_true(res_uvb1(1) > res_uvb0(1), "UV bump scaling (uvb 0 -> 1)", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_ccm89

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Calzetti et al. (Type 2)
    ! ------------------------------------------------------------------------
    subroutine test_calzetti()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(2) :: w_break
        real(sp), dimension(2) :: res_break
        real(sp), dimension(1) :: w_v
        real(sp), dimension(1) :: res_v

        call print_group("Calzetti (Type 2)")

        allocate(ctx)

        ! Test 3.1: 0.63um junction continuity
        w_break = [6300.000_sp, 6300.001_sp]
        res_break = compute_attenuation_curve(w_break, 2, settings, ctx)
        call assert_true(abs(res_break(1) - res_break(2)) <= LESS_SMALL_DELTA, &
                         "0.63um junction continuity", total_tests, total_failures)

        ! Test 3.2: Precision fix at 5500A
        w_v = [5500.0_sp]
        res_v = compute_attenuation_curve(w_v, 2, settings, ctx)
        call assert_float_equals(1.0_sp, res_v(1), 1.0e-3_sp, "Calzetti V-band normalization", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_calzetti

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Witt & Gordon (Type 3) and SMC (Type 5)
    ! ------------------------------------------------------------------------
    subroutine test_wg_smc()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: wavelengths
        real(sp), dimension(nlam) :: res

        call print_group("Witt & Gordon + SMC Tables (Types 3 & 5)")

        wavelengths = [1000.0_sp, 2000.0_sp, 3000.0_sp, 4000.0_sp, 5000.0_sp]
        call setup_mock_context(ctx, nlam)

        ! Test 4.1: W&G table mapping (index order)
        settings%wgp1 = 2
        settings%wgp2 = 3
        settings%wgp3 = 1
        res = compute_attenuation_curve(wavelengths, 3, settings, ctx)
        call assert_true(all(res == ctx%state%wgdust(:, 2, 3, 1)), "W&G table mapping", total_tests, total_failures)

        ! Test 4.2: SMC passthrough
        res = compute_attenuation_curve(wavelengths, 5, settings, ctx)
        call assert_true(all(res == ctx%state%g03smcextn), "SMC passthrough", total_tests, total_failures)

        call teardown_mock_context(ctx)

    end subroutine test_wg_smc

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Kriek & Conroy (Type 4)
    ! ------------------------------------------------------------------------
    subroutine test_kriek_conroy()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(3) :: wavelengths
        real(sp), dimension(3) :: res, expected
        real(sp), dimension(1) :: w_bump
        real(sp), dimension(1) :: res_bump, expected_no_bump
        real(sp), dimension(3) :: res_uvb0, res_uvb5

        call print_group("Kriek & Conroy (Type 4)")

        allocate(ctx)

        wavelengths = [1500.0_sp, 2175.0_sp, 5500.0_sp]

        ! Test 5.1: No-tilt baseline (delta=0)
        settings%dust_index = 0.0_sp
        res = compute_attenuation_curve(wavelengths, 4, settings, ctx)
        expected = (calc_calzetti_curve(wavelengths) * KC13_R_V_BASE + &
                    calc_drude_profile(wavelengths, 0.85_sp)) / KC13_R_V_BASE
        call assert_true(all(abs(res - expected) <= EPS), "No-tilt baseline (delta=0)", total_tests, total_failures)

        ! Test 5.2: No-bump index (delta=0.447)
        settings%dust_index = 0.447_sp
        w_bump = [2175.0_sp]
        res_bump = compute_attenuation_curve(w_bump, 4, settings, ctx)
        expected_no_bump = calc_calzetti_curve(w_bump) * (w_bump / V_BAND_ANGSTROMS)**settings%dust_index
        call assert_float_equals(expected_no_bump(1), res_bump(1), 5.0e-4_sp, &
                                 "No-bump index (delta=0.447)", total_tests, total_failures)

        ! Test 5.3: UVB independence (ghost variable)
        settings%dust_index = -0.2_sp
        settings%uvb = 0.0_sp
        res_uvb0 = compute_attenuation_curve(wavelengths, 4, settings, ctx)
        settings%uvb = 5.0_sp
        res_uvb5 = compute_attenuation_curve(wavelengths, 4, settings, ctx)
        call assert_true(all(res_uvb0 == res_uvb5), "UVB independence", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_kriek_conroy

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Reddy et al. (Type 6)
    ! ------------------------------------------------------------------------
    subroutine test_reddy()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(2) :: w_uv
        real(sp), dimension(2) :: res_uv
        real(sp), dimension(2) :: w_break
        real(sp), dimension(2) :: res_break

        call print_group("Reddy et al. (Type 6)")

        allocate(ctx)

        ! Test 6.1: UV extrapolation (lambda < 1500A constant)
        w_uv = [1400.0_sp, 1500.0_sp]
        res_uv = compute_attenuation_curve(w_uv, 6, settings, ctx)
        call assert_float_equals(res_uv(1), res_uv(2), 0.0_sp, "UV extrapolation clamp", total_tests, total_failures)

        ! Test 6.2: 6000A continuity
        w_break = [5999.999_sp, 6000.000_sp]
        res_break = compute_attenuation_curve(w_break, 6, settings, ctx)
        call assert_true(abs(res_break(1) - res_break(2)) <= SMALL_DELTA, "6000A continuity", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_reddy

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Robustness
    ! ------------------------------------------------------------------------
    subroutine test_robustness()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(:), allocatable :: w
        real(sp), dimension(:), allocatable :: res
        integer :: n
        integer :: i

        call print_group("Robustness")

        allocate(ctx)

        ! Test 7.1: Invalid ID returns NaN
        allocate(w(3))
        w = [1000.0_sp, 2000.0_sp, 3000.0_sp]
        res = compute_attenuation_curve(w, 99, settings, ctx)
        call assert_true(all(ieee_is_nan(res)), "Invalid ID returns NaN", total_tests, total_failures)
        deallocate(w, res)

        ! Test 7.2: Array shape conformance
        n = 1
        allocate(w(n), res(n))
        w = [5500.0_sp]
        res = compute_attenuation_curve(w, 0, settings, ctx)
        call assert_int_equals(size(w), size(res), "Array shape (size=1)", total_tests, total_failures)
        deallocate(w, res)

        n = 10
        allocate(w(n), res(n))
        w = [(real(i, sp), i=1, n)]
        res = compute_attenuation_curve(w, 0, settings, ctx)
        call assert_int_equals(size(w), size(res), "Array shape (size=10)", total_tests, total_failures)
        deallocate(w, res)

        n = 1000
        allocate(w(n), res(n))
        w = [(real(i, sp), i=1, n)]
        res = compute_attenuation_curve(w, 0, settings, ctx)
        call assert_int_equals(size(w), size(res), "Array shape (size=1000)", total_tests, total_failures)
        deallocate(w, res)

        deallocate(ctx)

    end subroutine test_robustness

    ! ------------------------------------------------------------------------
    ! Helper: mock context setup
    ! ------------------------------------------------------------------------
    subroutine setup_mock_context(ctx, nlam)
        type(fsps_context_t), allocatable, intent(out) :: ctx
        integer, intent(in) :: nlam
        integer :: i, j, k, l

        allocate(ctx)
        allocate(ctx%state%wgdust(nlam, 3, 4, 2))
        allocate(ctx%state%g03smcextn(nlam))

        do l = 1, 2
            do k = 1, 4
                do j = 1, 3
                    do i = 1, nlam
                        ctx%state%wgdust(i, j, k, l) = real(i + 100*j + 10000*k + 1000000*l, sp)
                    end do
                end do
            end do
        end do

        do i = 1, nlam
            ctx%state%g03smcextn(i) = real(i, sp)
        end do

    end subroutine setup_mock_context

    subroutine teardown_mock_context(ctx)
        type(fsps_context_t), allocatable, intent(inout) :: ctx

        if (associated(ctx%state%wgdust)) deallocate(ctx%state%wgdust)
        if (associated(ctx%state%g03smcextn)) deallocate(ctx%state%g03smcextn)
        deallocate(ctx)

    end subroutine teardown_mock_context

    ! ------------------------------------------------------------------------
    ! Helper: Calzetti curve (for expected values)
    ! ------------------------------------------------------------------------
    pure function calc_calzetti_curve(wavelengths) result(curve)
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), dimension(size(wavelengths)) :: curve
        real(sp), dimension(size(wavelengths)) :: wavenumbers, extinction_k

        wavenumbers = 1.0e4_sp / wavelengths
        extinction_k = 0.0_sp

        where (wavelengths > CALZ_LAM_BREAK .and. wavelengths <= CALZ_LAM_IR_MAX)
            extinction_k = CALZ_SCALE * (CALZ_OPT_COEFFS(0) + CALZ_OPT_COEFFS(1) * wavenumbers) + CALZ_R_V
        end where

        where (wavelengths >= CALZ_LAM_UV_MIN .and. wavelengths <= CALZ_LAM_BREAK)
            extinction_k = CALZ_R_V + CALZ_SCALE * ( &
                CALZ_UV_COEFFS(0) + wavenumbers * ( &
                    CALZ_UV_COEFFS(1) + wavenumbers * ( &
                        CALZ_UV_COEFFS(2) + wavenumbers * CALZ_UV_COEFFS(3) &
                    ) &
                ) &
            )
        end where

        curve = extinction_k / CALZ_R_V

    end function calc_calzetti_curve

    ! ------------------------------------------------------------------------
    ! Helper: Drude bump profile (for expected values)
    ! ------------------------------------------------------------------------
    pure function calc_drude_profile(wavelengths, amplitude) result(profile)
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), intent(in) :: amplitude
        real(sp), dimension(size(wavelengths)) :: profile

        profile = amplitude * (wavelengths * KC13_BUMP_WIDTH)**2 / &
                  ((wavelengths**2 - UV_BUMP_CENTER**2)**2 + (wavelengths * KC13_BUMP_WIDTH)**2)

    end function calc_drude_profile

end module test_fsps_dust_mod
