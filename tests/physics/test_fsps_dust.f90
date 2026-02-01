module test_fsps_dust_mod
    use fsps_types, only: sp, params, nemline, clight, gsig4pi, nteff_dagb, &
                          ntau_dagb, nagndust, msun, newton, rsun, yr2sc
    use fsps_context_types, only: fsps_context_t
    use fsps_dust
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: find_interval
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_float_equals, assert_true, assert_int_equals, assert_relative_error
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none

    ! Global failure counter
    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_dust_tests, total_failures, total_tests

    ! Local constants for expected-value calculations
    real(sp), parameter :: SAFE_FLOOR = tiny(0.0_sp)
    real(sp), parameter :: PI = acos(-1.0_sp)
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

    ! AGB / VW93 constants (copied from fsps_dust for test expectations)
    real(sp), parameter :: AGB_DELTA_C_RICH = 0.0025_sp
    real(sp), parameter :: AGB_DELTA_O_RICH = 0.01_sp
    real(sp), parameter :: AGB_KAPPA_C_RICH = 3200.0_sp
    real(sp), parameter :: AGB_KAPPA_O_RICH = 3000.0_sp
    real(sp), parameter :: AGB_RIN_FACTOR_C = 1.92E12_sp
    real(sp), parameter :: AGB_RIN_FACTOR_O = 4.74E12_sp
    real(sp), parameter :: AGB_PER_INTERCEPT = -2.07_sp
    real(sp), parameter :: AGB_PER_SLOPE_R   = 1.94_sp
    real(sp), parameter :: AGB_PER_SLOPE_M   = -0.9_sp
    real(sp), parameter :: AGB_VEXP_INTERCEPT = -13.5_sp
    real(sp), parameter :: AGB_VEXP_SLOPE     = 0.056_sp
    real(sp), parameter :: AGB_VEXP_MIN       = 3.0_sp
    real(sp), parameter :: AGB_VEXP_MAX       = 15.0_sp
    real(sp), parameter :: VW93_MDOT_LIMIT_ISO = 1.0e-4_sp
    real(sp), parameter :: VW93_PER_THRESH     = 500.0_sp
    real(sp), parameter :: VW93_MASS_THRESH    = 2.5_sp
    real(sp), parameter :: VW93_BASE_INTERCEPT = -11.4_sp
    real(sp), parameter :: VW93_BASE_SLOPE     = 0.0123_sp
    real(sp), parameter :: VW93_HIGH_SLOPE     = 0.0125_sp
    real(sp), parameter :: VW93_SUPERWIND_NORM = 1.93e3_sp
    real(sp), parameter :: AGB_DTG_VEL_NORM    = 225.0_sp
    real(sp), parameter :: AGB_DTG_LUM_NORM    = 1.0e4_sp
    real(sp), parameter :: AGB_DTG_LUM_EXP     = -0.6_sp

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
        call test_agn_dust_emission()
        call test_circumstellar_optical_depth()
        call test_agb_dust_screen()
        call test_energy_conservation()
        call test_differential_attenuation()
        call test_draine_li_energy_balance()
        call test_dust_self_absorption()
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
        call assert_float_equals(1.0_sp, res_v(1), 1.0e-3_sp, &
                                 "Calzetti V-band normalization", total_tests, total_failures)

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
    ! TEST SUITE: AGN Dust Emission
    ! ------------------------------------------------------------------------
    subroutine test_agn_dust_emission()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: wavelengths
        real(sp), dimension(nlam) :: spectrum, expected

        call print_group("AGN Dust Emission")

        call setup_physics_context(ctx, nlam)
        wavelengths = ctx%state%spec_lambda

        ctx%dust_type_val = 0
        settings%dust_index = 0.0_sp

        ! Test 1.1: Early exit when fagn=0
        settings%fagn = 0.0_sp
        spectrum = 100.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 0.0_sp, spectrum)
        call assert_true(all(abs(spectrum - 100.0_sp) <= EPS), "AGN early exit (fagn=0)", total_tests, total_failures)

        ! Test 1.2: Torus interpolation exact grid point
        settings%fagn = 1.0_sp
        settings%dust2 = 0.0_sp
        settings%agn_tau = 20.0_sp
        spectrum = 0.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 0.0_sp, spectrum)
        expected = 2.0_sp
        call assert_true(all(abs(spectrum - expected) <= EPS), "AGN grid exact (tau=20)", total_tests, total_failures)

        ! Test 1.3: Torus interpolation linear
        settings%agn_tau = 15.0_sp
        spectrum = 0.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 0.0_sp, spectrum)
        expected = 1.5_sp
        call assert_true(all(abs(spectrum - expected) <= EPS), "AGN grid linear (tau=15)", total_tests, total_failures)

        ! Test 1.4: Galaxy attenuation (power-law, dust2=1)
        settings%agn_tau = 10.0_sp
        settings%dust2 = 1.0_sp
        spectrum = 0.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 0.0_sp, spectrum)
        expected = exp(-1.0_sp)
        call assert_true(all(abs(spectrum - expected) <= EPS), "AGN host attenuation", total_tests, total_failures)

        ! Test 1.5: Type 3 exception (no dust2 scaling)
        ctx%dust_type_val = 3
        settings%dust2 = 50.0_sp
        spectrum = 0.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 0.0_sp, spectrum)
        expected = exp(-1.0_sp)
        call assert_true(all(abs(spectrum - expected) <= EPS), "AGN Type 3 exception", total_tests, total_failures)

        ! Test 1.6: Luminosity normalization
        ctx%dust_type_val = 0
        settings%dust2 = 0.0_sp
        settings%agn_tau = 10.0_sp
        settings%fagn = 0.1_sp
        spectrum = 0.0_sp
        call apply_agn_dust_emission(ctx, settings, wavelengths, 10.0_sp, spectrum)
        expected = 1.0e9_sp
        call assert_true(all(abs(spectrum - expected) <= EPS), "AGN luminosity normalization", total_tests, total_failures)

        call teardown_physics_context(ctx)

    end subroutine test_agn_dust_emission

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Circumstellar Optical Depth
    ! ------------------------------------------------------------------------
    subroutine test_circumstellar_optical_depth()
        type(fsps_context_t), allocatable :: ctx
        real(sp) :: tau_c, tau_o, tau_expected_c, tau_expected_o
        real(sp) :: m_act, log_l, log_g, log_mdot_iso
        real(sp) :: period_days, vexp_unclamped
        real(sp) :: tau_clamped, tau_unclamped

        call print_group("Circumstellar Optical Depth")

        allocate(ctx)
        ctx%use_isoc_mdot_val = 1

        m_act = 1.0_sp
        log_l = 3.0_sp
        log_g = 1.0_sp
        log_mdot_iso = -5.0_sp

        ! Test 2.1: Chemistry switching
        tau_c = compute_circumstellar_optical_depth(ctx, 1, m_act, log_l, log_g, log_mdot_iso)
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        tau_expected_c = calc_tau_expected(1, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .true.)
        tau_expected_o = calc_tau_expected(0, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .true.)
        call assert_float_equals(tau_expected_c, tau_c, 1.0e-8_sp, "C-rich kappa path", total_tests, total_failures)
        call assert_float_equals(tau_expected_o, tau_o, 1.0e-8_sp, "O-rich kappa path", total_tests, total_failures)

        ! Test 2.2: Expansion velocity clamping (low)
        log_g = 6.0_sp
        period_days = calc_period_days(m_act, log_g)
        vexp_unclamped = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        tau_clamped = calc_tau_expected(0, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .true.)
        tau_unclamped = calc_tau_expected(0, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .false.)
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        
        call assert_true(vexp_unclamped < AGB_VEXP_MIN, "Velocity unclamped < 3 km/s", total_tests, total_failures)
        call assert_float_equals(tau_clamped, tau_o, 1.0e-8_sp, "Velocity clamped low", total_tests, total_failures)
        call assert_true(abs(tau_unclamped - tau_o) > 0.0_sp, "Velocity clamp changes tau", total_tests, total_failures)

        ! Test 2.2b: Expansion velocity clamping (high) & Superwind
        ! FIX: Use log_g = -1.0 to generate a physical star large enough (P ~ 1500 days)
        !      to trigger the high velocity clamp (>15 km/s) AND the superwind branch (>500 days).
        log_g = -1.0_sp 
        
        period_days = calc_period_days(m_act, log_g)
        vexp_unclamped = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        
        tau_clamped = calc_tau_expected(0, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .true.)
        tau_unclamped = calc_tau_expected(0, m_act, log_l, log_g, 10.0_sp**log_mdot_iso, .false.)
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        
        call assert_true(vexp_unclamped > AGB_VEXP_MAX, "Velocity unclamped > 15 km/s", total_tests, total_failures)
        call assert_float_equals(tau_clamped, tau_o, 1.0e-8_sp, "Velocity clamped high", total_tests, total_failures)
        call assert_true(abs(tau_unclamped - tau_o) > 0.0_sp, "Velocity clamp changes tau (high)", total_tests, total_failures)

        ! Test 2.3: Mass loss logic (isochrone override + clamp)
        ! Reset log_g to standard value for isolation
        log_g = 1.0_sp
        log_mdot_iso = -5.0_sp
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        tau_expected_o = calc_tau_expected(0, m_act, log_l, log_g, 1.0e-5_sp, .true.)
        call assert_float_equals(tau_expected_o, tau_o, 1.0e-8_sp, "Isochrone mdot (1e-5)", total_tests, total_failures)

        log_mdot_iso = 0.0_sp
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        tau_expected_o = calc_tau_expected(0, m_act, log_l, log_g, VW93_MDOT_LIMIT_ISO, .true.)
        call assert_float_equals(tau_expected_o, tau_o, 1.0e-8_sp, "Isochrone mdot clamp (1e-4)", total_tests, total_failures)

        ! Test 2.4: Superwind transition
        ctx%use_isoc_mdot_val = 0
        ! FIX: Use log_g = -1.0 again to strictly test the superwind branch logic
        log_g = -1.0_sp
        period_days = calc_period_days(m_act, log_g)
        call assert_true(period_days > VW93_PER_THRESH, "Superwind period > 500", total_tests, total_failures)
        
        tau_o = compute_circumstellar_optical_depth(ctx, 0, m_act, log_l, log_g, log_mdot_iso)
        tau_expected_o = calc_tau_expected_superwind(0, m_act, log_l, log_g)
        call assert_float_equals(tau_expected_o, tau_o, 1.0e-8_sp, "Superwind branch", total_tests, total_failures)

        deallocate(ctx)

    end subroutine test_circumstellar_optical_depth

    ! ------------------------------------------------------------------------
    ! TEST SUITE: AGB Dust Screen
    ! ------------------------------------------------------------------------
    subroutine test_agb_dust_screen()
        type(fsps_context_t), allocatable :: ctx
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: spectrum, expected
        real(sp) :: mass_act, log_t, log_l, log_g, log_mdot, c_o_ratio
        real(sp) :: log_g_local, tau_1um

        call print_group("AGB Dust Screen")

        call setup_physics_context(ctx, nlam)

        ! Test 3.1: Identity weight
        spectrum = 100.0_sp
        call apply_agb_dust_screen(ctx, 0.0_sp, spectrum, 1.0_sp, 3.6_sp, 3.0_sp, 0.0_sp, 1.5_sp, -5.0_sp)
        call assert_true(all(abs(spectrum - 100.0_sp) <= EPS), "AGB weight=0 no change", total_tests, total_failures)

        ! Test 3.2: BaSTI log_g calculation
        ctx%state%isoc_type = 'bsti'
        ctx%use_isoc_mdot_val = 1
        spectrum = 100.0_sp
        mass_act = 1.0_sp
        log_t = 3.6_sp
        log_l = 3.0_sp
        log_g = 0.0_sp
        log_mdot = -5.0_sp
        c_o_ratio = 1.5_sp

        call set_flux_dagb_pattern(ctx, nlam)

        log_g_local = log10(gsig4pi * mass_act / (10.0_sp**log_l)) + 4.0_sp * log_t
        tau_1um = compute_circumstellar_optical_depth(ctx, 1, mass_act, log_l, log_g_local, log_mdot)
        expected = spectrum * calc_dusty_transfer(ctx, 1, log_t, tau_1um)

        call apply_agb_dust_screen(ctx, 1.0_sp, spectrum, mass_act, log_t, log_l, log_g, c_o_ratio, log_mdot)
        call assert_true(all(abs(spectrum - expected) <= 5.0e-6_sp), "BaSTI log_g fallback", total_tests, total_failures)

        ! Test 3.3: Extrapolation edge case
        spectrum = 100.0_sp
        call apply_agb_dust_screen(ctx, 1.0_sp, spectrum, 1.0_sp, 3.0_sp, 3.0_sp, 0.0_sp, 1.5_sp, -5.0_sp)
        call assert_true(all(ieee_is_finite(spectrum)) .and. all(spectrum > 0.0_sp), &
                         "AGB extrapolation safe", total_tests, total_failures)

        ! Test 3.4: Transfer function application (0.5)
        call set_flux_dagb_constant(ctx, nlam, 0.5_sp)
        spectrum = 100.0_sp
        call apply_agb_dust_screen(ctx, 1.0_sp, spectrum, 1.0_sp, 3.6_sp, 3.0_sp, 0.0_sp, 1.5_sp, -5.0_sp)
        call assert_true(all(abs(spectrum - 50.0_sp) <= LESS_SMALL_DELTA), "AGB transfer function", total_tests, total_failures)

        call teardown_physics_context(ctx)

    end subroutine test_agb_dust_screen

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Energy Conservation (AGN + AGB)
    ! ------------------------------------------------------------------------
    subroutine test_energy_conservation()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: spectrum, wavelengths
        real(sp), dimension(nlam) :: spectrum_out
        real(sp), dimension(nlam) :: freqs
        real(sp) :: lbol_in, lbol_out
        real(sp) :: norm_template

        call print_group("Energy Conservation")

        call setup_physics_context(ctx, nlam)
        wavelengths = ctx%state%spec_lambda
        freqs = clight / wavelengths

        ! Test 4.1: AGN energy conservation (template normalized)
        spectrum = 1.0_sp
        lbol_in = integrate_trapezoid_array(freqs, spectrum)

        norm_template = lbol_in
        ctx%state%agndust_spec(:, :) = 1.0_sp / norm_template

        settings%fagn = 0.1_sp
        settings%dust2 = 0.0_sp
        settings%agn_tau = 10.0_sp
        ctx%dust_type_val = 0

        call apply_agn_dust_emission(ctx, settings, wavelengths, log10(lbol_in), spectrum)
        lbol_out = integrate_trapezoid_array(freqs, spectrum)
        call assert_relative_error(lbol_in * 1.1_sp, lbol_out, 1.0e-6_sp, "AGN energy conservation", total_tests, total_failures)

        ! Test 4.2: AGB energy conservation (attenuation only)
        call set_flux_dagb_constant(ctx, nlam, 0.5_sp)
        spectrum_out = 1.0_sp
        lbol_in = integrate_trapezoid_array(freqs, spectrum_out)
        call apply_agb_dust_screen(ctx, 1.0_sp, spectrum_out, 1.0_sp, 3.6_sp, 3.0_sp, 0.0_sp, 1.5_sp, -5.0_sp)
        lbol_out = integrate_trapezoid_array(freqs, spectrum_out)
        call assert_true(lbol_out < lbol_in, "AGB screen reduces flux", total_tests, total_failures)

        call teardown_physics_context(ctx)

    end subroutine test_energy_conservation

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Differential Attenuation
    ! ------------------------------------------------------------------------
    subroutine test_differential_attenuation()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: spec_young, spec_old, spec_out
        real(sp), dimension(nemline) :: neb_young, neb_old, neb_out
        real(sp) :: dust_mass

        call print_group("Differential Attenuation")

        call setup_physics_context(ctx, nlam)
        ctx%dust_type_val = 0
        ctx%add_dust_emission_val = 0
        ctx%nebemlineinspec_val = 1

        settings%dust_index = 0.0_sp
        settings%dust1_index = 0.0_sp
        settings%frac_nodust = 0.0_sp
        settings%frac_obrun = 0.0_sp
        settings%dust3 = 0.0_sp

        neb_young = 0.0_sp
        neb_old = 0.0_sp

        spec_young = 0.0_sp
        spec_old = 100.0_sp
        settings%dust1 = 5.0_sp
        settings%dust2 = 1.0_sp
        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)
        call assert_true(all(abs(spec_out - 100.0_sp * exp(-1.0_sp)) <= EPS), &
                         "Old stars diffuse only", total_tests, total_failures)

        spec_young = 100.0_sp
        spec_old = 0.0_sp
        settings%dust1 = 1.0_sp
        settings%dust2 = 1.0_sp
        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)
        call assert_true(all(abs(spec_out - 100.0_sp * exp(-2.0_sp)) <= EPS), &
                         "Young stars birth+diffuse", total_tests, total_failures)

        spec_young = 100.0_sp
        spec_old = 0.0_sp
        settings%dust1 = 100.0_sp
        settings%dust2 = 0.0_sp
        settings%frac_obrun = 0.25_sp
        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)
        call assert_true(all(abs(spec_out - 25.0_sp) <= LESS_SMALL_DELTA), &
                         "OB runaways", total_tests, total_failures)

        spec_young = 0.0_sp
        spec_old = 100.0_sp
        settings%frac_obrun = 0.0_sp
        settings%dust2 = 100.0_sp
        settings%frac_nodust = 0.10_sp
        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)
        call assert_true(all(abs(spec_out - 10.0_sp) <= LESS_SMALL_DELTA), &
                         "Patchy ISM", total_tests, total_failures)

        settings%frac_nodust = 0.0_sp
        settings%dust2 = 0.0_sp
        settings%dust1 = 1.0_sp
        settings%dust1_index = -1.0_sp
        neb_young = 0.0_sp
        neb_old = 0.0_sp
        neb_young(1) = 1.0_sp
        ctx%state%nebem_line_pos = 5500.0_sp
        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)
        call assert_float_equals(exp(-1.0_sp), neb_out(1), 1.0e-6_sp, &
                                 "Nebular line attenuation", total_tests, total_failures)

        call teardown_physics_context(ctx)

    end subroutine test_differential_attenuation

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Draine & Li + Energy Balance
    ! ------------------------------------------------------------------------
    subroutine test_draine_li_energy_balance()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        integer, parameter :: nlam = 5
        real(sp), dimension(nlam) :: spec_young, spec_old, spec_out
        real(sp), dimension(nlam) :: freqs, dust_shape
        real(sp), dimension(nemline) :: neb_young, neb_old, neb_out
        real(sp) :: dust_mass, lbol_in, lbol_out
        real(sp) :: lbol_atten, lbol_abs, emission_norm
        real(sp) :: expected_mass
        real(sp), dimension(nlam) :: spec_atten

        call print_group("Draine & Li Energy Balance")

        call setup_physics_context(ctx, nlam)
        ctx%dust_type_val = 0
        ctx%add_dust_emission_val = 1
        ctx%nebemlineinspec_val = 1

        settings%dust_index = 0.0_sp
        settings%dust1_index = 0.0_sp
        settings%dust1 = 0.0_sp
        settings%dust2 = 1.0_sp
        settings%frac_nodust = 0.0_sp
        settings%frac_obrun = 0.0_sp
        settings%duste_gamma = 0.1_sp
        settings%duste_qpah = 1.0_sp
        settings%duste_umin = 0.5_sp

        spec_young = 100.0_sp
        spec_old = 0.0_sp
        neb_young = 0.0_sp
        neb_old = 0.0_sp

        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, neb_young, &
                                                 neb_old, spec_out, dust_mass, neb_out)

        freqs = clight / ctx%state%spec_lambda
        lbol_in = integrate_trapezoid_array(freqs, spec_young)
        lbol_out = integrate_trapezoid_array(freqs, spec_out)
        call assert_relative_error(lbol_in, lbol_out, 1.0e-6_sp, "Energy conservation", total_tests, total_failures)

        call interpolate_draine_li_dust_model(ctx, settings, dust_shape)
        call assert_true(all(abs(dust_shape - 1.9_sp) <= EPS), "Draine-Li gamma weighting", total_tests, total_failures)
        emission_norm = integrate_trapezoid_array(freqs, dust_shape)
        spec_atten = 100.0_sp * exp(-1.0_sp)
        lbol_atten = integrate_trapezoid_array(freqs, spec_atten)
        lbol_abs = lbol_in - lbol_atten
        expected_mass = 3.21e-3_sp / (4.0_sp * PI) * (lbol_abs / emission_norm)
        call assert_float_equals(expected_mass, dust_mass, 1.0e-6_sp, "Dust mass scaling", total_tests, total_failures)

        call teardown_physics_context(ctx)

    end subroutine test_draine_li_energy_balance

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Dust Self-Absorption
    ! ------------------------------------------------------------------------
    subroutine test_dust_self_absorption()
        real(sp), dimension(2) :: nu
        real(sp), dimension(2) :: shape, transmission, spec_final
        real(sp) :: lbol

        call print_group("Dust Self-Absorption")

        nu = [1.0_sp, 2.0_sp]
        shape = 1.0_sp

        ! Test 7.1: Transparent limit
        transmission = 1.0_sp
        call calculate_dust_self_absorption(nu, shape, transmission, 100.0_sp, spec_final)
        lbol = integrate_trapezoid_array(nu, spec_final)
        call assert_float_equals(100.0_sp, lbol, 1.0e-6_sp, "Transparent limit", total_tests, total_failures)

        ! Test 7.2: Opaque limit
        transmission = 0.5_sp
        call calculate_dust_self_absorption(nu, shape, transmission, 100.0_sp, spec_final)
        lbol = integrate_trapezoid_array(nu, spec_final)
        call assert_float_equals(100.0_sp, lbol, 1.0e-6_sp, "Opaque limit", total_tests, total_failures)

        ! Test 7.3: Differential self-absorption
        transmission = [0.1_sp, 1.0_sp]
        call calculate_dust_self_absorption(nu, shape, transmission, 100.0_sp, spec_final)
        call assert_float_equals(10.0_sp, spec_final(2) / spec_final(1), 1.0e-6_sp, &
                                 "Differential self-absorption", total_tests, total_failures)

        ! Test 8.2: Zero absorbed luminosity
        transmission = 1.0_sp
        call calculate_dust_self_absorption(nu, shape, transmission, 0.0_sp, spec_final)
        call assert_true(all(spec_final == 0.0_sp), "Zero absorbed luminosity", total_tests, total_failures)

    end subroutine test_dust_self_absorption

    ! ------------------------------------------------------------------------
    ! TEST SUITE: Robustness
    ! ------------------------------------------------------------------------
    subroutine test_robustness()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: settings
        real(sp), dimension(:), allocatable :: w
        real(sp), dimension(:), allocatable :: res
        real(sp), dimension(:), allocatable :: spec_young, spec_old, spec_out
        real(sp), dimension(:), allocatable :: neb_young, neb_old, neb_out
        real(sp) :: dust_mass
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

        deallocate(ctx)

        ! Test 8.1: No dust (output equals input, dust mass ~ 0)
        call setup_physics_context(ctx, 5)
        ctx%dust_type_val = 0
        ctx%add_dust_emission_val = 0
        ctx%nebemlineinspec_val = 1

        settings%dust_index = 0.0_sp
        settings%dust1 = 0.0_sp
        settings%dust2 = 0.0_sp
        settings%dust3 = 0.0_sp
        settings%frac_nodust = 0.0_sp
        settings%frac_obrun = 0.0_sp

        n = size(ctx%state%spec_lambda)
        allocate(spec_young(n), spec_old(n), spec_out(n))
        allocate(neb_young(nemline), neb_old(nemline), neb_out(nemline))
        spec_young = 100.0_sp
        spec_old = 0.0_sp
        neb_young = 0.0_sp
        neb_old = 0.0_sp

        call apply_dust_attenuation_and_emission(ctx, settings, spec_young, spec_old, &
                             neb_young, neb_old, spec_out, dust_mass, neb_out)
        call assert_true(all(abs(spec_out - 100.0_sp) <= EPS), "No dust spectrum unchanged", total_tests, total_failures)
        call assert_true(dust_mass <= SAFE_FLOOR, "No dust mass ~ 0", total_tests, total_failures)

        deallocate(spec_young, spec_old, spec_out, neb_young, neb_old, neb_out)
        call teardown_physics_context(ctx)

        ! Test 7.2: Array shape conformance
        allocate(ctx)
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

        deallocate(ctx)
        allocate(w(n), res(n))
        w = [(real(i, sp), i=1, n)]
        res = compute_attenuation_curve(w, 0, settings, ctx)
        call assert_int_equals(size(w), size(res), "Array shape (size=1000)", total_tests, total_failures)
        deallocate(w, res)

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
    ! Helper: mock context setup for full dust physics
    ! ------------------------------------------------------------------------
    subroutine setup_physics_context(ctx, nlam)
        type(fsps_context_t), allocatable, intent(out) :: ctx
        integer, intent(in) :: nlam
        integer :: i, j, k
        real(sp) :: teff_min, teff_max, tau_min, tau_max, dteff, dtau

        allocate(ctx)
        allocate(ctx%state%spec_lambda(nlam))
        allocate(ctx%state%wgdust(nlam, 3, 4, 2))
        allocate(ctx%state%g03smcextn(nlam))
        allocate(ctx%state%agndust_spec(nlam, nagndust))
        allocate(ctx%state%qpaharr(3))
        allocate(ctx%state%uminarr(3))
        allocate(ctx%state%dustem2_dustem(nlam, 3, 6))
        allocate(ctx%state%flux_dagb(nlam, 2, nteff_dagb, ntau_dagb))

        ctx%state%nqpah_dustem = 3
        ctx%state%numin_dustem = 3

        do i = 1, nlam
            ctx%state%spec_lambda(i) = 1000.0_sp * real(i, sp)
        end do

        ctx%state%nebem_line_pos = 5500.0_sp

        ctx%state%agndust_tau = [10.0_sp, 20.0_sp, 30.0_sp, 40.0_sp, 50.0_sp, 60.0_sp, 70.0_sp, 80.0_sp, 90.0_sp]
        ctx%state%agndust_spec(:, :) = 0.0_sp
        ctx%state%agndust_spec(:, 1) = 1.0_sp
        ctx%state%agndust_spec(:, 2) = 2.0_sp
        ctx%state%agndust_spec(:, 3) = 3.0_sp

        ctx%state%qpaharr = [1.0_sp, 2.0_sp, 3.0_sp]
        ctx%state%uminarr = [0.5_sp, 1.0_sp, 1.5_sp]

        ctx%state%dustem2_dustem(:, :, :) = 0.0_sp
        do j = 1, 3
            do k = 1, 3
                ctx%state%dustem2_dustem(:, j, 2*k - 1) = 1.0_sp
                ctx%state%dustem2_dustem(:, j, 2*k) = 10.0_sp
            end do
        end do

        ctx%state%wgdust(:, :, :, :) = 1.0_sp
        ctx%state%g03smcextn(:) = 1.0_sp

        teff_min = 3.5_sp
        teff_max = 4.0_sp
        dteff = (teff_max - teff_min) / real(nteff_dagb - 1, sp)
        do i = 1, nteff_dagb
            ctx%state%teff_dagb(1, i) = 10.0_sp**(teff_min + dteff * real(i - 1, sp))
            ctx%state%teff_dagb(2, i) = ctx%state%teff_dagb(1, i)
        end do

        tau_min = -2.0_sp
        tau_max = 1.0_sp
        dtau = (tau_max - tau_min) / real(ntau_dagb - 1, sp)
        do i = 1, ntau_dagb
            ctx%state%tau1_dagb(1, i) = tau_min + dtau * real(i - 1, sp)
            ctx%state%tau1_dagb(2, i) = ctx%state%tau1_dagb(1, i)
        end do

        ctx%state%flux_dagb(:, :, :, :) = 1.0_sp

    end subroutine setup_physics_context

    subroutine teardown_physics_context(ctx)
        type(fsps_context_t), allocatable, intent(inout) :: ctx

        if (associated(ctx%state%spec_lambda)) deallocate(ctx%state%spec_lambda)
        if (associated(ctx%state%wgdust)) deallocate(ctx%state%wgdust)
        if (associated(ctx%state%g03smcextn)) deallocate(ctx%state%g03smcextn)
        if (associated(ctx%state%agndust_spec)) deallocate(ctx%state%agndust_spec)
        if (associated(ctx%state%qpaharr)) deallocate(ctx%state%qpaharr)
        if (associated(ctx%state%uminarr)) deallocate(ctx%state%uminarr)
        if (associated(ctx%state%dustem2_dustem)) deallocate(ctx%state%dustem2_dustem)
        if (associated(ctx%state%flux_dagb)) deallocate(ctx%state%flux_dagb)
        deallocate(ctx)

    end subroutine teardown_physics_context

    subroutine set_flux_dagb_constant(ctx, nlam, value)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: nlam
        real(sp), intent(in) :: value

        ctx%state%flux_dagb(1:nlam, :, :, :) = value

    end subroutine set_flux_dagb_constant

    subroutine set_flux_dagb_pattern(ctx, nlam)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: nlam
        integer :: i, j, k, l

        do l = 1, ntau_dagb
            do k = 1, nteff_dagb
                do j = 1, 2
                    do i = 1, nlam
                        ctx%state%flux_dagb(i, j, k, l) = 0.1_sp * real(k, sp) + 0.01_sp * real(l, sp)
                    end do
                end do
            end do
        end do

    end subroutine set_flux_dagb_pattern

    pure function calc_period_days(m_act, log_g) result(period_days)
        real(sp), intent(in) :: m_act, log_g
        real(sp) :: period_days
        real(sp) :: radius_solar

        radius_solar = sqrt(m_act * msun * newton / (10.0_sp**log_g)) / rsun
        period_days = 10.0_sp**(AGB_PER_INTERCEPT + AGB_PER_SLOPE_R * log10(radius_solar) + &
                               AGB_PER_SLOPE_M * log10(m_act))
    end function calc_period_days

    pure function calc_tau_expected(c_rich_flag, m_act, log_l, log_g, mdot_sol_yr, clamp_velocity) result(tau)
        integer, intent(in) :: c_rich_flag
        real(sp), intent(in) :: m_act, log_l, log_g, mdot_sol_yr
        logical, intent(in) :: clamp_velocity
        real(sp) :: tau
        real(sp) :: period_days, velocity_exp
        real(sp) :: inner_radius_cm, dust_gas_ratio, kappa_eff

        if (c_rich_flag == 1) then
            kappa_eff = AGB_KAPPA_C_RICH
        else
            kappa_eff = AGB_KAPPA_O_RICH
        end if

        period_days = calc_period_days(m_act, log_g)
        velocity_exp = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        if (clamp_velocity) then
            velocity_exp = max(min(velocity_exp, AGB_VEXP_MAX), AGB_VEXP_MIN)
        end if

        if (c_rich_flag == 1) then
            inner_radius_cm = AGB_RIN_FACTOR_C * (10.0_sp**log_l)**0.5_sp
            dust_gas_ratio = AGB_DELTA_C_RICH
        else
            inner_radius_cm = AGB_RIN_FACTOR_O * (10.0_sp**log_l)**0.5_sp
            dust_gas_ratio = AGB_DELTA_O_RICH
        end if

        dust_gas_ratio = dust_gas_ratio * (velocity_exp**2 / AGB_DTG_VEL_NORM) * &
                         ((10.0_sp**log_l / AGB_DTG_LUM_NORM)**AGB_DTG_LUM_EXP)

        tau = kappa_eff * dust_gas_ratio * (mdot_sol_yr * msun / yr2sc) / &
              inner_radius_cm / (4.0_sp * PI) / (velocity_exp * 1.0e5_sp)

    end function calc_tau_expected

    pure function calc_tau_expected_superwind(c_rich_flag, m_act, log_l, log_g) result(tau)
        integer, intent(in) :: c_rich_flag
        real(sp), intent(in) :: m_act, log_l, log_g
        real(sp) :: tau
        real(sp) :: period_days, velocity_exp, mdot_sol_yr

        period_days = calc_period_days(m_act, log_g)
        velocity_exp = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        velocity_exp = max(min(velocity_exp, AGB_VEXP_MAX), AGB_VEXP_MIN)
        mdot_sol_yr = (10.0_sp**log_l) / velocity_exp * VW93_SUPERWIND_NORM * yr2sc / clight

        tau = calc_tau_expected(c_rich_flag, m_act, log_l, log_g, mdot_sol_yr, .true.)
    end function calc_tau_expected_superwind

    pure function calc_dusty_transfer(ctx, c_rich_flag, log_t, tau_1um) result(transfer)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in) :: c_rich_flag
        real(sp), intent(in) :: log_t, tau_1um
        real(sp), dimension(size(ctx%state%flux_dagb, 1)) :: transfer
        integer :: idx_teff, idx_tau, n_teff_grid, n_tau_grid
        real(sp) :: w_teff, w_tau

        n_teff_grid = size(ctx%state%teff_dagb, 2)
        n_tau_grid  = size(ctx%state%tau1_dagb, 2)

        idx_teff = find_interval(ctx%state%teff_dagb(c_rich_flag + 1, :), 10.0_sp**log_t)
        idx_teff = max(1, min(idx_teff, n_teff_grid - 1))

        w_teff = (10.0_sp**log_t - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff)) / &
                 (ctx%state%teff_dagb(c_rich_flag + 1, idx_teff + 1) - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff))
        w_teff = max(-1.0_sp, min(w_teff, 1.0_sp))

        idx_tau = find_interval(ctx%state%tau1_dagb(c_rich_flag + 1, :), log10(tau_1um))
        idx_tau = max(1, min(idx_tau, n_tau_grid - 1))

        w_tau = (log10(tau_1um) - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau)) / &
                (ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau + 1) - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau))
        w_tau = max(-1.0_sp, min(w_tau, 1.0_sp))

        transfer = &
            (1.0_sp - w_teff) * (1.0_sp - w_tau) * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff, idx_tau) + &
            w_teff            * (1.0_sp - w_tau) * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff + 1, idx_tau) + &
            w_teff            * w_tau            * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff + 1, idx_tau + 1) + &
            (1.0_sp - w_teff) * w_tau            * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff, idx_tau + 1)

    end function calc_dusty_transfer

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
