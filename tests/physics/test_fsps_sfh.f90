module test_fsps_sfh_mod
    use fsps_precision, only: WP
    use fsps_constants, only: SAFE_FLOOR
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params, sfhparams, compspout
    use fsps_sfh, only: compute_ssp_weights, get_sfh_properties_at_age, compute_sfh_statistics
    use fsps_interpolation, only: find_interval
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_float_equals, assert_true, assert_relative_error
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_sfh_tests, total_failures, total_tests

    real(WP), parameter :: EPS = 1.0e-6_wp
    real(WP), parameter :: REL_EPS = 1.0e-4_wp

contains

    subroutine run_fsps_sfh_tests()
        call print_minor_header("fsps_sfh")

        call test_mass_conservation_tau()
        call test_grid_invariance_delayed_tau()
        call test_youngest_edge_includes_zero_bin()

        call test_tabular_nonzero_weights()
        call test_analytic_vs_tabular_match()
        call test_simha_truncation_logic()

        call test_ssfr_surviving_mass_normalization()
        call test_mean_age_constant_sfh()
        call test_tabular_unit_indexing()

        call test_ssp_type_zero_properties()
        call test_tiny_tau_stability()
        call test_empty_segment_cycle()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_sfh_tests

    ! ------------------------------------------------------------------------
    ! Helper: setup minimal SFH context
    ! ------------------------------------------------------------------------
    subroutine setup_sfh_context(ctx, time_full, interpolation_type, tiny_logt)
        type(fsps_context_t), allocatable, intent(out) :: ctx
        real(WP), dimension(:), intent(in) :: time_full
        integer, intent(in), optional :: interpolation_type
        real(WP), intent(in), optional :: tiny_logt

        allocate(ctx)
        allocate(ctx%state%time_full(size(time_full)))
        ctx%state%time_full = time_full
        ctx%state%ntfull = size(time_full)
        ctx%state%ntabsfh = 0

        if (present(interpolation_type)) then
            ctx%interpolation_type_val = interpolation_type
        else
            ctx%interpolation_type_val = 1
        end if

        if (present(tiny_logt)) then
            ctx%tiny_logt_val = tiny_logt
        else
            ctx%tiny_logt_val = 0.0_wp
        end if
    end subroutine setup_sfh_context

    subroutine teardown_sfh_context(ctx)
        type(fsps_context_t), allocatable, intent(inout) :: ctx
        if (.not. allocated(ctx)) return
        if (associated(ctx%state%time_full)) deallocate(ctx%state%time_full)
        deallocate(ctx)
    end subroutine teardown_sfh_context

    subroutine set_sfh_tab(ctx, t_years, sfr)
        type(fsps_context_t), intent(inout) :: ctx
        real(WP), dimension(:), intent(in) :: t_years, sfr
        integer :: n

        n = size(t_years)
        ctx%state%sfh_tab = 0.0_wp
        ctx%state%ntabsfh = n
        ctx%state%sfh_tab(1, 1:n) = t_years
        ctx%state%sfh_tab(2, 1:n) = sfr
        ctx%state%sfh_tab(3, 1:n) = 0.0_wp
    end subroutine set_sfh_tab

    subroutine fill_log_grid(arr, log_min, log_max)
        real(WP), dimension(:), intent(out) :: arr
        real(WP), intent(in) :: log_min, log_max
        integer :: i, n
        n = size(arr)
        do i = 1, n
            arr(i) = log_min + (log_max - log_min) * real(i - 1, WP) / real(n - 1, WP)
        end do
    end subroutine fill_log_grid

    integer function find_exact_index(arr, value) result(idx)
        real(WP), dimension(:), intent(in) :: arr
        real(WP), intent(in) :: value
        real(WP) :: diff_min
        integer :: loc

        loc = minloc(abs(arr - value), dim=1)
        diff_min = abs(arr(loc) - value)
        if (diff_min > 1.0e-12_wp) then
            idx = loc
        else
            idx = loc
        end if
    end function find_exact_index

    ! ------------------------------------------------------------------------
    ! Group 1: Moment Method Integration Engine
    ! ------------------------------------------------------------------------
    subroutine test_mass_conservation_tau()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh
        real(WP), allocatable :: time_full(:), weights(:)
        integer :: n
        real(WP) :: expected, actual

        call print_group("Moment Method: Conservation of Mass (Tau)")

        n = 60
        allocate(time_full(n), weights(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        sfh%type = 1
        sfh%tau = 1.0e9_wp
        sfh%tage = 1.0e9_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx, sfh, 1, n - 1, weights)

        actual = sum(weights)
        expected = exp(sfh%tage / max(sfh%tau, SAFE_FLOOR)) - 1.0_wp

        call assert_relative_error(expected, actual, 5.0e-3_wp, &
                                   "Total mass matches analytic integral", total_tests, total_failures)

        deallocate(time_full, weights)
        call teardown_sfh_context(ctx)
    end subroutine test_mass_conservation_tau

    subroutine test_grid_invariance_delayed_tau()
        type(fsps_context_t), allocatable :: ctx_a, ctx_b
        type(sfhparams) :: sfh
        real(WP), allocatable :: grid_a(:), grid_b(:), weights_a(:), weights_b(:)
        integer :: n
        real(WP) :: mass_a, mass_b

        call print_group("Moment Method: Grid Invariance (Total Mass)")

        n = 40
        allocate(grid_a(n))
        call fill_log_grid(grid_a, 6.0_wp, 10.0_wp)

        ! Create Grid B by splitting one interval in Grid A
        allocate(grid_b(n + 1))
        grid_b(1:20) = grid_a(1:20)
        grid_b(21) = 0.5_wp * (grid_a(20) + grid_a(21)) ! Insert midpoint
        grid_b(22:n+1) = grid_a(21:n)

        call setup_sfh_context(ctx_a, grid_a, interpolation_type=1, tiny_logt=4.0_wp)
        call setup_sfh_context(ctx_b, grid_b, interpolation_type=1, tiny_logt=4.0_wp)

        allocate(weights_a(size(grid_a)))
        allocate(weights_b(size(grid_b)))

        sfh%type = 4 ! Delayed Tau
        sfh%tau = 2.0e9_wp
        sfh%tage = 6.0e9_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx_a, sfh, 1, size(grid_a) - 1, weights_a)
        call compute_ssp_weights(ctx_b, sfh, 1, size(grid_b) - 1, weights_b)

        ! Compare Total Mass formed.
        ! We cannot compare partial sums near index 20, because the weight distribution
        ! function W(t) changes shape locally when the grid spacing changes.
        ! However, the Total Mass (integral of SFR) must be conserved.
        mass_a = sum(weights_a)
        mass_b = sum(weights_b)

        call assert_relative_error(mass_a, mass_b, 1.0e-5_wp, &
                                   "Total mass invariant under grid refinement", total_tests, total_failures)

        deallocate(grid_a, grid_b, weights_a, weights_b)
        call teardown_sfh_context(ctx_a)
        call teardown_sfh_context(ctx_b)
    end subroutine test_grid_invariance_delayed_tau

    subroutine test_youngest_edge_includes_zero_bin()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh
        real(WP), allocatable :: time_full(:), w_min0(:), w_min1(:)
        integer :: n

        call print_group("Moment Method: Youngest Edge (idx_min == 0)")

        n = 30
        allocate(time_full(n), w_min0(n), w_min1(n))
        call fill_log_grid(time_full, 6.0_wp, 9.5_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        sfh%type = 0
        sfh%tau = 1.0_wp
        sfh%tage = 3.0e9_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx, sfh, 0, n - 1, w_min0)
        call compute_ssp_weights(ctx, sfh, 1, n - 1, w_min1)

        call assert_true(w_min0(1) > w_min1(1), &
                         "weights(1) includes [0, t1] when idx_min == 0", total_tests, total_failures)

        deallocate(time_full, w_min0, w_min1)
        call teardown_sfh_context(ctx)
    end subroutine test_youngest_edge_includes_zero_bin

    ! ------------------------------------------------------------------------
    ! Group 2: Analytic vs Tabular Dispatch
    ! ------------------------------------------------------------------------
    subroutine test_tabular_nonzero_weights()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh
        real(WP), allocatable :: time_full(:), weights(:)
        real(WP), allocatable :: t_tab(:), sfr_tab(:)
        integer :: n
        real(WP) :: total_mass

        call print_group("Tabular SFH: Non-zero Weight Regression")

        n = 50
        allocate(time_full(n), weights(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        allocate(t_tab(4), sfr_tab(4))
        t_tab = [0.0_wp, 1.0e9_wp, 1.001e9_wp, 2.0e9_wp]
        sfr_tab = [1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
        call set_sfh_tab(ctx, t_tab, sfr_tab)

        sfh%type = 3
        sfh%tau = 1.0_wp
        sfh%tage = 2.0e9_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx, sfh, 1, n - 1, weights)
        total_mass = sum(weights)

        call assert_true(any(abs(weights) > 0.0_wp), "Tabular weights are non-zero", total_tests, total_failures)
        call assert_relative_error(1.0e9_wp, total_mass, 1.0e-2_wp, &
                                   "Square-wave mass integrates to ~1 Gyr", total_tests, total_failures)

        deallocate(time_full, weights, t_tab, sfr_tab)
        call teardown_sfh_context(ctx)
    end subroutine test_tabular_nonzero_weights

    subroutine test_analytic_vs_tabular_match()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh_analytic, sfh_tabular
        real(WP), allocatable :: time_full(:), weights_a(:), weights_b(:)
        real(WP), allocatable :: t_tab(:), sfr_tab(:)
        integer :: n, ntab, i
        real(WP) :: tau, tmax, max_rel
        real(WP) :: denom

        call print_group("Analytic vs Tabular: High-Resolution Convergence")

        n = 80
        allocate(time_full(n), weights_a(n), weights_b(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        tau = 5.0e9_wp
        tmax = 1.0e10_wp

        ntab = 10001
        allocate(t_tab(ntab), sfr_tab(ntab))
        do i = 1, ntab
            t_tab(i) = tmax * real(i - 1, WP) / real(ntab - 1, WP)
            sfr_tab(i) = exp(t_tab(i) / tau) / tau
        end do
        call set_sfh_tab(ctx, t_tab, sfr_tab)

        sfh_analytic%type = 1
        sfh_analytic%tau = tau
        sfh_analytic%tage = tmax
        sfh_analytic%tq = 0.0_wp
        sfh_analytic%sf_trunc = 0.0_wp
        sfh_analytic%sf_slope = 0.0_wp
        sfh_analytic%t0 = 0.0_wp
        sfh_analytic%tb = 0.0_wp
        sfh_analytic%use_simha_limits = 0

        sfh_tabular = sfh_analytic
        sfh_tabular%type = 3

        call compute_ssp_weights(ctx, sfh_analytic, 1, n - 1, weights_a)
        call compute_ssp_weights(ctx, sfh_tabular, 1, n - 1, weights_b)

        max_rel = 0.0_wp
        do i = 1, n
            denom = max(1.0_wp, abs(weights_a(i)))
            max_rel = max(max_rel, abs(weights_a(i) - weights_b(i)) / denom)
        end do

        call assert_true(max_rel < 1.0e-3_wp, "Tabular matches analytic weights", total_tests, total_failures)

        deallocate(time_full, weights_a, weights_b, t_tab, sfr_tab)
        call teardown_sfh_context(ctx)
    end subroutine test_analytic_vs_tabular_match

    subroutine test_simha_truncation_logic()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:)
        real(WP) :: mass_frac, sfr_norm, frac_linear

        call print_group("Simha: Linear-Exponential Truncation")

        allocate(time_full(30))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        pset%sfh = 5
        pset%tau = 1.0_wp
        pset%sf_start = 0.0_wp
        pset%sf_trunc = 2.0_wp
        pset%sf_slope = -1.0_wp

        call get_sfh_properties_at_age(ctx, pset, 4.0_wp, mass_frac, sfr_norm, frac_linear)

        call assert_float_equals(0.0_wp, sfr_norm, 1.0e-12_wp, "SFR truncated at t > t_zero", total_tests, total_failures)
        call assert_true(frac_linear > 0.0_wp, "Linear mass fraction positive", total_tests, total_failures)

        deallocate(time_full)
        call teardown_sfh_context(ctx)
    end subroutine test_simha_truncation_logic

    ! ------------------------------------------------------------------------
    ! Group 3: Statistics & Physical Normalization
    ! ------------------------------------------------------------------------
    subroutine test_ssfr_surviving_mass_normalization()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout) :: model
        real(WP), allocatable :: time_full(:)
        real(WP) :: ssfr_log(3), mean_age
        real(WP) :: expected

        call print_group("sSFR Normalization: Surviving Mass")

        allocate(time_full(10))
        call fill_log_grid(time_full, 6.0_wp, 8.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        pset%sfh = 1
        pset%const = 1.0_wp
        pset%sf_start = 0.0_wp
        pset%tau = 1.0_wp

        model%age = 0.0_wp
        model%mass_csp = 0.6_wp

        call compute_sfh_statistics(ctx, pset, model, ssfr_log, mean_age)

        expected = log10(1.0_wp / model%mass_csp)
        call assert_float_equals(expected, ssfr_log(1), 1.0e-6_wp, &
                                 "sSFR uses surviving mass (1 Myr window)", total_tests, total_failures)

        deallocate(time_full)
        call teardown_sfh_context(ctx)
    end subroutine test_ssfr_surviving_mass_normalization

    subroutine test_mean_age_constant_sfh()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout) :: model
        real(WP), allocatable :: time_full(:)
        real(WP) :: ssfr_log(3), mean_age

        call print_group("Mean Age: Constant SFH")

        allocate(time_full(20))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        pset%sfh = 1
        pset%const = 1.0_wp
        pset%sf_start = 0.0_wp

        model%age = 10.0_wp
        model%mass_csp = 1.0_wp

        call compute_sfh_statistics(ctx, pset, model, ssfr_log, mean_age)

        call assert_float_equals(5.0_wp, mean_age, 1.0e-6_wp, &
                                 "Mass-weighted mean age is 5 Gyr", total_tests, total_failures)

        deallocate(time_full)
        call teardown_sfh_context(ctx)
    end subroutine test_mean_age_constant_sfh

    subroutine test_tabular_unit_indexing()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:)
        real(WP) :: mass_frac, sfr_norm, frac_linear
        real(WP), allocatable :: t_tab(:), sfr_tab(:)

        call print_group("Units: Tabular Time Indexing in Years")

        allocate(time_full(10))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        allocate(t_tab(3), sfr_tab(3))
        t_tab = [1.0e8_wp, 1.0e9_wp, 2.0e9_wp]
        sfr_tab = [0.1_wp, 0.5_wp, 0.9_wp]
        call set_sfh_tab(ctx, t_tab, sfr_tab)

        pset%sfh = 3

        call get_sfh_properties_at_age(ctx, pset, 1.0_wp, mass_frac, sfr_norm, frac_linear)

        call assert_float_equals(0.5_wp, sfr_norm, 1.0e-12_wp, &
                                 "Tabular lookup uses GYR_TO_YR conversion", total_tests, total_failures)

        deallocate(time_full, t_tab, sfr_tab)
        call teardown_sfh_context(ctx)
    end subroutine test_tabular_unit_indexing

    ! ------------------------------------------------------------------------
    ! Group 4: Edge Cases & Stability
    ! ------------------------------------------------------------------------
    subroutine test_ssp_type_zero_properties()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:)
        real(WP) :: mass_frac, sfr_norm, frac_linear

        call print_group("Edge Case: SSP Type 0 Properties")

        allocate(time_full(10))
        call fill_log_grid(time_full, 6.0_wp, 8.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        pset%sfh = 0

        call get_sfh_properties_at_age(ctx, pset, 1.0_wp, mass_frac, sfr_norm, frac_linear)

        call assert_float_equals(1.0_wp, mass_frac, EPS, "SSP mass_frac=1", total_tests, total_failures)
        call assert_float_equals(0.0_wp, sfr_norm, EPS, "SSP sfr_norm=0", total_tests, total_failures)

        deallocate(time_full)
        call teardown_sfh_context(ctx)
    end subroutine test_ssp_type_zero_properties

    subroutine test_tiny_tau_stability()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh
        real(WP), allocatable :: time_full(:), weights(:)
        integer :: n

        call print_group("Edge Case: Tiny Tau Stability")

        n = 25
        allocate(time_full(n), weights(n))
        call fill_log_grid(time_full, 6.0_wp, 9.0_wp)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        sfh%type = 1
        sfh%tau = 1.0e4_wp
        sfh%tage = 1.0e6_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx, sfh, 1, n - 1, weights)

        call assert_true(all(ieee_is_finite(weights)), "No NaN/Inf for tiny tau", total_tests, total_failures)

        deallocate(time_full, weights)
        call teardown_sfh_context(ctx)
    end subroutine test_tiny_tau_stability

    subroutine test_empty_segment_cycle()
        type(fsps_context_t), allocatable :: ctx
        type(sfhparams) :: sfh
        real(WP), allocatable :: time_full(:), weights(:)
        integer :: n

        call print_group("Edge Case: Empty Segment Cycle")

        n = 15
        allocate(time_full(n), weights(n))
        call fill_log_grid(time_full, 6.0_wp, 9.0_wp)
        time_full(7) = time_full(6)
        call setup_sfh_context(ctx, time_full, interpolation_type=1, tiny_logt=4.0_wp)

        sfh%type = 0
        sfh%tau = 1.0_wp
        sfh%tage = 1.0e9_wp
        sfh%tq = 0.0_wp
        sfh%sf_trunc = 0.0_wp
        sfh%sf_slope = 0.0_wp
        sfh%t0 = 0.0_wp
        sfh%tb = 0.0_wp
        sfh%use_simha_limits = 0

        call compute_ssp_weights(ctx, sfh, 1, n - 1, weights)

        call assert_true(all(ieee_is_finite(weights)), "Empty segment skipped safely", total_tests, total_failures)

        deallocate(time_full, weights)
        call teardown_sfh_context(ctx)
    end subroutine test_empty_segment_cycle

end module test_fsps_sfh_mod
