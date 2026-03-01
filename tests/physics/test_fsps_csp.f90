module test_fsps_csp_mod
    use fsps_precision, only: WP
    use fsps_constants, only: SAFE_FLOOR, NEMLINE, C_LIGHT
    use fsps_types, only: params, sfhparams, compspout
    use fsps_context_types, only: fsps_context_t
    use fsps_csp, only: compute_sfh_weights, convert_sfhparams, integrate_csp_step, apply_dust_physics, &
                        apply_post_processing, compute_csp_scenario, csp_buffer_t, init_csp_buffer, free_csp_buffer
    use fsps_interpolation, only: find_interval
    use fsps_cosmology, only: get_igm_transmission
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, &
                              assert_float_equals, assert_true, assert_relative_error, assert_int_equals
    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    real(WP), parameter :: EPS = 1.0e-6_wp
    real(WP), parameter :: REL_EPS = 1.0e-3_wp

    public :: run_fsps_csp_tests, total_failures, total_tests

contains

    subroutine run_fsps_csp_tests()
        call print_minor_header("fsps_csp")

        call test_convert_sfhparams_units()
        call test_single_ssp_dirac()
        call test_exponential_tau_ratio()
        call test_delayed_tau_peak()
        call test_burst_addition()
        call test_tabular_zmix_interpolation()

        call test_all_old_limit()
        call test_all_young_limit()
        call test_split_boundary()
        call test_mass_sum_consistency()

        call test_dust_screen_logic()
        call test_igm_absorption_toggle()
        call test_smoothing_conserves_flux()

        call test_snapshot_vs_history()
        call test_mass_normalization_surviving()
        call test_nebular_precalculation()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_csp_tests

    ! ---------------------------------------------------------------------
    ! Helpers
    ! ---------------------------------------------------------------------
    subroutine fill_log_grid(arr, log_min, log_max)
        real(WP), dimension(:), intent(out) :: arr
        real(WP), intent(in) :: log_min, log_max
        integer :: i, n

        n = size(arr)
        do i = 1, n
            arr(i) = log_min + (log_max - log_min) * real(i - 1, WP) / real(n - 1, WP)
        end do
    end subroutine fill_log_grid

    integer function nearest_index(arr, value) result(idx)
        real(WP), dimension(:), intent(in) :: arr
        real(WP), intent(in) :: value

        idx = minloc(abs(arr - value), dim=1)
    end function nearest_index

    subroutine setup_basic_context(ctx, time_full, spec_lambda, nz, n_bands)
        type(fsps_context_t), allocatable, intent(out) :: ctx
        real(WP), dimension(:), intent(in) :: time_full
        real(WP), dimension(:), intent(in) :: spec_lambda
        integer, intent(in) :: nz, n_bands
        integer :: i

        allocate(ctx)

        allocate(ctx%state%time_full(size(time_full)))
        ctx%state%time_full = time_full
        ctx%state%ntfull = size(time_full)
        ctx%state%nspec = size(spec_lambda)
        ctx%state%nz = nz
        ctx%state%nbands = n_bands
        ctx%state%nindx = 0
        ctx%state%check_sps_setup = 1

        ctx%interpolation_type_val = 1
        ctx%tiny_logt_val = 4.0_wp

        allocate(ctx%state%spec_lambda(size(spec_lambda)))
        allocate(ctx%state%spec_nu(size(spec_lambda)))
        ctx%state%spec_lambda = spec_lambda
        ctx%state%spec_nu = C_LIGHT / spec_lambda

        allocate(ctx%state%zlegend(nz))
        do i = 1, nz
            ctx%state%zlegend(i) = 0.01_wp + 0.01_wp * real(i - 1, WP)
        end do

        if (n_bands > 0) then
            allocate(ctx%state%bands(size(spec_lambda), n_bands))
            ctx%state%bands = 1.0_wp
        end if

        ctx%compute_vega_mags_val = 0
        ctx%compute_light_ages_val = 0
        ctx%redshift_colors_val = 0
        ctx%add_dust_emission_val = 0
        ctx%dust_type_val = 0
        ctx%add_agn_dust_val = 0
        ctx%nebemlineinspec_val = 0
        ctx%add_neb_continuum_val = 0
        ctx%add_neb_emission_val = 0
        ctx%add_igm_absorption_val = 0
        ctx%smooth_velocity_val = 1
        ctx%smoothspec_fast_val = 0

        ctx%state%nebem_line_pos = 5500.0_wp

        !$acc enter data copyin(ctx)
        !$acc enter data copyin(ctx%state)
        !$acc enter data copyin(ctx%state%time_full, ctx%state%spec_lambda, ctx%state%spec_nu, ctx%state%zlegend)
        !$acc enter data attach(ctx%state%time_full)
        !$acc enter data attach(ctx%state%spec_lambda)
        !$acc enter data attach(ctx%state%spec_nu)
        !$acc enter data attach(ctx%state%zlegend)
        if (associated(ctx%state%bands)) then
            !$acc enter data copyin(ctx%state%bands)
            !$acc enter data attach(ctx%state%bands)
        end if
    end subroutine setup_basic_context

    subroutine teardown_basic_context(ctx)
        type(fsps_context_t), allocatable, intent(inout) :: ctx

        if (.not. allocated(ctx)) return

        if (associated(ctx%state%time_full)) then
            !$acc exit data delete(ctx%state%time_full)
        end if
        if (associated(ctx%state%spec_lambda)) then
            !$acc exit data delete(ctx%state%spec_lambda)
        end if
        if (associated(ctx%state%spec_nu)) then
            !$acc exit data delete(ctx%state%spec_nu)
        end if
        if (associated(ctx%state%bands)) then
            !$acc exit data delete(ctx%state%bands)
        end if
        if (associated(ctx%state%zlegend)) then
            !$acc exit data delete(ctx%state%zlegend)
        end if
        if (associated(ctx%state%neb_res_min)) then
            !$acc exit data delete(ctx%state%neb_res_min)
        end if
        if (associated(ctx%state%gaussnebarr)) then
            !$acc exit data delete(ctx%state%gaussnebarr)
        end if
        !$acc exit data delete(ctx%state)
        !$acc exit data delete(ctx)

        if (associated(ctx%state%time_full)) deallocate(ctx%state%time_full)
        if (associated(ctx%state%spec_lambda)) deallocate(ctx%state%spec_lambda)
        if (associated(ctx%state%spec_nu)) deallocate(ctx%state%spec_nu)
        if (associated(ctx%state%bands)) deallocate(ctx%state%bands)
        if (associated(ctx%state%zlegend)) deallocate(ctx%state%zlegend)
        if (associated(ctx%state%indexdefined)) deallocate(ctx%state%indexdefined)
        if (associated(ctx%state%neb_res_min)) deallocate(ctx%state%neb_res_min)
        if (associated(ctx%state%gaussnebarr)) deallocate(ctx%state%gaussnebarr)

        deallocate(ctx)
    end subroutine teardown_basic_context

    subroutine map_csp_buffer_to_device(buf)
        type(csp_buffer_t), intent(inout) :: buf

        !$acc enter data copyin(buf)
        !$acc enter data copyin(buf%ssp_weights, buf%spec_young, buf%spec_old)
        !$acc enter data copyin(buf%emlin_young, buf%emlin_old)
        !$acc enter data attach(buf%ssp_weights)
        !$acc enter data attach(buf%spec_young)
        !$acc enter data attach(buf%spec_old)
        !$acc enter data attach(buf%emlin_young)
        !$acc enter data attach(buf%emlin_old)
    end subroutine map_csp_buffer_to_device

    subroutine unmap_csp_buffer_from_device(buf)
        type(csp_buffer_t), intent(inout) :: buf

        !$acc exit data delete(buf%ssp_weights, buf%spec_young, buf%spec_old)
        !$acc exit data delete(buf%emlin_young, buf%emlin_old)
        !$acc exit data delete(buf)
    end subroutine unmap_csp_buffer_from_device

    subroutine setup_nebular_line_context(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: i

        ctx%add_neb_emission_val = 1
        ctx%add_neb_continuum_val = 0
        ctx%nebemlineinspec_val = 1
        ctx%setup_nebular_gaussians_val = 0
        ctx%smooth_velocity_val = 0

        if (.not. associated(ctx%state%neb_res_min)) then
            allocate(ctx%state%neb_res_min(NEMLINE))
        end if
        if (.not. associated(ctx%state%gaussnebarr)) then
            allocate(ctx%state%gaussnebarr(ctx%state%nspec, NEMLINE))
        end if

        ctx%state%neb_res_min = 1.0_wp
        ctx%state%nebem_line_pos = 6563.0_wp

        do i = 1, size(ctx%state%nebem_logz)
            ctx%state%nebem_logz(i) = -2.0_wp + 0.2_wp * real(i - 1, WP)
        end do
        do i = 1, size(ctx%state%nebem_logu)
            ctx%state%nebem_logu(i) = -4.0_wp + 0.5_wp * real(i - 1, WP)
        end do
        do i = 1, size(ctx%state%nebem_age)
            ctx%state%nebem_age(i) = 6.0_wp + 0.4_wp * real(i - 1, WP)
        end do

        ctx%state%nebem_line = -99.0_wp
        ctx%state%nebem_line(1, :, :, :) = 5.0_wp

        if (associated(ctx%state%neb_res_min)) then
            !$acc enter data copyin(ctx%state%neb_res_min)
            !$acc enter data attach(ctx%state%neb_res_min)
        end if
        if (associated(ctx%state%gaussnebarr)) then
            !$acc enter data copyin(ctx%state%gaussnebarr)
            !$acc enter data attach(ctx%state%gaussnebarr)
        end if
        !$acc update device(ctx%add_neb_emission_val, ctx%add_neb_continuum_val, ctx%nebemlineinspec_val)
        !$acc update device(ctx%setup_nebular_gaussians_val, ctx%smooth_velocity_val)
        !$acc update device(ctx%state%nebem_line_pos, ctx%state%nebem_logz, ctx%state%nebem_logu)
        !$acc update device(ctx%state%nebem_age, ctx%state%nebem_line)
    end subroutine setup_nebular_line_context

    ! ---------------------------------------------------------------------
    ! Group 1: SFH Weight Logic
    ! ---------------------------------------------------------------------
    subroutine test_convert_sfhparams_units()
        type(params) :: pset
        type(sfhparams) :: sfh
        real(WP) :: tage

        call print_group("SFH Params: Unit Conversion")

        pset%sfh = 1
        pset%tau = 2.0_wp
        pset%tburst = 3.0_wp
        pset%sf_trunc = 1.0_wp
        pset%sf_start = 0.5_wp
        pset%sf_slope = 0.1_wp
        tage = 5.0_wp

        call convert_sfhparams(pset, tage, sfh)

        call assert_float_equals(2.0e9_wp, sfh%tau, EPS, &
                                 "tau in years", total_tests, total_failures)
        call assert_float_equals((tage - pset%sf_start) * 1.0e9_wp, sfh%tage, EPS, &
                                 "tage offset conversion", total_tests, total_failures)
        call assert_float_equals((tage - pset%tburst) * 1.0e9_wp, sfh%tb, EPS, &
                                 "tburst lookback", total_tests, total_failures)
    end subroutine test_convert_sfhparams_units

    subroutine test_single_ssp_dirac()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:), weights(:,:)
        real(WP), allocatable :: spec_lambda(:)
        integer :: n, idx_max, idx_tage
        real(WP) :: log_tage

        call print_group("SFH Weights: Single SSP (Dirac)")

        n = 80
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(2))
        spec_lambda = [5500.0_wp, 5600.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(weights(n, 1))

        pset%sfh = 0
        log_tage = log10(5.0_wp * 1.0e9_wp)

        call compute_sfh_weights(ctx, pset, 5.0_wp, 1, weights)

        call assert_float_equals(1.0_wp, sum(weights(:,1)), EPS, &
                                 "SSP weights sum to 1", total_tests, total_failures)
        call assert_true(count(weights(:,1) > 1.0e-8_wp) <= 2, &
                         "SSP weights localized to one/two bins", total_tests, total_failures)

        idx_max = maxloc(weights(:,1), dim=1)
        idx_tage = find_interval(time_full, log_tage)
        call assert_true(idx_max == idx_tage .or. idx_max == idx_tage + 1, &
                         "SSP peak at target age bin", total_tests, total_failures)

        deallocate(time_full, weights, spec_lambda)
        call teardown_basic_context(ctx)
    end subroutine test_single_ssp_dirac

    subroutine test_exponential_tau_ratio()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:), weights(:,:)
        real(WP), allocatable :: spec_lambda(:)
        integer :: n, i1, i2
        real(WP) :: ratio

        call print_group("SFH Weights: Exponential Tau Ratio")

        n = 120
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(2))
        spec_lambda = [5500.0_wp, 5600.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(weights(n, 1))

        pset%sfh = 1
        pset%tau = 1.0_wp

        call compute_sfh_weights(ctx, pset, 10.0_wp, 1, weights)

        call assert_float_equals(1.0_wp, sum(weights(:,1)), EPS, &
                                 "Tau weights sum to 1", total_tests, total_failures)

        i1 = nearest_index(time_full, log10(1.0e9_wp))
        i2 = nearest_index(time_full, log10(2.0e9_wp))

        call assert_true(weights(i1,1) > SAFE_FLOOR .and. weights(i2,1) > SAFE_FLOOR, &
                         "Non-zero tau weights", total_tests, total_failures)

        ratio = weights(i1,1) / weights(i2,1)
        call assert_relative_error(0.5_wp * exp(-1.0_wp), ratio, 0.2_wp, &
                                   "Tau ratio ~ 0.5 * e^-1", total_tests, total_failures)

        deallocate(time_full, weights, spec_lambda)
        call teardown_basic_context(ctx)
    end subroutine test_exponential_tau_ratio

    subroutine test_delayed_tau_peak()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:), weights(:,:)
        real(WP), allocatable :: spec_lambda(:)
        integer :: n, i_peak

        call print_group("SFH Weights: Delayed Tau Peak")

        n = 100
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(2))
        spec_lambda = [5500.0_wp, 5600.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(weights(n, 1))

        pset%sfh = 4
        pset%tau = 2.0_wp

        call compute_sfh_weights(ctx, pset, 10.0_wp, 1, weights)

        i_peak = nearest_index(time_full, log10(8.0e9_wp))

        call assert_true(weights(i_peak,1) > weights(1,1), &
                         "Delayed tau: Peak > Start", total_tests, total_failures)
        call assert_true(weights(i_peak,1) > weights(n,1), &
                         "Delayed tau: Peak > End", total_tests, total_failures)

        deallocate(time_full, weights, spec_lambda)
        call teardown_basic_context(ctx)
    end subroutine test_delayed_tau_peak

    subroutine test_burst_addition()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:), weights_base(:,:), weights_burst(:,:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP) :: fburst, excess_mass
        integer :: n, i_burst

        call print_group("SFH Weights: Burst Addition")

        n = 100
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(2))
        spec_lambda = [5500.0_wp, 5600.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(weights_base(n, 1), weights_burst(n, 1))

        pset%sfh = 1
        pset%tau = 1.0_wp
        pset%tburst = 9.0_wp
        fburst = 0.1_wp

        pset%fburst = 0.0_wp
        call compute_sfh_weights(ctx, pset, 10.0_wp, 1, weights_base)

        pset%fburst = fburst
        call compute_sfh_weights(ctx, pset, 10.0_wp, 1, weights_burst)

        call assert_float_equals(1.0_wp, sum(weights_burst(:,1)), EPS, &
                                 "Burst weights sum to 1", total_tests, total_failures)

        excess_mass = sum(weights_burst(:,1) - (1.0_wp - fburst) * weights_base(:,1))
        call assert_relative_error(fburst, excess_mass, 1.0e-3_wp, &
                                   "Burst excess mass ~ fburst", total_tests, total_failures)

        i_burst = nearest_index(time_full, log10(1.0e9_wp))
        call assert_true(weights_burst(i_burst,1) > weights_base(i_burst,1), &
                         "Burst creates localized spike", total_tests, total_failures)

        deallocate(time_full, weights_base, weights_burst, spec_lambda)
        call teardown_basic_context(ctx)
    end subroutine test_burst_addition

    subroutine test_tabular_zmix_interpolation()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        real(WP), allocatable :: time_full(:), weights(:,:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP) :: mass_z1, mass_z2, ratio, dz, zbin
        integer :: n

        call print_group("SFH Weights: Tabular Z-Mix Interpolation")

        n = 60
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 9.5_wp)

        allocate(spec_lambda(2))
        spec_lambda = [5500.0_wp, 5600.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 2, 1)

        ctx%state%zlegend(1) = 0.01_wp
        ctx%state%zlegend(2) = 0.03_wp

        ctx%state%ntabsfh = 2
        ctx%state%sfh_tab(:, 1) = [0.0_wp, 1.0_wp, 0.019_wp]
        ctx%state%sfh_tab(:, 2) = [1.0e9_wp, 1.0_wp, 0.019_wp]

        allocate(weights(n, 2))

        pset%sfh = 3
        call compute_sfh_weights(ctx, pset, 1.0_wp, 2, weights)

        mass_z1 = sum(weights(:,1))
        mass_z2 = sum(weights(:,2))
        ratio = mass_z2 / max(mass_z1 + mass_z2, SAFE_FLOOR)

        zbin = 0.019_wp
        dz = (log10(zbin) - log10(ctx%state%zlegend(1))) / &
             (log10(ctx%state%zlegend(2)) - log10(ctx%state%zlegend(1)))

        call assert_relative_error(dz, ratio, 2.0e-2_wp, &
                                   "Metallicity split matches log interpolation", total_tests, total_failures)

        deallocate(time_full, weights, spec_lambda)
        call teardown_basic_context(ctx)
    end subroutine test_tabular_zmix_interpolation

    ! ---------------------------------------------------------------------
    ! Group 2: Integration & Young/Old Split
    ! ---------------------------------------------------------------------
    subroutine test_all_old_limit()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(csp_buffer_t) :: buf
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: ssp_grid(:,:,:), emlin_grid(:,:,:), mass_ssp(:,:), ssp_lum(:,:)
        real(WP) :: mass_csp, lbol_csp
        integer :: n

        call print_group("CSP Integrator: All Old Limit")

        n = 40
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(1))
        spec_lambda = [5500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(ssp_grid(1, n, 1))
        allocate(emlin_grid(NEMLINE, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(ssp_lum(n, 1))

        ssp_grid = 1.0_wp
        emlin_grid = 0.0_wp
        mass_ssp = 1.0_wp
        ssp_lum = 1.0_wp

        call init_csp_buffer(buf, 1, n, 1)
        call map_csp_buffer_to_device(buf)

        pset%sfh = 0
        pset%dust_tesc = 3.0_wp

        !$acc data copyin(ssp_grid, emlin_grid)
        call integrate_csp_step(ctx, pset, 10.0_wp, 1, ssp_grid, emlin_grid, mass_ssp, ssp_lum, buf, mass_csp, lbol_csp)
        !$acc end data
        call unmap_csp_buffer_from_device(buf)

        call assert_true(all(abs(buf%spec_young) <= 1.0e-8_wp), "Young component ~0", total_tests, total_failures)
        call assert_true(sum(buf%spec_old) > 0.0_wp, "Old component carries flux", total_tests, total_failures)

        call free_csp_buffer(buf)
        deallocate(time_full, spec_lambda, ssp_grid, emlin_grid, mass_ssp, ssp_lum)
        call teardown_basic_context(ctx)
    end subroutine test_all_old_limit

    subroutine test_all_young_limit()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(csp_buffer_t) :: buf
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: ssp_grid(:,:,:), emlin_grid(:,:,:), mass_ssp(:,:), ssp_lum(:,:)
        real(WP) :: mass_csp, lbol_csp
        integer :: n

        call print_group("CSP Integrator: All Young Limit")

        n = 40
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(1))
        spec_lambda = [5500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(ssp_grid(1, n, 1))
        allocate(emlin_grid(NEMLINE, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(ssp_lum(n, 1))

        ssp_grid = 1.0_wp
        emlin_grid = 0.0_wp
        mass_ssp = 1.0_wp
        ssp_lum = 1.0_wp

        call init_csp_buffer(buf, 1, n, 1)
        call map_csp_buffer_to_device(buf)

        pset%sfh = 1
        pset%tau = 2.0_wp
        pset%dust_tesc = 10.15_wp

        !$acc data copyin(ssp_grid, emlin_grid)
        call integrate_csp_step(ctx, pset, 10.0_wp, 1, ssp_grid, emlin_grid, mass_ssp, ssp_lum, buf, mass_csp, lbol_csp)
        !$acc end data
        call unmap_csp_buffer_from_device(buf)

        call assert_true(all(abs(buf%spec_old) <= 1.0e-8_wp), "Old component ~0", total_tests, total_failures)
        call assert_true(sum(buf%spec_young) > 0.0_wp, "Young component carries flux", total_tests, total_failures)

        call free_csp_buffer(buf)
        deallocate(time_full, spec_lambda, ssp_grid, emlin_grid, mass_ssp, ssp_lum)
        call teardown_basic_context(ctx)
    end subroutine test_all_young_limit

    subroutine test_split_boundary()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(csp_buffer_t) :: buf
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: ssp_grid(:,:,:), emlin_grid(:,:,:), mass_ssp(:,:), ssp_lum(:,:)
        real(WP) :: mass_csp, lbol_csp, expected_young, expected_old
        integer :: n, k

        call print_group("CSP Integrator: Split Boundary")

        n = 50
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(1))
        spec_lambda = [5500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(ssp_grid(1, n, 1))
        allocate(emlin_grid(NEMLINE, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(ssp_lum(n, 1))

        ssp_grid = 1.0_wp
        emlin_grid = 0.0_wp
        mass_ssp = 1.0_wp
        ssp_lum = 1.0_wp

        call init_csp_buffer(buf, 1, n, 1)
        call map_csp_buffer_to_device(buf)

        pset%sfh = 1
        pset%tau = 2.0_wp
        pset%const = 1.0_wp
        pset%dust_tesc = 7.0_wp

        !$acc data copyin(ssp_grid, emlin_grid)
        call integrate_csp_step(ctx, pset, 10.0_wp, 1, ssp_grid, emlin_grid, mass_ssp, ssp_lum, buf, mass_csp, lbol_csp)
        !$acc end data
        call unmap_csp_buffer_from_device(buf)

        k = find_interval(time_full, 7.0_wp)
        expected_young = sum(buf%ssp_weights(1:k, 1))
        expected_old = sum(buf%ssp_weights(k+1:n, 1))

        call assert_relative_error(expected_young, buf%spec_young(1), REL_EPS, &
                                   "Young sum matches weights", total_tests, total_failures)
        call assert_relative_error(expected_old, buf%spec_old(1), REL_EPS, &
                                   "Old sum matches weights", total_tests, total_failures)

        call free_csp_buffer(buf)
        deallocate(time_full, spec_lambda, ssp_grid, emlin_grid, mass_ssp, ssp_lum)
        call teardown_basic_context(ctx)
    end subroutine test_split_boundary

    subroutine test_mass_sum_consistency()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(csp_buffer_t) :: buf
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: ssp_grid(:,:,:), emlin_grid(:,:,:), mass_ssp(:,:), ssp_lum(:,:)
        real(WP) :: mass_csp, lbol_csp
        integer :: n

        call print_group("CSP Integrator: Mass Sum Consistency")

        n = 60
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(1))
        spec_lambda = [5500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(ssp_grid(1, n, 1))
        allocate(emlin_grid(NEMLINE, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(ssp_lum(n, 1))

        ssp_grid = 1.0_wp
        emlin_grid = 0.0_wp
        mass_ssp = 1.0_wp
        ssp_lum = 1.0_wp

        call init_csp_buffer(buf, 1, n, 1)
        call map_csp_buffer_to_device(buf)

        pset%sfh = 1
        pset%tau = 1.0_wp
        pset%dust_tesc = 7.0_wp

        !$acc data copyin(ssp_grid, emlin_grid)
        call integrate_csp_step(ctx, pset, 10.0_wp, 1, ssp_grid, emlin_grid, mass_ssp, ssp_lum, buf, mass_csp, lbol_csp)
        !$acc end data
        call unmap_csp_buffer_from_device(buf)

        call assert_relative_error(sum(buf%ssp_weights(:,1)), mass_csp, REL_EPS, &
                                   "Mass matches sum of weights", total_tests, total_failures)
        call assert_relative_error(1.0_wp, mass_csp, 1.0e-3_wp, &
                                   "Mass ~1 for normalized weights", total_tests, total_failures)

        call free_csp_buffer(buf)
        deallocate(time_full, spec_lambda, ssp_grid, emlin_grid, mass_ssp, ssp_lum)
        call teardown_basic_context(ctx)
    end subroutine test_mass_sum_consistency

    ! ---------------------------------------------------------------------
    ! Group 3: Physics & Post-Processing
    ! ---------------------------------------------------------------------
    subroutine test_dust_screen_logic()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(csp_buffer_t) :: buf
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: spec_total(:), emlin_total(:)
        real(WP) :: mdust

        call print_group("Dust: Birth Cloud vs Diffuse Screen")

        allocate(time_full(2))
        time_full = [6.0_wp, 7.0_wp]

        allocate(spec_lambda(3))
        spec_lambda = [5500.0_wp, 5500.0_wp, 5500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        call init_csp_buffer(buf, size(spec_lambda), 2, 1)

        buf%spec_young = 1.0_wp
        buf%spec_old = 1.0_wp
        buf%emlin_young = 0.0_wp
        buf%emlin_old = 0.0_wp

        call map_csp_buffer_to_device(buf)

        allocate(spec_total(size(spec_lambda)))
        allocate(emlin_total(NEMLINE))

        pset%dust1 = 1.0_wp
        pset%dust2 = 0.0_wp
        pset%dust3 = 0.0_wp
        pset%dust1_index = -1.0_wp
        pset%frac_obrun = 0.0_wp
        pset%frac_nodust = 0.0_wp
        pset%uvb = 1.0_wp
        pset%wgp1 = 1
        pset%wgp2 = 1
        pset%wgp3 = 1

        !$acc data copyin(pset) copy(spec_total, emlin_total)
        call apply_dust_physics(ctx, pset, buf, spec_total, emlin_total, mdust)
        !$acc end data

        call assert_true(spec_total(1) > 1.0_wp .and. spec_total(1) < 2.0_wp, &
                         "Young attenuated, old unattenuated", total_tests, total_failures)
        call assert_relative_error(exp(-1.0_wp), spec_total(1) - 1.0_wp, 1.0e-3_wp, &
                                   "Young component attenuated", total_tests, total_failures)

        deallocate(time_full, spec_lambda, spec_total, emlin_total)
        call unmap_csp_buffer_from_device(buf)
        call free_csp_buffer(buf)
        call teardown_basic_context(ctx)
    end subroutine test_dust_screen_logic

    subroutine test_igm_absorption_toggle()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout) :: result_no, result_yes
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: spec_no(:), spec_yes(:), emlines(:), igm(:)
        real(WP) :: mass_csp, lbol_csp, mdust
        integer :: i_blue, i_red

        call print_group("IGM: Absorption Toggle")

        allocate(time_full(3))
        time_full = [6.0_wp, 8.0_wp, 10.0_wp]

        allocate(spec_lambda(4))
        spec_lambda = [900.0_wp, 1100.0_wp, 1300.0_wp, 1500.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(spec_no(size(spec_lambda)))
        allocate(spec_yes(size(spec_lambda)))
        allocate(emlines(NEMLINE))
        spec_no = 1.0_wp
        spec_yes = 1.0_wp
        emlines = 0.0_wp

        allocate(result_no%mags(ctx%state%nbands))
        allocate(result_no%indx(ctx%state%nindx))
        allocate(result_yes%mags(ctx%state%nbands))
        allocate(result_yes%indx(ctx%state%nindx))

        pset%sfh = 0
        pset%tage = 5.0_wp
        pset%zred = 3.0_wp
        pset%igm_factor = 1.0_wp

        mass_csp = 1.0_wp
        lbol_csp = 0.0_wp
        mdust = 0.0_wp

        !$acc data copyin(pset) copy(spec_no, emlines)
        call apply_post_processing(ctx, pset, pset%tage, mass_csp, lbol_csp, mdust, spec_no, emlines, result=result_no)
        !$acc end data

        igm = get_igm_transmission(spec_lambda, pset%zred, pset%igm_factor)
        !$acc data copyin(pset, igm) copy(spec_yes, emlines)
        call apply_post_processing(ctx, pset, pset%tage, mass_csp, lbol_csp, mdust, spec_yes, emlines, igm, result_yes)
        !$acc end data

        i_blue = 1
        i_red = 4

        call assert_true(result_yes%spec(i_blue) < result_no%spec(i_blue), &
                         "Blueward flux attenuated", total_tests, total_failures)
        call assert_relative_error(result_no%spec(i_red), result_yes%spec(i_red), 1.0e-6_wp, &
                                   "Redward flux unchanged", total_tests, total_failures)

        deallocate(time_full, spec_lambda, spec_no, spec_yes, emlines, igm)
        call teardown_basic_context(ctx)
    end subroutine test_igm_absorption_toggle

    subroutine test_smoothing_conserves_flux()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout) :: result
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: spec(:), emlines(:)
        real(WP) :: mass_csp, lbol_csp, mdust
        real(WP) :: area_before, area_after, peak_before, peak_after
        integer :: n

        call print_group("Smoothing: Flux Conservation")

        n = 120
        allocate(time_full(3))
        time_full = [6.0_wp, 8.0_wp, 10.0_wp]

        allocate(spec_lambda(n))
        call fill_log_grid(spec_lambda, 3.6_wp, 3.85_wp)
        spec_lambda = 10.0_wp**spec_lambda

        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(spec(n))
        allocate(emlines(NEMLINE))
        spec = 0.0_wp
        spec(n/2) = 1.0_wp
        emlines = 0.0_wp

        allocate(result%mags(ctx%state%nbands))
        allocate(result%indx(ctx%state%nindx))

        pset%sfh = 0
        pset%tage = 5.0_wp
        pset%sigma_smooth = 200.0_wp
        pset%min_wave_smooth = minval(spec_lambda)
        pset%max_wave_smooth = maxval(spec_lambda)

        ctx%smooth_velocity_val = 1

        area_before = sum(spec)
        peak_before = maxval(spec)

        mass_csp = 1.0_wp
        lbol_csp = 0.0_wp
        mdust = 0.0_wp

        !$acc data copyin(pset) copy(spec, emlines)
        call apply_post_processing(ctx, pset, pset%tage, mass_csp, lbol_csp, mdust, spec, emlines, result=result)
        !$acc end data

        area_after = sum(result%spec)
        peak_after = maxval(result%spec)

        call assert_true(peak_after < peak_before, "Smoothing lowers peak", total_tests, total_failures)
        call assert_relative_error(area_before, area_after, 1.0e-2_wp, "Smoothing conserves flux", total_tests, total_failures)

        deallocate(time_full, spec_lambda, spec, emlines)
        call teardown_basic_context(ctx)
    end subroutine test_smoothing_conserves_flux

    ! ---------------------------------------------------------------------
    ! Group 4: Driver Logic & Output Modes
    ! ---------------------------------------------------------------------
    subroutine test_snapshot_vs_history()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout), allocatable :: results(:)
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: tspec_ssp(:,:,:), mass_ssp(:,:), lbol_ssp(:,:)
        integer :: n, nspec

        call print_group("Driver: Snapshot vs History")

        n = 6
        nspec = 5
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(nspec))
        spec_lambda = [800.0_wp, 1200.0_wp, 3000.0_wp, 5500.0_wp, 9000.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(tspec_ssp(nspec, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(lbol_ssp(n, 1))

        tspec_ssp = 1.0_wp
        mass_ssp = 1.0_wp
        lbol_ssp = 0.0_wp

        pset%sfh = 0
        pset%tage = 5.0_wp

        call compute_csp_scenario(ctx, pset, 1, tspec_ssp, mass_ssp, lbol_ssp, results)
        call assert_int_equals(1, size(results), "Snapshot returns single output", total_tests, total_failures)
        call assert_float_equals(log10(5.0_wp * 1.0e9_wp), results(1)%age, 1.0e-5_wp, &
                                 "Snapshot age matches", total_tests, total_failures)

        deallocate(results)

        pset%tage = 0.0_wp
        call compute_csp_scenario(ctx, pset, 1, tspec_ssp, mass_ssp, lbol_ssp, results)
        call assert_int_equals(n, size(results), "History returns ntfull outputs", total_tests, total_failures)
        call assert_float_equals(time_full(1), results(1)%age, 1.0e-5_wp, &
                                 "History age(1) matches grid", total_tests, total_failures)
        call assert_float_equals(time_full(n), results(n)%age, 1.0e-5_wp, &
                                 "History age(end) matches grid", total_tests, total_failures)

        deallocate(results)
        deallocate(time_full, spec_lambda, tspec_ssp, mass_ssp, lbol_ssp)
        call teardown_basic_context(ctx)
    end subroutine test_snapshot_vs_history

    subroutine test_mass_normalization_surviving()
        type(fsps_context_t), allocatable :: ctx
        type(params) :: pset
        type(compspout), allocatable :: results(:)
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: tspec_ssp(:,:,:), mass_ssp(:,:), lbol_ssp(:,:)
        integer :: n, nspec

        call print_group("Driver: Mass Normalization (Surviving)")

        n = 8
        nspec = 5
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.15_wp)

        allocate(spec_lambda(nspec))
        spec_lambda = [800.0_wp, 1200.0_wp, 3000.0_wp, 5500.0_wp, 9000.0_wp]
        call setup_basic_context(ctx, time_full, spec_lambda, 1, 1)

        allocate(tspec_ssp(nspec, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(lbol_ssp(n, 1))

        tspec_ssp = 1.0_wp
        mass_ssp = 1.0_wp
        lbol_ssp = 0.0_wp

        pset%sfh = 1
        pset%tau = 2.0_wp
        pset%tage = 13.0_wp

        call compute_csp_scenario(ctx, pset, 1, tspec_ssp, mass_ssp, lbol_ssp, results)

        call assert_true(results(1)%mass_csp < 1.0_wp, "Surviving mass < 1", total_tests, total_failures)
        call assert_relative_error(1.0_wp, results(1)%mformed, 1.0e-6_wp, &
                                   "Formed mass normalized to 1", total_tests, total_failures)

        deallocate(results)
        deallocate(time_full, spec_lambda, tspec_ssp, mass_ssp, lbol_ssp)
        call teardown_basic_context(ctx)
    end subroutine test_mass_normalization_surviving

    subroutine test_nebular_precalculation()
        type(fsps_context_t), allocatable :: ctx_no, ctx_yes
        type(params) :: pset
        type(compspout), allocatable :: res_no(:), res_yes(:)
        real(WP), allocatable :: time_full(:)
        real(WP), allocatable :: spec_lambda(:)
        real(WP), allocatable :: tspec_ssp(:,:,:), mass_ssp(:,:), lbol_ssp(:,:)
        integer :: n, nspec, line_idx

        call print_group("Driver: Nebular Pre-Calculation")

        n = 6
        nspec = 7
        allocate(time_full(n))
        call fill_log_grid(time_full, 6.0_wp, 10.0_wp)

        allocate(spec_lambda(nspec))
        spec_lambda = [800.0_wp, 6000.0_wp, 6400.0_wp, 6563.0_wp, 6700.0_wp, 7000.0_wp, 9000.0_wp]

        call setup_basic_context(ctx_no, time_full, spec_lambda, 1, 1)
        call setup_basic_context(ctx_yes, time_full, spec_lambda, 1, 1)

        call setup_nebular_line_context(ctx_yes)

        allocate(tspec_ssp(nspec, n, 1))
        allocate(mass_ssp(n, 1))
        allocate(lbol_ssp(n, 1))

        tspec_ssp = 1.0_wp
        mass_ssp = 1.0_wp
        lbol_ssp = 0.0_wp

        pset%sfh = 0
        pset%tage = 5.0_wp
        pset%gas_logz = -1.0_wp
        pset%gas_logu = -2.0_wp
        pset%sigma_smooth = 1.0_wp
        pset%frac_obrun = 0.0_wp
        ctx_yes%state%whlylim = 2

        call compute_csp_scenario(ctx_no, pset, 1, tspec_ssp, mass_ssp, lbol_ssp, res_no)
        call compute_csp_scenario(ctx_yes, pset, 1, tspec_ssp, mass_ssp, lbol_ssp, res_yes)

        line_idx = nearest_index(spec_lambda, 6563.0_wp)
        
        call assert_true(res_yes(1)%spec(line_idx) > res_no(1)%spec(line_idx), &
                         "Nebular lines increase flux", total_tests, total_failures)

        deallocate(res_no, res_yes)
        deallocate(time_full, spec_lambda, tspec_ssp, mass_ssp, lbol_ssp)
        call teardown_basic_context(ctx_no)
        call teardown_basic_context(ctx_yes)
    end subroutine test_nebular_precalculation

end module test_fsps_csp_mod
