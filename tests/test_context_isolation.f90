program test_context_isolation
    !> Integration tests for contexts that share a setup cache entry.
    !>
    !> Contexts with identical library selections share one read-only setup cache (stellar
    !> libraries, isochrones, tables). Anything a context writes after setup, or that depends on a
    !> context-specific setting, must not leak into another context through that cache:
    !>   1. setup-time nebular line kernels built from nebular_smooth_init,
    !>   2. per-call nebular line kernels rebuilt from each context's sigma_smooth,
    !>   3. cache reference counting across create/destroy.
    !> Requires the FSPS data tree (FSPS_DATA_HOME), like the regression tests.

    use fsps_precision, only: WP
    use fsps_api, only: fsps_create, fsps_setup, fsps_destroy, fsps_set_param_int, fsps_set_param_float
    use fsps_context_types, only: fsps_context_t
    use test_utils_mod, only: print_major_header, print_group, print_summary_line, assert_true, &
                              assert_int_equals
    implicit none

    integer, parameter :: ZMET = 10
    integer :: total_tests = 0
    integer :: total_failures = 0

    call print_major_header("Context Isolation (shared setup cache)")

    call test_setup_kernels_follow_own_smoothing()
    call test_per_call_kernels_are_private()
    call test_cache_refcount_bookkeeping()

    call print_summary_line("Context Isolation", total_tests - total_failures, total_tests)
    if (total_failures > 0) stop 1

contains

    !> Configures a context for nebular line emission and runs setup.
    subroutine setup_nebular_context(ctx, setup_gaussians, smooth_init)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: setup_gaussians
        real(WP), intent(in) :: smooth_init
        integer :: status

        call fsps_create(ctx)
        call fsps_set_param_int(ctx, 'add_neb_emission', 1, status)
        call fsps_set_param_int(ctx, 'setup_nebular_gaussians', setup_gaussians, status)
        call fsps_set_param_float(ctx, 'nebular_smooth_init', smooth_init, status)
        call fsps_setup(ctx, ZMET, 'mist', 'miles', 'DL07')
    end subroutine setup_nebular_context

    !> With setup_nebular_gaussians = 1 the line kernels are built once at setup from
    !> nebular_smooth_init. A second context with a different nebular_smooth_init must neither
    !> overwrite the first context's kernels nor inherit them.
    subroutine test_setup_kernels_follow_own_smoothing()
        type(fsps_context_t) :: ctx_a, ctx_b
        real(WP), allocatable :: kernels_a(:, :)

        call print_group("Setup-time line kernels follow each context's nebular_smooth_init")

        call setup_nebular_context(ctx_a, 1, 100.0_wp)
        kernels_a = ctx_a%state%gaussnebarr

        call setup_nebular_context(ctx_b, 1, 500.0_wp)

        call assert_true(all(ctx_a%state%gaussnebarr == kernels_a), &
                         "context A kernels unchanged by context B setup", total_tests, total_failures)
        call assert_true(any(ctx_b%state%gaussnebarr /= kernels_a), &
                         "context B kernels built from its own nebular_smooth_init", total_tests, total_failures)

        call fsps_destroy(ctx_b)
        call fsps_destroy(ctx_a)
    end subroutine test_setup_kernels_follow_own_smoothing

    !> With setup_nebular_gaussians = 0 (the default) every nebular call rebuilds the line kernels
    !> from that context's sigma_smooth, so the kernel array must be owned by the context; the
    !> read-only spectral library must still be shared.
    subroutine test_per_call_kernels_are_private()
        type(fsps_context_t) :: ctx_c, ctx_d

        call print_group("Per-call line kernels are private to each context")

        call setup_nebular_context(ctx_c, 0, 100.0_wp)
        call setup_nebular_context(ctx_d, 0, 100.0_wp)

        call assert_true(.not. associated(ctx_c%state%gaussnebarr, ctx_d%state%gaussnebarr), &
                         "per-call kernel arrays are distinct", total_tests, total_failures)
        call assert_true(associated(ctx_c%state%speclib, ctx_d%state%speclib), &
                         "read-only spectral library is still shared", total_tests, total_failures)

        call fsps_destroy(ctx_d)
        call assert_true(associated(ctx_c%state%gaussnebarr), &
                         "kernels survive destruction of the other context", total_tests, total_failures)
        call fsps_destroy(ctx_c)
    end subroutine test_per_call_kernels_are_private

    !> Contexts with the same key share one cache entry whose reference count tracks the number of
    !> live contexts; destroying one must leave the shared data intact for the others.
    subroutine test_cache_refcount_bookkeeping()
        type(fsps_context_t) :: ctx_e, ctx_f, ctx_g

        call print_group("Setup cache reference counting")

        call setup_nebular_context(ctx_e, 0, 100.0_wp)
        call setup_nebular_context(ctx_f, 0, 100.0_wp)
        call setup_nebular_context(ctx_g, 0, 100.0_wp)

        call assert_true(associated(ctx_e%setup_cache, ctx_g%setup_cache), &
                         "identical selections share one cache entry", total_tests, total_failures)
        call assert_int_equals(3, ctx_e%setup_cache%refcount, "refcount = 3 live contexts", &
                               total_tests, total_failures)

        call fsps_destroy(ctx_g)
        call assert_int_equals(2, ctx_e%setup_cache%refcount, "refcount = 2 after one destroy", &
                               total_tests, total_failures)
        call fsps_destroy(ctx_f)
        call assert_int_equals(1, ctx_e%setup_cache%refcount, "refcount = 1 after two destroys", &
                               total_tests, total_failures)
        call assert_true(associated(ctx_e%state%speclib), "surviving context keeps shared data", &
                         total_tests, total_failures)

        call fsps_destroy(ctx_e)
    end subroutine test_cache_refcount_bookkeeping

end program test_context_isolation
