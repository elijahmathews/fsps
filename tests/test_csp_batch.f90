program test_csp_batch
    !> The batched CSP (compute_csp_batch) against compute_csp_scenario, model by model: spectra,
    !> magnitudes, mass, L_bol, dust mass and SFR, for SFH types 0/1/4/5, birth-cloud and diffuse
    !> dust (dust types 0, 1, 2, 4), dust emission (Draine & Li parameters varied per model), and
    !> redshifts 0, 0.5 (with IGM absorption) and 1.2. Requires the FSPS data tree, like the
    !> regression tests.

    use fsps_precision, only: WP, KP
    use fsps_types, only: PARAMS, COMPSPOUT
    use fsps_api, only: fsps_create, fsps_setup, fsps_destroy
    use fsps_csp, only: compute_csp_scenario
    use fsps_csp_batch, only: compute_csp_batch
    use fsps_ssp, only: generate_ssp_grid
    use fsps_context, only: fsps_context_update_ssp_basis, fsps_context_move_to_device, &
                            fsps_context_set_pset, fsps_context_get_pset, fsps_context_prepare_pset, &
                            fsps_context_set_param_int, fsps_context_set_param_float
    use fsps_context_types, only: fsps_context_t
    use test_utils_mod, only: tap_finish, print_major_header, print_group, print_summary_line, assert_true
    implicit none

    integer, parameter :: N_MOD = 12
    integer :: total_tests = 0
    integer :: total_failures = 0
    ! Spectra and magnitudes: the batched path sums in a different order (and integrates the
    ! filters over each band's support), so they agree to the rounding of the kernel precision
    ! FP32 kernels: 100 eps(KP), as elsewhere; OpenACC multicore runs the FP32 reductions as
    ! plain per-core sums (~4e-6), the GPU as trees (~3e-7)
    real(WP), parameter :: SPEC_RTOL = merge(1.0e-10_wp, 100.0_wp * real(epsilon(1.0_KP), WP), KP == WP)
    real(WP), parameter :: MAG_TOL = merge(1.0e-9_wp, 1.0e-5_wp, KP == WP)
    type(fsps_context_t), allocatable :: ctx

    call print_major_header("Batched CSP vs compute_csp_scenario")

    allocate(ctx)
    call fsps_create(ctx)
    ctx%pset%zmet = 10
    call fsps_setup(ctx, 10, isoc_type_in='mist', spec_type_in='miles', dust_type_in='DL07')
    call fsps_context_move_to_device(ctx)
    call set_int('zmet', 10)
    call set_int('imf_type', 1)
    call set_int('add_dust_emission', 1)
    call build_basis()

    call check_batch(2, 0.0_wp, 0, "Calzetti, z = 0")
    call check_batch(0, 0.5_wp, 1, "power law, z = 0.5, IGM")
    call check_batch(1, 0.0_wp, 0, "Cardelli et al., z = 0")
    call check_batch(4, 1.2_wp, 0, "Kriek & Conroy, z = 1.2")

    call fsps_destroy(ctx)
    call print_summary_line("Batched CSP", total_tests - total_failures, total_tests)
    call tap_finish()
    if (total_failures > 0) stop 1

contains

    !> SSP basis for the default (Chabrier, solar) parameters; the batch varies only CSP parameters.
    subroutine build_basis()
        type(PARAMS) :: pset
        real(WP), allocatable :: mass_ssp(:), lbol_ssp(:), spec_ssp(:, :)

        call fsps_context_prepare_pset(ctx)
        call fsps_context_get_pset(ctx, pset)
        allocate(mass_ssp(ctx%state%ntfull), lbol_ssp(ctx%state%ntfull), &
                 spec_ssp(ctx%state%nspec, ctx%state%ntfull))
        call fsps_context_set_pset(ctx, pset)
        call generate_ssp_grid(ctx, pset, mass_ssp, lbol_ssp, spec_ssp)
        call fsps_context_update_ssp_basis(ctx, reshape(spec_ssp, [size(spec_ssp, 1), size(spec_ssp, 2), 1]), &
                                           reshape(mass_ssp, [size(mass_ssp), 1]), &
                                           reshape(lbol_ssp, [size(lbol_ssp), 1]), 1)
    end subroutine build_basis

    !> Model m of the batch: SFH and dust parameters varied over their usual ranges.
    subroutine set_model(m, zred)
        integer, intent(in) :: m
        real(WP), intent(in) :: zred
        integer, parameter :: sfh(N_MOD) = [1, 4, 5, 0, 1, 4, 5, 1, 4, 1, 5, 4]
        real(WP), parameter :: tage(N_MOD) = [10.0_wp, 3.0_wp, 7.0_wp, 0.05_wp, 1.0_wp, 12.0_wp, &
                                              5.0_wp, 0.3_wp, 8.0_wp, 13.0_wp, 2.0_wp, 6.0_wp]
        real(WP), parameter :: tau(N_MOD) = [1.0_wp, 0.5_wp, 2.0_wp, 1.0_wp, 0.3_wp, 5.0_wp, &
                                             1.5_wp, 0.1_wp, 3.0_wp, 0.7_wp, 0.4_wp, 1.0_wp]
        real(WP), parameter :: dust1(N_MOD) = [0.0_wp, 0.6_wp, 1.0_wp, 0.3_wp, 0.0_wp, 0.8_wp, &
                                               0.2_wp, 1.5_wp, 0.0_wp, 0.4_wp, 0.9_wp, 0.1_wp]
        real(WP), parameter :: dust2(N_MOD) = [0.3_wp, 0.0_wp, 0.5_wp, 1.0_wp, 0.05_wp, 0.3_wp, &
                                               0.8_wp, 2.0_wp, 0.0_wp, 0.2_wp, 0.6_wp, 0.01_wp]
        real(WP), parameter :: dindex(N_MOD) = [-0.7_wp, -0.3_wp, -1.2_wp, 0.0_wp, -0.5_wp, -0.9_wp, &
                                                0.2_wp, -0.7_wp, -1.0_wp, -0.4_wp, -0.6_wp, -0.8_wp]
        real(WP), parameter :: umin(N_MOD) = [1.0_wp, 5.0_wp, 15.0_wp, 0.5_wp, 2.0_wp, 1.0_wp, &
                                              8.0_wp, 20.0_wp, 1.0_wp, 3.0_wp, 0.2_wp, 10.0_wp]
        real(WP), parameter :: qpah(N_MOD) = [3.5_wp, 0.5_wp, 4.5_wp, 2.0_wp, 1.0_wp, 3.5_wp, &
                                              0.1_wp, 3.0_wp, 2.5_wp, 4.0_wp, 1.5_wp, 3.5_wp]
        real(WP), parameter :: gam(N_MOD) = [0.01_wp, 0.1_wp, 0.0_wp, 0.3_wp, 0.01_wp, 0.05_wp, &
                                             0.5_wp, 0.01_wp, 0.2_wp, 0.0_wp, 0.02_wp, 1.0_wp]

        call set_int('sfh', sfh(m))
        call set_float('tage', tage(m))
        call set_float('tau', tau(m))
        call set_float('const', merge(0.2_wp, 0.0_wp, mod(m, 3) == 0))
        call set_float('fburst', merge(0.1_wp, 0.0_wp, mod(m, 4) == 1 .and. sfh(m) == 1))
        call set_float('tburst', 0.5_wp * tage(m))
        call set_float('sf_start', merge(0.5_wp, 0.0_wp, m == 6))
        call set_float('sf_trunc', merge(0.7_wp * tage(m), 0.0_wp, sfh(m) == 5))
        call set_float('sf_slope', merge(-1.0_wp, 0.0_wp, sfh(m) == 5))
        call set_float('dust1', dust1(m))
        call set_float('dust2', dust2(m))
        call set_float('dust_index', dindex(m))
        call set_float('dust1_index', merge(-1.0_wp, -0.7_wp, mod(m, 2) == 0))
        call set_float('dust_tesc', merge(7.5_wp, 7.0_wp, mod(m, 5) == 0))
        call set_float('frac_obrun', merge(0.2_wp, 0.0_wp, mod(m, 3) == 1))
        call set_float('frac_nodust', merge(0.1_wp, 0.0_wp, mod(m, 4) == 2))
        call set_float('mwr', merge(4.0_wp, 3.1_wp, mod(m, 2) == 1))
        call set_float('uvb', merge(0.5_wp, 1.0_wp, mod(m, 3) == 2))
        call set_float('duste_umin', umin(m))
        call set_float('duste_qpah', qpah(m))
        call set_float('duste_gamma', gam(m))
        call set_float('zred', zred)
    end subroutine set_model

    subroutine check_batch(dust_type, zred, igm, label)
        integer, intent(in) :: dust_type, igm
        real(WP), intent(in) :: zred
        character(len=*), intent(in) :: label
        type(PARAMS) :: psets(N_MOD)
        type(COMPSPOUT), allocatable :: results(:)
        real(WP), allocatable :: spec(:, :), mags(:, :), ref_spec(:, :), ref_mags(:, :)
        real(WP) :: mass(N_MOD), lbol(N_MOD), mdust(N_MOD), sfr(N_MOD)
        real(WP) :: ref_mass(N_MOD), ref_lbol(N_MOD), ref_mdust(N_MOD), ref_sfr(N_MOD)
        real(WP) :: d_spec, d_mag, d_scal, peak
        integer :: m, status, nspec, nbands, i
        logical :: bright(size(ctx%state%bands, 2))

        call print_group("Batched CSP: " // label)
        nspec = ctx%state%nspec
        nbands = ctx%state%nbands
        allocate(spec(nspec, N_MOD), mags(nbands, N_MOD), ref_spec(nspec, N_MOD), ref_mags(nbands, N_MOD))
        call set_int('dust_type', dust_type)
        call set_int('add_igm_absorption', igm)

        do m = 1, N_MOD
            call set_model(m, zred)
            call fsps_context_prepare_pset(ctx)
            call fsps_context_get_pset(ctx, psets(m))
            call compute_csp_scenario(ctx, psets(m), 1, results, status)
            if (status /= 0) then
                write(*, *) 'ERROR: compute_csp_scenario returned status ', status
                stop 1
            end if
            ref_spec(:, m) = results(1)%spec(1:nspec)
            ref_mags(:, m) = results(1)%mags(1:nbands)
            ref_mass(m) = results(1)%mass_csp
            ref_lbol(m) = results(1)%lbol_csp
            ref_mdust(m) = results(1)%mdust
            ref_sfr(m) = results(1)%sfr
        end do

        call compute_csp_batch(ctx, psets, mags, status, spec=spec, mass=mass, lbol=lbol, mdust=mdust, sfr=sfr)
        call assert_true(status == 0, "status 0", total_tests, total_failures)
        if (status /= 0) return

        d_spec = 0.0_wp
        d_mag = 0.0_wp
        do m = 1, N_MOD
            peak = maxval(ref_spec(:, m))
            do i = 1, nspec
                if (ref_spec(i, m) > 1.0e-6_wp * peak) &
                    d_spec = max(d_spec, abs(spec(i, m) / ref_spec(i, m) - 1.0_wp))
            end do
            ! Bands within 30 mag of the brightest (fainter ones carry numerically zero flux)
            bright = ref_mags(:, m) == ref_mags(:, m)
            where (bright) bright = ref_mags(:, m) < minval(ref_mags(:, m), mask=bright) + 30.0_wp
            d_mag = max(d_mag, maxval(abs(mags(:, m) - ref_mags(:, m)), mask=bright))
            if (any(bright .neqv. (mags(:, m) == mags(:, m) .and. bright))) d_mag = huge(1.0_wp)
        end do
        d_scal = max(maxval(abs(mass / ref_mass - 1.0_wp)), maxval(abs(lbol - ref_lbol)), &
                     maxval(abs(sfr - ref_sfr) / max(abs(ref_sfr), tiny(1.0_wp))))
        write(*, '(A,4ES10.2)') "# max difference (spectrum, mag, mass/lbol/sfr, mdust): ", d_spec, d_mag, d_scal, &
            maxval(abs(mdust / ref_mdust - 1.0_wp))
        call assert_true(d_spec <= SPEC_RTOL, "spectra", total_tests, total_failures)
        call assert_true(d_mag <= MAG_TOL, "magnitudes", total_tests, total_failures)
        call assert_true(d_scal <= 1.0e-12_wp, "mass, L_bol, SFR", total_tests, total_failures)
        call assert_true(maxval(abs(mdust / ref_mdust - 1.0_wp)) <= SPEC_RTOL, "dust mass", &
                         total_tests, total_failures)
    end subroutine check_batch

    subroutine set_int(key, value)
        character(len=*), intent(in) :: key
        integer, intent(in) :: value
        integer :: status

        call fsps_context_set_param_int(ctx, key, value, status)
        if (status /= 0) then
            write(*, *) 'ERROR: cannot set ', key
            stop 1
        end if
    end subroutine set_int

    subroutine set_float(key, value)
        character(len=*), intent(in) :: key
        real(WP), intent(in) :: value
        integer :: status

        call fsps_context_set_param_float(ctx, key, value, status)
        if (status /= 0) then
            write(*, *) 'ERROR: cannot set ', key
            stop 1
        end if
    end subroutine set_float

end program test_csp_batch
