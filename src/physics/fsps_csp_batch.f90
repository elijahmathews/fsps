module fsps_csp_batch
    !> @brief Batched CSP evaluation: many parameter sets against resident SSP bases.
    !>
    !> @details
    !> The fitting workload: N models that share the libraries and the resident SSP bases (one
    !> metallicity, or all of them for per-model logzsol), and differ in everything else that
    !> acts in the CSP: star formation history (parameters or a weight matrix), metallicity,
    !> redshift and IGM, dust, dust emission, nebular and AGN parameters. Per model the work is
    !> the same as `compute_csp_scenario` with `tage > 0`; here each stage runs once for the
    !> whole batch:
    !>
    !> 1. SFH weights per model in FP64 (on the device, one thread per model; on the host below
    !>    16 models), split into young and old ages by each model's `dust_tesc`, cast once to
    !>    kernel precision: W_young, W_old (ages x N).
    !> 2. Age sums: young = basis . W_young, old = basis . W_old (matrix products with K = number
    !>    of ages). With the SSP bases of all metallicities resident (fsps_context_update_ssp_basis
    !>    with nz bases), each model's basis is interpolated in log Z at its logzsol between the
    !>    two bracketing metallicities (as compute_interpolated_ssp: linear in the spectra, the
    !>    masses and log L_bol, not clamped); mass, L_bol and the nebular emission (Q, EUV) use the
    !>    interpolated SSPs. This mode uses the portable kernels (also in cuBLAS builds).
    !> 3. Nebular emission (add_neb_emission = 1), as apply_nebular_emission does to each young SSP
    !>    before the age sum, but summed over ages per model: the EUV rows of the age sums are
    !>    redone with the escaping fraction (frac_obrun) for the nebular ages; the nebular
    !>    continuum (Q(t) times the CLOUDY grid interpolated in gas_logz, gas_logu and age, in
    !>    FP64 log space) is summed with the SFH weights per (pixel, model); line luminosities per
    !>    (line, model), then their profiles added over each pixel's overlapping lines.
    !> 4. Attenuation, the birth-cloud and diffuse dust of `fsps_dust`, and the absorbed light.
    !> 5. Dust emission (energy balance): per-model reductions over pixels, then the
    !>    Draine & Li template, interpolated per model from the resident grid, is added.
    !> 6. IGM transmission (per model when the redshifts or igm_factor differ), then AGN dust
    !>    emission (fagn, agn_tau per model), attenuated by the model's diffuse dust, as
    !>    apply_agn_dust_emission.
    !> 7. Velocity smoothing (sigma_smooth, shared by the batch; smooth_velocity = 1 with the
    !>    default smoothspec_fast = 1), the sparse operator of apply_smoothing applied on the device.
    !> 8. Magnitudes: per (band, model) filter integrals over each band's support, compensated
    !>    at FP32 kernel precision (each model resampled to the observed frame at its redshift);
    !>    optionally spectral indices.
    !>
    !> Every stage runs on the device with OpenACC (-Daccel=gpu) or OpenMP target offload
    !> (-Daccel=omp_gpu), both directive sets on the same loops; the context must be on the
    !> device (moved there on the first call if it is not). The two matrix products (age sums,
    !> band fluxes) use cuBLAS when built with -Dcublas=true (NVIDIA only; a ceiling for the
    !> portable loops).
    !>
    !> Redshift and igm_factor may differ per model: each model's IGM transmission is computed on
    !> the device and its spectrum resampled to the (fixed) observed-frame grid for the shared
    !> filters, with its own distance term.
    !>
    !> Instead of SFH parameters, a precomputed SSP-weight matrix (ages x N, ssp_weights) can be
    !> given (non-parametric SFHs): the weights of each SSP age (mass formed in that bin, as
    !> compute_csp_scenario's csp_weights; spectra, mass and L_bol are per that mass), used as
    !> they are; tage is then not used, and SFR is not defined (NaN).
    !>
    !> Not supported (status 2): smoothing other than in velocity with smoothspec_fast = 1,
    !> redshifted colors, light-weighted ages, and models with tage <= 0 (all-ages output). Models
    !> must share the smoothing parameters (sigma_smooth, min/max_wave_smooth).
    !>
    !> Mass, L_bol and SFR follow `compute_csp_scenario` in FP64; spectra, magnitudes, indices,
    !> line luminosities and dust masses agree with it to the rounding of the kernel precision
    !> (summation order differs).

    use fsps_precision, only: WP, KP, kahan_add
    use fsps_constants, only: SAFE_FLOOR, PI, NEMLINE, NEBNAGE, NEBNZ, NEBNIP, H_PLANCK, L_SOL
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t, fsps_context_sync_device_settings
    use fsps_context, only: fsps_context_move_to_device
    use fsps_csp, only: sfh_weights_core
    use fsps_sfh, only: get_sfh_properties_at_age
    use fsps_dust, only: compute_attenuation_curve_point, one_minus_trans
    use fsps_interpolation, only: find_interval
    use fsps_cosmology, only: igm_tau
    use fsps_omp_check, only: omp_check_enabled, omp_check_present, omp_check_context
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fsps_gas, only: get_grid_indices_weights, prepare_line_gaussians
    use fsps_spectral_indices, only: compute_spectral_indices_columns
    use fsps_batch_photometry, only: compute_magnitudes_columns, prepare_batch_grid, smoothing_on_device, &
                                     prepare_smoothing_kernel, smooth_columns
#ifdef FSPS_CUBLAS
    use fsps_batch_photometry, only: blas, blas_init
    use cublas
#endif

    implicit none
    private

    public :: compute_csp_batch


    real(WP), parameter :: V_BAND_ANGSTROMS = 5500.0_wp   ! as in fsps_dust

    !> Nebular emission data of a batch (nebular_setup)
    type :: nebular_setup_t
        integer :: nt_neb = 0         !< ages with nebular emission (1..nt_neb)
        integer :: whlylim = 0        !< pixels below the Lyman limit
        logical :: cont = .false., lines_in_spec = .false., xrb = .false.
        integer,  allocatable :: ia(:)          !< per age: CLOUDY age interval and weight
        real(WP), allocatable :: wa(:)
        real(WP), allocatable :: dwa(:)         !< per age: wa(t) - wa(t-1) within an interval
        logical,  allocatable :: newr(:)        !< per age: the step dwa differs from the previous one
        real(WP), allocatable :: cy(:,:), co(:,:) !< (ages, N): SFH weight * Q, young and old (continuum)
        real(WP), allocatable :: q0(:,:)        !< per age and basis: ionizing photons (frac_obrun = 0)
        integer,  allocatable :: iz(:), iu(:)   !< per model: gas_logz and gas_logu intervals and weights
        real(WP), allocatable :: wz(:), wu(:)
        real(WP), allocatable :: fesc(:)        !< per model: escaping EUV fraction (clamped frac_obrun)
        real(WP), allocatable :: fabs(:)        !< per model: 1 - frac_obrun (absorbed fraction in Q)
        real(KP), allocatable :: lum(:,:)       !< line luminosities (lines, 2N): young, then old
        integer,  allocatable :: lptr(:), lid(:) !< per pixel (CSR): lines whose profile reaches it
        integer,  allocatable :: llo(:)         !< per line: pixel interval of its wavelength
    end type nebular_setup_t

contains

    !> @brief Evaluates a batch of single-age CSPs.
    !>
    !> @param[inout] ctx     Context with the SSP basis set (fsps_context_update_ssp_basis) and,
    !>                       in OpenACC builds, moved to the device.
    !> @param[in]    psets   Parameter sets, one per model (prepared: fsps_context_prepare_pset).
    !> @param[out]   mags    Magnitudes (nbands, N); bands not requested in psets(1)%mag_compute
    !>                       are NaN, as in compute_magnitudes.
    !> @param[out]   status  0 on success; 1 no setup or basis; 2 unsupported option; 3 bad shape.
    !> @param[out]   spec    Optional spectra (nspec, N), L_sun/Hz per solar mass formed.
    !> @param[out]   mass, lbol, mdust, sfr  Optional per-model scalars as in compspout.
    !> @param[out]   emlines Optional emission-line luminosities (NEMLINE, N), as compspout%emlines
    !>                       (zero without nebular emission).
    !> @param[out]   indx    Optional spectral indices (nindx, N), as compspout%indx.
    !> @param[in]    ssp_weights Optional SSP weights (ntfull, N) instead of the SFH parameters.
    subroutine compute_csp_batch(ctx, psets, mags, status, spec, mass, lbol, mdust, sfr, emlines, indx, ssp_weights)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: psets(:)
        real(WP), intent(out)               :: mags(:,:)
        integer, intent(out)                :: status
        real(WP), intent(out), optional     :: spec(:,:)
        real(WP), intent(out), optional     :: mass(:), lbol(:), mdust(:), sfr(:)
        real(WP), intent(out), optional     :: emlines(:,:), indx(:,:)
        real(WP), intent(in), optional      :: ssp_weights(:,:)

        integer :: n_mod, nt, nspec, nbands, m, k
        real(WP) :: zred, dust_age_log, w, mass_frac, sfr_norm, frac_linear
        logical  :: igm, any_emis, neb, agn, smooth, zmode, zvar, given_w, host_w
        !> Below this batch size the SFH weights are computed on the host: on the device each
        !> model's SFH integral is one serial thread (~0.3 ms on the A100, the same on the P100),
        !> on the host ~20 us per model. (Hardware note: the crossover depends on the host core
        !> and on the GPU's single-thread FP64 speed; measured on the A100 node.)
        integer, parameter :: HOST_SFH_MAX = 16
        real(WP), allocatable :: wh(:,:), t1h(:), t2h(:)
        real(WP), allocatable :: wext(:,:), zcol(:)
        integer  :: nzb
        integer,  allocatable :: zlo(:), zhi(:)
        real(WP), allocatable :: dz(:)

        ! Per-model parameters and results on the host
        real(KP), allocatable :: w2(:,:)   ! (ages, 2N): young weights, then old
        real(WP), allocatable :: weights(:,:,:), tmp1(:,:), tmp2(:,:), lum_linear(:), msum(:), lsum(:)
        integer,  allocatable :: it_tesc(:)
        integer,  allocatable :: iq(:), iu(:)
        real(WP), allocatable :: wq(:), wu(:), gam(:)
        logical,  allocatable :: emis(:)
        real(WP), allocatable :: s_abs(:), s_norm(:), s_esc(:), s_line(:)

        ! Nebular emission and AGN (per age, per model, per pixel line lists)
        type(nebular_setup_t) :: nb
        integer,  allocatable :: ita(:)
        real(WP), allocatable :: wta(:)
        logical,  allocatable :: agn_on(:)
        real(KP), allocatable :: eml(:,:)

        ! Work arrays (pixels x models) at kernel precision
        real(KP), allocatable :: yo(:,:), total(:,:), absorbed(:,:), tdiff(:,:)

        ! Device builds: the batch runs on the device, so a context not yet there is moved
#if (defined(_OPENACC) && !defined(ACC_MULTICORE)) || defined(OMP_TARGET)
        if (.not. ctx%on_device) call fsps_context_move_to_device(ctx)
#endif
        call fsps_context_sync_device_settings(ctx)
        status = 0
        n_mod = size(psets)
        nt = ctx%state%ntfull
        nspec = ctx%state%nspec
        nbands = ctx%state%nbands
        if (n_mod == 0) return

        ! ------------------------------------------------------------------
        ! 0. Checks
        ! ------------------------------------------------------------------
        if (ctx%state%check_sps_setup == 0 .or. .not. allocated(ctx%state%ssp_basis_spec)) then
            status = 1; return
        end if
        if (size(mags, 1) < nbands .or. size(mags, 2) < n_mod) then
            status = 3; return
        end if
        if (present(spec)) then
            if (size(spec, 1) < nspec .or. size(spec, 2) < n_mod) then
                status = 3; return
            end if
        end if
        if (present(emlines)) then
            if (size(emlines, 1) < NEMLINE .or. size(emlines, 2) < n_mod) then
                status = 3; return
            end if
        end if
        if (present(indx)) then
            if (size(indx, 1) < ctx%state%nindx .or. size(indx, 2) < n_mod) then
                status = 3; return
            end if
        end if
        if (ctx%redshift_colors_val /= 0 .or. ctx%compute_light_ages_val /= 0) then
            status = 2; return
        end if
        given_w = present(ssp_weights)
        if (given_w) then
            if (size(ssp_weights, 1) < nt .or. size(ssp_weights, 2) < n_mod) then
                status = 3; return
            end if
        end if
        zred = psets(1)%zred
        zvar = .false.
        do m = 1, n_mod
            if ((psets(m)%tage <= 0.0_wp .and. .not. given_w) .or. psets(m)%sigma_smooth /= psets(1)%sigma_smooth .or. &
                psets(m)%min_wave_smooth /= psets(1)%min_wave_smooth .or. &
                psets(m)%max_wave_smooth /= psets(1)%max_wave_smooth .or. &
                psets(m)%uvb < 0.0_wp .or. psets(m)%wgp1 < 1 .or. psets(m)%wgp2 < 1) then
                status = 2; return
            end if
            if (psets(m)%zred /= zred .or. psets(m)%igm_factor /= psets(1)%igm_factor) zvar = .true.
        end do
        allocate(zcol(n_mod))
        zcol = psets(1:n_mod)%zred
        host_w = .false.
#if (defined(_OPENACC) && !defined(ACC_MULTICORE)) || defined(OMP_TARGET)
        host_w = (.not. given_w) .and. n_mod < HOST_SFH_MAX
#endif
        if (given_w) then
            allocate(wext(nt, n_mod))
            wext = ssp_weights(1:nt, 1:n_mod)
        else if (host_w) then
            ! Few models: the SFH weights on the host (as compute_csp_scenario does for one output),
            ! then used as a given weight matrix
            allocate(wext(nt, n_mod), wh(nt, 1), t1h(nt), t2h(nt))
            do m = 1, n_mod
                call sfh_weights_core(ctx, psets(m), psets(m)%tage, 1, wh, t1h, t2h)
                wext(:, m) = wh(:, 1)
            end do
        else
            allocate(wext(1, 1))
        end if

        ! Smoothing as apply_smoothing (which does nothing for sigma <= SAFE_FLOOR)
        smooth = psets(1)%sigma_smooth > SAFE_FLOOR
        if (smooth .and. .not. smoothing_on_device(ctx)) then
            status = 2; return
        end if

        ! Per-model metallicity: a basis per metallicity of the context
        nzb = size(ctx%state%ssp_basis_spec, 3)
        zmode = nzb > 1
        if (zmode .and. nzb /= ctx%state%nz) then
            status = 3; return
        end if
        allocate(zlo(n_mod), zhi(n_mod), dz(n_mod))
        zlo = 1; zhi = 1; dz = 0.0_wp
        if (zmode) then
            do m = 1, n_mod
                call metallicity_interval(ctx, psets(m)%logzsol, zlo(m), dz(m))
                zhi(m) = zlo(m) + 1
            end do
        end if

        ! ------------------------------------------------------------------
        ! 1. Host: per-model indices (birth-cloud age, dust emission grid, nebular grid, AGN
        !    template) and SFR
        ! ------------------------------------------------------------------
        allocate(w2(nt, 2 * n_mod), lum_linear(nt), it_tesc(n_mod), msum(n_mod), lsum(n_mod))
        ! SFH scratch: sfh_weights_core uses ntfull entries (sized from the ages, not from the
        ! context's CSP workspace, which a context used only through the batch may not have)
        allocate(weights(nt, 1, n_mod), tmp1(nt, n_mod), tmp2(nt, n_mod))
        allocate(iq(n_mod), iu(n_mod), wq(n_mod), wu(n_mod), gam(n_mod), emis(n_mod))
        allocate(ita(n_mod), wta(n_mod), agn_on(n_mod))
        lum_linear = 10.0_wp**ctx%state%ssp_basis_lbol(1:nt, 1)
        do m = 1, n_mod
            if (psets(m)%dust_tesc > SAFE_FLOOR) then
                dust_age_log = psets(m)%dust_tesc
            else
                dust_age_log = 7.0_wp
            end if
            it_tesc(m) = max(1, min(find_interval(ctx%state%time_full, dust_age_log), nt))
            if (present(sfr) .and. given_w) then
                sfr(m) = ieee_value(1.0_wp, ieee_quiet_nan)
            else if (present(sfr)) then
                call get_sfh_properties_at_age(ctx, psets(m), psets(m)%tage, mass_frac, sfr_norm, frac_linear)
                if (mass_frac > SAFE_FLOOR) sfr_norm = sfr_norm / mass_frac
                sfr(m) = sfr_norm
            end if

            emis(m) = ctx%add_dust_emission_val == 1 .and. &
                      (psets(m)%dust1 > SAFE_FLOOR .or. psets(m)%dust2 > SAFE_FLOOR)
            call dl07_indices(ctx, psets(m), iq(m), wq(m), iu(m), wu(m), gam(m))
            call agn_indices(ctx, psets(m), agn_on(m), ita(m), wta(m))
        end do
        any_emis = any(emis)
        agn = any(agn_on)
        neb = ctx%add_neb_emission_val == 1
        call nebular_setup(ctx, psets, neb, nzb, nb)
        if (neb) then
            allocate(eml(NEMLINE, n_mod))
        else
            allocate(eml(1, 1))
        end if

        ! Grid-only data for this redshift (cached in ctx%state%batch, on the device too)
        ! (with per-model redshifts: the z = 0 data, for the redshift-independent parts)
        if (zvar) then
            call prepare_batch_grid(ctx, 0.0_wp, psets(1)%igm_factor)
            igm = ctx%add_igm_absorption_val == 1 .and. any(zcol > SAFE_FLOOR)
        else
            call prepare_batch_grid(ctx, zred, psets(1)%igm_factor)
            igm = ctx%state%batch%igm
        end if
        if (smooth) call prepare_smoothing_kernel(ctx, psets(1)%sigma_smooth, psets(1)%min_wave_smooth, &
                                                  psets(1)%max_wave_smooth)

        ! ------------------------------------------------------------------
        ! 2-7. Device: age sums, nebular emission, attenuation, dust emission, IGM, AGN, smoothing
        ! ------------------------------------------------------------------
        allocate(yo(nspec, 2 * n_mod), total(nspec, n_mod), absorbed(nspec, n_mod))
        ! Diffuse-dust transmission kept from the attenuation for the AGN emission
        if (agn) then
            allocate(tdiff(nspec, n_mod))
        else
            allocate(tdiff(1, 1))
        end if
        allocate(s_abs(n_mod), s_norm(n_mod), s_esc(n_mod), s_line(n_mod))

        !$acc data present(ctx) copyin(psets, iq, iu, wq, wu, gam, emis, it_tesc, lum_linear, ita, wta, agn_on) &
        !$acc      copyin(zlo, zhi, dz, wext, zcol) &
        !$acc      copyin(nb%ia, nb%wa, nb%dwa, nb%newr, nb%iz, nb%wz, nb%iu, nb%wu, nb%fesc, nb%fabs, nb%lptr) &
        !$acc      copyin(nb%lid, nb%llo) create(w2, weights, tmp1, tmp2, yo, total, absorbed, nb%q0, nb%lum, s_line) &
        !$acc      create(nb%cy, nb%co, tdiff) &
        !$acc      copyout(s_abs, s_norm, s_esc, msum, lsum, eml)
#ifdef OMP_TARGET
        !$omp target data map(to: psets, iq, iu, wq, wu, gam, emis, it_tesc, lum_linear, ita, wta, agn_on) &
        !$omp map(to: zlo, zhi, dz, wext, zcol) &
        !$omp map(alloc: w2, weights, tmp1, tmp2, yo, total, absorbed, s_line, tdiff) &
        !$omp map(from: s_abs, s_norm, s_esc, msum, lsum, eml)
#endif
#ifdef OMP_TARGET
        !$omp target data map(to: nb%ia, nb%wa, nb%dwa, nb%newr, nb%iz, nb%wz, nb%iu, nb%wu, nb%fesc, nb%fabs) &
        !$omp map(to: nb%lptr, nb%lid, nb%llo) map(alloc: nb%q0, nb%lum, nb%cy, nb%co)
#endif

#ifdef OMP_TARGET
        if (omp_check_enabled()) call check_batch_maps()
#endif

        ! SFH weights (FP64), mass and L_bol per model, one thread per model; the weights are split
        ! into young and old ages at the model's dust_tesc and cast once to kernel precision
        !$acc parallel loop gang vector present(ctx, psets, weights, tmp1, tmp2, w2, it_tesc, lum_linear, msum, lsum, &
        !$acc zlo, zhi, dz, wext) private(k, w)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(k, w)
#endif
        do m = 1, n_mod
            if (given_w .or. host_w) then
                !$acc loop seq
                do k = 1, nt
                    weights(k, 1, m) = wext(k, m)
                end do
            else
                call sfh_weights_core(ctx, psets(m), psets(m)%tage, 1, weights(:, :, m), tmp1(:, m), tmp2(:, m))
            end if
            msum(m) = 0.0_wp
            lsum(m) = 0.0_wp
            !$acc loop seq
            do k = 1, nt
                w = weights(k, 1, m)
                if (w > SAFE_FLOOR) then
                    if (zmode) then
                        ! The interpolated SSP's mass and L_bol (compute_interpolated_ssp)
                        msum(m) = msum(m) + w * ((1.0_wp - dz(m)) * ctx%state%ssp_basis_mass(k, zlo(m)) + &
                                                 dz(m) * ctx%state%ssp_basis_mass(k, zhi(m)))
                        lsum(m) = lsum(m) + w * 10.0_wp**((1.0_wp - dz(m)) * ctx%state%ssp_basis_lbol(k, zlo(m)) + &
                                                          dz(m) * ctx%state%ssp_basis_lbol(k, zhi(m)))
                    else
                        msum(m) = msum(m) + w * ctx%state%ssp_basis_mass(k, 1)
                        lsum(m) = lsum(m) + w * lum_linear(k)
                    end if
                else
                    w = 0.0_wp
                end if
                if (k <= it_tesc(m)) then
                    w2(k, m) = real(w, KP); w2(k, n_mod + m) = 0.0_KP
                else
                    w2(k, m) = 0.0_KP; w2(k, n_mod + m) = real(w, KP)
                end if
            end do
        end do

        if (zmode) then
            call age_sums_z(nspec, nt, n_mod, ctx%state%ssp_basis_spec, w2, it_tesc, zlo, zhi, dz, yo)
        else
            call age_sums(nspec, nt, n_mod, ctx%state%ssp_basis_spec(:, :, 1), w2, it_tesc, yo)
        end if
        if (neb) call add_nebular_batch(ctx, n_mod, nspec, nt, nzb, ctx%state%ssp_basis_spec, w2, it_tesc, &
                                        zlo, zhi, dz, nb, yo, eml)
        if (neb) call attenuate_lines(ctx, psets, n_mod, nb%lum, nb%llo, eml, s_line)
        call attenuate_batch(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, agn_on, yo, total, absorbed, tdiff)
        if (any_emis) then
            call dust_emission_reductions(n_mod, nspec, emis, ctx%state%batch%wt, yo, absorbed, s_abs, s_norm, s_esc)
            ! Lines not in the spectrum: their absorbed part (intrinsic - attenuated) is added, as
            ! in apply_dust_attenuation_and_emission
            if (neb .and. ctx%nebemlineinspec_val == 0) then
                !$acc parallel loop present(s_abs, s_line)
#ifdef OMP_TARGET
                !$omp target teams distribute parallel do
#endif
                do m = 1, n_mod
                    s_abs(m) = s_abs(m) + s_line(m)
                end do
            end if
            call add_dust_emission(n_mod, nspec, emis, yo, s_abs, s_esc, total)
        end if

        if (igm .and. zvar) then
            call apply_igm_z(ctx, psets, n_mod, nspec, absorbed, total)   ! absorbed: free by now
        else if (igm) then
            call apply_igm(n_mod, nspec, ctx%state%batch%igm_t, total)
        end if
        if (agn) call add_agn_batch(ctx, psets, n_mod, nspec, agn_on, ita, wta, lsum, tdiff, total)
        if (smooth) call smooth_columns(ctx, total, n_mod)

        ! ------------------------------------------------------------------
        ! 8. Magnitudes and indices
        ! ------------------------------------------------------------------
        if (zvar .and. allocated(psets(1)%mag_compute)) then
            call compute_magnitudes_columns(ctx, zred, psets(1)%igm_factor, total, n_mod, mags, psets(1)%mag_compute, &
                                            zred_col=zcol)
        else if (zvar) then
            call compute_magnitudes_columns(ctx, zred, psets(1)%igm_factor, total, n_mod, mags, zred_col=zcol)
        else if (allocated(psets(1)%mag_compute)) then
            call compute_magnitudes_columns(ctx, zred, psets(1)%igm_factor, total, n_mod, mags, psets(1)%mag_compute)
        else
            call compute_magnitudes_columns(ctx, zred, psets(1)%igm_factor, total, n_mod, mags)
        end if
        if (present(indx)) call compute_spectral_indices_columns(ctx, ctx%state%spec_lambda, total, n_mod, indx)

        if (present(spec)) then
            !$acc update self(total)
#ifdef OMP_TARGET
            !$omp target update from(total)
#endif
            do m = 1, n_mod
                spec(1:nspec, m) = real(total(:, m), WP)
            end do
        end if
#ifdef OMP_TARGET
        !$omp end target data
#endif
#ifdef OMP_TARGET
        !$omp end target data
        ! LLVM flang (22) leaves the data of psets (derived type with allocatable components)
        ! mapped after the region; a later host allocation at the same address would then be
        ! taken as present and not copied. Delete the stale entry (no-op where it is gone).
        !$omp target exit data map(delete: psets)
#endif
        !$acc end data

        if (present(emlines)) then
            do m = 1, n_mod
                if (neb) then
                    emlines(1:NEMLINE, m) = real(eml(:, m), WP)
                else
                    emlines(1:NEMLINE, m) = 0.0_wp
                end if
            end do
        end if
        do m = 1, n_mod
            if (present(mass)) mass(m) = msum(m)
            if (present(lbol)) then
                if (lsum(m) > SAFE_FLOOR) then
                    lbol(m) = log10(lsum(m))
                else
                    lbol(m) = 0.0_wp
                end if
            end if
            if (present(mdust)) then
                if (emis(m) .and. s_norm(m) > SAFE_FLOOR) then
                    mdust(m) = 3.21e-3_wp / (4.0_wp * PI) * (s_abs(m) / s_norm(m))
                else
                    mdust(m) = SAFE_FLOOR
                end if
            end if
        end do

#ifdef OMP_TARGET
    contains

        !> FSPS_OMP_CHECK_PRE: every array the batch kernels use is on the device, and the
        !> context components read in device routines are attached there (fsps_omp_check)
        subroutine check_batch_maps()
            character(len=*), parameter :: W = 'compute_csp_batch'

            call omp_check_context(ctx, W)
            call omp_check_present(loc(psets), 'psets', W, size(psets))
            call omp_check_present(loc(iq), 'iq', W, size(iq)); call omp_check_present(loc(iu), 'iu', W, size(iu))
            call omp_check_present(loc(wq), 'wq', W, size(wq)); call omp_check_present(loc(wu), 'wu', W, size(wu))
            call omp_check_present(loc(gam), 'gam', W, size(gam)); call omp_check_present(loc(emis), 'emis', W, size(emis))
            call omp_check_present(loc(it_tesc), 'it_tesc', W, size(it_tesc)); call omp_check_present(loc(lum_linear), 'lum_linear', W, size(lum_linear))
            call omp_check_present(loc(ita), 'ita', W, size(ita)); call omp_check_present(loc(wta), 'wta', W, size(wta))
            call omp_check_present(loc(agn_on), 'agn_on', W, size(agn_on))
            call omp_check_present(loc(zlo), 'zlo', W, size(zlo)); call omp_check_present(loc(zhi), 'zhi', W, size(zhi))
            call omp_check_present(loc(dz), 'dz', W, size(dz)); call omp_check_present(loc(wext), 'wext', W, size(wext))
            call omp_check_present(loc(zcol), 'zcol', W, size(zcol)); call omp_check_present(loc(w2), 'w2', W, size(w2))
            call omp_check_present(loc(weights), 'weights', W, size(weights)); call omp_check_present(loc(tmp1), 'tmp1', W, size(tmp1))
            call omp_check_present(loc(tmp2), 'tmp2', W, size(tmp2)); call omp_check_present(loc(yo), 'yo', W, size(yo))
            call omp_check_present(loc(total), 'total', W, size(total)); call omp_check_present(loc(absorbed), 'absorbed', W, size(absorbed))
            call omp_check_present(loc(s_line), 's_line', W, size(s_line)); call omp_check_present(loc(tdiff), 'tdiff', W, size(tdiff))
            call omp_check_present(loc(s_abs), 's_abs', W, size(s_abs)); call omp_check_present(loc(s_norm), 's_norm', W, size(s_norm))
            call omp_check_present(loc(s_esc), 's_esc', W, size(s_esc)); call omp_check_present(loc(msum), 'msum', W, size(msum))
            call omp_check_present(loc(lsum), 'lsum', W, size(lsum)); call omp_check_present(loc(eml), 'eml', W, size(eml))
            if (allocated(nb%ia)) then
                call omp_check_present(loc(nb%ia), 'nb%ia', W, size(nb%ia)); call omp_check_present(loc(nb%wa), 'nb%wa', W, size(nb%wa))
                call omp_check_present(loc(nb%dwa), 'nb%dwa', W, size(nb%dwa)); call omp_check_present(loc(nb%newr), 'nb%newr', W, size(nb%newr))
                call omp_check_present(loc(nb%iz), 'nb%iz', W, size(nb%iz)); call omp_check_present(loc(nb%wz), 'nb%wz', W, size(nb%wz))
                call omp_check_present(loc(nb%iu), 'nb%iu', W, size(nb%iu)); call omp_check_present(loc(nb%wu), 'nb%wu', W, size(nb%wu))
                call omp_check_present(loc(nb%fesc), 'nb%fesc', W, size(nb%fesc)); call omp_check_present(loc(nb%fabs), 'nb%fabs', W, size(nb%fabs))
                call omp_check_present(loc(nb%lptr), 'nb%lptr', W, size(nb%lptr)); call omp_check_present(loc(nb%lid), 'nb%lid', W, size(nb%lid))
                call omp_check_present(loc(nb%llo), 'nb%llo', W, size(nb%llo)); call omp_check_present(loc(nb%q0), 'nb%q0', W, size(nb%q0))
                call omp_check_present(loc(nb%lum), 'nb%lum', W, size(nb%lum)); call omp_check_present(loc(nb%cy), 'nb%cy', W, size(nb%cy))
                call omp_check_present(loc(nb%co), 'nb%co', W, size(nb%co))
            end if
            call omp_check_present(loc(ctx%state%ssp_basis_spec), 'ctx%state%ssp_basis_spec', W, size(ctx%state%ssp_basis_spec))
            if (allocated(ctx%state%batch%wt)) then
                call omp_check_present(loc(ctx%state%batch%wt), 'ctx%state%batch%wt', W, size(ctx%state%batch%wt))
                call omp_check_present(loc(ctx%state%batch%igm_t), 'ctx%state%batch%igm_t', W, size(ctx%state%batch%igm_t))
                call omp_check_present(loc(ctx%state%batch%i_lo), 'ctx%state%batch%i_lo', W, size(ctx%state%batch%i_lo))
                call omp_check_present(loc(ctx%state%batch%t_obs), 'ctx%state%batch%t_obs', W, size(ctx%state%batch%t_obs))
                call omp_check_present(loc(ctx%state%batch%wl), 'ctx%state%batch%wl', W, size(ctx%state%batch%wl))
                call omp_check_present(loc(ctx%state%batch%wtl), 'ctx%state%batch%wtl', W, size(ctx%state%batch%wtl))
                call omp_check_present(loc(ctx%state%batch%bands_t), 'ctx%state%batch%bands_t', W, size(ctx%state%batch%bands_t))
                call omp_check_present(loc(ctx%state%batch%j_lo), 'ctx%state%batch%j_lo', W, size(ctx%state%batch%j_lo))
                call omp_check_present(loc(ctx%state%batch%j_hi), 'ctx%state%batch%j_hi', W, size(ctx%state%batch%j_hi))
            end if
            if (smooth) then
                call omp_check_present(loc(ctx%state%batch%sm_lo), 'ctx%state%batch%sm_lo', W, size(ctx%state%batch%sm_lo))
                call omp_check_present(loc(ctx%state%batch%sm_hi), 'ctx%state%batch%sm_hi', W, size(ctx%state%batch%sm_hi))
                call omp_check_present(loc(ctx%state%batch%sm_off), 'ctx%state%batch%sm_off', W, size(ctx%state%batch%sm_off))
                call omp_check_present(loc(ctx%state%batch%sm_coef), 'ctx%state%batch%sm_coef', W, size(ctx%state%batch%sm_coef))
            end if
        end subroutine check_batch_maps
#endif
    end subroutine compute_csp_batch

    !> Host-side nebular data of a batch: young ages and their grid interval, per-model grid
    !> indices and escape fractions, and per-pixel lists of the lines whose profile reaches it.
    subroutine nebular_setup(ctx, psets, neb, nzb, nb)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: psets(:)
        logical, intent(in)                 :: neb
        integer, intent(in)                 :: nzb
        type(nebular_setup_t), intent(out)  :: nb

        real(WP), parameter :: LINE_PROFILE_CUT = 1.0e-14_wp   ! as add_lines_to_spectrum (CPU)
        integer  :: n_mod, nspec, t, m, i, j, jc, n, p
        integer, allocatable :: j_lo(:), j_hi(:), cnt(:)
        real(WP) :: cut

        n_mod = size(psets)
        nspec = ctx%state%nspec
        nb%nt_neb = 0
        if (neb) nb%nt_neb = max(0, min(find_interval(ctx%state%time_full, ctx%state%nebem_age(NEBNAGE)), &
                                        ctx%state%ntfull))
        nb%cont = neb .and. ctx%add_neb_continuum_val == 1
        nb%lines_in_spec = neb .and. ctx%nebemlineinspec_val == 1
        nb%xrb = (ctx%state%isoc_type == 'bpss') .and. (ctx%add_xrb_emission_val == 1)
        nb%whlylim = ctx%state%whlylim
        n = max(1, nb%nt_neb)
        allocate(nb%ia(n), nb%wa(n), nb%q0(n, nzb))
        allocate(nb%iz(n_mod), nb%wz(n_mod), nb%iu(n_mod), nb%wu(n_mod), nb%fesc(n_mod), nb%fabs(n_mod))
        nb%ia = 1; nb%wa = 0.0_wp
        do t = 1, nb%nt_neb
            call get_grid_indices_weights(ctx%state%nebem_age, ctx%state%time_full(t), NEBNAGE, nb%ia(t), nb%wa(t))
        end do
        ! Steps of wa within each interval (constant on the uniform log-age grid, zero where wa is
        ! clamped), and where a new step starts
        allocate(nb%dwa(n), nb%newr(n))
        nb%dwa = 0.0_wp
        nb%newr = .false.
        do t = 2, nb%nt_neb
            if (nb%ia(t) /= nb%ia(t - 1)) cycle
            nb%dwa(t) = nb%wa(t) - nb%wa(t - 1)
            nb%newr(t) = nb%ia(t - 1) /= nb%ia(max(1, t - 2)) .or. t == 2 .or. nb%dwa(t) /= nb%dwa(t - 1)
        end do
        allocate(nb%cy(n, n_mod), nb%co(n, n_mod))
        do m = 1, n_mod
            nb%iz(m) = 1; nb%wz(m) = 0.0_wp; nb%iu(m) = 1; nb%wu(m) = 0.0_wp
            if (neb) then
                call get_grid_indices_weights(ctx%state%nebem_logz, psets(m)%gas_logz, NEBNZ, nb%iz(m), nb%wz(m))
                call get_grid_indices_weights(ctx%state%nebem_logu, psets(m)%gas_logu, NEBNIP, nb%iu(m), nb%wu(m))
            end if
            nb%fesc(m) = max(0.0_wp, min(psets(m)%frac_obrun, 1.0_wp))
            nb%fabs(m) = 1.0_wp - psets(m)%frac_obrun
        end do
        if (neb) then
            allocate(nb%lum(NEMLINE, 2 * n_mod))
        else
            allocate(nb%lum(1, 1))
        end if
        ! Pixel interval of each line wavelength (the search of apply_dust_attenuation_and_emission)
        allocate(nb%llo(NEMLINE))
        do i = 1, NEMLINE
            jc = 1
            j = nspec + 1
            do while (j - jc > 1)
                p = (jc + j) / 2
                if (ctx%state%nebem_line_pos(i) >= ctx%state%spec_lambda(p)) then
                    jc = p
                else
                    j = p
                end if
            end do
            nb%llo(i) = max(1, min(jc, nspec - 1))
        end do

        ! Line windows (as the CPU add_lines_to_spectrum: from the pixel nearest the line centre
        ! until the profile falls below LINE_PROFILE_CUT of its peak), then per pixel the lines
        ! whose window contains it, in line order (CSR)
        allocate(nb%lptr(nspec + 1))
        nb%lptr = 1
        if (.not. nb%lines_in_spec) then
            allocate(nb%lid(1)); nb%lid = 1
            return
        end if
        call prepare_line_gaussians(ctx, psets(1))
        allocate(j_lo(NEMLINE), j_hi(NEMLINE), cnt(nspec))
        cnt = 0
        do i = 1, NEMLINE
            j_lo(i) = 1; j_hi(i) = 0
            jc = max(1, min(find_interval(ctx%state%spec_lambda(1:nspec), ctx%state%nebem_line_pos(i)), nspec))
            j = min(jc + 1, nspec)
            cut = LINE_PROFILE_CUT * max(real(ctx%state%gaussnebarr(jc, i), WP), real(ctx%state%gaussnebarr(j, i), WP))
            if (cut <= 0.0_wp) cycle
            j_hi(i) = j
            j_lo(i) = jc
            do while (j_lo(i) > 1)
                if (ctx%state%gaussnebarr(j_lo(i) - 1, i) <= cut) exit
                j_lo(i) = j_lo(i) - 1
            end do
            do while (j_hi(i) < nspec)
                if (ctx%state%gaussnebarr(j_hi(i) + 1, i) <= cut) exit
                j_hi(i) = j_hi(i) + 1
            end do
            cnt(j_lo(i):j_hi(i)) = cnt(j_lo(i):j_hi(i)) + 1
        end do
        do j = 1, nspec
            nb%lptr(j + 1) = nb%lptr(j) + cnt(j)
        end do
        allocate(nb%lid(max(1, nb%lptr(nspec + 1) - 1)))
        cnt = 0
        do i = 1, NEMLINE
            do j = j_lo(i), j_hi(i)
                p = nb%lptr(j) + cnt(j)
                nb%lid(p) = i
                cnt(j) = cnt(j) + 1
            end do
        end do
    end subroutine nebular_setup

    !> AGN template interval for pset (apply_agn_dust_emission); on = AGN emission is added.
    subroutine agn_indices(ctx, pset, on, ita, wta)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: pset
        logical, intent(out)             :: on
        integer, intent(out)             :: ita
        real(WP), intent(out)            :: wta
        integer :: n

        on = ctx%add_agn_dust_val == 1 .and. pset%fagn > SAFE_FLOOR .and. pset%fagn > tiny(0.0_wp)
        ita = 1; wta = 0.0_wp
        if (.not. on) return
        n = size(ctx%state%agndust_tau)
        ita = max(1, min(find_interval(ctx%state%agndust_tau, pset%agn_tau), n - 1))
        wta = (pset%agn_tau - ctx%state%agndust_tau(ita)) / (ctx%state%agndust_tau(ita + 1) - ctx%state%agndust_tau(ita))
        wta = max(0.0_wp, min(wta, 1.0_wp))
    end subroutine agn_indices

    !> Interpolation indices of the Draine & Li grid (interpolate_draine_li_dust_model).
    subroutine dl07_indices(ctx, pset, iq, wq, iu, wu, gam)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: pset
        integer, intent(out)             :: iq, iu
        real(WP), intent(out)            :: wq, wu, gam

        iq = 1; iu = 1; wq = 0.0_wp; wu = 0.0_wp; gam = 0.0_wp
        if (ctx%add_dust_emission_val /= 1 .or. .not. associated(ctx%state%qpaharr)) return
        iq = max(1, min(find_interval(ctx%state%qpaharr, pset%duste_qpah), ctx%state%nqpah_dustem - 1))
        wq = (pset%duste_qpah - ctx%state%qpaharr(iq)) / (ctx%state%qpaharr(iq + 1) - ctx%state%qpaharr(iq))
        wq = max(0.0_wp, min(wq, 1.0_wp))
        iu = max(1, min(find_interval(ctx%state%uminarr, pset%duste_umin), ctx%state%numin_dustem - 1))
        wu = (pset%duste_umin - ctx%state%uminarr(iu)) / (ctx%state%uminarr(iu + 1) - ctx%state%uminarr(iu))
        wu = max(0.0_wp, min(wu, 1.0_wp))
        gam = max(0.0_wp, min(pset%duste_gamma, 1.0_wp))
    end subroutine dl07_indices

    !> Draine & Li emission template at one pixel (as interpolate_draine_li_dust_model).
    function dl07_shape(ctx, i, iq, wq, iu, wu, gam) result(shape)
        !$acc routine seq
#ifdef OMP_TARGET
        !$omp declare target
#endif
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)              :: i, iq, iu
        real(WP), intent(in)             :: wq, wu, gam
        real(KP) :: shape
        real(WP) :: lo, hi

        lo = (1-wq)*(1-wu) * ctx%state%dustem2_dustem(i, iq,   2*iu - 1) + &
             wq    *(1-wu) * ctx%state%dustem2_dustem(i, iq+1, 2*iu - 1) + &
             wq    *wu     * ctx%state%dustem2_dustem(i, iq+1, 2*(iu+1) - 1) + &
             (1-wq)*wu     * ctx%state%dustem2_dustem(i, iq,   2*(iu+1) - 1)
        hi = (1-wq)*(1-wu) * ctx%state%dustem2_dustem(i, iq,   2*iu) + &
             wq    *(1-wu) * ctx%state%dustem2_dustem(i, iq+1, 2*iu) + &
             wq    *wu     * ctx%state%dustem2_dustem(i, iq+1, 2*(iu+1)) + &
             (1-wq)*wu     * ctx%state%dustem2_dustem(i, iq,   2*(iu+1))
        shape = max(real((1.0_wp - gam) * lo + gam * hi, KP), tiny(1.0_KP))
    end function dl07_shape

    !> Diffuse-dust transmission of model pset at pixel i (and its optical depth).
    subroutine diffuse_transmission(ctx, pset, i, curve, tau_d, trans)
        !$acc routine seq
#ifdef OMP_TARGET
        !$omp declare target
#endif
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: pset
        integer, intent(in)              :: i
        real(WP), intent(out)            :: curve, tau_d, trans

        curve = compute_attenuation_curve_point(ctx%state%spec_lambda(i), i, ctx%dust_type_val, pset, ctx)
        if (ctx%dust_type_val == 3) then
            tau_d = curve
        else
            tau_d = pset%dust2 * curve
        end if
        trans = exp(-tau_d)
    end subroutine diffuse_transmission

    !> Age sums (young and old) and attenuation; absorbed light density for dust emission.

    !> Nebular emission of every model, added to its young and old age sums (yo), as
    !> apply_nebular_emission modifies each SSP of age t <= nt_neb before the age sum of
    !> compute_csp_scenario: (1) the EUV escapes with fraction fesc, (2) the continuum
    !> 10**log(t; Z, U) * Q(t), (3) the lines, whose luminosities (also returned in eml) are summed
    !> over ages first and then spread over their profiles.
    subroutine add_nebular_batch(ctx, n_mod, nspec, nt, nzb, basis, w2, it_tesc, zlo, zhi, dz, nb, yo, eml)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)              :: n_mod, nspec, nt, nzb, it_tesc(:), zlo(:), zhi(:)
        real(KP), intent(in), contiguous :: basis(:,:,:), w2(:,:)
        real(WP), intent(in)             :: dz(:)
        type(nebular_setup_t), intent(inout) :: nb
        real(KP), intent(inout)          :: yo(:,:)
        real(KP), intent(out)            :: eml(:,:)

        ! The kernels take the components as plain arrays: device code never refers to the
        ! derived type itself, which is mapped only through its components (mapping it whole and
        ! by component, or implicitly through a dummy, conflicts in LLVM's offload runtime)
        call nebular_kernels(ctx, n_mod, nspec, nt, nzb, basis, w2, it_tesc, zlo, zhi, dz, &
                             nb%nt_neb, nb%whlylim, nb%cont, nb%lines_in_spec, nb%xrb, nb%ia, nb%wa, nb%dwa, &
                             nb%newr, nb%q0, nb%iz, nb%wz, nb%iu, nb%wu, nb%fesc, nb%fabs, nb%cy, nb%co, nb%lum, &
                             nb%lptr, nb%lid, yo, eml)
    end subroutine add_nebular_batch

    !> The kernels of add_nebular_batch, on the components of its nebular_setup_t.
    subroutine nebular_kernels(ctx, n_mod, nspec, nt, nzb, basis, w2, it_tesc, zlo, zhi, dz, &
                               nb_nt_neb, nb_whlylim, nb_cont, nb_lines_in_spec, nb_xrb, nb_ia, nb_wa, nb_dwa, &
                               nb_newr, nb_q0, nb_iz, nb_wz, nb_iu, nb_wu, nb_fesc, nb_fabs, nb_cy, nb_co, nb_lum, &
                               nb_lptr, nb_lid, yo, eml)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)              :: n_mod, nspec, nt, nzb, it_tesc(:), zlo(:), zhi(:)
        real(KP), intent(in), contiguous :: basis(:,:,:), w2(:,:)
        real(WP), intent(in)             :: dz(:)
        integer, intent(in)              :: nb_nt_neb, nb_whlylim
        logical, intent(in)              :: nb_cont, nb_lines_in_spec, nb_xrb
        integer, intent(in), contiguous  :: nb_ia(:), nb_iz(:), nb_iu(:), nb_lptr(:), nb_lid(:)
        real(WP), intent(in), contiguous :: nb_wa(:), nb_dwa(:), nb_wz(:), nb_wu(:), nb_fesc(:), nb_fabs(:)
        logical, intent(in), contiguous  :: nb_newr(:)
        real(WP), intent(inout), contiguous :: nb_q0(:,:), nb_cy(:,:), nb_co(:,:)
        real(KP), intent(inout), contiguous :: nb_lum(:,:)
        real(KP), intent(inout)          :: yo(:,:)
        real(KP), intent(out)            :: eml(:,:)

        integer  :: m, i, j, k, t, nt_neb, whlylim, a, a_prev, p, iz, iu
        real(KP) :: s, c, wt, bk, y, o, wk, v, dzk
        real(WP) :: q, lg, ra, rb, dl, f, r, sy, so, c00, c01, c10, c11
        real(WP), parameter :: LN10 = log(10.0_wp)
        logical  :: xrb

        nt_neb = nb_nt_neb
        whlylim = nb_whlylim
        xrb = nb_xrb
        if (nt_neb < 1) then
            !$acc parallel loop collapse(2) present(eml)
#ifdef OMP_TARGET
            !$omp target teams distribute parallel do collapse(2)
#endif
            do m = 1, n_mod
                do i = 1, NEMLINE
                    eml(i, m) = 0.0_KP
                end do
            end do
            return
        end if

        ! 1. Ionizing photons of each young basis spectrum (the photon-number integral of
        !    process_ionizing_radiation, without the absorbed fraction)
        !    (linear in the spectrum, so interpolated SSPs take the interpolated Q)
        !$acc parallel loop gang collapse(2) present(ctx, basis, nb_q0) private(s, c)
#ifdef OMP_TARGET
        !$omp target teams distribute collapse(2) private(s, c)
#endif
        do iz = 1, nzb
        do t = 1, nt_neb
            s = 0.0_KP
            c = 0.0_KP
            if (whlylim >= 2) then
                !$acc loop vector reduction(+:s) private(wt, v)
#ifdef OMP_TARGET
                !$omp parallel do reduction(+:s) private(wt, v)
#endif
                do i = 1, whlylim - 1
                    wt = real(0.5_wp * abs(ctx%state%spec_nu(i+1) - ctx%state%spec_nu(i)), KP)
                    v = wt * (basis(i, t, iz) / real(ctx%state%spec_nu(i), KP) + &
                              basis(i+1, t, iz) / real(ctx%state%spec_nu(i+1), KP))
#if defined(_OPENACC) || defined(OMP_TARGET)
                    s = s + v
#else
                    if (KP == WP) then
                        s = s + v
                    else
                        call kahan_add(s, c, v)
                    end if
#endif
                end do
                nb_q0(t, iz) = real(s, WP) / H_PLANCK * L_SOL
            else
                nb_q0(t, iz) = 0.0_wp
            end if
        end do
        end do

        ! 2. EUV rows: the age sums again, the nebular ages' spectra times the escaping fraction
        if (whlylim > 0) then
            !$acc parallel loop gang vector collapse(2) present(basis, w2, it_tesc, nb_fesc, yo, zlo, zhi, dz) &
            !$acc private(y, o, bk, k, dzk)
#ifdef OMP_TARGET
            !$omp target teams distribute parallel do collapse(2) private(y, o, bk, k, dzk)
#endif
            do m = 1, n_mod
                do i = 1, whlylim
                    y = 0.0_KP
                    o = 0.0_KP
                    dzk = real(dz(m), KP)
                    !$acc loop seq
                    do k = 1, nt
                        if (nzb > 1) then
                            bk = (1.0_KP - dzk) * basis(i, k, zlo(m)) + dzk * basis(i, k, zhi(m))
                        else
                            bk = basis(i, k, 1)
                        end if
                        if (k <= nt_neb) bk = bk * real(nb_fesc(m), KP)
                        if (k <= it_tesc(m)) then
                            y = y + w2(k, m) * bk
                        else
                            o = o + w2(k, n_mod + m) * bk
                        end if
                    end do
                    yo(i, m) = y
                    yo(i, n_mod + m) = o
                end do
            end do
        end if

        ! 3. Continuum: per (pixel, model), sum over the nebular ages t of C(t) 10**lg(t), with
        !    C = SFH weight * Q (per model, young or old; FP64) and lg the log grid interpolated
        !    in Z, U (FP64) and age. Within one CLOUDY age interval [ra, rb], lg = ra + wa (rb - ra)
        !    and on the uniform log-age grid wa advances by a constant step, so 10**lg is a
        !    geometric progression: per interval two FP64 exponentials (the first age's value and
        !    the ratio), per age one multiply and one add, all in FP64 (each interval's sum is cast
        !    once). The FP64 exponentials are the cost of this kernel. Hardware note: measured on
        !    the A100 (FP64 at half the FP32 rate); on GPUs with slow FP64 (1/32-1/64 of FP32, most
        !    consumer cards) these exponentials weigh much more, and the next step there would be
        !    the progression's ratio in FP32 with the argument reduced in FP64.
        if (nb_cont) then
            !$acc parallel loop gang vector collapse(2) present(w2, it_tesc, nb_q0, nb_fabs, nb_cy, nb_co, zlo, zhi, dz)
#ifdef OMP_TARGET
            !$omp target teams distribute parallel do collapse(2)
#endif
            do m = 1, n_mod
                do t = 1, nt_neb
                    q = abs(((1.0_wp - dz(m)) * nb_q0(t, zlo(m)) + dz(m) * nb_q0(t, zhi(m))) * nb_fabs(m))
                    if (t <= it_tesc(m)) then
                        wk = w2(t, m)
                    else
                        wk = w2(t, n_mod + m)
                    end if
                    nb_cy(t, m) = 0.0_wp
                    nb_co(t, m) = 0.0_wp
                    if (wk /= 0.0_KP .and. q > SAFE_FLOOR) then
                        if (t <= it_tesc(m)) then
                            nb_cy(t, m) = real(wk, WP) * q
                        else
                            nb_co(t, m) = real(wk, WP) * q
                        end if
                    end if
                end do
            end do

            !$acc parallel loop gang vector collapse(2) present(ctx, nb_ia, nb_wa, nb_dwa, nb_newr, nb_iz, nb_wz, &
            !$acc nb_iu, nb_wu, nb_cy, nb_co, yo) &
            !$acc private(y, o, sy, so, ra, rb, dl, f, r, a, a_prev, t, iz, iu, c00, c01, c10, c11)
#ifdef OMP_TARGET
            !$omp target teams distribute parallel do collapse(2) private(y, o, sy, so, ra, rb, dl, f, r, &
            !$omp& a, a_prev, t, iz, iu, c00, c01, c10, c11)
#endif
            do m = 1, n_mod
                do j = 1, nspec
                    y = 0.0_KP
                    o = 0.0_KP
                    sy = 0.0_wp
                    so = 0.0_wp
                    a_prev = -1
                    ra = 0.0_wp; rb = 0.0_wp
                    dl = 0.0_wp; f = 0.0_wp; r = 1.0_wp
                    iz = nb_iz(m)
                    iu = nb_iu(m)
                    ! As interpolate_zu_slice_point, written out: with nvfortran 25.1 (A100) a call
                    ! with the 4-d assumed-shape grid per interval and thread cost more than the
                    ! arithmetic (descriptor set-up in local memory); 53 -> 41 ms at N = 10000
                    c00 = (1.0_wp - nb_wz(m)) * (1.0_wp - nb_wu(m))
                    c01 = (1.0_wp - nb_wz(m)) * (         nb_wu(m))
                    c10 = (         nb_wz(m)) * (1.0_wp - nb_wu(m))
                    c11 = (         nb_wz(m)) * (         nb_wu(m))
                    !$acc loop seq
                    do t = 1, nt_neb
                        a = nb_ia(t)
                        if (a /= a_prev) then
                            ! Close the previous interval
                            y = y + real(sy, KP)
                            o = o + real(so, KP)
                            sy = 0.0_wp
                            so = 0.0_wp
                            ! The interval only advances: its lower end is the previous upper end
                            if (a == a_prev + 1) then
                                ra = rb
                            else
                                if (xrb) then
                                    ra = c00 * ctx%state%xnebem_cont(j, iz, a, iu) + &
                                         c01 * ctx%state%xnebem_cont(j, iz, a, iu + 1) + &
                                         c10 * ctx%state%xnebem_cont(j, iz + 1, a, iu) + &
                                         c11 * ctx%state%xnebem_cont(j, iz + 1, a, iu + 1)
                                else
                                    ra = c00 * ctx%state%nebem_cont(j, iz, a, iu) + &
                                         c01 * ctx%state%nebem_cont(j, iz, a, iu + 1) + &
                                         c10 * ctx%state%nebem_cont(j, iz + 1, a, iu) + &
                                         c11 * ctx%state%nebem_cont(j, iz + 1, a, iu + 1)
                                end if
                            end if
                            if (xrb) then
                                rb = c00 * ctx%state%xnebem_cont(j, iz, a + 1, iu) + &
                                     c01 * ctx%state%xnebem_cont(j, iz, a + 1, iu + 1) + &
                                     c10 * ctx%state%xnebem_cont(j, iz + 1, a + 1, iu) + &
                                     c11 * ctx%state%xnebem_cont(j, iz + 1, a + 1, iu + 1)
                            else
                                rb = c00 * ctx%state%nebem_cont(j, iz, a + 1, iu) + &
                                     c01 * ctx%state%nebem_cont(j, iz, a + 1, iu + 1) + &
                                     c10 * ctx%state%nebem_cont(j, iz + 1, a + 1, iu) + &
                                     c11 * ctx%state%nebem_cont(j, iz + 1, a + 1, iu + 1)
                            end if
                            a_prev = a
                            dl = (rb - ra) * LN10
                            f = exp(ra * LN10 + nb_wa(t) * dl)
                        else
                            if (nb_newr(t)) r = exp(nb_dwa(t) * dl)
                            f = f * r
                        end if
                        sy = sy + nb_cy(t, m) * f
                        so = so + nb_co(t, m) * f
                    end do
                    y = y + real(sy, KP)
                    o = o + real(so, KP)
                    yo(j, m) = yo(j, m) + y
                    yo(j, n_mod + m) = yo(j, n_mod + m) + o
                end do
            end do
        end if

        ! 4. Line luminosities per (line, model), young and old: sum over ages of the weight times
        !    10**(log grid interpolated) * Q
        !$acc parallel loop gang vector collapse(2) present(ctx, w2, it_tesc, nb_ia, nb_wa, nb_q0, nb_iz, nb_wz, &
        !$acc nb_iu, nb_wu, nb_fabs, nb_lum, eml, zlo, zhi, dz) private(y, o, wk, q, lg, ra, rb, a, t)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2) private(y, o, wk, q, lg, ra, rb, a, t)
#endif
        do m = 1, n_mod
            do i = 1, NEMLINE
                y = 0.0_KP
                o = 0.0_KP
                !$acc loop seq
                do t = 1, nt_neb
                    if (t <= it_tesc(m)) then
                        wk = w2(t, m)
                    else
                        wk = w2(t, n_mod + m)
                    end if
                    if (wk == 0.0_KP) cycle
                    q = abs(((1.0_wp - dz(m)) * nb_q0(t, zlo(m)) + dz(m) * nb_q0(t, zhi(m))) * nb_fabs(m))
                    if (q <= SAFE_FLOOR) cycle
                    a = nb_ia(t)
                    if (xrb) then
                        ra = zu_point(ctx%state%xnebem_line, i, a, nb_iz(m), nb_iu(m), nb_wz(m), nb_wu(m))
                        rb = zu_point(ctx%state%xnebem_line, i, a + 1, nb_iz(m), nb_iu(m), nb_wz(m), nb_wu(m))
                    else
                        ra = zu_point(ctx%state%nebem_line, i, a, nb_iz(m), nb_iu(m), nb_wz(m), nb_wu(m))
                        rb = zu_point(ctx%state%nebem_line, i, a + 1, nb_iz(m), nb_iu(m), nb_wz(m), nb_wu(m))
                    end if
                    lg = (1.0_wp - nb_wa(t)) * ra + nb_wa(t) * rb
                    if (t <= it_tesc(m)) then
                        y = y + wk * real((10.0_wp**lg) * q, KP)
                    else
                        o = o + wk * real((10.0_wp**lg) * q, KP)
                    end if
                end do
                nb_lum(i, m) = y
                nb_lum(i, n_mod + m) = o
                eml(i, m) = y + o
            end do
        end do

        ! 5. Line profiles: per pixel, the lines whose profile reaches it
        if (nb_lines_in_spec) then
            !$acc parallel loop gang vector collapse(2) present(ctx, nb_lptr, nb_lid, nb_lum, yo) private(s, p, i)
#ifdef OMP_TARGET
            !$omp target teams distribute parallel do collapse(2) private(s, p, i)
#endif
            do m = 1, 2 * n_mod
                do j = 1, nspec
                    s = 0.0_KP
                    !$acc loop seq
                    do p = nb_lptr(j), nb_lptr(j + 1) - 1
                        i = nb_lid(p)
                        s = s + ctx%state%gaussnebarr(j, i) * nb_lum(i, m)
                    end do
                    yo(j, m) = yo(j, m) + s
                end do
            end do
        end if
    end subroutine nebular_kernels

    !> Line luminosities after dust (apply_dust_attenuation_and_emission): young lines through
    !> the birth cloud, all through the diffuse dust at the line wavelength (the transmission
    !> interpolated between pixels, at KP as there). eml = attenuated luminosities; s_line =
    !> intrinsic - attenuated, per model.
    subroutine attenuate_lines(ctx, psets, n_mod, nb_lum, nb_llo, eml, s_line)
        type(fsps_context_t), intent(in)  :: ctx
        type(params), intent(in)          :: psets(:)
        integer, intent(in)               :: n_mod
        real(KP), intent(in), contiguous  :: nb_lum(:,:)   ! line luminosities (nebular_setup_t)
        integer, intent(in), contiguous   :: nb_llo(:)     ! line pixel intervals
        real(KP), intent(out)             :: eml(:,:)
        real(WP), intent(out)             :: s_line(:)

        integer  :: m, i, lo
        real(WP) :: si, sa, t_lo, t_hi, curve, tau_d, trans, td, birth, ay, ao, fob, fnd, pos
        real(KP) :: ly, lo_

        !$acc parallel loop gang present(ctx, psets, nb_lum, nb_llo, eml, s_line) private(si, sa, fob, fnd)
#ifdef OMP_TARGET
        !$omp target teams distribute private(si, sa, fob, fnd)
#endif
        do m = 1, n_mod
            si = 0.0_wp
            sa = 0.0_wp
            fob = psets(m)%frac_obrun
            fnd = psets(m)%frac_nodust
            !$acc loop vector reduction(+:si, sa) private(lo, t_lo, t_hi, curve, tau_d, trans, td, birth, ay, ao, &
            !$acc ly, lo_, pos)
#ifdef OMP_TARGET
            !$omp parallel do reduction(+:si, sa) private(lo, t_lo, t_hi, curve, tau_d, trans, td, birth, &
            !$omp& ay, ao, ly, lo_, pos)
#endif
            do i = 1, NEMLINE
                pos = ctx%state%nebem_line_pos(i)
                if (1.0_wp - fnd <= SAFE_FLOOR) then
                    td = 1.0_wp
                else
                    lo = nb_llo(i)
                    call diffuse_transmission(ctx, psets(m), lo, curve, tau_d, trans)
                    t_lo = real(real(trans, KP), WP)
                    call diffuse_transmission(ctx, psets(m), lo + 1, curve, tau_d, trans)
                    t_hi = real(real(trans, KP), WP)
                    td = t_lo + (t_hi - t_lo) / (ctx%state%spec_lambda(lo + 1) - ctx%state%spec_lambda(lo)) * &
                         (pos - ctx%state%spec_lambda(lo))
                end if
                if (td /= td) td = 1.0_wp
                if (psets(m)%dust1 <= SAFE_FLOOR .or. 1.0_wp - fob <= SAFE_FLOOR) then
                    birth = 1.0_wp
                else
                    birth = exp(-psets(m)%dust1 * (pos / V_BAND_ANGSTROMS)**psets(m)%dust1_index)
                end if
                ly = nb_lum(i, m)
                lo_ = nb_lum(i, n_mod + m)
                ay = ly * (birth * (1.0_wp - fob) + fob)
                ao = lo_
                if (1.0_wp - fnd > SAFE_FLOOR) then
                    ay = ay * (td * (1.0_wp - fnd) + fnd)
                    ao = ao * (td * (1.0_wp - fnd) + fnd)
                end if
                si = si + (ly + lo_)
                sa = sa + (ay + ao)
                eml(i, m) = real(ay, KP) + real(ao, KP)
            end do
            s_line(m) = si - sa
        end do
    end subroutine attenuate_lines

    !> CLOUDY grid value at (i, age a) bilinear in (Z, U), as interpolate_zu_slice_point.
    pure function zu_point(grid, i, a, iz, iu, wz, wu) result(val)
        !$acc routine seq
#ifdef OMP_TARGET
        !$omp declare target
#endif
        real(WP), intent(in) :: grid(:,:,:,:)
        integer, intent(in)  :: i, a, iz, iu
        real(WP), intent(in) :: wz, wu
        real(WP) :: val

        val = (1.0_wp - wz) * (1.0_wp - wu) * grid(i, iz,     a, iu    ) + &
              (1.0_wp - wz) * (         wu) * grid(i, iz,     a, iu + 1) + &
              (         wz) * (1.0_wp - wu) * grid(i, iz + 1, a, iu    ) + &
              (         wz) * (         wu) * grid(i, iz + 1, a, iu + 1)
    end function zu_point

    !> AGN dust emission (apply_agn_dust_emission): fagn L_bol times the agn_tau template,
    !> attenuated by the model's diffuse dust; after the IGM transmission, as there.
    subroutine add_agn_batch(ctx, psets, n_mod, nspec, agn_on, ita, wta, lsum, tdiff, total)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, ita(:)
        logical, intent(in)              :: agn_on(:)
        real(WP), intent(in)             :: wta(:), lsum(:)
        real(KP), intent(in)             :: tdiff(:,:)   ! the diffuse transmission (attenuate_batch)
        real(KP), intent(inout)          :: total(:,:)

        integer  :: m, i
        real(WP) :: lagn, tmpl, lb

        !$acc parallel loop gang vector collapse(2) present(ctx, psets, agn_on, ita, wta, lsum, tdiff, total) &
        !$acc private(lagn, tmpl, lb)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2) private(lagn, tmpl, lb)
#endif
        do m = 1, n_mod
            do i = 1, nspec
                if (agn_on(m)) then
                    ! L_bol through its logarithm, as compute_csp_scenario stores it
                    if (lsum(m) > SAFE_FLOOR) then
                        lb = log10(lsum(m))
                    else
                        lb = 0.0_wp
                    end if
                    lagn = (10.0_wp**lb) * psets(m)%fagn
                    tmpl = (1.0_wp - wta(m)) * ctx%state%agndust_spec(i, ita(m)) + &
                           wta(m) * ctx%state%agndust_spec(i, ita(m) + 1)
                    ! exp(-dust2 curve) (exp(-curve) for dust_type 3): the transmission of attenuate_batch,
                    ! kept at KP as fsps_dust keeps it
                    tmpl = tmpl * real(tdiff(i, m), WP)
                    total(i, m) = total(i, m) + lagn * tmpl
                end if
            end do
        end do
    end subroutine add_agn_batch

    !> Interval and weight of log Z = logzsol in the context's metallicities, as
    !> compute_interpolated_ssp (the weight is not clamped: outside the grid it extrapolates).
    subroutine metallicity_interval(ctx, logzsol, zlo, dz)
        type(fsps_context_t), intent(in) :: ctx
        real(WP), intent(in)             :: logzsol
        integer, intent(out)             :: zlo
        real(WP), intent(out)            :: dz
        integer :: nz

        nz = ctx%state%nz
        zlo = max(1, min(find_interval(log10(ctx%state%zlegend(1:nz) / ctx%state%zsol), logzsol), nz - 1))
        dz = (logzsol - log10(ctx%state%zlegend(zlo) / ctx%state%zsol)) / &
             (log10(ctx%state%zlegend(zlo + 1) / ctx%state%zsol) - log10(ctx%state%zlegend(zlo) / ctx%state%zsol))
    end subroutine metallicity_interval

    !> Age sums with per-model metallicity: each model's basis is (1 - dz) basis(:, :, zlo) +
    !> dz basis(:, :, zhi), formed at kernel precision term by term (portable kernel; also used in
    !> cuBLAS builds). One pass over the ages (young or old per model), 4 models per thread.
    subroutine age_sums_z(nspec, nt, n_mod, basis, w2, it_tesc, zlo, zhi, dz, yo)
        integer, intent(in)              :: nspec, nt, n_mod, it_tesc(:), zlo(:), zhi(:)
        real(KP), intent(in), contiguous :: basis(:,:,:), w2(:,:)
        real(WP), intent(in)             :: dz(:)
        real(KP), intent(out)            :: yo(:,:)
        integer, parameter :: MT = 4
        integer  :: mq, i, k, l, mm(MT)
        real(KP) :: bk, y(MT), o(MT), d(MT)

        !$acc parallel loop gang vector collapse(2) present(basis, w2, it_tesc, zlo, zhi, dz, yo) &
        !$acc private(k, l, bk, mm, y, o, d)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2) private(k, l, bk, mm, y, o, d)
#endif
        do mq = 1, (n_mod + MT - 1) / MT
            do i = 1, nspec
                !$acc loop seq
                do l = 1, MT
                    mm(l) = min((mq - 1) * MT + l, n_mod)
                    d(l) = real(dz(mm(l)), KP)
                    y(l) = 0.0_KP
                    o(l) = 0.0_KP
                end do
                !$acc loop seq
                do k = 1, nt
                    !$acc loop seq
                    do l = 1, MT
                        bk = (1.0_KP - d(l)) * basis(i, k, zlo(mm(l))) + d(l) * basis(i, k, zhi(mm(l)))
                        if (k <= it_tesc(mm(l))) then
                            y(l) = y(l) + w2(k, mm(l)) * bk
                        else
                            o(l) = o(l) + w2(k, n_mod + mm(l)) * bk
                        end if
                    end do
                end do
                !$acc loop seq
                do l = 1, MT
                    yo(i, mm(l)) = y(l)
                    yo(i, n_mod + mm(l)) = o(l)
                end do
            end do
        end do
    end subroutine age_sums_z

    !> IGM transmission per model (redshift and igm_factor per model), as compute_igm_transmission:
    !> the optical depth at KP, capped below the wavelength of its maximum, then exp(-tau factor).
    !> tau is work space (pixels x models).
    subroutine apply_igm_z(ctx, psets, n_mod, nspec, tau, total)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec
        real(KP), intent(inout)          :: tau(:,:), total(:,:)
        integer  :: m, j, jmax
        real(WP) :: tmax

        !$acc parallel loop gang vector collapse(2) present(ctx, psets, tau)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2)
#endif
        do m = 1, n_mod
            do j = 1, nspec
                if (psets(m)%zred > SAFE_FLOOR) tau(j, m) = real(igm_tau(ctx%state%spec_lambda(j), 1.0_wp + psets(m)%zred), KP)
            end do
        end do
        ! The short-wavelength cap, one thread per model
        !$acc parallel loop gang vector present(psets, tau) private(tmax, jmax, j)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do private(tmax, jmax, j)
#endif
        do m = 1, n_mod
            if (psets(m)%zred > SAFE_FLOOR) then
                tmax = 0.0_wp
                jmax = 1
                !$acc loop seq
                do j = 1, nspec
                    if (tau(j, m) > tmax) then
                        tmax = tau(j, m)
                        jmax = j
                    end if
                end do
                if (jmax > 1) then
                    !$acc loop seq
                    do j = 1, jmax
                        tau(j, m) = real(tmax, KP)
                    end do
                end if
            end if
        end do
        !$acc parallel loop gang vector collapse(2) present(psets, tau, total)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2)
#endif
        do m = 1, n_mod
            do j = 1, nspec
                if (psets(m)%zred > SAFE_FLOOR) &
                    total(j, m) = total(j, m) * real(exp(-real(tau(j, m), WP) * psets(m)%igm_factor), KP)
            end do
        end do
    end subroutine apply_igm_z

    !> Common IGM transmission.
    subroutine apply_igm(n_mod, nspec, igm_t, total)
        integer, intent(in)     :: n_mod, nspec
        real(KP), intent(in)    :: igm_t(:)
        real(KP), intent(inout) :: total(:,:)
        integer :: m, i

        !$acc parallel loop collapse(2) present(total, igm_t)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2)
#endif
        do m = 1, n_mod
            do i = 1, nspec
                total(i, m) = total(i, m) * igm_t(i)
            end do
        end do
    end subroutine apply_igm

    !> Age sums of all models: yo = basis . w2 (pixels x ages . ages x 2N), plain sums over the ages
    !> at kernel precision as in compute_csp_scenario.
    subroutine age_sums(nspec, nt, n_mod, basis, w2, it_tesc, yo)
        integer, intent(in)              :: nspec, nt, n_mod, it_tesc(:)
        real(KP), intent(in), contiguous :: basis(:,:), w2(:,:)
        real(KP), intent(out)            :: yo(:,:)
#ifdef FSPS_CUBLAS
        integer :: istat

        call blas_init()
        !$acc host_data use_device(basis, w2, yo)
#ifdef KERNEL_SINGLE_PRECISION
        istat = cublasSgemm_v2(blas, CUBLAS_OP_N, CUBLAS_OP_N, nspec, 2 * n_mod, nt, 1.0_KP, basis, nspec, &
                               w2, nt, 0.0_KP, yo, nspec)
#else
        istat = cublasDgemm_v2(blas, CUBLAS_OP_N, CUBLAS_OP_N, nspec, 2 * n_mod, nt, 1.0_KP, basis, nspec, &
                               w2, nt, 0.0_KP, yo, nspec)
#endif
        !$acc end host_data
        if (istat /= 0) error stop 'compute_csp_batch: cuBLAS GEMM failed'
#else
        ! Register tile: each thread does PT pixels (strided by npt, so each load is coalesced
        ! across the threads) x MT models, PT + 2 MT loads per PT * MT multiply-adds.
        ! Hardware note: MT = 8, PT = 4 were measured on the A100 (nvfortran 25.1, FP32: 64
        ! accumulators per thread, no spills); retune on other GPUs, where registers per thread
        ! and the load / FMA balance differ. Each age is young or old for a given model: one
        ! running sum per (pixel, model), stored as the young sum at the model's last young age
        ! (it_tesc) and restarted for the old ages. Same sums, in the same order, as two separate
        ! products. OpenACC and OpenMP target (-Daccel=omp_gpu) run the same loop.
        integer, parameter :: MT = 8, PT = 4
        integer  :: mq, ip, npt, i, k, l, p, mm(MT), it(MT), ii(PT)
        real(KP) :: bk(PT), wk(MT), acc(PT, MT), yy(PT, MT)

        npt = (nspec + PT - 1) / PT
        !$acc parallel loop gang vector collapse(2) present(basis, w2, it_tesc, yo) &
        !$acc private(mm, it, ii, bk, wk, acc, yy, k, l, p, i)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2) private(mm, it, ii, bk, wk, acc, yy, k, l, p, i)
#endif
        do mq = 1, (n_mod + MT - 1) / MT
            do ip = 1, npt
                !$acc loop seq
                do l = 1, MT
                    ! Models of this thread (clamped at the end; duplicates are written twice, same value)
                    mm(l) = min((mq - 1) * MT + l, n_mod)
                    it(l) = it_tesc(mm(l))
                end do
                !$acc loop seq
                do p = 1, PT
                    ii(p) = min(ip + (p - 1) * npt, nspec)
                end do
                acc = 0.0_KP
                yy = 0.0_KP
                !$acc loop seq
                do k = 1, nt
                    !$acc loop seq
                    do p = 1, PT
                        bk(p) = basis(ii(p), k)
                    end do
                    !$acc loop seq
                    do l = 1, MT
                        wk(l) = w2(k, mm(l)) + w2(k, n_mod + mm(l))   ! one of the two is zero
                    end do
                    !$acc loop seq
                    do l = 1, MT
                        !$acc loop seq
                        do p = 1, PT
                            acc(p, l) = acc(p, l) + wk(l) * bk(p)
                        end do
                        if (k == it(l)) then
                            !$acc loop seq
                            do p = 1, PT
                                yy(p, l) = acc(p, l)
                                acc(p, l) = 0.0_KP
                            end do
                        end if
                    end do
                end do
                !$acc loop seq
                do l = 1, MT
                    !$acc loop seq
                    do p = 1, PT
                        yo(ii(p), mm(l)) = yy(p, l)
                        yo(ii(p), n_mod + mm(l)) = acc(p, l)
                    end do
                end do
            end do
        end do
#endif
    end subroutine age_sums

    !> Attenuation and absorbed light per (pixel, model). For models with dust emission the age
    !> sums are then replaced by what the emission stages need: yo(:, m) = diffuse transmission
    !> (stored at KP, as fsps_dust stores it), yo(:, N + m) = Draine & Li template.
    subroutine attenuate_batch(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, agn_on, yo, total, absorbed, tdiff)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, iq(:), iu(:)
        real(WP), intent(in)             :: wq(:), wu(:), gam(:)
        logical, intent(in)              :: emis(:), agn_on(:)
        real(KP), intent(inout)          :: yo(:,:)   ! in: young (1:N), old (N+1:2N)
        real(KP), intent(out)            :: total(:,:), absorbed(:,:)
        real(KP), intent(inout)          :: tdiff(:,:)   ! diffuse transmission of AGN models

        integer  :: m, i
        real(KP) :: young, old
        real(WP) :: curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, frac_nodust
        real(WP) :: one_m_tb, one_m_to, one_m_td, a_young, spec_sum

        !$acc parallel loop gang vector collapse(2) present(ctx, psets, iq, iu, wq, wu, gam, yo, emis, total, absorbed) &
        !$acc present(agn_on, tdiff) &
        !$acc private(young, old, curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, &
        !$acc         frac_nodust, one_m_tb, one_m_to, one_m_td, a_young, spec_sum)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2) private(young, old, curve, tau_b, tau_o, &
        !$omp& tau_d, trans_birth, trans_old, trans_d, frac_obrun, frac_nodust, one_m_tb, one_m_to, &
        !$omp& one_m_td, a_young, spec_sum)
#endif
        do m = 1, n_mod
            do i = 1, nspec
                young = yo(i, m)
                old = yo(i, n_mod + m)

                frac_obrun = psets(m)%frac_obrun
                frac_nodust = psets(m)%frac_nodust
                call diffuse_transmission(ctx, psets(m), i, curve, tau_d, trans_d)
                if (agn_on(m)) tdiff(i, m) = real(trans_d, KP)
                if (psets(m)%dust1 <= SAFE_FLOOR .or. 1.0_wp - frac_obrun <= SAFE_FLOOR) then
                    tau_b = 0.0_wp
                    trans_birth = 1.0_wp
                else
                    tau_b = psets(m)%dust1 * (ctx%state%spec_lambda(i) / V_BAND_ANGSTROMS)**psets(m)%dust1_index
                    trans_birth = exp(-tau_b)
                end if
                if (abs(psets(m)%dust3) <= SAFE_FLOOR) then
                    tau_o = 0.0_wp
                    trans_old = 1.0_wp
                else
                    tau_o = psets(m)%dust3 * curve
                    trans_old = exp(-tau_o)
                end if

                spec_sum = (young * trans_birth * (1.0_wp - frac_obrun) + young * frac_obrun) + old * trans_old
                if (1.0_wp - frac_nodust <= SAFE_FLOOR) then
                    total(i, m) = real(spec_sum, KP)
                else
                    total(i, m) = real(spec_sum * (trans_d * (1.0_wp - frac_nodust) + frac_nodust), KP)
                end if

                if (emis(m)) then
                    one_m_tb = one_minus_trans(tau_b, trans_birth)
                    one_m_to = one_minus_trans(tau_o, trans_old)
                    if (1.0_wp - frac_nodust <= SAFE_FLOOR) then
                        one_m_td = 0.0_wp
                    else
                        one_m_td = (1.0_wp - frac_nodust) * one_minus_trans(tau_d, trans_d)
                    end if
                    a_young = trans_birth * (1.0_wp - frac_obrun) + frac_obrun
                    absorbed(i, m) = young * real((1.0_wp - frac_obrun) * one_m_tb + a_young * one_m_td, KP) + &
                                     old * real(one_m_to + trans_old * one_m_td, KP)
                    yo(i, m) = real(trans_d, KP)
                    yo(i, n_mod + m) = dl07_shape(ctx, i, iq(m), wq(m), iu(m), wu(m), gam(m))
                else
                    absorbed(i, m) = 0.0_KP
                end if
            end do
        end do
    end subroutine attenuate_batch

    !> Per model: absorbed luminosity, template normalization and escaped template fraction.
    subroutine dust_emission_reductions(n_mod, nspec, emis, wt, yo, absorbed, s_abs, s_norm, s_esc)
        integer, intent(in)              :: n_mod, nspec
        logical, intent(in)              :: emis(:)
        real(KP), intent(in)             :: wt(:), yo(:,:), absorbed(:,:)
        real(WP), intent(out)            :: s_abs(:), s_norm(:), s_esc(:)

        integer  :: m, i
        real(KP) :: sa, sn, se, ca, cn, ce

        !$acc parallel loop gang present(emis, wt, yo, absorbed, s_abs, s_norm, s_esc) private(sa, sn, se, ca, cn, ce)
#ifdef OMP_TARGET
        !$omp target teams distribute private(sa, sn, se, ca, cn, ce)
#endif
        do m = 1, n_mod
            sa = 0.0_KP; sn = 0.0_KP; se = 0.0_KP
            ca = 0.0_KP; cn = 0.0_KP; ce = 0.0_KP
            if (emis(m)) then
                !$acc loop vector reduction(+:sa, sn, se)
#ifdef OMP_TARGET
                !$omp parallel do reduction(+:sa, sn, se)
#endif
                do i = 1, nspec - 1
#if defined(_OPENACC) || defined(OMP_TARGET)
                    sa = sa + wt(i) * (absorbed(i, m) + absorbed(i + 1, m))
                    sn = sn + wt(i) * (yo(i, n_mod + m) + yo(i + 1, n_mod + m))
                    se = se + wt(i) * (yo(i, n_mod + m) * yo(i, m) + yo(i + 1, n_mod + m) * yo(i + 1, m))
#else
                    call kadd(sa, ca, wt(i) * (absorbed(i, m) + absorbed(i + 1, m)))
                    call kadd(sn, cn, wt(i) * (yo(i, n_mod + m) + yo(i + 1, n_mod + m)))
                    call kadd(se, ce, wt(i) * (yo(i, n_mod + m) * yo(i, m) + yo(i + 1, n_mod + m) * yo(i + 1, m)))
#endif
                end do
            end if
            s_abs(m) = real(sa, WP)
            s_norm(m) = real(sn, WP)
            s_esc(m) = real(se, WP)
        end do
    end subroutine dust_emission_reductions

    !> Adds the dust emission, scaled to the absorbed luminosity (energy balance).
    subroutine add_dust_emission(n_mod, nspec, emis, yo, s_abs, s_esc, total)
        integer, intent(in)              :: n_mod, nspec
        logical, intent(in)              :: emis(:)
        real(KP), intent(in)             :: yo(:,:)
        real(WP), intent(in)             :: s_abs(:), s_esc(:)
        real(KP), intent(inout)          :: total(:,:)

        integer  :: m, i

        !$acc parallel loop gang vector collapse(2) present(emis, yo, s_abs, s_esc, total)
#ifdef OMP_TARGET
        !$omp target teams distribute parallel do collapse(2)
#endif
        do m = 1, n_mod
            do i = 1, nspec
                if (emis(m) .and. s_esc(m) > SAFE_FLOOR) then
                    ! As dust_emission_final in fsps_dust: shape * transmission (KP), then * norm (WP)
                    total(i, m) = total(i, m) + real(yo(i, n_mod + m) * yo(i, m) * (s_abs(m) / s_esc(m)), KP)
                end if
            end do
        end do
    end subroutine add_dust_emission




    !> Compensated add (CPU builds at FP32 kernel precision).
    subroutine kadd(s, c, x)
        real(KP), intent(inout) :: s, c
        real(KP), intent(in)    :: x
        real(KP) :: t, y
        if (KP == WP) then
            s = s + x
        else
            y = x - c
            t = s + y
            c = (t - s) - y
            s = t
        end if
    end subroutine kadd

end module fsps_csp_batch
