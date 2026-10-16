module fsps_csp_batch
    !> @brief Batched CSP evaluation: many parameter sets against one resident SSP basis.
    !>
    !> @details
    !> The fitting workload: N models that share the libraries, the SSP basis (one metallicity)
    !> and the redshift, and differ in their star formation history and dust parameters. Per
    !> model the work is the same as `compute_csp_scenario` with `tage > 0`; here each stage
    !> runs once for the whole batch:
    !>
    !> 1. SFH weights per model on the host (FP64), split into young and old ages by each
    !>    model's `dust_tesc`, cast once to kernel precision: W_young, W_old (ages x N).
    !> 2. Age sums, attenuation and the absorbed luminosity density in one kernel over
    !>    (pixel, model): young = basis . W_young, old = basis . W_old (matrix products with
    !>    K = number of ages), then the birth-cloud and diffuse attenuation of `fsps_dust`.
    !> 3. Dust emission (energy balance): per-model reductions over pixels, then the
    !>    Draine & Li template, interpolated per model from the resident grid, is added.
    !> 4. Magnitudes: per (band, model) filter integrals over each band's support, compensated
    !>    at FP32 kernel precision.
    !>
    !> The two matrix products (age sums, band fluxes) use cuBLAS when built with -Dcublas=true
    !> (NVIDIA only; a ceiling for the portable OpenACC loops, which are the default).
    !>
    !> Not supported (status 2): nebular emission, AGN, velocity smoothing, redshifted colors,
    !> light-weighted ages, and models with tage <= 0 (all-ages output). Models must share
    !> zred and igm_factor. Spectral indices are not computed.
    !>
    !> Mass, L_bol and SFR follow `compute_csp_scenario` in FP64; spectra, magnitudes and dust
    !> masses agree with it to FP32 rounding with FP32 kernels (summation order differs).

    use fsps_precision, only: WP, KP
    use fsps_constants, only: SAFE_FLOOR, C_LIGHT, PI, ABS_MAG_ZEROPOINT_LOG
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t, fsps_context_sync_device_settings
    use fsps_csp, only: sfh_weights_core
    use fsps_sfh, only: get_sfh_properties_at_age
    use fsps_dust, only: compute_attenuation_curve_point, one_minus_trans
    use fsps_cosmology, only: compute_igm_transmission
    use fsps_interpolation, only: find_interval, interpolate_linear
    use fsps_special_functions, only: mag_from_flux
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#ifdef FSPS_CUBLAS
    use cublas
    use openacc, only: acc_get_cuda_stream, acc_async_sync
#endif

    implicit none
    private

    public :: compute_csp_batch

    !> Pixels per block of the split-K band product (cuBLAS path): each block is a plain FP32 sum
    !> of at most 256 positive terms; the blocks are added in FP64
    integer, parameter :: BAND_BLOCK = 256
#ifdef FSPS_CUBLAS
    type(cublasHandle), save :: blas
    logical, save :: blas_ready = .false.
#endif

    real(WP), parameter :: AB_ZEROPOINT = 48.60_wp
    real(WP), parameter :: V_BAND_ANGSTROMS = 5500.0_wp   ! as in fsps_dust
    integer, parameter  :: IDX_V_BAND = 1

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
    subroutine compute_csp_batch(ctx, psets, mags, status, spec, mass, lbol, mdust, sfr)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in)            :: psets(:)
        real(WP), intent(out)               :: mags(:,:)
        integer, intent(out)                :: status
        real(WP), intent(out), optional     :: spec(:,:)
        real(WP), intent(out), optional     :: mass(:), lbol(:), mdust(:), sfr(:)

        integer :: n_mod, nt, nspec, nbands, m, k, i, j, b, ilo, ihi
        real(WP) :: zred, dust_age_log, w, mass_frac, sfr_norm, frac_linear
        real(WP) :: z_factor, target, dm, dist_term, flux
        logical  :: igm, any_emis

        ! Per-model parameters and results on the host
        real(KP), allocatable :: w2(:,:)   ! (ages, 2N): young weights, then old
        real(WP), allocatable :: weights(:,:,:), tmp1(:,:), tmp2(:,:), lum_linear(:), msum(:), lsum(:)
        integer,  allocatable :: it_tesc(:)
        integer,  allocatable :: iq(:), iu(:), i_lo(:), j_lo(:), j_hi(:), do_band(:)
        real(WP), allocatable :: wq(:), wu(:), gam(:), t_obs(:)
        logical,  allocatable :: emis(:)
        real(WP), allocatable :: s_abs(:), s_norm(:), s_esc(:)

        ! Grid-only factors (FP64 formed, KP stored)
        real(KP), allocatable :: wt(:), wtl(:), wl(:), bands_t(:,:), igm_t(:)

        ! Work arrays (pixels x models) at kernel precision
        real(KP), allocatable :: yo(:,:), total(:,:), absorbed(:,:), fobs_t(:,:)
        real(WP), allocatable :: fband(:,:)
        integer :: npad

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
        if (ctx%add_neb_emission_val /= 0 .or. ctx%redshift_colors_val /= 0 .or. &
            ctx%compute_light_ages_val /= 0) then
            status = 2; return
        end if
        zred = psets(1)%zred
        do m = 1, n_mod
            if (psets(m)%tage <= 0.0_wp .or. psets(m)%sigma_smooth > 0.0_wp .or. &
                psets(m)%zred /= zred .or. psets(m)%igm_factor /= psets(1)%igm_factor .or. &
                psets(m)%uvb < 0.0_wp .or. psets(m)%wgp1 < 1 .or. psets(m)%wgp2 < 1 .or. &
                (ctx%add_agn_dust_val == 1 .and. psets(m)%fagn > SAFE_FLOOR)) then
                status = 2; return
            end if
        end do

        ! ------------------------------------------------------------------
        ! 1. Host: per-model indices (birth-cloud age, dust emission grid) and SFR
        ! ------------------------------------------------------------------
        allocate(w2(nt, 2 * n_mod), lum_linear(nt), it_tesc(n_mod), msum(n_mod), lsum(n_mod))
        allocate(weights(nt, 1, n_mod), tmp1(size(ctx%state%sfh_w_tmp1, 1), n_mod), &
                 tmp2(size(ctx%state%sfh_w_tmp2, 1), n_mod))
        allocate(iq(n_mod), iu(n_mod), wq(n_mod), wu(n_mod), gam(n_mod), emis(n_mod))
        lum_linear = 10.0_wp**ctx%state%ssp_basis_lbol(1:nt, 1)
        do m = 1, n_mod
            if (psets(m)%dust_tesc > SAFE_FLOOR) then
                dust_age_log = psets(m)%dust_tesc
            else
                dust_age_log = 7.0_wp
            end if
            it_tesc(m) = max(1, min(find_interval(ctx%state%time_full, dust_age_log), nt))
            if (present(sfr)) then
                call get_sfh_properties_at_age(ctx, psets(m), psets(m)%tage, mass_frac, sfr_norm, frac_linear)
                if (mass_frac > SAFE_FLOOR) sfr_norm = sfr_norm / mass_frac
                sfr(m) = sfr_norm
            end if

            emis(m) = ctx%add_dust_emission_val == 1 .and. &
                      (psets(m)%dust1 > SAFE_FLOOR .or. psets(m)%dust2 > SAFE_FLOOR)
            call dl07_indices(ctx, psets(m), iq(m), wq(m), iu(m), wu(m), gam(m))
        end do
        any_emis = any(emis)

        ! Grid-only trapezoid weights in frequency (formed in FP64, rule 4)
        allocate(wt(nspec))
        do i = 1, nspec - 1
            wt(i) = real(0.5_wp * abs(C_LIGHT / ctx%state%spec_lambda(i+1) - C_LIGHT / ctx%state%spec_lambda(i)), KP)
        end do
        wt(nspec) = 0.0_KP

        ! IGM transmission (common redshift; computed on the device in OpenACC builds, below)
        igm = ctx%add_igm_absorption_val == 1 .and. zred > SAFE_FLOOR
        allocate(igm_t(nspec))

        ! ------------------------------------------------------------------
        ! 2-3. Device: age sums, attenuation, dust emission
        ! ------------------------------------------------------------------
        ! The absorbed-light array later holds the observed spectrum / lambda for the band product,
        ! padded with zero rows to whole pixel blocks on the cuBLAS path
        npad = nspec
#ifdef FSPS_CUBLAS
        npad = BAND_BLOCK * ((nspec + BAND_BLOCK - 1) / BAND_BLOCK)
#endif
        allocate(yo(nspec, 2 * n_mod), total(nspec, n_mod), absorbed(npad, n_mod))
        allocate(s_abs(n_mod), s_norm(n_mod), s_esc(n_mod))

        !$acc data present(ctx) copyin(psets, iq, iu, wq, wu, gam, emis, wt, it_tesc, lum_linear) &
        !$acc      create(w2, weights, tmp1, tmp2, yo, total, absorbed, igm_t) copyout(s_abs, s_norm, s_esc, msum, lsum)
        if (igm) call compute_igm_transmission(ctx%state%spec_lambda, zred, psets(1)%igm_factor, igm_t)

        ! SFH weights (FP64), mass and L_bol per model, one thread per model; the weights are split
        ! into young and old ages at the model's dust_tesc and cast once to kernel precision
        !$acc parallel loop gang vector present(ctx, psets, weights, tmp1, tmp2, w2, it_tesc, lum_linear, msum, lsum) &
        !$acc private(k, w)
        do m = 1, n_mod
            call sfh_weights_core(ctx, psets(m), psets(m)%tage, 1, weights(:, :, m), tmp1(:, m), tmp2(:, m))
            msum(m) = 0.0_wp
            lsum(m) = 0.0_wp
            !$acc loop seq
            do k = 1, nt
                w = weights(k, 1, m)
                if (w > SAFE_FLOOR) then
                    msum(m) = msum(m) + w * ctx%state%ssp_basis_mass(k, 1)
                    lsum(m) = lsum(m) + w * lum_linear(k)
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

        call age_sums(nspec, nt, n_mod, ctx%state%ssp_basis_spec(:, :, 1), w2, it_tesc, yo)
        call attenuate_batch(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, yo, total, absorbed)
        if (any_emis) then
            call dust_emission_reductions(n_mod, nspec, emis, wt, yo, absorbed, s_abs, s_norm, s_esc)
            call add_dust_emission(n_mod, nspec, emis, yo, s_abs, s_esc, total)
        end if

        ! IGM (common transmission)
        if (igm) then
            !$acc parallel loop collapse(2) present(total, igm_t)
            do m = 1, n_mod
                do i = 1, nspec
                    total(i, m) = total(i, m) * igm_t(i)
                end do
            end do
        end if

        ! ------------------------------------------------------------------
        ! 4. Magnitudes
        ! ------------------------------------------------------------------
        ! Observed-frame interpolation (prepare_observed_spectrum): per pixel, the rest-frame
        ! interval and position; then flux / lambda into the absorbed-light work array
        allocate(i_lo(nspec), t_obs(nspec), wl(nspec))
        if (abs(zred) > SAFE_FLOOR) then
            z_factor = 1.0_wp + zred
            do j = 1, nspec
                target = ctx%state%spec_lambda(j) / z_factor
                i_lo(j) = max(1, min(find_interval(ctx%state%spec_lambda, target), nspec - 1))
                t_obs(j) = (target - ctx%state%spec_lambda(i_lo(j))) / &
                           (ctx%state%spec_lambda(i_lo(j)+1) - ctx%state%spec_lambda(i_lo(j)))
            end do
            dm = interpolate_linear(ctx%state%cosmospl(:,1), ctx%state%cosmospl(:,3), zred)
            if (dm /= dm .or. dm <= SAFE_FLOOR) then
                dist_term = 0.0_wp
            else
                dist_term = 5.0_wp * log10(dm / 10.0_wp) + mag_from_flux(z_factor)
            end if
        else
            i_lo = [(j, j = 1, nspec)]
            t_obs = 0.0_wp
            dist_term = 0.0_wp
        end if
        do j = 1, nspec
            if (ctx%state%spec_lambda(j) > tiny(0.0_wp)) then
                wl(j) = real(1.0_wp / ctx%state%spec_lambda(j), KP)
            else
                wl(j) = 0.0_KP
            end if
        end do

        ! Bands: transposed at kernel precision, the support of each, trapezoid weights in lambda
        allocate(bands_t(nbands, nspec), j_lo(nbands), j_hi(nbands), do_band(nbands))
        do_band = 1
        if (allocated(psets(1)%mag_compute)) then
            b = min(nbands, size(psets(1)%mag_compute))
            do_band(1:b) = psets(1)%mag_compute(1:b)
            if (ctx%compute_vega_mags_val == 1) do_band(IDX_V_BAND) = 1
        end if
        do b = 1, nbands
            bands_t(b, :) = real(ctx%state%bands(1:nspec, b), KP)
            ilo = nspec; ihi = 1
            do j = 1, nspec
                if (ctx%state%bands(j, b) /= 0.0_wp) then
                    ilo = min(ilo, j); ihi = max(ihi, j)
                end if
            end do
            ! Trapezoid intervals touching the support: j = ilo-1 .. ihi
            j_lo(b) = max(1, ilo - 1)
            j_hi(b) = min(nspec - 1, ihi)
            if (ilo > ihi) then
                j_lo(b) = 1; j_hi(b) = 0
            end if
        end do
        allocate(wtl(nspec))
        do j = 1, nspec - 1
            wtl(j) = real(0.5_wp * abs(ctx%state%spec_lambda(j + 1) - ctx%state%spec_lambda(j)), KP)
        end do
        wtl(nspec) = 0.0_KP

        ! The absorbed-light array is free again: it holds the observed spectrum / lambda
        allocate(fband(nbands, n_mod), fobs_t(n_mod, nspec))
        !$acc data copyin(i_lo, t_obs, wl, bands_t, j_lo, j_hi, do_band, wtl) copyout(fband) create(fobs_t)
#ifdef FSPS_CUBLAS
        call observed_frame(n_mod, nspec, npad, i_lo, t_obs, wl, total, absorbed)
        call band_fluxes_gemm(ctx, n_mod, nspec, npad, nbands, wtl, absorbed, fband)
#else
        call observed_frame_t(n_mod, nspec, i_lo, t_obs, wl, total, fobs_t)
        call band_fluxes(n_mod, nbands, bands_t, j_lo, j_hi, do_band, wtl, fobs_t, fband)
#endif
        !$acc end data

        if (present(spec)) then
            !$acc update self(total)
            do m = 1, n_mod
                spec(1:nspec, m) = real(total(:, m), WP)
            end do
        end if
        !$acc end data

        do m = 1, n_mod
            if (present(mass)) mass(m) = msum(m)
            if (present(lbol)) then
                if (lsum(m) > SAFE_FLOOR) then
                    lbol(m) = log10(lsum(m))
                else
                    lbol(m) = 0.0_wp
                end if
            end if
            mags(:, m) = ieee_value(1.0_wp, ieee_quiet_nan)
            do b = 1, nbands
                if (do_band(b) == 0) cycle
                flux = fband(b, m)
                mags(b, m) = mag_from_flux(flux) - AB_ZEROPOINT - 2.5_wp * ABS_MAG_ZEROPOINT_LOG + dist_term
            end do
            if (ctx%compute_vega_mags_val == 1) then
                do b = 2, nbands
                    if (mags(IDX_V_BAND, m) == mags(IDX_V_BAND, m)) then
                        mags(b, m) = mags(b, m) - ctx%state%magvega(b) + ctx%state%magvega(IDX_V_BAND)
                    else
                        mags(b, m) = ieee_value(1.0_wp, ieee_quiet_nan)
                    end if
                end do
            end if
            if (present(mdust)) then
                if (emis(m) .and. s_norm(m) > SAFE_FLOOR) then
                    mdust(m) = 3.21e-3_wp / (4.0_wp * PI) * (s_abs(m) / s_norm(m))
                else
                    mdust(m) = SAFE_FLOOR
                end if
            end if
        end do
    end subroutine compute_csp_batch

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
        ! Each age is young or old for a given model, so one pass with one product per age; each
        ! thread does 4 models of one pixel, reusing every basis value it loads. Same sums, in the
        ! same order, as the two separate products.
        integer, parameter :: MT = 4
        integer  :: mq, i, k, m1, m2, m3, m4
        real(KP) :: bk, y1, y2, y3, y4, o1, o2, o3, o4

        !$acc parallel loop gang vector collapse(2) present(basis, w2, it_tesc, yo) &
        !$acc private(k, bk, m1, m2, m3, m4, y1, y2, y3, y4, o1, o2, o3, o4)
        do mq = 1, (n_mod + MT - 1) / MT
            do i = 1, nspec
                ! Models of this thread (clamped at the end; duplicates are written twice, same value)
                m1 = min((mq - 1) * MT + 1, n_mod)
                m2 = min((mq - 1) * MT + 2, n_mod)
                m3 = min((mq - 1) * MT + 3, n_mod)
                m4 = min((mq - 1) * MT + 4, n_mod)
                y1 = 0.0_KP; y2 = 0.0_KP; y3 = 0.0_KP; y4 = 0.0_KP
                o1 = 0.0_KP; o2 = 0.0_KP; o3 = 0.0_KP; o4 = 0.0_KP
                !$acc loop seq
                do k = 1, nt
                    bk = basis(i, k)
                    if (k <= it_tesc(m1)) then
                        y1 = y1 + w2(k, m1) * bk
                    else
                        o1 = o1 + w2(k, n_mod + m1) * bk
                    end if
                    if (k <= it_tesc(m2)) then
                        y2 = y2 + w2(k, m2) * bk
                    else
                        o2 = o2 + w2(k, n_mod + m2) * bk
                    end if
                    if (k <= it_tesc(m3)) then
                        y3 = y3 + w2(k, m3) * bk
                    else
                        o3 = o3 + w2(k, n_mod + m3) * bk
                    end if
                    if (k <= it_tesc(m4)) then
                        y4 = y4 + w2(k, m4) * bk
                    else
                        o4 = o4 + w2(k, n_mod + m4) * bk
                    end if
                end do
                yo(i, m1) = y1; yo(i, n_mod + m1) = o1
                yo(i, m2) = y2; yo(i, n_mod + m2) = o2
                yo(i, m3) = y3; yo(i, n_mod + m3) = o3
                yo(i, m4) = y4; yo(i, n_mod + m4) = o4
            end do
        end do
#endif
    end subroutine age_sums

    !> Attenuation and absorbed light per (pixel, model). For models with dust emission the age
    !> sums are then replaced by what the emission stages need: yo(:, m) = diffuse transmission
    !> (stored at KP, as fsps_dust stores it), yo(:, N + m) = Draine & Li template.
    subroutine attenuate_batch(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, yo, total, absorbed)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, iq(:), iu(:)
        real(WP), intent(in)             :: wq(:), wu(:), gam(:)
        logical, intent(in)              :: emis(:)
        real(KP), intent(inout)          :: yo(:,:)   ! in: young (1:N), old (N+1:2N)
        real(KP), intent(out)            :: total(:,:), absorbed(:,:)

        integer  :: m, i
        real(KP) :: young, old
        real(WP) :: curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, frac_nodust
        real(WP) :: one_m_tb, one_m_to, one_m_td, a_young, spec_sum

        !$acc parallel loop gang vector collapse(2) present(ctx, psets, iq, iu, wq, wu, gam, yo, emis, total, absorbed) &
        !$acc private(young, old, curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, &
        !$acc         frac_nodust, one_m_tb, one_m_to, one_m_td, a_young, spec_sum)
        do m = 1, n_mod
            do i = 1, nspec
                young = yo(i, m)
                old = yo(i, n_mod + m)

                frac_obrun = psets(m)%frac_obrun
                frac_nodust = psets(m)%frac_nodust
                call diffuse_transmission(ctx, psets(m), i, curve, tau_d, trans_d)
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
        do m = 1, n_mod
            sa = 0.0_KP; sn = 0.0_KP; se = 0.0_KP
            ca = 0.0_KP; cn = 0.0_KP; ce = 0.0_KP
            if (emis(m)) then
                !$acc loop vector reduction(+:sa, sn, se)
                do i = 1, nspec - 1
#ifdef _OPENACC
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
        do m = 1, n_mod
            do i = 1, nspec
                if (emis(m) .and. s_esc(m) > SAFE_FLOOR) then
                    ! As dust_emission_final in fsps_dust: shape * transmission (KP), then * norm (WP)
                    total(i, m) = total(i, m) + real(yo(i, n_mod + m) * yo(i, m) * (s_abs(m) / s_esc(m)), KP)
                end if
            end do
        end do
    end subroutine add_dust_emission

    !> Observed-frame spectrum / lambda (prepare_observed_spectrum), zero in the padding rows.
    subroutine observed_frame(n_mod, nspec, npad, i_lo, t_obs, wl, total, fobs)
        integer, intent(in)     :: n_mod, nspec, npad, i_lo(:)
        real(WP), intent(in)    :: t_obs(:)
        real(KP), intent(in)    :: wl(:), total(:,:)
        real(KP), intent(out)   :: fobs(:,:)

        integer  :: m, j, il
        real(WP) :: v

        !$acc parallel loop gang vector collapse(2) present(i_lo, t_obs, wl, total, fobs) private(il, v)
        do m = 1, n_mod
            do j = 1, npad
                if (j <= nspec) then
                    il = i_lo(j)
                    if (t_obs(j) == 0.0_wp) then   ! rest frame (i_lo = j: no right neighbour at j = nspec)
                        v = real(total(il, m), WP)
                    else
                        v = real(total(il, m), WP) + t_obs(j) * (real(total(il + 1, m), WP) - real(total(il, m), WP))
                    end if
                    if (v /= v) v = 0.0_wp
                    fobs(j, m) = real(max(v, 0.0_wp), KP) * wl(j)
                else
                    fobs(j, m) = 0.0_KP
                end if
            end do
        end do
    end subroutine observed_frame

    !> observed_frame, written transposed (models, pixels) for band_fluxes.
    subroutine observed_frame_t(n_mod, nspec, i_lo, t_obs, wl, total, fobs_t)
        integer, intent(in)     :: n_mod, nspec, i_lo(:)
        real(WP), intent(in)    :: t_obs(:)
        real(KP), intent(in)    :: wl(:), total(:,:)
        real(KP), intent(out)   :: fobs_t(:,:)

        integer  :: m, j, il
        real(WP) :: v

        !$acc parallel loop gang vector collapse(2) present(i_lo, t_obs, wl, total, fobs_t) private(il, v)
        do j = 1, nspec
            do m = 1, n_mod
                il = i_lo(j)
                if (t_obs(j) == 0.0_wp) then
                    v = real(total(il, m), WP)
                else
                    v = real(total(il, m), WP) + t_obs(j) * (real(total(il + 1, m), WP) - real(total(il, m), WP))
                end if
                if (v /= v) v = 0.0_wp
                fobs_t(m, j) = real(max(v, 0.0_wp), KP) * wl(j)
            end do
        end do
    end subroutine observed_frame_t

#ifdef FSPS_CUBLAS
    !> Band fluxes as a split-K matrix product: per block of BAND_BLOCK pixels,
    !> partial = C^T . fobs with C(j, b) = (dl(j-1) + dl(j)) / 2 * band(j, b) (the trapezoid rule
    !> regrouped by pixel, formed in FP64); the blocks are summed in FP64.
    subroutine band_fluxes_gemm(ctx, n_mod, nspec, npad, nbands, wtl, fobs, fband)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)     :: n_mod, nspec, npad, nbands
        real(KP), intent(in)    :: wtl(:), fobs(:,:)
        real(WP), intent(out)   :: fband(:,:)

        real(KP), allocatable :: cmat(:,:), partial(:,:,:)
        integer  :: j, b, m, q, nkb, istat
        real(WP) :: sw

        nkb = npad / BAND_BLOCK
        allocate(cmat(npad, nbands), partial(nbands, n_mod, nkb))
        do b = 1, nbands
            do j = 1, npad
                if (j > nspec) then
                    cmat(j, b) = 0.0_KP
                else
                    sw = real(wtl(j), WP)
                    if (j > 1) sw = sw + real(wtl(j - 1), WP)
                    cmat(j, b) = real(sw * ctx%state%bands(j, b), KP)
                end if
            end do
        end do
        call blas_init()
        !$acc data copyin(cmat) create(partial) present(fobs, fband)
        !$acc host_data use_device(cmat, fobs, partial)
        do q = 1, nkb
#ifdef KERNEL_SINGLE_PRECISION
            istat = cublasSgemm_v2(blas, CUBLAS_OP_T, CUBLAS_OP_N, nbands, n_mod, BAND_BLOCK, 1.0_KP, &
                                   cmat((q - 1) * BAND_BLOCK + 1, 1), npad, fobs((q - 1) * BAND_BLOCK + 1, 1), npad, &
                                   0.0_KP, partial(1, 1, q), nbands)
#else
            istat = cublasDgemm_v2(blas, CUBLAS_OP_T, CUBLAS_OP_N, nbands, n_mod, BAND_BLOCK, 1.0_KP, &
                                   cmat((q - 1) * BAND_BLOCK + 1, 1), npad, fobs((q - 1) * BAND_BLOCK + 1, 1), npad, &
                                   0.0_KP, partial(1, 1, q), nbands)
#endif
            if (istat /= 0) exit
        end do
        !$acc end host_data
        if (istat /= 0) error stop 'compute_csp_batch: cuBLAS GEMM failed'
        !$acc parallel loop gang vector collapse(2) present(partial, fband) private(sw, q)
        do m = 1, n_mod
            do b = 1, nbands
                sw = 0.0_wp
                !$acc loop seq
                do q = 1, nkb
                    sw = sw + real(partial(b, m, q), WP)
                end do
                fband(b, m) = sw
            end do
        end do
        !$acc end data
    end subroutine band_fluxes_gemm

    subroutine blas_init()
        integer :: istat
        if (blas_ready) return
        istat = cublasCreate(blas)
        istat = cublasSetStream(blas, acc_get_cuda_stream(acc_async_sync))
        blas_ready = .true.
    end subroutine blas_init
#endif

    !> Band fluxes: the trapezoid filter integral of fobs over each band's support, compensated at
    !> FP32 kernel precision.
    subroutine band_fluxes(n_mod, nbands, bands_t, j_lo, j_hi, do_band, wt, fobs_t, fband)
        integer, intent(in)     :: n_mod, nbands, j_lo(:), j_hi(:), do_band(:)
        !> Observed spectrum / lambda, transposed (models, pixels): the threads of a gang read
        !> consecutive models
        real(KP), intent(in)    :: bands_t(:,:), wt(:), fobs_t(:,:)
        real(WP), intent(out)   :: fband(:,:)

        !> Models per gang: the threads of a gang share one band, hence one trip count (no
        !> divergence)
        integer, parameter :: MBLK = 128
        integer  :: m, j, b, mb, ml
        real(KP) :: s, c, term, tsum

        !$acc parallel loop gang collapse(2) vector_length(MBLK) present(bands_t, j_lo, j_hi, do_band, wt, fobs_t, fband)
        do b = 1, nbands
            do mb = 1, (n_mod + MBLK - 1) / MBLK
                !$acc loop vector private(m, s, c, term, tsum, j)
                do ml = 1, MBLK
                    m = (mb - 1) * MBLK + ml
                    if (m <= n_mod) then
                        s = 0.0_KP
                        c = 0.0_KP
                        if (do_band(b) /= 0) then
                            !$acc loop seq
                            do j = j_lo(b), j_hi(b)
                                if (KP == WP) then
                                    s = s + wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1))
                                else
                                    term = wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1)) - c
                                    tsum = s + term
                                    c = (tsum - s) - term
                                    s = tsum
                                end if
                            end do
                        end if
                        fband(b, m) = real(s, WP)
                    end if
                end do
            end do
        end do
    end subroutine band_fluxes

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
