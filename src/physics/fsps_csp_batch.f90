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
    use fsps_csp, only: compute_sfh_weights
    use fsps_sfh, only: get_sfh_properties_at_age
    use fsps_dust, only: compute_attenuation_curve_point, one_minus_trans
    use fsps_cosmology, only: compute_igm_transmission
    use fsps_interpolation, only: find_interval, interpolate_linear
    use fsps_special_functions, only: mag_from_flux
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none
    private

    public :: compute_csp_batch

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

        integer :: n_mod, nt, nspec, nbands, m, k, i, j, b, it_tesc, ilo, ihi
        real(WP) :: zred, dust_age_log, w, msum, lsum, mass_frac, sfr_norm, frac_linear
        real(WP) :: z_factor, target, dm, dist_term, flux
        logical  :: igm, any_emis

        ! Per-model parameters and results on the host
        real(KP), allocatable :: w_young(:,:), w_old(:,:)
        real(WP), allocatable :: weights(:,:,:), lum_linear(:)
        integer,  allocatable :: iq(:), iu(:), i_lo(:), j_lo(:), j_hi(:), do_band(:)
        real(WP), allocatable :: wq(:), wu(:), gam(:), t_obs(:)
        logical,  allocatable :: emis(:)
        real(WP), allocatable :: s_abs(:), s_norm(:), s_esc(:)

        ! Grid-only factors (FP64 formed, KP stored)
        real(KP), allocatable :: wt(:), wtl(:), wl(:), bands_t(:,:), igm_t(:)

        ! Work arrays (pixels x models) at kernel precision
        real(KP), allocatable :: total(:,:), absorbed(:,:), fband(:,:)

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
        ! 1. Host: SFH weights, mass, L_bol, SFR (FP64); per-model dust emission indices
        ! ------------------------------------------------------------------
        allocate(w_young(nt, n_mod), w_old(nt, n_mod), weights(nt, 1, 1), lum_linear(nt))
        allocate(iq(n_mod), iu(n_mod), wq(n_mod), wu(n_mod), gam(n_mod), emis(n_mod))
        lum_linear = 10.0_wp**ctx%state%ssp_basis_lbol(1:nt, 1)
        do m = 1, n_mod
            call compute_sfh_weights(ctx, psets(m), psets(m)%tage, 1, 1, weights)
            if (psets(m)%dust_tesc > SAFE_FLOOR) then
                dust_age_log = psets(m)%dust_tesc
            else
                dust_age_log = 7.0_wp
            end if
            it_tesc = max(1, min(find_interval(ctx%state%time_full, dust_age_log), nt))
            msum = 0.0_wp
            lsum = 0.0_wp
            do k = 1, nt
                w = weights(k, 1, 1)
                if (w > SAFE_FLOOR) then
                    msum = msum + w * ctx%state%ssp_basis_mass(k, 1)
                    lsum = lsum + w * lum_linear(k)
                else
                    w = 0.0_wp
                end if
                if (k <= it_tesc) then
                    w_young(k, m) = real(w, KP); w_old(k, m) = 0.0_KP
                else
                    w_young(k, m) = 0.0_KP; w_old(k, m) = real(w, KP)
                end if
            end do
            if (present(mass)) mass(m) = msum
            if (present(lbol)) then
                if (lsum > SAFE_FLOOR) then
                    lbol(m) = log10(lsum)
                else
                    lbol(m) = 0.0_wp
                end if
            end if
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
        allocate(total(nspec, n_mod), absorbed(nspec, n_mod))
        allocate(s_abs(n_mod), s_norm(n_mod), s_esc(n_mod))

        !$acc data present(ctx) copyin(psets, w_young, w_old, iq, iu, wq, wu, gam, emis, wt) &
        !$acc      create(total, absorbed, igm_t) copyout(s_abs, s_norm, s_esc)
        if (igm) call compute_igm_transmission(ctx%state%spec_lambda, zred, psets(1)%igm_factor, igm_t)

        call attenuate_batch(ctx, psets, n_mod, nspec, nt, w_young, w_old, emis, total, absorbed)
        if (any_emis) then
            call dust_emission_reductions(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, wt, &
                                          absorbed, s_abs, s_norm, s_esc)
            call add_dust_emission(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, s_abs, s_esc, total)
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
        allocate(fband(nbands, n_mod))
        !$acc data copyin(i_lo, t_obs, wl, bands_t, j_lo, j_hi, do_band, wtl) copyout(fband)
        call band_fluxes(n_mod, nspec, nbands, i_lo, t_obs, wl, bands_t, j_lo, j_hi, do_band, wtl, &
                         total, absorbed, fband)
        !$acc end data

        if (present(spec)) then
            !$acc update self(total)
            do m = 1, n_mod
                spec(1:nspec, m) = real(total(:, m), WP)
            end do
        end if
        !$acc end data

        do m = 1, n_mod
            mags(:, m) = ieee_value(1.0_wp, ieee_quiet_nan)
            do b = 1, nbands
                if (do_band(b) == 0) cycle
                flux = real(fband(b, m), WP)
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
    subroutine attenuate_batch(ctx, psets, n_mod, nspec, nt, w_young, w_old, emis, total, absorbed)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, nt
        real(KP), intent(in)             :: w_young(:,:), w_old(:,:)
        logical, intent(in)              :: emis(:)
        real(KP), intent(out)            :: total(:,:), absorbed(:,:)

        integer  :: m, i, k
        real(KP) :: young, old, bk
        real(WP) :: curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, frac_nodust
        real(WP) :: one_m_tb, one_m_to, one_m_td, a_young, spec_sum

        !$acc parallel loop gang vector collapse(2) present(ctx, psets, w_young, w_old, emis, total, absorbed) &
        !$acc private(young, old, bk, curve, tau_b, tau_o, tau_d, trans_birth, trans_old, trans_d, frac_obrun, &
        !$acc         frac_nodust, one_m_tb, one_m_to, one_m_td, a_young, spec_sum, k)
        do m = 1, n_mod
            do i = 1, nspec
                young = 0.0_KP
                old = 0.0_KP
                !$acc loop seq
                do k = 1, nt
                    bk = ctx%state%ssp_basis_spec(i, k, 1)
                    young = young + w_young(k, m) * bk
                    old = old + w_old(k, m) * bk
                end do

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
                else
                    absorbed(i, m) = 0.0_KP
                end if
            end do
        end do
    end subroutine attenuate_batch

    !> Per model: absorbed luminosity, template normalization and escaped template fraction.
    subroutine dust_emission_reductions(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, wt, &
                                        absorbed, s_abs, s_norm, s_esc)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, iq(:), iu(:)
        real(WP), intent(in)             :: wq(:), wu(:), gam(:)
        logical, intent(in)              :: emis(:)
        real(KP), intent(in)             :: wt(:), absorbed(:,:)
        real(WP), intent(out)            :: s_abs(:), s_norm(:), s_esc(:)

        integer  :: m, i
        real(KP) :: sa, sn, se, ca, cn, ce, sh0, sh1, t0, t1
        real(WP) :: curve, tau_d, trans0, trans1

        !$acc parallel loop gang present(ctx, psets, iq, iu, wq, wu, gam, emis, wt, absorbed, s_abs, s_norm, s_esc) &
        !$acc private(sa, sn, se, ca, cn, ce)
        do m = 1, n_mod
            sa = 0.0_KP; sn = 0.0_KP; se = 0.0_KP
            ca = 0.0_KP; cn = 0.0_KP; ce = 0.0_KP
            if (emis(m)) then
                !$acc loop vector reduction(+:sa, sn, se) private(sh0, sh1, t0, t1, curve, tau_d, trans0, trans1)
                do i = 1, nspec - 1
                    sh0 = dl07_shape(ctx, i, iq(m), wq(m), iu(m), wu(m), gam(m))
                    sh1 = dl07_shape(ctx, i + 1, iq(m), wq(m), iu(m), wu(m), gam(m))
                    call diffuse_transmission(ctx, psets(m), i, curve, tau_d, trans0)
                    call diffuse_transmission(ctx, psets(m), i + 1, curve, tau_d, trans1)
                    t0 = real(trans0, KP)
                    t1 = real(trans1, KP)
#ifdef _OPENACC
                    sa = sa + wt(i) * (absorbed(i, m) + absorbed(i + 1, m))
                    sn = sn + wt(i) * (sh0 + sh1)
                    se = se + wt(i) * (sh0 * t0 + sh1 * t1)
#else
                    call kadd(sa, ca, wt(i) * (absorbed(i, m) + absorbed(i + 1, m)))
                    call kadd(sn, cn, wt(i) * (sh0 + sh1))
                    call kadd(se, ce, wt(i) * (sh0 * t0 + sh1 * t1))
#endif
                end do
            end if
            s_abs(m) = real(sa, WP)
            s_norm(m) = real(sn, WP)
            s_esc(m) = real(se, WP)
        end do
    end subroutine dust_emission_reductions

    !> Adds the dust emission, scaled to the absorbed luminosity (energy balance).
    subroutine add_dust_emission(ctx, psets, n_mod, nspec, iq, iu, wq, wu, gam, emis, s_abs, s_esc, total)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in)         :: psets(:)
        integer, intent(in)              :: n_mod, nspec, iq(:), iu(:)
        real(WP), intent(in)             :: wq(:), wu(:), gam(:), s_abs(:), s_esc(:)
        logical, intent(in)              :: emis(:)
        real(KP), intent(inout)          :: total(:,:)

        integer  :: m, i
        real(WP) :: curve, tau_d, trans, norm
        real(KP) :: em

        !$acc parallel loop gang vector collapse(2) &
        !$acc present(ctx, psets, iq, iu, wq, wu, gam, emis, s_abs, s_esc, total) private(curve, tau_d, trans, norm, em)
        do m = 1, n_mod
            do i = 1, nspec
                if (emis(m) .and. s_esc(m) > SAFE_FLOOR) then
                    norm = s_abs(m) / s_esc(m)
                    call diffuse_transmission(ctx, psets(m), i, curve, tau_d, trans)
                    ! As dust_emission_final in fsps_dust: shape * transmission (KP), then * norm (WP)
                    em = real(dl07_shape(ctx, i, iq(m), wq(m), iu(m), wu(m), gam(m)) * real(trans, KP) * norm, KP)
                    total(i, m) = total(i, m) + em
                end if
            end do
        end do
    end subroutine add_dust_emission

    !> Band fluxes: observed-frame spectrum / lambda, then the trapezoid filter integral over each
    !> band's support, compensated at FP32 kernel precision.
    subroutine band_fluxes(n_mod, nspec, nbands, i_lo, t_obs, wl, bands_t, j_lo, j_hi, do_band, wt, &
                           total, fobs, fband)
        integer, intent(in)     :: n_mod, nspec, nbands, i_lo(:), j_lo(:), j_hi(:), do_band(:)
        real(WP), intent(in)    :: t_obs(:)
        real(KP), intent(in)    :: wl(:), bands_t(:,:), wt(:), total(:,:)
        real(KP), intent(inout) :: fobs(:,:)
        real(KP), intent(out)   :: fband(:,:)

        integer  :: m, j, b, il
        real(WP) :: v
        real(KP) :: s, c, term, tsum

        !$acc parallel loop gang vector collapse(2) present(i_lo, t_obs, wl, total, fobs) private(il, v)
        do m = 1, n_mod
            do j = 1, nspec
                il = i_lo(j)
                v = real(total(il, m), WP) + t_obs(j) * (real(total(il + 1, m), WP) - real(total(il, m), WP))
                if (v /= v) v = 0.0_wp
                fobs(j, m) = real(max(v, 0.0_wp), KP) * wl(j)
            end do
        end do

        !$acc parallel loop gang vector collapse(2) present(bands_t, j_lo, j_hi, do_band, wt, fobs, fband) &
        !$acc private(s, c, term, tsum, j)
        do m = 1, n_mod
            do b = 1, nbands
                s = 0.0_KP
                c = 0.0_KP
                if (do_band(b) /= 0) then
                    !$acc loop seq
                    do j = j_lo(b), j_hi(b)
                        if (KP == WP) then
                            s = s + wt(j) * (fobs(j, m) * bands_t(b, j) + fobs(j + 1, m) * bands_t(b, j + 1))
                        else
                            term = wt(j) * (fobs(j, m) * bands_t(b, j) + fobs(j + 1, m) * bands_t(b, j + 1)) - c
                            tsum = s + term
                            c = (tsum - s) - term
                            s = tsum
                        end if
                    end do
                end if
                fband(b, m) = s
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
