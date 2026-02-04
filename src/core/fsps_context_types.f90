MODULE FSPS_CONTEXT_TYPES
  USE fsps_constants, ONLY: SP, &
       NDIM_LOGT, NDIM_LOGG, NDIM_WMB_LOGT, NDIM_WMB_LOGG, &
       N_AGB_CAR, NDIM_PAGB, NDIM_WR, NTAU_DAGB, NTEFF_DAGB, &
       NEMLINE, NEBNZ, NEBNAGE, NEBNIP, NAGNDUST, NTABMAX
  USE fsps_types, ONLY: PARAMS, COMPSPOUT, TLSF, OBSDAT
  USE fsps_cache, ONLY: fsps_setup_cache_t
  IMPLICIT NONE

  TYPE :: fsps_context_state_t
     REAL(SP) :: zsol = 0.0
     REAL(SP) :: zsol_spec = 0.0
    CHARACTER(LEN=64) :: isoc_type = ''
    CHARACTER(LEN=64) :: spec_type = ''
     INTEGER :: nt = 0
     INTEGER :: nz = 0
     INTEGER :: nspec = 0
     INTEGER :: nzinit = 0
     INTEGER :: nbands = 0
     INTEGER :: nindx = 0
     INTEGER :: ntfull = 0
     INTEGER :: nspec_xrb = 0
     INTEGER :: nt_xrb = 0
     INTEGER :: nz_xrb = 0
     INTEGER :: check_sps_setup = 0
     REAL(SP) :: tuniv = 0.0
     INTEGER :: whlam5000 = 0
     INTEGER :: whlylim = 0
     REAL(SP) :: zpow2 = 1.0
     INTEGER, DIMENSION(6) :: mwdindex = 0
     REAL(SP), DIMENSION(500,3) :: cosmospl = 0.0
     INTEGER :: ntabsfh = 0
     REAL(SP), DIMENSION(3,NTABMAX) :: sfh_tab = 0.0
     REAL(SP), DIMENSION(3) :: imf_alpha = 1.3
     REAL(SP) :: imf_vdmc = 0.08
     REAL(SP) :: imf_mdave = 0.5
     INTEGER :: n_user_imf = 0
     REAL(SP), DIMENSION(3,100) :: imf_user_alpha = 0.0
     REAL(SP) :: salp_ind = 2.35
     REAL(SP) :: imf_lower_limit = 0.08
     REAL(SP) :: imf_upper_limit = 120.0
     REAL(SP) :: imf_lower_bound = 0.0
     REAL(SP) :: mlim_bh = 40.0
     REAL(SP) :: mlim_ns = 8.5
     CHARACTER(30) :: alt_filter_file = ''
    REAL(SP), POINTER :: indexdefined(:,:) => NULL()
    REAL(SP), POINTER :: wgdust(:,:,:,:) => NULL()
    REAL(SP), POINTER :: g03smcextn(:) => NULL()
    REAL(SP), POINTER :: bands(:,:) => NULL()
    REAL(SP), POINTER :: magsun(:) => NULL()
    REAL(SP), POINTER :: magvega(:) => NULL()
    REAL(SP), POINTER :: filter_leff(:) => NULL()
    REAL(SP), POINTER :: vega_spec(:) => NULL()
    REAL(SP), POINTER :: sun_spec(:) => NULL()
    REAL(SP), POINTER :: spec_lambda(:) => NULL()
    REAL(SP), POINTER :: spec_nu(:) => NULL()
    REAL(SP), POINTER :: spec_res(:) => NULL()
     REAL(SP), DIMENSION(NDIM_LOGT) :: speclib_logt = 0.0
     REAL(SP), DIMENSION(NDIM_LOGG) :: speclib_logg = 0.0
    REAL(KIND(1.0)), POINTER :: speclib(:,:,:,:) => NULL()
     REAL(SP), DIMENSION(NDIM_WMB_LOGT) :: wmb_logt = 0.0
     REAL(SP), DIMENSION(NDIM_WMB_LOGG) :: wmb_logg = 0.0
    REAL(KIND(1.0)), POINTER :: wmb_spec(:,:,:,:) => NULL()
    REAL(SP), POINTER :: agb_spec_o(:,:) => NULL()
    REAL(SP), POINTER :: agb_logt_o(:,:) => NULL()
    REAL(SP), POINTER :: agb_spec_c(:,:) => NULL()
    REAL(SP), POINTER :: agb_logt_c(:) => NULL()
     REAL(SP), DIMENSION(N_AGB_CAR) :: agb_logt_car = 0.0
    REAL(SP), POINTER :: agb_spec_car(:,:) => NULL()
    REAL(SP), POINTER :: pagb_spec(:,:,:) => NULL()
     REAL(SP), DIMENSION(NDIM_PAGB) :: pagb_logt = 0.0
    REAL(SP), POINTER :: wrn_spec(:,:,:) => NULL()
    REAL(SP), POINTER :: wrc_spec(:,:,:) => NULL()
     REAL(SP), DIMENSION(NDIM_WR) :: wrn_logt = 0.0
     REAL(SP), DIMENSION(NDIM_WR) :: wrc_logt = 0.0
     INTEGER :: ndim_dustem = 0
     INTEGER :: numin_dustem = 0
     INTEGER :: nqpah_dustem = 0
     CHARACTER(6) :: str_dustem = 'DL07'
    REAL(SP), POINTER :: qpaharr(:) => NULL()
    REAL(SP), POINTER :: uminarr(:) => NULL()
    REAL(SP), POINTER :: lambda_dustem(:) => NULL()
    REAL(SP), POINTER :: dustem_dustem(:,:) => NULL()
    REAL(SP), POINTER :: dustem2_dustem(:,:,:) => NULL()
    REAL(SP), POINTER :: flux_dagb(:,:,:,:) => NULL()
     REAL(SP), DIMENSION(2,NTAU_DAGB) :: tau1_dagb = 0.0
     REAL(SP), DIMENSION(2,NTEFF_DAGB) :: teff_dagb = 0.0
     REAL(SP), DIMENSION(NEMLINE) :: nebem_line_pos = 0.0
     REAL(SP), DIMENSION(NEMLINE,NEBNZ,NEBNAGE,NEBNIP) :: nebem_line = 0.0
     REAL(SP), DIMENSION(NEMLINE,NEBNZ,NEBNAGE,NEBNIP) :: xnebem_line = 0.0
    REAL(SP), POINTER :: nebem_cont(:,:,:,:) => NULL()
    REAL(SP), POINTER :: xnebem_cont(:,:,:,:) => NULL()
     REAL(SP), DIMENSION(NEBNZ) :: nebem_logz = 0.0
     REAL(SP), DIMENSION(NEBNAGE) :: nebem_age = 0.0
     REAL(SP), DIMENSION(NEBNIP) :: nebem_logu = 0.0
    REAL(SP), POINTER :: neb_res_min(:) => NULL()
    REAL(SP), POINTER :: gaussnebarr(:,:) => NULL()
     REAL(SP), DIMENSION(NAGNDUST) :: agndust_tau = 0.0
    REAL(SP), POINTER :: agndust_spec(:,:) => NULL()
    REAL(SP), POINTER :: mact_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: logl_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: logt_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: logg_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: ffco_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: phase_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: mini_isoc(:,:,:) => NULL()
    REAL(SP), POINTER :: lmdot_isoc(:,:,:) => NULL()
    INTEGER, POINTER :: nmass_isoc(:,:) => NULL()
    REAL(SP), POINTER :: timestep_isoc(:,:) => NULL()
    REAL(SP), POINTER :: zlegend(:) => NULL()
    REAL(SP), POINTER :: zlegendinit(:) => NULL()
     REAL(SP), ALLOCATABLE :: spec_ssp_zz(:,:,:)
     REAL(SP), ALLOCATABLE :: mass_ssp_zz(:,:)
     REAL(SP), ALLOCATABLE :: lbol_ssp_zz(:,:)
    REAL(SP), POINTER :: time_full(:) => NULL()
     REAL(SP), ALLOCATABLE :: weight_ssp(:,:)
     REAL(SP), ALLOCATABLE :: spec_young(:)
     REAL(SP), ALLOCATABLE :: spec_old(:)
    REAL(SP), POINTER :: bpass_spec_ssp(:,:,:) => NULL()
    REAL(SP), POINTER :: bpass_mass_ssp(:,:) => NULL()
    REAL(SP), POINTER :: lam_xrb(:) => NULL()
    REAL(SP), POINTER :: spec_xrb(:,:,:) => NULL()
    REAL(SP), POINTER :: ages_xrb(:) => NULL()
    REAL(SP), POINTER :: zmet_xrb(:) => NULL()
     TYPE(TLSF) :: lsfinfo
     TYPE(OBSDAT) :: powell_data
     TYPE(OBSDAT) :: sedfit_data
  END TYPE fsps_context_state_t

  TYPE :: fsps_context_t
     LOGICAL :: initialized = .FALSE.
     INTEGER :: zin = 0
     CHARACTER(LEN=64) :: isoc_type_name = ''
     CHARACTER(LEN=64) :: spec_type_name = ''
     CHARACTER(LEN=64) :: dust_type_name = ''
     CHARACTER(LEN=250) :: sps_home = ''
     CHARACTER(LEN=250) :: data_home = ''
     CHARACTER(LEN=250) :: output_home = ''
      TYPE(fsps_setup_cache_t), POINTER :: setup_cache => NULL()
     TYPE(fsps_context_state_t) :: state
     REAL(SP) :: om0_val = 0.0
     REAL(SP) :: ol0_val = 0.0
     REAL(SP) :: H0_val = 0.0
     REAL(SP) :: tiny_logt_val = 0.0
     REAL(SP) :: imf_upper_limit_val = 0.0
     REAL(SP) :: imf_lower_limit_val = 0.0
     REAL(SP) :: logt_wmb_hot_val = 0.0
     REAL(SP) :: nebular_smooth_init_val = 0.0
     INTEGER :: imf_type_val = 0
     INTEGER :: tpagb_norm_type_val = 0
     INTEGER :: pzcon_val = 0
     INTEGER :: interpolation_type_val = 0
     INTEGER :: add_agb_dust_model_val = 0
     INTEGER :: dust_type_val = 0
     INTEGER :: add_dust_emission_val = 0
     INTEGER :: compute_vega_mags_val = 0
     INTEGER :: vactoair_flag_val = 0
     INTEGER :: add_agn_dust_val = 0
     INTEGER :: use_wr_spectra_val = 0
     INTEGER :: add_neb_emission_val = 0
     INTEGER :: add_neb_continuum_val = 0
     INTEGER :: cloudy_dust_val = 0
     INTEGER :: add_igm_absorption_val = 0
     INTEGER :: nebemlineinspec_val = 0
     INTEGER :: add_xrb_emission_val = 0
     INTEGER :: add_stellar_remnants_val = 0
     INTEGER :: smooth_velocity_val = 0
     INTEGER :: smooth_lsf_val = 0
     INTEGER :: smoothspec_fast_val = 0
     INTEGER :: redshift_colors_val = 0
     INTEGER :: compute_light_ages_val = 0
     INTEGER :: use_isoc_mdot_val = 0
     INTEGER :: setup_nebular_gaussians_val = 0
     TYPE(PARAMS) :: pset
  END TYPE fsps_context_t

CONTAINS

  SUBROUTINE fsps_context_state_destroy(state)
    TYPE(fsps_context_state_t), INTENT(INOUT) :: state

    IF (ASSOCIATED(state%indexdefined)) NULLIFY(state%indexdefined)
    IF (ASSOCIATED(state%wgdust)) NULLIFY(state%wgdust)
    IF (ASSOCIATED(state%g03smcextn)) NULLIFY(state%g03smcextn)
    IF (ASSOCIATED(state%bands)) NULLIFY(state%bands)
    IF (ASSOCIATED(state%magsun)) NULLIFY(state%magsun)
    IF (ASSOCIATED(state%magvega)) NULLIFY(state%magvega)
    IF (ASSOCIATED(state%filter_leff)) NULLIFY(state%filter_leff)
    IF (ASSOCIATED(state%vega_spec)) NULLIFY(state%vega_spec)
    IF (ASSOCIATED(state%sun_spec)) NULLIFY(state%sun_spec)
    IF (ASSOCIATED(state%spec_lambda)) NULLIFY(state%spec_lambda)
    IF (ASSOCIATED(state%spec_nu)) NULLIFY(state%spec_nu)
    IF (ASSOCIATED(state%spec_res)) NULLIFY(state%spec_res)
    IF (ASSOCIATED(state%speclib)) NULLIFY(state%speclib)
    IF (ASSOCIATED(state%wmb_spec)) NULLIFY(state%wmb_spec)
    IF (ASSOCIATED(state%agb_spec_o)) NULLIFY(state%agb_spec_o)
    IF (ASSOCIATED(state%agb_logt_o)) NULLIFY(state%agb_logt_o)
    IF (ASSOCIATED(state%agb_spec_c)) NULLIFY(state%agb_spec_c)
    IF (ASSOCIATED(state%agb_logt_c)) NULLIFY(state%agb_logt_c)
    IF (ASSOCIATED(state%agb_spec_car)) NULLIFY(state%agb_spec_car)
    IF (ASSOCIATED(state%pagb_spec)) NULLIFY(state%pagb_spec)
    IF (ASSOCIATED(state%wrn_spec)) NULLIFY(state%wrn_spec)
    IF (ASSOCIATED(state%wrc_spec)) NULLIFY(state%wrc_spec)
    IF (ASSOCIATED(state%qpaharr)) NULLIFY(state%qpaharr)
    IF (ASSOCIATED(state%uminarr)) NULLIFY(state%uminarr)
    IF (ASSOCIATED(state%lambda_dustem)) NULLIFY(state%lambda_dustem)
    IF (ASSOCIATED(state%dustem_dustem)) NULLIFY(state%dustem_dustem)
    IF (ASSOCIATED(state%dustem2_dustem)) NULLIFY(state%dustem2_dustem)
    IF (ASSOCIATED(state%flux_dagb)) NULLIFY(state%flux_dagb)
    IF (ASSOCIATED(state%nebem_cont)) NULLIFY(state%nebem_cont)
    IF (ASSOCIATED(state%xnebem_cont)) NULLIFY(state%xnebem_cont)
    IF (ASSOCIATED(state%neb_res_min)) NULLIFY(state%neb_res_min)
    IF (ASSOCIATED(state%gaussnebarr)) NULLIFY(state%gaussnebarr)
    IF (ASSOCIATED(state%agndust_spec)) NULLIFY(state%agndust_spec)
    IF (ASSOCIATED(state%mact_isoc)) NULLIFY(state%mact_isoc)
    IF (ASSOCIATED(state%logl_isoc)) NULLIFY(state%logl_isoc)
    IF (ASSOCIATED(state%logt_isoc)) NULLIFY(state%logt_isoc)
    IF (ASSOCIATED(state%logg_isoc)) NULLIFY(state%logg_isoc)
    IF (ASSOCIATED(state%ffco_isoc)) NULLIFY(state%ffco_isoc)
    IF (ASSOCIATED(state%phase_isoc)) NULLIFY(state%phase_isoc)
    IF (ASSOCIATED(state%mini_isoc)) NULLIFY(state%mini_isoc)
    IF (ASSOCIATED(state%lmdot_isoc)) NULLIFY(state%lmdot_isoc)
    IF (ASSOCIATED(state%nmass_isoc)) NULLIFY(state%nmass_isoc)
    IF (ASSOCIATED(state%timestep_isoc)) NULLIFY(state%timestep_isoc)
    IF (ASSOCIATED(state%zlegend)) NULLIFY(state%zlegend)
    IF (ASSOCIATED(state%zlegendinit)) NULLIFY(state%zlegendinit)
    IF (ALLOCATED(state%spec_ssp_zz)) DEALLOCATE(state%spec_ssp_zz)
    IF (ALLOCATED(state%mass_ssp_zz)) DEALLOCATE(state%mass_ssp_zz)
    IF (ALLOCATED(state%lbol_ssp_zz)) DEALLOCATE(state%lbol_ssp_zz)
    IF (ASSOCIATED(state%time_full)) NULLIFY(state%time_full)
    IF (ALLOCATED(state%weight_ssp)) DEALLOCATE(state%weight_ssp)
    IF (ALLOCATED(state%spec_young)) DEALLOCATE(state%spec_young)
    IF (ALLOCATED(state%spec_old)) DEALLOCATE(state%spec_old)
    IF (ASSOCIATED(state%bpass_spec_ssp)) NULLIFY(state%bpass_spec_ssp)
    IF (ASSOCIATED(state%bpass_mass_ssp)) NULLIFY(state%bpass_mass_ssp)
    IF (ASSOCIATED(state%lam_xrb)) NULLIFY(state%lam_xrb)
    IF (ASSOCIATED(state%spec_xrb)) NULLIFY(state%spec_xrb)
    IF (ASSOCIATED(state%ages_xrb)) NULLIFY(state%ages_xrb)
    IF (ASSOCIATED(state%zmet_xrb)) NULLIFY(state%zmet_xrb)
    IF (ALLOCATED(state%lsfinfo%lsf)) DEALLOCATE(state%lsfinfo%lsf)
    IF (ALLOCATED(state%powell_data%mags)) DEALLOCATE(state%powell_data%mags)
    IF (ALLOCATED(state%powell_data%magerr)) DEALLOCATE(state%powell_data%magerr)
    IF (ALLOCATED(state%powell_data%spec)) DEALLOCATE(state%powell_data%spec)
    IF (ALLOCATED(state%powell_data%specerr)) DEALLOCATE(state%powell_data%specerr)
    IF (ALLOCATED(state%sedfit_data%mags)) DEALLOCATE(state%sedfit_data%mags)
    IF (ALLOCATED(state%sedfit_data%magerr)) DEALLOCATE(state%sedfit_data%magerr)
    IF (ALLOCATED(state%sedfit_data%spec)) DEALLOCATE(state%sedfit_data%spec)
    IF (ALLOCATED(state%sedfit_data%specerr)) DEALLOCATE(state%sedfit_data%specerr)

    state%nt = 0
    state%nz = 0
    state%nspec = 0
    state%nzinit = 0
    state%nbands = 0
    state%nindx = 0
    state%ntfull = 0
    state%nspec_xrb = 0
    state%nt_xrb = 0
    state%nz_xrb = 0
    state%check_sps_setup = 0
    state%tuniv = 0.0
    state%whlam5000 = 0
    state%whlylim = 0
    state%zsol = 0.0
    state%zsol_spec = 0.0
    state%zpow2 = 1.0
    state%ntabsfh = 0
    state%lsfinfo%minlam = 0.0
    state%lsfinfo%maxlam = 0.0
    state%str_dustem = 'DL07'
    state%imf_alpha = 1.3
    state%imf_vdmc = 0.08
    state%imf_mdave = 0.5
    state%n_user_imf = 0
    state%imf_user_alpha = 0.0
    state%salp_ind = 2.35
    state%imf_lower_limit = 0.08
    state%imf_upper_limit = 120.0
    state%imf_lower_bound = 0.0
    state%mlim_bh = 40.0
    state%mlim_ns = 8.5
    state%alt_filter_file = ''
    state%powell_data%zred = 0.0
    state%powell_data%logsmass = 0.0
    state%sedfit_data%zred = 0.0
    state%sedfit_data%logsmass = 0.0
  END SUBROUTINE fsps_context_state_destroy

END MODULE FSPS_CONTEXT_TYPES
