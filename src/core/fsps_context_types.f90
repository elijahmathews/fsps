module fsps_context_types
    !> @brief
    !> Core derived types that define FSPS contexts and shared state.
    !>
    !> @details
    !> Provides the `fsps_context_t` and `fsps_context_state_t` structures along
    !> with lifecycle helpers for releasing associated memory.

    use fsps_precision, only: WP
    use fsps_constants, only: &
        NDIM_LOGT, NDIM_LOGG, NDIM_WMB_LOGT, NDIM_WMB_LOGG, &
        N_AGB_CAR, NDIM_PAGB, NDIM_WR, NTAU_DAGB, NTEFF_DAGB, &
        NEMLINE, NEBNZ, NEBNAGE, NEBNIP, NAGNDUST, NTABMAX
    use fsps_types, only: PARAMS, COMPSPOUT, TLSF, OBSDAT
    use fsps_cache, only: fsps_setup_cache_t

    implicit none
    private

    public :: fsps_context_state_t
    public :: fsps_context_t
    public :: fsps_context_state_destroy

    ! ---------------------------------------------------------------------
    ! Module constants
    ! ---------------------------------------------------------------------

    type :: fsps_context_state_t
        real(WP) :: zsol = 0.0
        real(WP) :: zsol_spec = 0.0
        character(len=64) :: isoc_type = ''
        character(len=64) :: spec_type = ''
        integer :: nt = 0
        integer :: nz = 0
        integer :: nspec = 0
        integer :: nzinit = 0
        integer :: nbands = 0
        integer :: nindx = 0
        integer :: ntfull = 0
        integer :: nspec_xrb = 0
        integer :: nt_xrb = 0
        integer :: nz_xrb = 0
        integer :: check_sps_setup = 0
        real(WP) :: tuniv = 0.0
        integer :: whlam5000 = 0
        integer :: whlylim = 0
        real(WP) :: zpow2 = 1.0
        integer, dimension(6) :: mwdindex = 0
        real(WP), dimension(500, 3) :: cosmospl = 0.0
        integer :: ntabsfh = 0
        real(WP), dimension(3, NTABMAX) :: sfh_tab = 0.0
        real(WP), dimension(3) :: imf_alpha = 1.3
        real(WP) :: imf_vdmc = 0.08
        real(WP) :: imf_mdave = 0.5
        integer :: n_user_imf = 0
        real(WP), dimension(3, 100) :: imf_user_alpha = 0.0
        real(WP) :: salp_ind = 2.35
        real(WP) :: imf_lower_limit = 0.08
        real(WP) :: imf_upper_limit = 120.0
        real(WP) :: imf_lower_bound = 0.0
        real(WP) :: mlim_bh = 40.0
        real(WP) :: mlim_ns = 8.5
        character(len=30) :: alt_filter_file = ''
        real(WP), pointer :: indexdefined(:, :) => null()
        real(WP), pointer :: wgdust(:, :, :, :) => null()
        real(WP), pointer :: g03smcextn(:) => null()
        real(WP), pointer :: bands(:, :) => null()
        real(WP), pointer :: magsun(:) => null()
        real(WP), pointer :: magvega(:) => null()
        real(WP), pointer :: filter_leff(:) => null()
        real(WP), pointer :: vega_spec(:) => null()
        real(WP), pointer :: sun_spec(:) => null()
        real(WP), pointer :: spec_lambda(:) => null()
        real(WP), pointer :: spec_nu(:) => null()
        real(WP), pointer :: spec_res(:) => null()
        real(WP), dimension(NDIM_LOGT) :: speclib_logt = 0.0
        real(WP), dimension(NDIM_LOGG) :: speclib_logg = 0.0
        real(kind(1.0)), pointer :: speclib(:, :, :, :) => null()
        real(WP), dimension(NDIM_WMB_LOGT) :: wmb_logt = 0.0
        real(WP), dimension(NDIM_WMB_LOGG) :: wmb_logg = 0.0
        real(kind(1.0)), pointer :: wmb_spec(:, :, :, :) => null()
        real(WP), pointer :: agb_spec_o(:, :) => null()
        real(WP), pointer :: agb_logt_o(:, :) => null()
        real(WP), pointer :: agb_spec_c(:, :) => null()
        real(WP), pointer :: agb_logt_c(:) => null()
        real(WP), dimension(N_AGB_CAR) :: agb_logt_car = 0.0
        real(WP), pointer :: agb_spec_car(:, :) => null()
        real(WP), pointer :: pagb_spec(:, :, :) => null()
        real(WP), dimension(NDIM_PAGB) :: pagb_logt = 0.0
        real(WP), pointer :: wrn_spec(:, :, :) => null()
        real(WP), pointer :: wrc_spec(:, :, :) => null()
        real(WP), dimension(NDIM_WR) :: wrn_logt = 0.0
        real(WP), dimension(NDIM_WR) :: wrc_logt = 0.0
        integer :: ndim_dustem = 0
        integer :: numin_dustem = 0
        integer :: nqpah_dustem = 0
        character(len=6) :: str_dustem = 'DL07'
        real(WP), pointer :: qpaharr(:) => null()
        real(WP), pointer :: uminarr(:) => null()
        real(WP), pointer :: lambda_dustem(:) => null()
        real(WP), pointer :: dustem_dustem(:, :) => null()
        real(WP), pointer :: dustem2_dustem(:, :, :) => null()
        real(WP), pointer :: flux_dagb(:, :, :, :) => null()
        real(WP), dimension(2, NTAU_DAGB) :: tau1_dagb = 0.0
        real(WP), dimension(2, NTEFF_DAGB) :: teff_dagb = 0.0
        real(WP), dimension(NEMLINE) :: nebem_line_pos = 0.0
        real(WP), dimension(NEMLINE, NEBNZ, NEBNAGE, NEBNIP) :: nebem_line = 0.0
        real(WP), dimension(NEMLINE, NEBNZ, NEBNAGE, NEBNIP) :: xnebem_line = 0.0
        real(WP), pointer :: nebem_cont(:, :, :, :) => null()
        real(WP), pointer :: xnebem_cont(:, :, :, :) => null()
        real(WP), dimension(NEBNZ) :: nebem_logz = 0.0
        real(WP), dimension(NEBNAGE) :: nebem_age = 0.0
        real(WP), dimension(NEBNIP) :: nebem_logu = 0.0
        real(WP), pointer :: neb_res_min(:) => null()
        real(WP), pointer :: gaussnebarr(:, :) => null()
        real(WP), dimension(NAGNDUST) :: agndust_tau = 0.0
        real(WP), pointer :: agndust_spec(:, :) => null()
        real(WP), pointer :: mact_isoc(:, :, :) => null()
        real(WP), pointer :: logl_isoc(:, :, :) => null()
        real(WP), pointer :: logt_isoc(:, :, :) => null()
        real(WP), pointer :: logg_isoc(:, :, :) => null()
        real(WP), pointer :: ffco_isoc(:, :, :) => null()
        real(WP), pointer :: phase_isoc(:, :, :) => null()
        real(WP), pointer :: mini_isoc(:, :, :) => null()
        real(WP), pointer :: lmdot_isoc(:, :, :) => null()
        integer, pointer :: nmass_isoc(:, :) => null()
        real(WP), pointer :: timestep_isoc(:, :) => null()
        real(WP), pointer :: zlegend(:) => null()
        real(WP), pointer :: zlegendinit(:) => null()
        real(WP), allocatable :: spec_ssp_zz(:, :, :)
        real(WP), allocatable :: mass_ssp_zz(:, :)
        real(WP), allocatable :: lbol_ssp_zz(:, :)
        real(WP), pointer :: time_full(:) => null()
        real(WP), allocatable :: weight_ssp(:, :)
        real(WP), allocatable :: spec_young(:)
        real(WP), allocatable :: spec_old(:)
        real(WP), allocatable :: ssp_temp_grid(:, :)
        integer, allocatable :: ssp_active_idx(:)
        real(WP), allocatable :: ssp_active_w(:)
        integer :: ssp_temp_nspec = 0
        integer :: ssp_temp_nstars = 0
        real(WP), pointer :: bpass_spec_ssp(:, :, :) => null()
        real(WP), pointer :: bpass_mass_ssp(:, :) => null()
        real(WP), pointer :: lam_xrb(:) => null()
        real(WP), pointer :: spec_xrb(:, :, :) => null()
        real(WP), pointer :: ages_xrb(:) => null()
        real(WP), pointer :: zmet_xrb(:) => null()
        type(TLSF) :: lsfinfo
        type(OBSDAT) :: powell_data
        type(OBSDAT) :: sedfit_data

        ! --- Persistent CSP Workspace ---
        real(WP), allocatable :: ssp_basis_spec(:, :, :)
        real(WP), allocatable :: ssp_basis_mass(:, :)
        real(WP), allocatable :: ssp_basis_lbol(:, :)
        real(WP), allocatable :: csp_ssp_grid(:,:,:)
        real(WP), allocatable :: csp_emlin_grid(:,:,:)
        real(WP), allocatable :: csp_ssp_lum_linear(:,:)
        real(WP), allocatable :: csp_igm_transmission(:)
        real(WP), allocatable :: csp_spec_final(:)
        real(WP), allocatable :: csp_emlin_final(:)

        ! --- Replacements for csp_buffer_t ---
        real(WP), allocatable :: csp_weights(:,:)
        real(WP), allocatable :: csp_emlin_young(:)
        real(WP), allocatable :: csp_emlin_old(:)
    end type fsps_context_state_t

    type :: fsps_context_t
        logical :: initialized = .false.
        integer :: zin = 0
        character(len=64) :: isoc_type_name = ''
        character(len=64) :: spec_type_name = ''
        character(len=64) :: dust_type_name = ''
        character(len=250) :: sps_home = ''
        character(len=250) :: data_home = ''
        character(len=250) :: output_home = ''
        type(fsps_setup_cache_t), pointer :: setup_cache => null()
        type(fsps_context_state_t) :: state
        real(WP) :: om0_val = 0.0
        real(WP) :: ol0_val = 0.0
        real(WP) :: H0_val = 0.0
        real(WP) :: tiny_logt_val = 0.0
        real(WP) :: imf_upper_limit_val = 0.0
        real(WP) :: imf_lower_limit_val = 0.0
        real(WP) :: logt_wmb_hot_val = 0.0
        real(WP) :: nebular_smooth_init_val = 0.0
        integer :: imf_type_val = 0
        integer :: tpagb_norm_type_val = 0
        integer :: pzcon_val = 0
        integer :: interpolation_type_val = 0
        integer :: add_agb_dust_model_val = 0
        integer :: dust_type_val = 0
        integer :: add_dust_emission_val = 0
        integer :: compute_vega_mags_val = 0
        integer :: vactoair_flag_val = 0
        integer :: add_agn_dust_val = 0
        integer :: use_wr_spectra_val = 0
        integer :: add_neb_emission_val = 0
        integer :: add_neb_continuum_val = 0
        integer :: cloudy_dust_val = 0
        integer :: add_igm_absorption_val = 0
        integer :: nebemlineinspec_val = 0
        integer :: add_xrb_emission_val = 0
        integer :: add_stellar_remnants_val = 0
        integer :: smooth_velocity_val = 0
        integer :: smooth_lsf_val = 0
        integer :: smoothspec_fast_val = 0
        integer :: redshift_colors_val = 0
        integer :: compute_light_ages_val = 0
        integer :: use_isoc_mdot_val = 0
        integer :: setup_nebular_gaussians_val = 0
        type(PARAMS) :: pset
    end type fsps_context_t

contains

    !> @brief Release pointer members in a context state object.
    !> @param[inout] state State to detach from shared cache arrays.
    subroutine fsps_context_state_destroy(state)
        type(fsps_context_state_t), intent(inout) :: state

        if (associated(state%indexdefined)) nullify (state%indexdefined)
        if (associated(state%wgdust)) nullify (state%wgdust)
        if (associated(state%g03smcextn)) nullify (state%g03smcextn)
        if (associated(state%bands)) nullify (state%bands)
        if (associated(state%magsun)) nullify (state%magsun)
        if (associated(state%magvega)) nullify (state%magvega)
        if (associated(state%filter_leff)) nullify (state%filter_leff)
        if (associated(state%vega_spec)) nullify (state%vega_spec)
        if (associated(state%sun_spec)) nullify (state%sun_spec)
        if (associated(state%spec_lambda)) nullify (state%spec_lambda)
        if (associated(state%spec_nu)) nullify (state%spec_nu)
        if (associated(state%spec_res)) nullify (state%spec_res)
        if (associated(state%speclib)) nullify (state%speclib)
        if (associated(state%wmb_spec)) nullify (state%wmb_spec)
        if (associated(state%agb_spec_o)) nullify (state%agb_spec_o)
        if (associated(state%agb_logt_o)) nullify (state%agb_logt_o)
        if (associated(state%agb_spec_c)) nullify (state%agb_spec_c)
        if (associated(state%agb_logt_c)) nullify (state%agb_logt_c)
        if (associated(state%agb_spec_car)) nullify (state%agb_spec_car)
        if (associated(state%pagb_spec)) nullify (state%pagb_spec)
        if (associated(state%wrn_spec)) nullify (state%wrn_spec)
        if (associated(state%wrc_spec)) nullify (state%wrc_spec)
        if (associated(state%qpaharr)) nullify (state%qpaharr)
        if (associated(state%uminarr)) nullify (state%uminarr)
        if (associated(state%lambda_dustem)) nullify (state%lambda_dustem)
        if (associated(state%dustem_dustem)) nullify (state%dustem_dustem)
        if (associated(state%dustem2_dustem)) nullify (state%dustem2_dustem)
        if (associated(state%flux_dagb)) nullify (state%flux_dagb)
        if (associated(state%nebem_cont)) nullify (state%nebem_cont)
        if (associated(state%xnebem_cont)) nullify (state%xnebem_cont)
        if (associated(state%neb_res_min)) nullify (state%neb_res_min)
        if (associated(state%gaussnebarr)) nullify (state%gaussnebarr)
        if (associated(state%agndust_spec)) nullify (state%agndust_spec)
        if (associated(state%mact_isoc)) nullify (state%mact_isoc)
        if (associated(state%logl_isoc)) nullify (state%logl_isoc)
        if (associated(state%logt_isoc)) nullify (state%logt_isoc)
        if (associated(state%logg_isoc)) nullify (state%logg_isoc)
        if (associated(state%ffco_isoc)) nullify (state%ffco_isoc)
        if (associated(state%phase_isoc)) nullify (state%phase_isoc)
        if (associated(state%mini_isoc)) nullify (state%mini_isoc)
        if (associated(state%lmdot_isoc)) nullify (state%lmdot_isoc)
        if (associated(state%nmass_isoc)) nullify (state%nmass_isoc)
        if (associated(state%timestep_isoc)) nullify (state%timestep_isoc)
        if (associated(state%zlegend)) nullify (state%zlegend)
        if (associated(state%zlegendinit)) nullify (state%zlegendinit)
        if (allocated(state%spec_ssp_zz)) deallocate (state%spec_ssp_zz)
        if (allocated(state%mass_ssp_zz)) deallocate (state%mass_ssp_zz)
        if (allocated(state%lbol_ssp_zz)) deallocate (state%lbol_ssp_zz)
        if (associated(state%time_full)) nullify (state%time_full)
        if (allocated(state%weight_ssp)) deallocate (state%weight_ssp)
        if (allocated(state%spec_young)) deallocate (state%spec_young)
        if (allocated(state%spec_old)) deallocate (state%spec_old)
        if (allocated(state%ssp_temp_grid)) deallocate (state%ssp_temp_grid)
        if (allocated(state%ssp_active_idx)) deallocate (state%ssp_active_idx)
        if (allocated(state%ssp_active_w)) deallocate (state%ssp_active_w)
        if (associated(state%bpass_spec_ssp)) nullify (state%bpass_spec_ssp)
        if (associated(state%bpass_mass_ssp)) nullify (state%bpass_mass_ssp)
        if (associated(state%lam_xrb)) nullify (state%lam_xrb)
        if (associated(state%spec_xrb)) nullify (state%spec_xrb)
        if (associated(state%ages_xrb)) nullify (state%ages_xrb)
        if (associated(state%zmet_xrb)) nullify (state%zmet_xrb)
        if (allocated(state%lsfinfo%lsf)) deallocate (state%lsfinfo%lsf)
        if (allocated(state%powell_data%mags)) deallocate (state%powell_data%mags)
        if (allocated(state%powell_data%magerr)) deallocate (state%powell_data%magerr)
        if (allocated(state%powell_data%spec)) deallocate (state%powell_data%spec)
        if (allocated(state%powell_data%specerr)) deallocate (state%powell_data%specerr)
        if (allocated(state%sedfit_data%mags)) deallocate (state%sedfit_data%mags)
        if (allocated(state%sedfit_data%magerr)) deallocate (state%sedfit_data%magerr)
        if (allocated(state%sedfit_data%spec)) deallocate (state%sedfit_data%spec)
        if (allocated(state%sedfit_data%specerr)) deallocate (state%sedfit_data%specerr)

        ! --- Persistent CSP Workspace ---
        if (allocated(state%ssp_basis_spec)) deallocate (state%ssp_basis_spec)
        if (allocated(state%ssp_basis_mass)) deallocate (state%ssp_basis_mass)
        if (allocated(state%ssp_basis_lbol)) deallocate (state%ssp_basis_lbol)
        if (allocated(state%csp_ssp_grid)) deallocate (state%csp_ssp_grid)
        if (allocated(state%csp_emlin_grid)) deallocate (state%csp_emlin_grid)
        if (allocated(state%csp_ssp_lum_linear)) deallocate (state%csp_ssp_lum_linear)
        if (allocated(state%csp_igm_transmission)) deallocate (state%csp_igm_transmission)
        if (allocated(state%csp_spec_final)) deallocate (state%csp_spec_final)
        if (allocated(state%csp_emlin_final)) deallocate (state%csp_emlin_final)

        ! --- Replacements for csp_buffer_t ---
        if (allocated(state%csp_weights)) deallocate (state%csp_weights)
        if (allocated(state%csp_emlin_young)) deallocate (state%csp_emlin_young)
        if (allocated(state%csp_emlin_old)) deallocate (state%csp_emlin_old)

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
        state%ssp_temp_nspec = 0
        state%ssp_temp_nstars = 0
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
        
    end subroutine fsps_context_state_destroy

end module fsps_context_types
