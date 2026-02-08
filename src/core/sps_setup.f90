MODULE sps_setup_utils

   USE fsps_context_types, ONLY: fsps_context_t

  INTERFACE
       SUBROUTINE SPS_SETUP(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
          USE fsps_context_types, ONLY: fsps_context_t
          INTEGER, INTENT(in) :: zin
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          CHARACTER(LEN=*), INTENT(in), OPTIONAL :: isoc_type_in
          CHARACTER(LEN=*), INTENT(in), OPTIONAL :: spec_type_in
          CHARACTER(LEN=*), INTENT(in), OPTIONAL :: dust_type_in
       END SUBROUTINE SPS_SETUP
  END INTERFACE

  PUBLIC :: SPS_SETUP, fsps_resolve_paths, sps_takedown

CONTAINS

  LOGICAL FUNCTION fsps_data_exists(path)
    CHARACTER(LEN=*), INTENT(IN) :: path
    LOGICAL :: ok

    ok = .FALSE.
    INQUIRE(FILE=TRIM(path)//'/allfilters.dat', EXIST=ok)
    IF (.NOT. ok) INQUIRE(FILE=TRIM(path)//'/FILTER_LIST', EXIST=ok)
    fsps_data_exists = ok
  END FUNCTION fsps_data_exists

  SUBROUTINE fsps_resolve_paths(ctx)
    USE fsps_context_types, ONLY: fsps_context_t
    CHARACTER(250) :: env, candidate
    LOGICAL :: system_prefix
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx

    ctx%data_home = ''
    ctx%output_home = ''

    CALL getenv('SPS_HOME', ctx%sps_home)
    IF (LEN_TRIM(ctx%sps_home) > 0) THEN
        candidate = TRIM(ctx%sps_home)
        IF (.NOT. fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = ''
    ENDIF

    CALL getenv('FSPS_DATA_HOME', env)
    IF (LEN_TRIM(env) > 0) THEN
       candidate = TRIM(env)
       IF (fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = candidate
    ENDIF

    IF (LEN_TRIM(ctx%sps_home) == 0) THEN
       CALL getenv('XDG_DATA_HOME', env)
       IF (LEN_TRIM(env) > 0) THEN
          candidate = TRIM(env)//'/fsps'
          IF (fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = candidate
       ENDIF
    ENDIF

    IF (LEN_TRIM(ctx%sps_home) == 0) THEN
       CALL getenv('HOME', env)
       IF (LEN_TRIM(env) > 0) THEN
          candidate = TRIM(env)//'/.local/share/fsps'
          IF (fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = candidate
       ENDIF
    ENDIF

    IF (LEN_TRIM(ctx%sps_home) == 0) THEN
       candidate = '/usr/share/fsps'
       IF (fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = candidate
    ENDIF

    IF (LEN_TRIM(ctx%sps_home) == 0) THEN
       candidate = '/usr/local/share/fsps'
       IF (fsps_data_exists(TRIM(candidate)//'/data')) ctx%sps_home = candidate
    ENDIF

    IF (LEN_TRIM(ctx%sps_home) == 0) THEN
       WRITE(*,*) 'SPS_SETUP ERROR: FSPS data path not found. Set FSPS_DATA_HOME or SPS_HOME.'
       STOP
    ENDIF

    ctx%data_home = TRIM(ctx%sps_home)//'/data'

    CALL getenv('FSPS_OUTPUT_HOME', env)
    IF (LEN_TRIM(env) > 0) THEN
       ctx%output_home = TRIM(env)
    ELSE
       system_prefix = .FALSE.
       IF (LEN_TRIM(ctx%sps_home) >= 5) THEN
          IF (ctx%sps_home(1:5) == '/usr/') system_prefix = .TRUE.
       ENDIF
       IF (LEN_TRIM(ctx%sps_home) >= 10) THEN
          IF (ctx%sps_home(1:10) == '/usr/local') system_prefix = .TRUE.
       ENDIF
       IF (system_prefix) THEN
          CALL getenv('HOME', env)
          IF (LEN_TRIM(env) > 0) THEN
             ctx%output_home = TRIM(env)//'/.local/share/fsps'
          ELSE
             ctx%output_home = '.'
          ENDIF
       ELSE
          ctx%output_home = TRIM(ctx%sps_home)
       ENDIF
    ENDIF
  END SUBROUTINE fsps_resolve_paths

  SUBROUTINE sps_takedown(ctx)
    USE fsps_cache, ONLY: fsps_cache_release_setup
    USE fsps_context_types, ONLY: fsps_context_t, fsps_context_state_destroy
    IMPLICIT NONE
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx

    IF (ASSOCIATED(ctx%setup_cache)) THEN
       CALL fsps_cache_release_setup(ctx%setup_cache)
       NULLIFY(ctx%setup_cache)
    ENDIF

    CALL fsps_context_state_destroy(ctx%state)
  END SUBROUTINE sps_takedown

END MODULE sps_setup_utils

SUBROUTINE SPS_SETUP(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)

  !read in isochrones and spectral libraries for all metallicities.
  !read in band-pass info and the spectrum for Vega.
  !Arrays are stored in a common block defined in sps_vars.f90

  !If zin=-1 then all metallicities are read in, otherwise only the
  !metallicity corresponding to zin in the look-up table zlegend.dat
  !is read.  Specifying only the metallicity of interest results
  !in a much faster setup.

   USE fsps_precision, ONLY: WP
   USE fsps_constants, ONLY: VERBOSE, TIME_RES_INCR, NM, NDIM_LOGT, NDIM_LOGG, &
            N_AGB_O, N_AGB_C, N_AGB_CAR, NDIM_PAGB, NDIM_WR, NDIM_WMB_LOGT, NDIM_WMB_LOGG, &
            NTAU_DAGB, NTEFF_DAGB, NEMLINE, NLAM_NEBCONT, NEBNZ, NEBNAGE, NEBNIP, &
            NAGNDUST, NAGNDUST_SPEC, SAFE_FLOOR, C_LIGHT, PI, L_SOL
   USE fsps_cache, ONLY: fsps_setup_cache_t, fsps_cache_get_setup
   USE fsps_context_types, ONLY: fsps_context_t
      USE sps_setup_utils, ONLY: sps_takedown, fsps_resolve_paths
      USE fsps_io, ONLY: load_zlegend_file, load_wavelength_grid, load_spectral_resolution, &
                     read_isochrone_database, read_spectral_binary, read_bpass_data, &
                     load_dust_emission_table, load_nebular_grid, load_filter_definitions, &
                     load_standard_sed, load_index_definitions, load_wr_spectra, &
                     load_agb_spectra, load_post_agb_spectra, load_attenuation_curves, &
                     load_wmbasic_spectra, load_lsf_data, load_agn_dust_models, &
                     apply_legacy_filter_norm
   USE fsps_cosmology, ONLY: get_universe_age, get_luminosity_distance
   USE fsps_interpolation, ONLY: find_interval, interpolate_linear
   USE fsps_integration, ONLY: integrate_trapezoid_array
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  IMPLICIT NONE
  TYPE(fsps_context_t), INTENT(INOUT) :: ctx
  INTEGER, INTENT(in) :: zin
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: isoc_type_in
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: spec_type_in
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: dust_type_in

   INTEGER :: stat=1,i,j,jj,k,i1,i2
   INTEGER :: z,zmin,zmax,nlam
   REAL(WP) :: d1,x,a,zero=0.0_wp,d,dz,dlam
  CHARACTER(LEN=512) :: cache_key
   CHARACTER(LEN=250) :: SPS_HOME
  LOGICAL :: cache_new
  TYPE(fsps_setup_cache_t), POINTER :: cache
  
  CHARACTER(5), ALLOCATABLE :: zz_str_xrb(:)
   REAL(WP), ALLOCATABLE :: tspec(:)
   REAL(WP), DIMENSION(10000) :: lambda_dagb=0.,fluxin_dagb=0.
   REAL(WP), DIMENSION(NLAM_NEBCONT) :: readlambneb=0.,readcontneb=0.
   REAL(WP), DIMENSION(NAGNDUST_SPEC)           :: agndust_lam=0.
   REAL(WP), DIMENSION(NAGNDUST_SPEC,NAGNDUST)  :: agndust_specinit=0.
   REAL(WP), ALLOCATABLE :: speclibinit(:,:,:,:)
   REAL(WP), ALLOCATABLE :: tspec_xrb(:)

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  ASSOCIATE( &
     zsol => ctx%state%zsol, zsol_spec => ctx%state%zsol_spec, &
     isoc_type => ctx%state%isoc_type, spec_type => ctx%state%spec_type, &
     nt => ctx%state%nt, nz => ctx%state%nz, nspec => ctx%state%nspec, &
     nzinit => ctx%state%nzinit, nbands => ctx%state%nbands, &
     nindx => ctx%state%nindx, ntfull => ctx%state%ntfull, &
     nspec_xrb => ctx%state%nspec_xrb, nt_xrb => ctx%state%nt_xrb, &
     nz_xrb => ctx%state%nz_xrb, check_sps_setup => ctx%state%check_sps_setup, &
     tuniv => ctx%state%tuniv, whlam5000 => ctx%state%whlam5000, &
     whlylim => ctx%state%whlylim, zpow2 => ctx%state%zpow2, &
     mwdindex => ctx%state%mwdindex, cosmospl => ctx%state%cosmospl, &
     ntabsfh => ctx%state%ntabsfh, sfh_tab => ctx%state%sfh_tab, &
     imf_alpha => ctx%state%imf_alpha, imf_vdmc => ctx%state%imf_vdmc, &
     imf_mdave => ctx%state%imf_mdave, n_user_imf => ctx%state%n_user_imf, &
     imf_user_alpha => ctx%state%imf_user_alpha, salp_ind => ctx%state%salp_ind, &
     imf_lower_limit => ctx%state%imf_lower_limit, &
     imf_upper_limit => ctx%state%imf_upper_limit, &
     imf_lower_bound => ctx%state%imf_lower_bound, mlim_bh => ctx%state%mlim_bh, &
     mlim_ns => ctx%state%mlim_ns, alt_filter_file => ctx%state%alt_filter_file, &
     indexdefined => ctx%state%indexdefined, wgdust => ctx%state%wgdust, &
     g03smcextn => ctx%state%g03smcextn, bands => ctx%state%bands, &
     magsun => ctx%state%magsun, magvega => ctx%state%magvega, &
     filter_leff => ctx%state%filter_leff, vega_spec => ctx%state%vega_spec, &
     sun_spec => ctx%state%sun_spec, spec_lambda => ctx%state%spec_lambda, &
     spec_nu => ctx%state%spec_nu, spec_res => ctx%state%spec_res, &
     speclib_logt => ctx%state%speclib_logt, speclib_logg => ctx%state%speclib_logg, &
     speclib => ctx%state%speclib, wmb_logt => ctx%state%wmb_logt, &
     wmb_logg => ctx%state%wmb_logg, wmb_spec => ctx%state%wmb_spec, &
     agb_spec_o => ctx%state%agb_spec_o, agb_logt_o => ctx%state%agb_logt_o, &
     agb_spec_c => ctx%state%agb_spec_c, agb_logt_c => ctx%state%agb_logt_c, &
     agb_logt_car => ctx%state%agb_logt_car, agb_spec_car => ctx%state%agb_spec_car, &
     pagb_spec => ctx%state%pagb_spec, pagb_logt => ctx%state%pagb_logt, &
     wrn_spec => ctx%state%wrn_spec, wrc_spec => ctx%state%wrc_spec, &
     wrn_logt => ctx%state%wrn_logt, wrc_logt => ctx%state%wrc_logt, &
     ndim_dustem => ctx%state%ndim_dustem, numin_dustem => ctx%state%numin_dustem, &
     nqpah_dustem => ctx%state%nqpah_dustem, str_dustem => ctx%state%str_dustem, &
     qpaharr => ctx%state%qpaharr, uminarr => ctx%state%uminarr, &
     lambda_dustem => ctx%state%lambda_dustem, dustem_dustem => ctx%state%dustem_dustem, &
     dustem2_dustem => ctx%state%dustem2_dustem, flux_dagb => ctx%state%flux_dagb, &
     tau1_dagb => ctx%state%tau1_dagb, teff_dagb => ctx%state%teff_dagb, &
     nebem_line_pos => ctx%state%nebem_line_pos, nebem_line => ctx%state%nebem_line, &
     xnebem_line => ctx%state%xnebem_line, nebem_cont => ctx%state%nebem_cont, &
     xnebem_cont => ctx%state%xnebem_cont, nebem_logz => ctx%state%nebem_logz, &
     nebem_age => ctx%state%nebem_age, nebem_logu => ctx%state%nebem_logu, &
     neb_res_min => ctx%state%neb_res_min, gaussnebarr => ctx%state%gaussnebarr, &
     agndust_tau => ctx%state%agndust_tau, agndust_spec => ctx%state%agndust_spec, &
     mact_isoc => ctx%state%mact_isoc, logl_isoc => ctx%state%logl_isoc, &
     logt_isoc => ctx%state%logt_isoc, logg_isoc => ctx%state%logg_isoc, &
     ffco_isoc => ctx%state%ffco_isoc, phase_isoc => ctx%state%phase_isoc, &
     mini_isoc => ctx%state%mini_isoc, lmdot_isoc => ctx%state%lmdot_isoc, &
     nmass_isoc => ctx%state%nmass_isoc, timestep_isoc => ctx%state%timestep_isoc, &
     zlegend => ctx%state%zlegend, zlegendinit => ctx%state%zlegendinit, &
     spec_ssp_zz => ctx%state%spec_ssp_zz, mass_ssp_zz => ctx%state%mass_ssp_zz, &
     lbol_ssp_zz => ctx%state%lbol_ssp_zz, time_full => ctx%state%time_full, &
     weight_ssp => ctx%state%weight_ssp, spec_young => ctx%state%spec_young, &
     spec_old => ctx%state%spec_old, bpass_spec_ssp => ctx%state%bpass_spec_ssp, &
     bpass_mass_ssp => ctx%state%bpass_mass_ssp, lam_xrb => ctx%state%lam_xrb, &
   spec_xrb => ctx%state%spec_xrb, ages_xrb => ctx%state%ages_xrb, &
   zmet_xrb => ctx%state%zmet_xrb, lsfinfo => ctx%state%lsfinfo, &
   cloudy_dust => ctx%cloudy_dust_val, smooth_velocity => ctx%smooth_velocity_val, &
   smooth_lsf => ctx%smooth_lsf_val, setup_nebular_gaussians => ctx%setup_nebular_gaussians_val, &
   add_neb_emission => ctx%add_neb_emission_val, add_neb_continuum => ctx%add_neb_continuum_val, &
   add_dust_emission => ctx%add_dust_emission_val, add_agn_dust => ctx%add_agn_dust_val, &
   add_xrb_emission => ctx%add_xrb_emission_val, add_agb_dust_model => ctx%add_agb_dust_model_val, &
   use_wr_spectra => ctx%use_wr_spectra_val, &
   powell_data => ctx%state%powell_data, sedfit_data => ctx%state%sedfit_data )

  CALL SPS_TAKEDOWN(ctx)

  IF (VERBOSE.EQ.1) THEN
     WRITE(*,*)
     WRITE(*,*) '    Setting up SPS...'
  ENDIF

  ! Initialize Library Variables
  IF (PRESENT(isoc_type_in)) THEN
     IF (LEN_TRIM(isoc_type_in) > 0) THEN
        isoc_type = isoc_type_in
     ELSE
        isoc_type = 'mist'
     END IF
  ELSE
     isoc_type = 'mist'
  END IF
  
  IF (PRESENT(spec_type_in)) THEN
     IF (LEN_TRIM(spec_type_in) > 0) THEN
        spec_type = spec_type_in
     ELSE
        spec_type = 'miles'
     END IF
  ELSE
     spec_type = 'miles'
  END IF

  IF (TRIM(isoc_type) == 'bpss' .AND. TRIM(spec_type) /= 'bpass') THEN
     WRITE(*,*) 'Notice: BPASS isochrones selected; forcing spec_type="bpass"'
     spec_type = 'bpass'
  END IF

  ! Set isochrone dimensions
  SELECT CASE (isoc_type)
  CASE ('mist')
     zsol = 0.0142
     nt=107
     nz=12
  CASE ('pdva')
     zsol = 0.019
     nt=94
     nz=22
  CASE ('prsc')
     zsol = 0.01524
     nt=93
     nz=15
  CASE ('bsti')
     zsol = 0.020
     nt=94
     nz=10
  CASE ('gnva')
     zsol = 0.020
     nt=51
     nz=5
  CASE ('bpss')
     zsol = 0.020
     nt=43
     nz=12
  CASE DEFAULT
      WRITE(*,*) 'SPS_SETUP ERROR: Unknown isoc_type: ', isoc_type
      STOP
  END SELECT
  
  ! Set spectral library dimensions
  SELECT CASE (spec_type)
  CASE ('miles')
     zsol_spec = 0.019
     nzinit=5
     nspec=5994
  CASE ('basel')
     zsol_spec = 0.020
     nzinit=6
     nspec=1963
  CASE ('bpass')
     zsol_spec = 0.020
     nzinit=1
     nspec=15000
  CASE DEFAULT
     ! Check for C3K
     IF (spec_type(1:3).EQ.'c3k') THEN
        zsol_spec = 0.0134
        nzinit=11
        nspec=11149
     ELSE
        WRITE(*,*) 'SPS_SETUP ERROR: Unknown spec_type: ', spec_type
        STOP
     END IF
  END SELECT

  ! Set other dimensions
  nbands = 159
  nindx = 30
  ntfull = TIME_RES_INCR * nt
  
  ! Set XRB dimensions
  nspec_xrb=15000
  nt_xrb=10
  nz_xrb=11

  ! Set Dust Emission Model
  IF (PRESENT(dust_type_in)) THEN
     IF (dust_type_in == 'themis' .OR. dust_type_in == 'THEMIS') THEN
        str_dustem = 'THEMIS'
        ndim_dustem = 576
        numin_dustem = 37
        nqpah_dustem = 11
     ELSE
        ! Default to DL07
        str_dustem = 'DL07'
        ndim_dustem = 1001
        numin_dustem = 22
        nqpah_dustem = 7
     END IF
  ELSE
     str_dustem = 'DL07'
     ndim_dustem = 1001
     numin_dustem = 22
     nqpah_dustem = 7
  END IF

  IF (zin.GT.nz) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zin GT nz', zin,nz
     STOP
  ENDIF

  ! Resolve paths and attach or create shared caches
   CALL fsps_resolve_paths(ctx)
   SPS_HOME = ctx%sps_home

  WRITE(cache_key, '(A,"|",A,"|",A,"|",A,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",A)') &
      TRIM(ctx%sps_home), TRIM(isoc_type), TRIM(spec_type), TRIM(str_dustem), &
       zin, smooth_velocity, setup_nebular_gaussians, add_neb_emission, add_neb_continuum, &
       add_dust_emission, add_agn_dust, add_xrb_emission, add_agb_dust_model, TRIM(alt_filter_file)

  CALL fsps_cache_get_setup(cache_key, ctx%setup_cache, cache_new)
  cache => ctx%setup_cache

  IF (cache_new) THEN
     cache%isoc_type = TRIM(isoc_type)
     cache%spec_type = TRIM(spec_type)
     cache%dust_type = TRIM(str_dustem)
     cache%alt_filter_file = alt_filter_file
     cache%nz = nz
     cache%nt = nt
     cache%nspec = nspec
     cache%nzinit = nzinit
     cache%nbands = nbands
     cache%nindx = nindx
     cache%ntfull = ntfull
     cache%nspec_xrb = nspec_xrb
     cache%nt_xrb = nt_xrb
     cache%nz_xrb = nz_xrb
     cache%smooth_velocity = smooth_velocity
     cache%setup_nebular_gaussians = setup_nebular_gaussians
     cache%add_neb_emission = add_neb_emission
     cache%add_neb_continuum = add_neb_continuum
     cache%add_dust_emission = add_dust_emission
     cache%add_agn_dust = add_agn_dust
     cache%add_xrb_emission = add_xrb_emission
     cache%add_agb_dust_model = add_agb_dust_model
     cache%use_wr_spectra = use_wr_spectra

     ALLOCATE(cache%indexdefined(7,nindx))
     ALLOCATE(cache%wgdust(nspec,18,6,2))
     ALLOCATE(cache%g03smcextn(nspec))
     ALLOCATE(cache%bands(nspec,nbands))
     ALLOCATE(cache%magsun(nbands),cache%magvega(nbands),cache%filter_leff(nbands))
     ALLOCATE(cache%vega_spec(nspec),cache%sun_spec(nspec))
     ALLOCATE(cache%spec_lambda(nspec),cache%spec_nu(nspec))
     ALLOCATE(cache%spec_res(nspec))
     ALLOCATE(cache%speclib(nspec,nz,NDIM_LOGT,NDIM_LOGG))
     ALLOCATE(cache%wmb_spec(nspec,nz,NDIM_WMB_LOGT,NDIM_WMB_LOGG))
     ALLOCATE(cache%agb_spec_o(nspec,N_AGB_O))
     ALLOCATE(cache%agb_logt_o(nz,N_AGB_O))
     ALLOCATE(cache%agb_spec_c(nspec,N_AGB_C))
     ALLOCATE(cache%agb_logt_c(N_AGB_C))
     ALLOCATE(cache%agb_spec_car(nspec,N_AGB_CAR))
     ALLOCATE(cache%pagb_spec(nspec,NDIM_PAGB,2))
     ALLOCATE(cache%wrn_spec(nspec,NDIM_WR,nz),cache%wrc_spec(nspec,NDIM_WR,nz))

     ! Dust models
     ALLOCATE(cache%qpaharr(nqpah_dustem))
     ALLOCATE(cache%uminarr(numin_dustem))
     ALLOCATE(cache%lambda_dustem(ndim_dustem))
     ALLOCATE(cache%dustem_dustem(ndim_dustem,numin_dustem*2))
     ALLOCATE(cache%dustem2_dustem(nspec,nqpah_dustem,numin_dustem*2))

     ALLOCATE(cache%flux_dagb(nspec,2,NTEFF_DAGB,NTAU_DAGB))
     ALLOCATE(cache%nebem_cont(nspec,NEBNZ,NEBNAGE,NEBNIP),cache%xnebem_cont(nspec,NEBNZ,NEBNAGE,NEBNIP))
     ALLOCATE(cache%neb_res_min(nspec))
     ALLOCATE(cache%gaussnebarr(nspec,NEMLINE))
     ALLOCATE(cache%agndust_spec(nspec,NAGNDUST))
     ALLOCATE(cache%mact_isoc(nz,nt,NM),cache%logl_isoc(nz,nt,NM),cache%logt_isoc(nz,nt,NM),cache%logg_isoc(nz,nt,NM))
     ALLOCATE(cache%ffco_isoc(nz,nt,NM),cache%phase_isoc(nz,nt,NM),cache%mini_isoc(nz,nt,NM),cache%lmdot_isoc(nz,nt,NM))
     ALLOCATE(cache%nmass_isoc(nz,nt))
     ALLOCATE(cache%timestep_isoc(nz,nt))
     ALLOCATE(cache%zlegend(nz))
     ALLOCATE(cache%zlegendinit(nzinit))
     ALLOCATE(cache%time_full(ntfull))
     ALLOCATE(cache%bpass_spec_ssp(nspec,nt,nz))
     ALLOCATE(cache%bpass_mass_ssp(nt,nz))
     ALLOCATE(cache%lam_xrb(nspec_xrb))
     ALLOCATE(cache%spec_xrb(nspec,nt_xrb,nz_xrb))
     ALLOCATE(cache%ages_xrb(nt_xrb))
     ALLOCATE(cache%zmet_xrb(nz_xrb))
  ENDIF

   ctx%state%indexdefined => cache%indexdefined
   ctx%state%wgdust => cache%wgdust
   ctx%state%g03smcextn => cache%g03smcextn
   ctx%state%bands => cache%bands
   ctx%state%magsun => cache%magsun
   ctx%state%magvega => cache%magvega
   ctx%state%filter_leff => cache%filter_leff
   ctx%state%vega_spec => cache%vega_spec
   ctx%state%sun_spec => cache%sun_spec
   ctx%state%spec_lambda => cache%spec_lambda
   ctx%state%spec_nu => cache%spec_nu
   ctx%state%spec_res => cache%spec_res
   ctx%state%speclib => cache%speclib
   ctx%state%wmb_spec => cache%wmb_spec
   ctx%state%agb_spec_o => cache%agb_spec_o
   ctx%state%agb_logt_o => cache%agb_logt_o
   ctx%state%agb_spec_c => cache%agb_spec_c
   ctx%state%agb_logt_c => cache%agb_logt_c
   ctx%state%agb_spec_car => cache%agb_spec_car
   ctx%state%pagb_spec => cache%pagb_spec
   ctx%state%wrn_spec => cache%wrn_spec
   ctx%state%wrc_spec => cache%wrc_spec
   ctx%state%qpaharr => cache%qpaharr
   ctx%state%uminarr => cache%uminarr
   ctx%state%lambda_dustem => cache%lambda_dustem
   ctx%state%dustem_dustem => cache%dustem_dustem
   ctx%state%dustem2_dustem => cache%dustem2_dustem
   ctx%state%flux_dagb => cache%flux_dagb
   ctx%state%nebem_cont => cache%nebem_cont
   ctx%state%xnebem_cont => cache%xnebem_cont
   ctx%state%neb_res_min => cache%neb_res_min
   ctx%state%gaussnebarr => cache%gaussnebarr
   ctx%state%agndust_spec => cache%agndust_spec
   ctx%state%mact_isoc => cache%mact_isoc
   ctx%state%logl_isoc => cache%logl_isoc
   ctx%state%logt_isoc => cache%logt_isoc
   ctx%state%logg_isoc => cache%logg_isoc
   ctx%state%ffco_isoc => cache%ffco_isoc
   ctx%state%phase_isoc => cache%phase_isoc
   ctx%state%mini_isoc => cache%mini_isoc
   ctx%state%lmdot_isoc => cache%lmdot_isoc
   ctx%state%nmass_isoc => cache%nmass_isoc
   ctx%state%timestep_isoc => cache%timestep_isoc
   ctx%state%zlegend => cache%zlegend
   ctx%state%zlegendinit => cache%zlegendinit
   ctx%state%bpass_spec_ssp => cache%bpass_spec_ssp
   ctx%state%bpass_mass_ssp => cache%bpass_mass_ssp
   ctx%state%lam_xrb => cache%lam_xrb
   ctx%state%spec_xrb => cache%spec_xrb
   ctx%state%ages_xrb => cache%ages_xrb
   ctx%state%zmet_xrb => cache%zmet_xrb
   ctx%state%time_full => cache%time_full

  ! Allocate per-context arrays
  ALLOCATE(ctx%state%spec_ssp_zz(nspec,ntfull,nz))
  ALLOCATE(ctx%state%mass_ssp_zz(ntfull,nz),ctx%state%lbol_ssp_zz(ntfull,nz))
  ALLOCATE(ctx%state%weight_ssp(ntfull,nz))
  ALLOCATE(ctx%state%spec_young(nspec),ctx%state%spec_old(nspec))
  ALLOCATE(ctx%state%lsfinfo%lsf(nspec))
  
  ! Allocate local arrays
  ALLOCATE(tspec(nspec))
  ALLOCATE(tspec_xrb(nspec_xrb))
  ALLOCATE(zz_str_xrb(nz_xrb))

  ASSOCIATE( &
     indexdefined => ctx%state%indexdefined, wgdust => ctx%state%wgdust, &
     g03smcextn => ctx%state%g03smcextn, bands => ctx%state%bands, &
     magsun => ctx%state%magsun, magvega => ctx%state%magvega, &
     filter_leff => ctx%state%filter_leff, vega_spec => ctx%state%vega_spec, &
     sun_spec => ctx%state%sun_spec, spec_lambda => ctx%state%spec_lambda, &
     spec_nu => ctx%state%spec_nu, spec_res => ctx%state%spec_res, &
     speclib => ctx%state%speclib, wmb_spec => ctx%state%wmb_spec, &
     agb_spec_o => ctx%state%agb_spec_o, agb_logt_o => ctx%state%agb_logt_o, &
     agb_spec_c => ctx%state%agb_spec_c, agb_logt_c => ctx%state%agb_logt_c, &
     agb_spec_car => ctx%state%agb_spec_car, pagb_spec => ctx%state%pagb_spec, &
     wrn_spec => ctx%state%wrn_spec, wrc_spec => ctx%state%wrc_spec, &
     qpaharr => ctx%state%qpaharr, uminarr => ctx%state%uminarr, &
     lambda_dustem => ctx%state%lambda_dustem, dustem_dustem => ctx%state%dustem_dustem, &
     dustem2_dustem => ctx%state%dustem2_dustem, flux_dagb => ctx%state%flux_dagb, &
     nebem_cont => ctx%state%nebem_cont, xnebem_cont => ctx%state%xnebem_cont, &
     neb_res_min => ctx%state%neb_res_min, gaussnebarr => ctx%state%gaussnebarr, &
     agndust_spec => ctx%state%agndust_spec, mact_isoc => ctx%state%mact_isoc, &
     logl_isoc => ctx%state%logl_isoc, logt_isoc => ctx%state%logt_isoc, &
     logg_isoc => ctx%state%logg_isoc, ffco_isoc => ctx%state%ffco_isoc, &
     phase_isoc => ctx%state%phase_isoc, mini_isoc => ctx%state%mini_isoc, &
     lmdot_isoc => ctx%state%lmdot_isoc, nmass_isoc => ctx%state%nmass_isoc, &
     timestep_isoc => ctx%state%timestep_isoc, zlegend => ctx%state%zlegend, &
     zlegendinit => ctx%state%zlegendinit, spec_ssp_zz => ctx%state%spec_ssp_zz, &
     mass_ssp_zz => ctx%state%mass_ssp_zz, lbol_ssp_zz => ctx%state%lbol_ssp_zz, &
     time_full => ctx%state%time_full, weight_ssp => ctx%state%weight_ssp, &
     spec_young => ctx%state%spec_young, spec_old => ctx%state%spec_old, &
     bpass_spec_ssp => ctx%state%bpass_spec_ssp, bpass_mass_ssp => ctx%state%bpass_mass_ssp, &
     lam_xrb => ctx%state%lam_xrb, spec_xrb => ctx%state%spec_xrb, &
     ages_xrb => ctx%state%ages_xrb, zmet_xrb => ctx%state%zmet_xrb )

  ! Initialize new arrays to 0.0 or default values
  IF (cache_new) THEN
     bands = 0.0
     magsun = 0.0
     magvega = 0.0
     filter_leff = 0.0
     vega_spec = 0.0
     sun_spec = 0.0
     spec_lambda = 0.0
     spec_nu = 0.0
     spec_res = 0.0
     speclib = 0.0
     wmb_spec = 0.0
     agb_spec_o = 0.0
     agb_logt_o = 0.0
     agb_spec_c = 0.0
     agb_logt_c = 0.0
     agb_spec_car = 0.0
     pagb_spec = 0.0
     wrn_spec = 0.0
     wrc_spec = 0.0

     ! Dust initialization
     dustem2_dustem = 0.0

     IF (TRIM(str_dustem) == 'THEMIS') THEN
        qpaharr = (/0.02,0.06,0.10,0.14,0.17,0.20,0.24,0.28,0.32,0.36,0.40/)/2.2*100
        uminarr = (/0.1,0.12,0.15,0.17,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,1.0,&
          1.2,1.5,1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,&
          12.0, 15.0, 17.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 80.0/)
     ELSE
        ! DL07
        qpaharr = (/0.47,1.12,1.77,2.50,3.19,3.90,4.58/)
        uminarr = (/0.1,0.15,0.2,0.3,0.4,0.5,0.7,0.8,1.0,1.2,1.5,2.0,&
          2.5,3.0,4.0,5.0,7.0,8.0,12.0,15.0,20.0,25.0/)
     END IF

     lambda_dustem = 0.0
     dustem_dustem = 0.0

     flux_dagb = 0.0
     nebem_cont = 0.0
     xnebem_cont = 0.0
     neb_res_min = 0.0
     gaussnebarr = 0.0
     agndust_spec = 0.0
     mact_isoc = 0.0
     logl_isoc = 0.0
     logt_isoc = 0.0
     logg_isoc = 0.0
     ffco_isoc = 0.0
     phase_isoc = 0.0
     mini_isoc = 0.0
     lmdot_isoc = 0.0
     nmass_isoc = 0
     timestep_isoc = 0.0
     zlegend = -99.0
     zlegendinit = -99.0
     time_full = 0.0
     bpass_spec_ssp = 0.0
     bpass_mass_ssp = 0.0
     lam_xrb = 0.0
     spec_xrb = 0.0
     ages_xrb = 0.0
     zmet_xrb = 0.0
  ENDIF

  spec_ssp_zz = 0.0
  mass_ssp_zz = 0.0
  lbol_ssp_zz = 0.0
  weight_ssp = 0.0
  spec_young = 0.0
  spec_old = 0.0
  lsfinfo%lsf = 0.0
  
  tspec = 0.0
  tspec_xrb = 0.0
  zz_str_xrb = ''
  IF (cache_new) THEN
     indexdefined = 0.0
     wgdust = 0.0
     g03smcextn = 0.0
  ENDIF
  mwdindex = 0
  sfh_tab = 0.0
  ntabsfh = 0

  !clean out all the common block arrays
  IF (cache_new) THEN
     mini_isoc     = 0.
     mact_isoc     = 0.
     logl_isoc     = 0.
     logt_isoc     = 0.
     logg_isoc     = 0.
     ffco_isoc     = 0.
     phase_isoc    = 0.
     nmass_isoc    = 0
     timestep_isoc = 0.
     spec_lambda   = 0.
     vega_spec     = 0.
     sun_spec      = 0.
     speclib       = 0.
     speclib_logg  = 0.
     speclib_logt  = 0.
     agb_spec_o    = 0.
     agb_logt_o    = 0.
     agb_spec_c    = 0.
     agb_logt_c    = 0.
  ENDIF


  !----------------------------------------------------------------!
  !--------------Confirm that variables are properly set-----------!
  !----------------------------------------------------------------!

  IF (cache_new) THEN

  !----------------------------------------------------------------!
  !----------------Read in metallicity values----------------------!
  !----------------------------------------------------------------!

  !units are simply metal fraction by mass (e.g. Z=0.0190 for Zsun)
  call load_zlegend_file(ctx, isoc_type, .false.)

  IF (zin.LE.0) THEN
     zmin = 1
     zmax = nz
  ELSE
     zmin = zin
     zmax = zin
  ENDIF

  !----------------------------------------------------------------!
  !---------------------Read in BPASS SSPs-------------------------!
  !----------------------------------------------------------------!

  IF (isoc_type.EQ.'bpss') THEN

     IF (TIME_RES_INCR.NE.1) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: cannot have TIME_RES_INCR>1 w/ BPASS models'
        STOP
     ENDIF

     call read_bpass_data(ctx)

  ENDIF

  !----------------------------------------------------------------!
  !-----------------Read in spectral libraries---------------------!
  !----------------------------------------------------------------!

  IF (isoc_type.NE.'bpss') THEN

  !read in wavelength array, spectral resolution, and spectral metallicity grid
  call load_wavelength_grid(ctx, spec_type)
  call load_spectral_resolution(ctx, spec_type)
  call load_zlegend_file(ctx, spec_type, .true.)

  !read in primary logg and logt arrays
  !NB: these are the same for all spectral libraries
   OPEN(91,FILE=TRIM(SPS_HOME)//'/data/spectra/BaSeL3.1/basel_logt.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,NDIM_LOGT
     READ(91,*) speclib_logt(i)
  ENDDO
  CLOSE(91)
   OPEN(91,FILE=TRIM(SPS_HOME)//'/data/spectra/BaSeL3.1/basel_logg.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,NDIM_LOGG
     READ(91,*) speclib_logg(i)
  ENDDO
  CLOSE(91)

   ALLOCATE(speclibinit(nspec,nzinit,NDIM_LOGT,NDIM_LOGG))
   speclibinit = 0.0

  !read in each metallicity
  DO z=1,nzinit
     call read_spectral_binary(ctx, spec_type, z, speclibinit(:,z,:,:))
  ENDDO

  !interpolate the input spectral library to the isochrone grid
  !notice that we're interpolating at fixed Z/Zsol even in cases
  !where the isochrones and spectra might have different Zsol. This might
  !in fact be the best thing to do.  Either way, its not ideal.
  DO z=1,nz

   i1 = MIN(MAX(find_interval(LOG10(zlegendinit/zsol_spec),&
      LOG10(zlegend(z)/zsol)),1),nzinit-1)
     dz = (LOG10(zlegend(z)/zsol)-LOG10(zlegendinit(i1)/zsol_spec)) / &
          (LOG10(zlegendinit(i1+1)/zsol_spec)-LOG10(zlegendinit(i1)/zsol_spec))
     dz = MIN(MAX(dz,0.0),1.0) !no extrapolation!

   speclib(:,z,:,:) = REAL((1-dz)*LOG10(speclibinit(:,i1,:,:)+SAFE_FLOOR) + &
      dz*LOG10(speclibinit(:,i1+1,:,:)+SAFE_FLOOR), KIND(speclib))
     speclib(:,z,:,:) = 10**speclib(:,z,:,:)

  ENDDO

  DEALLOCATE(speclibinit)

  !--------------Read WMBasic Grid from JJ Eldridge----------------;
  call load_wmbasic_spectra(ctx)

  !-----------Read in TP-AGB Library from Lancon & Wood------------;
  call load_agb_spectra(ctx)

  !------------read in post-AGB spectra from Rauch 2003------------;
  call load_post_agb_spectra(ctx)

  !--------------read in WR spectra from Smith et al.--------------;
  call load_wr_spectra(ctx)

  !----------------------------------------------------------------!
  !--------------------Read in isochrones--------------------------!
  !----------------------------------------------------------------!

  !read in all metallicities
  DO z=zmin,zmax
     call read_isochrone_database(ctx, isoc_type, z)
  ENDDO

  !this is necessary because Geneva does not extend below 1.0 Msun
  !see imf_weight.f90 for details
  IF (isoc_type.EQ.'gnva') THEN
     imf_lower_bound = MINVAL(mini_isoc(zmin,1,1:nmass_isoc(zmin,1)))*0.99
  ELSE
     imf_lower_bound = imf_lower_limit
  ENDIF

  ENDIF

  !----------------------------------------------------------------!
  !--------Read in dust emission spectra from Draine & Li----------!
  !----------------------------------------------------------------!

  call load_dust_emission_table(ctx, str_dustem)

  !----------------------------------------------------------------!
  !-------------Read in circumstellar AGB dust models--------------!
  !----------------------------------------------------------------!

  !O-rich spectra
   OPEN(99,FILE=TRIM(SPS_HOME)//'/data/dust/dusty/Orich_dusty.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: error opening dusty models'
     STOP
  ENDIF

  !number of wavelength points in the AGB grid
  READ(99,*) nlam
  !read in the wavelength grid
  READ(99,*,IOSTAT=stat) lambda_dagb(1:nlam)

  lambda_dagb(1:nlam) = lambda_dagb(1:nlam)

  DO i=1,NTEFF_DAGB
     DO j=1,NTAU_DAGB
        READ(99,*,IOSTAT=stat) teff_dagb(1,i), tau1_dagb(1,j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
      jj = MAX(find_interval(spec_lambda,lambda_dagb(1)), 1)
      flux_dagb(jj:,1,i,j) = interpolate_linear(lambda_dagb(1:nlam),&
             fluxin_dagb(1:nlam),spec_lambda(jj:))
     ENDDO
  ENDDO
  CLOSE(99)

  !C-rich spectra
   OPEN(99,FILE=TRIM(SPS_HOME)//'/data/dust/dusty/Crich_dusty.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: error opening dusty models'
     STOP
  ENDIF

  !number of wavelength points in the AGB grid
  READ(99,*) nlam
  !read in the wavelength grid
  READ(99,*,IOSTAT=stat) lambda_dagb(1:nlam)

  DO i=1,NTEFF_DAGB
     DO j=1,NTAU_DAGB
        READ(99,*,IOSTAT=stat) teff_dagb(2,i), tau1_dagb(2,j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
      jj = MAX(find_interval(spec_lambda,lambda_dagb(1)), 1)
      flux_dagb(jj:,2,i,j) = interpolate_linear(lambda_dagb(1:nlam),&
             fluxin_dagb(1:nlam),spec_lambda(jj:))
     ENDDO
  ENDDO
  CLOSE(99)

  !----------------------------------------------------------------!
  !--------------------Set up AGN dust model-----------------------!
  !----------------------------------------------------------------!

  !models from Nenkova et al. 2008

   OPEN(99,FILE=TRIM(SPS_HOME)//'/data/dust/Nenkova08_y010_torusg_n10_q2.0.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: error opening AGN dust models'
     STOP
  ENDIF

  !burn the header
  DO i=1,3
     READ(99,*)
  ENDDO

  !read in the optical depths
  READ(99,*) agndust_tau

  DO i=1,NAGNDUST_SPEC
     READ(99,*) agndust_lam(i),agndust_specinit(i,:)
  ENDDO

   i1 = MAX(find_interval(spec_lambda,agndust_lam(1)), 1)
   i2 = MAX(find_interval(spec_lambda,agndust_lam(NAGNDUST_SPEC)), 1)
  DO i=1,NAGNDUST
   agndust_spec(i1:i2,i) = 10**interpolate_linear(LOG10(agndust_lam),&
          LOG10(agndust_specinit(:,i)+SAFE_FLOOR),LOG10(spec_lambda(i1:i2)))-SAFE_FLOOR
  ENDDO


  !----------------------------------------------------------------!
  !----------------Set up nebular emission arrays------------------!
  !----------------------------------------------------------------!

  IF (isoc_type.EQ.'mist'.OR.isoc_type.EQ.'pdva'.OR.&
     isoc_type.EQ.'prsc'.OR.isoc_type.EQ.'bpss') THEN

     call load_nebular_grid(ctx, isoc_type, cloudy_dust.EQ.1)

     !define the minimum resolution of the emission lines
     !based on the resolution of the spectral library
     DO i=1,NEMLINE
      j = MIN(MAX(find_interval(spec_lambda,nebem_line_pos(i)),1),nspec-1)
        neb_res_min(i) = spec_lambda(j+1)-spec_lambda(j)
     ENDDO

     !set up a "master" array of normalized Gaussians
     !this makes the code much faster
     IF (setup_nebular_gaussians.EQ.1) THEN
        DO i=1,NEMLINE
           IF (smooth_velocity.EQ.1) THEN
              !smoothing variable is km/s
              dlam = nebem_line_pos(i)*ctx%nebular_smooth_init_val/C_LIGHT*1E13
           ELSE
              !smoothing variable is A
              dlam = ctx%nebular_smooth_init_val
           ENDIF
           !broaden the line to at least the resolution element
           !of the spectrum (x2).
           dlam = MAX(dlam,neb_res_min(i)*2)
           gaussnebarr(:,i) = 1/SQRT(2*PI)/dlam*&
                EXP(-(spec_lambda-nebem_line_pos(i))**2/2/dlam**2)  / &
                C_LIGHT*nebem_line_pos(i)**2
        ENDDO
     ENDIF

  ENDIF

  !----------------------------------------------------------------!
  !------------------Set up X-ray nebular --------------------!
  !----------------------------------------------------------------!

  IF (isoc_type.EQ.'bpss') THEN
      !read in nebular continuum arrays.  Units are Lsun/Hz/Q
      IF (cloudy_dust.EQ.1) THEN
         OPEN(99,FILE=TRIM(SPS_HOME)//'/data/nebular/ZAU_WX_WD_'//TRIM(isoc_type)//'.cont',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ELSE
         OPEN(99,FILE=TRIM(SPS_HOME)//'/data/nebular/ZAU_WX_ND_'//TRIM(isoc_type)//'.cont',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ENDIF
      IF (stat.NE.0) THEN
         WRITE(*,*) 'SPS_SETUP ERROR: nebular cont file cannot be opened. '
         STOP
      ENDIF
      !burn the header
      READ(99,*)
      !read the wavelength array
      READ(99,*) readlambneb
      DO i=1,NEBNZ
         DO j=1,NEBNAGE
            DO k=1,NEBNIP
               READ(99,*,iostat=stat) nebem_logz(i),nebem_age(j),nebem_logu(k)
               READ(99,*,iostat=stat) readcontneb
               !interpolate onto the main wavelength grid
               !some values in the table are 0.0, set a floor of 1E-95
                  xnebem_cont(:,i,j,k) = interpolate_linear(readlambneb,&
                     LOG10(readcontneb+10**(-95.d0)),spec_lambda)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(99)

      !read in nebular emission line luminosities.  Units are Lsun/Q
      IF (cloudy_dust.EQ.1) THEN
         OPEN(99,FILE=TRIM(SPS_HOME)//'/data/nebular/ZAU_WX_WD_'//TRIM(isoc_type)//'.lines',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ELSE
         OPEN(99,FILE=TRIM(SPS_HOME)//'/data/nebular/ZAU_WX_ND_'//TRIM(isoc_type)//'.lines',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ENDIF
      IF (stat.NE.0) THEN
         WRITE(*,*) 'SPS_SETUP ERROR: nebular line file cannot be opened. Only available for Padova or MIST isochrones.'
         STOP
      ENDIF
      !burn the header
      READ(99,*)
      !read the wavelength array
      READ(99,*) nebem_line_pos
      DO i=1,NEBNZ
         DO j=1,NEBNAGE
            DO k=1,NEBNIP
               READ(99,*,iostat=stat) nebem_logz(i),nebem_age(j),nebem_logu(k)
               READ(99,*,iostat=stat) xnebem_line(:,i,j,k)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(99)

      !convert the nebem_age array to log(age), and log the emission arrays
      nebem_age  = LOG10(nebem_age)
      xnebem_line = LOG10(xnebem_line+10**(-95.d0))
   ENDIF


  !----------------------------------------------------------------!
  !------------------Set up X-ray binary arrays--------------------!
  !----------------------------------------------------------------!

   OPEN(98,FILE=TRIM(SPS_HOME)//'/data/spectra/xrb/xsp.lambda',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: xsp.lambda cannot be opened'
     STOP
  ENDIF

  DO i=1,nspec_xrb
     READ(98,*) lam_xrb(i)
  ENDDO

  CLOSE(98)

  ages_xrb = LOG10((/1.0,2.0,3.0,4.0,5.0,8.0,10.0,12.6,16.0,20.0/))+6.0
  zmet_xrb = (/-1.3,-1.0,-0.8,-0.7,-0.5,-0.4,-0.3,-0.2,+0.0,+0.2,+0.3/)

  zz_str_xrb = (/'-1.30','-1.00','-0.80','-0.70','-0.50','-0.40','-0.30','-0.20','+0.00','+0.20','+0.30'/)

  DO j=1,nz_xrb
   OPEN(98,FILE=TRIM(SPS_HOME)//'/data/spectra/xrb/xsp_feh'//zz_str_xrb(j)&
          //'.spec',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: xsp_feh'//zz_str_xrb(j)//&
             '.spec cannot be opened'
        STOP
     ENDIF
     DO i=1,nt_xrb
        READ(98,*) tspec_xrb
        !interpolate to the main wavelength array
      spec_xrb(:,i,j) = MAX(interpolate_linear(lam_xrb,tspec_xrb,spec_lambda),SAFE_FLOOR)
     ENDDO
     CLOSE(98)
  ENDDO

  !convert to Lsun/Hz/Msun
  spec_xrb = spec_xrb * L_SOL

  !----------------------------------------------------------------!
  !-------------------Set up magnitude info------------------------!
  !----------------------------------------------------------------!

  !read in and set up band-pass filters
  IF (TRIM(alt_filter_file).EQ.'' .OR. TRIM(alt_filter_file).EQ.'allfilters.dat') THEN
     call load_filter_definitions(ctx)
  ELSE
     call load_filter_definitions(ctx, TRIM(alt_filter_file))
  ENDIF

  !read in Vega/Sun SEDs and compute mag zero-points
  call load_standard_sed(ctx)

  !put Sun magnitudes in the Vega system if keyword is set
  IF (ctx%compute_vega_mags_val.EQ.1) THEN
     DO i=1,nbands
        IF (magsun(i).NE.99.0) magsun(i) = (magsun(i)-magsun(1)) - (magvega(i)-magvega(1)) + magsun(1)
     ENDDO
  ENDIF

  ! Apply Legacy Normalization HERE (After Zero Points are calculated)
  ! only execute this loop for the standard filter list
  IF (TRIM(alt_filter_file).EQ.'' .OR. TRIM(alt_filter_file).EQ.'allfilters.dat') THEN
     call apply_legacy_filter_norm(ctx)
  ENDIF

  !compute the effective wavelength of each filter
  !NB: These are sometimes referred to as "pivot" wavelengths
  ! in the literature.  See Bessell & Murphy 2012 A.2.1 for details'
  DO i=1,nbands
   filter_leff(i) = integrate_trapezoid_array(spec_lambda,spec_lambda*bands(:,i))
     d = integrate_trapezoid_array(spec_lambda,bands(:,i)/spec_lambda)
     IF (ieee_is_nan(filter_leff(i)) .OR. ieee_is_nan(d) .OR. d.LE.SAFE_FLOOR) THEN
        filter_leff(i) = 0.0
     ELSE
        filter_leff(i) = filter_leff(i) / d
        filter_leff(i) = SQRT(filter_leff(i))
     ENDIF
  ENDDO

  !----------------------------------------------------------------!
  !---------------Set up extinction curve indices------------------!
  !----------------------------------------------------------------!

  !----------------------------------------------------------------!
  !----------Set up Witt & Gordon 2000 attenuation curves----------!
  !----------------------------------------------------------------!

  call load_attenuation_curves(ctx)

  call load_agn_dust_models(ctx)

  ENDIF

  IF (zin.LE.0) THEN
     zmin = 1
     zmax = nz
  ELSE
     zmin = zin
     zmax = zin
  ENDIF

  !these are the breakpoints for the CCM89 MW parameterization
  DO j=1,nspec
     x = 1E4/spec_lambda(j)
     IF (x.GT.12.) mwdindex(6)=j
     IF (x.GE.8.)  mwdindex(5)=j
     IF (x.GE.5.9) mwdindex(4)=j
     IF (x.GE.3.3) mwdindex(3)=j
     IF (x.GE.1.1) mwdindex(2)=j
     IF (x.GE.0.1) mwdindex(1)=j
  ENDDO

  !----------------------------------------------------------------!
  !--------------set up the redshift-age-DL relations--------------!
  !----------------------------------------------------------------!

  DO i=1,500
     a = (i-1)/499.*(1-1/1001.)+1/1001.
     cosmospl(i,1) = 1/a-1  !redshift
   cosmospl(i,2) = get_universe_age(ctx, cosmospl(i,1))   ! Tuniv in Gyr
   cosmospl(i,3) = get_luminosity_distance(ctx, cosmospl(i,1)) ! Lum Dist in pc
  ENDDO

  !set Tuniv
   tuniv = get_universe_age(ctx, zero)

  !----------------------------------------------------------------!
  !-----------------read in index definitions----------------------!
  !----------------------------------------------------------------!

  IF (cache_new) THEN
     call load_index_definitions(ctx)
  ENDIF

  !----------------------------------------------------------------!
  !-----------------set up expanded time array---------------------!
  !----------------------------------------------------------------!

  IF (cache_new) THEN
     IF (isoc_type.NE.'bpss') THEN

        DO i=1,ntfull
           IF (MOD(i-1,TIME_RES_INCR).EQ.0) THEN
              time_full(i) = timestep_isoc(zmin,(i-1)/TIME_RES_INCR+1)
           ELSE
              IF ((i-1)/TIME_RES_INCR+2.LT.nt) THEN
                 d1 = (timestep_isoc(zmin,(i-1)/TIME_RES_INCR+2)-&
                      timestep_isoc(zmin,(i-1)/TIME_RES_INCR+1))/TIME_RES_INCR
              ENDIF
              time_full(i) = timestep_isoc(zmin,(i-1)/TIME_RES_INCR+1)+d1
              time_full(i) = time_full(i-1)+d1
           ENDIF
        ENDDO

     ENDIF
  ENDIF

  !----------------------------------------------------------------!
  !------------------------set up the LSF--------------------------!
  !----------------------------------------------------------------!

  IF (smooth_lsf.EQ.1) THEN
     call load_lsf_data(ctx)
  ENDIF


  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

   whlam5000 = find_interval(spec_lambda,5000.d0)
   whlylim   = find_interval(spec_lambda,912.d0)
  !define the frequency array
  IF (cache_new) THEN
     spec_nu   = C_LIGHT / spec_lambda
  ENDIF

  !set flag indicating that sps_setup has been run, initializing
  !important common block vars/arrays
  check_sps_setup = 1

  IF (VERBOSE.EQ.1) THEN
     WRITE(*,*) '      ...done'
     WRITE(*,*)
  ENDIF

   END ASSOCIATE

   END ASSOCIATE

END SUBROUTINE SPS_SETUP
