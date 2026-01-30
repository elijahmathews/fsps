MODULE FSPS_CONTEXT
            USE fsps_types, ONLY: SP, PARAMS, COMPSPOUT
  USE fsps_context_types, ONLY: fsps_context_t, fsps_context_state_destroy
  USE sps_utils
  IMPLICIT NONE

   INTEGER, PARAMETER :: FSPS_ERR_UNKNOWN_INT_PARAM = 101
   INTEGER, PARAMETER :: FSPS_ERR_UNKNOWN_FLOAT_PARAM = 102
   INTEGER, PARAMETER :: FSPS_ERR_UNKNOWN_STRING_PARAM = 103

CONTAINS

  SUBROUTINE fsps_context_create(ctx)
      TYPE(fsps_context_t), INTENT(OUT) :: ctx
    ctx%initialized = .FALSE.
    ctx%zin = 0
      ctx%isoc_type_name = ''
      ctx%spec_type_name = ''
      ctx%dust_type_name = ''
      ctx%sps_home = ''
      ctx%data_home = ''
      ctx%output_home = ''
         ctx%om0_val = 0.27
         ctx%ol0_val = 0.73
         ctx%H0_val = 72.0
         ctx%tpagb_norm_type_val = 2
         ctx%pzcon_val = 0
         ctx%interpolation_type_val = 0
         ctx%tiny_logt_val = 0.0
         ctx%compute_light_ages_val = 0
         ctx%add_dust_emission_val = 1
         ctx%add_agn_dust_val = 1
         ctx%add_agb_dust_model_val = 1
         ctx%use_wr_spectra_val = 1
         ctx%logt_wmb_hot_val = 0.0
         ctx%add_neb_emission_val = 0
         ctx%add_neb_continuum_val = 1
         ctx%cloudy_dust_val = 0
         ctx%add_igm_absorption_val = 0
         ctx%add_xrb_emission_val = 0
         ctx%add_stellar_remnants_val = 1
         ctx%smoothspec_fast_val = 1
         ctx%smooth_velocity_val = 1
         ctx%smooth_lsf_val = 0
         ctx%dust_type_val = 0
         ctx%imf_type_val = 2
         ctx%compute_vega_mags_val = 0
         ctx%vactoair_flag_val = 0
         ctx%redshift_colors_val = 0
         ctx%use_isoc_mdot_val = 0
         ctx%setup_nebular_gaussians_val = 0
         ctx%nebular_smooth_init_val = 100.0
         ctx%nebemlineinspec_val = 1
         ctx%imf_lower_limit_val = 0.08
         ctx%imf_upper_limit_val = 120.0
  END SUBROUTINE fsps_context_create

  SUBROUTINE fsps_context_setup(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    INTEGER, INTENT(IN) :: zin
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: isoc_type_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: spec_type_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dust_type_in

    IF (PRESENT(isoc_type_in) .AND. PRESENT(spec_type_in) .AND. PRESENT(dust_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
       ctx%isoc_type_name = TRIM(isoc_type_in)
       ctx%spec_type_name = TRIM(spec_type_in)
       ctx%dust_type_name = TRIM(dust_type_in)
    ELSE IF (PRESENT(isoc_type_in) .AND. PRESENT(spec_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, isoc_type_in, spec_type_in)
       ctx%isoc_type_name = TRIM(isoc_type_in)
       ctx%spec_type_name = TRIM(spec_type_in)
    ELSE IF (PRESENT(isoc_type_in) .AND. PRESENT(dust_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, isoc_type_in, dust_type_in=dust_type_in)
       ctx%isoc_type_name = TRIM(isoc_type_in)
       ctx%dust_type_name = TRIM(dust_type_in)
    ELSE IF (PRESENT(spec_type_in) .AND. PRESENT(dust_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, spec_type_in=spec_type_in, dust_type_in=dust_type_in)
       ctx%spec_type_name = TRIM(spec_type_in)
       ctx%dust_type_name = TRIM(dust_type_in)
    ELSE IF (PRESENT(isoc_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, isoc_type_in)
       ctx%isoc_type_name = TRIM(isoc_type_in)
    ELSE IF (PRESENT(spec_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, spec_type_in=spec_type_in)
       ctx%spec_type_name = TRIM(spec_type_in)
    ELSE IF (PRESENT(dust_type_in)) THEN
      CALL SPS_SETUP(ctx, zin, dust_type_in=dust_type_in)
       ctx%dust_type_name = TRIM(dust_type_in)
    ELSE
      CALL SPS_SETUP(ctx, zin)
    END IF

      ctx%zin = zin
      ctx%initialized = .TRUE.
      IF (LEN_TRIM(ctx%isoc_type_name) == 0) THEN
         ctx%isoc_type_name = TRIM(ctx%state%isoc_type)
      END IF
      IF (LEN_TRIM(ctx%spec_type_name) == 0) THEN
         ctx%spec_type_name = TRIM(ctx%state%spec_type)
      END IF
      IF (LEN_TRIM(ctx%dust_type_name) == 0) THEN
         ctx%dust_type_name = TRIM(ctx%state%str_dustem)
      END IF
      CALL fsps_context_prepare_pset(ctx)
  END SUBROUTINE fsps_context_setup

  SUBROUTINE fsps_context_ensure_setup(ctx)
     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
     LOGICAL :: need_setup
     CHARACTER(LEN=64) :: isoc_name
     CHARACTER(LEN=64) :: spec_name
     CHARACTER(LEN=64) :: dust_name

     need_setup = .FALSE.
     isoc_name = TRIM(ctx%isoc_type_name)
     spec_name = TRIM(ctx%spec_type_name)
     dust_name = TRIM(ctx%dust_type_name)

     IF (LEN_TRIM(isoc_name) /= 0 .AND. TRIM(ctx%state%isoc_type) /= isoc_name) THEN
        need_setup = .TRUE.
     END IF

     IF (LEN_TRIM(spec_name) /= 0 .AND. TRIM(ctx%state%spec_type) /= spec_name) THEN
        need_setup = .TRUE.
     END IF

     IF (.NOT. need_setup) THEN
        IF (LEN_TRIM(isoc_name) == 0) THEN
           ctx%isoc_type_name = TRIM(ctx%state%isoc_type)
        END IF
        IF (LEN_TRIM(spec_name) == 0) THEN
           ctx%spec_type_name = TRIM(ctx%state%spec_type)
        END IF
        IF (LEN_TRIM(dust_name) == 0) THEN
           ctx%dust_type_name = TRIM(ctx%state%str_dustem)
        END IF
        RETURN
     END IF

     IF (LEN_TRIM(isoc_name) == 0 .AND. LEN_TRIM(spec_name) == 0 .AND. LEN_TRIM(dust_name) == 0) THEN
      CALL SPS_SETUP(ctx, ctx%zin)
     ELSE IF (LEN_TRIM(spec_name) == 0 .AND. LEN_TRIM(dust_name) == 0) THEN
      CALL SPS_SETUP(ctx, ctx%zin, isoc_name)
     ELSE IF (LEN_TRIM(dust_name) == 0) THEN
      CALL SPS_SETUP(ctx, ctx%zin, isoc_name, spec_name)
     ELSE IF (LEN_TRIM(spec_name) == 0) THEN
      CALL SPS_SETUP(ctx, ctx%zin, isoc_name, dust_type_in=dust_name)
     ELSE
      CALL SPS_SETUP(ctx, ctx%zin, isoc_name, spec_name, dust_name)
     END IF

   ctx%isoc_type_name = TRIM(ctx%state%isoc_type)
   ctx%spec_type_name = TRIM(ctx%state%spec_type)
   ctx%dust_type_name = TRIM(ctx%state%str_dustem)
  END SUBROUTINE fsps_context_ensure_setup

  SUBROUTINE fsps_context_destroy(ctx)
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    IF (ctx%initialized) THEN
          CALL SPS_TAKEDOWN(ctx)
       ctx%initialized = .FALSE.
    END IF
      CALL fsps_context_state_destroy(ctx%state)
  END SUBROUTINE fsps_context_destroy


  SUBROUTINE fsps_context_set_pset(ctx, pset_in)
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    TYPE(PARAMS), INTENT(IN) :: pset_in
    ctx%pset = pset_in
  END SUBROUTINE fsps_context_set_pset

  SUBROUTINE fsps_context_get_pset(ctx, pset_out)
    TYPE(fsps_context_t), INTENT(IN) :: ctx
    TYPE(PARAMS), INTENT(OUT) :: pset_out
    pset_out = ctx%pset
  END SUBROUTINE fsps_context_get_pset

  SUBROUTINE fsps_context_set_param_int(ctx, key, value, status)
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    CHARACTER(LEN=*), INTENT(IN) :: key
    INTEGER, INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: status

    status = 0
    SELECT CASE (TRIM(key))
    CASE ('sfh')
       ctx%pset%sfh = value
    CASE ('zmet')
       ctx%pset%zmet = value
    CASE ('wgp1')
       ctx%pset%wgp1 = value
    CASE ('wgp2')
       ctx%pset%wgp2 = value
    CASE ('wgp3')
       ctx%pset%wgp3 = value
    CASE ('evtype')
       ctx%pset%evtype = value
    CASE ('imf_type')
       ctx%imf_type_val = value
    CASE ('tpagb_norm_type')
       ctx%tpagb_norm_type_val = value
    CASE ('pzcon')
       ctx%pzcon_val = value
    CASE ('interpolation_type')
       ctx%interpolation_type_val = value
    CASE ('add_agb_dust_model')
       ctx%add_agb_dust_model_val = value
    CASE ('dust_type')
       ctx%dust_type_val = value
    CASE ('add_dust_emission')
       ctx%add_dust_emission_val = value
    CASE ('compute_vega_mags')
       ctx%compute_vega_mags_val = value
    CASE ('vactoair_flag')
       ctx%vactoair_flag_val = value
    CASE ('add_agn_dust')
       ctx%add_agn_dust_val = value
    CASE ('use_wr_spectra')
       ctx%use_wr_spectra_val = value
    CASE ('add_neb_emission')
       ctx%add_neb_emission_val = value
    CASE ('add_neb_continuum')
       ctx%add_neb_continuum_val = value
    CASE ('cloudy_dust')
       ctx%cloudy_dust_val = value
    CASE ('add_igm_absorption')
       ctx%add_igm_absorption_val = value
    CASE ('nebemlineinspec')
       ctx%nebemlineinspec_val = value
    CASE ('add_xrb_emission')
       ctx%add_xrb_emission_val = value
    CASE ('add_stellar_remnants')
       ctx%add_stellar_remnants_val = value
    CASE ('smooth_velocity')
       ctx%smooth_velocity_val = value
    CASE ('smooth_lsf')
       ctx%smooth_lsf_val = value
    CASE ('smoothspec_fast')
       ctx%smoothspec_fast_val = value
    CASE ('redshift_colors')
       ctx%redshift_colors_val = value
    CASE ('compute_light_ages')
       ctx%compute_light_ages_val = value
    CASE ('use_isoc_mdot')
       ctx%use_isoc_mdot_val = value
    CASE ('setup_nebular_gaussians')
       ctx%setup_nebular_gaussians_val = value
    CASE DEFAULT
       status = FSPS_ERR_UNKNOWN_INT_PARAM
    END SELECT
  END SUBROUTINE fsps_context_set_param_int

  SUBROUTINE fsps_context_set_param_float(ctx, key, value, status)
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    CHARACTER(LEN=*), INTENT(IN) :: key
    REAL(SP), INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: status

    status = 0
    SELECT CASE (TRIM(key))
    CASE ('om0')
       ctx%om0_val = value
    CASE ('ol0')
       ctx%ol0_val = value
    CASE ('H0')
       ctx%H0_val = value
    CASE ('tiny_logt')
       ctx%tiny_logt_val = value
    CASE ('imf_upper_limit')
       ctx%imf_upper_limit_val = value
    CASE ('imf_lower_limit')
       ctx%imf_lower_limit_val = value
    CASE ('logt_wmb_hot')
       ctx%logt_wmb_hot_val = value
    CASE ('nebular_smooth_init')
       ctx%nebular_smooth_init_val = value
    CASE ('imf1')
       ctx%pset%imf1 = value
    CASE ('imf2')
       ctx%pset%imf2 = value
    CASE ('imf3')
       ctx%pset%imf3 = value
    CASE ('vdmc')
       ctx%pset%vdmc = value
    CASE ('mdave')
       ctx%pset%mdave = value
    CASE ('dell')
       ctx%pset%dell = value
    CASE ('delt')
       ctx%pset%delt = value
    CASE ('sbss')
       ctx%pset%sbss = value
    CASE ('fbhb')
       ctx%pset%fbhb = value
    CASE ('pagb')
       ctx%pset%pagb = value
    CASE ('agb_dust')
       ctx%pset%agb_dust = value
    CASE ('redgb')
       ctx%pset%redgb = value
    CASE ('agb')
       ctx%pset%agb = value
    CASE ('masscut')
       ctx%pset%masscut = value
    CASE ('fcstar')
       ctx%pset%fcstar = value
    CASE ('frac_xrb')
       ctx%pset%frac_xrb = value
    CASE ('logzsol')
       ctx%pset%logzsol = value
    CASE ('tau')
       ctx%pset%tau = value
    CASE ('const')
       ctx%pset%const = value
    CASE ('tage')
       ctx%pset%tage = value
    CASE ('fburst')
       ctx%pset%fburst = value
    CASE ('tburst')
       ctx%pset%tburst = value
    CASE ('dust1')
       ctx%pset%dust1 = value
    CASE ('dust2')
       ctx%pset%dust2 = value
    CASE ('dust3')
       ctx%pset%dust3 = value
    CASE ('zred')
       ctx%pset%zred = value
    CASE ('pmetals')
       ctx%pset%pmetals = value
    CASE ('dust_clumps')
       ctx%pset%dust_clumps = value
    CASE ('frac_nodust')
       ctx%pset%frac_nodust = value
    CASE ('dust_index')
       ctx%pset%dust_index = value
    CASE ('dust_tesc')
       ctx%pset%dust_tesc = value
    CASE ('frac_obrun')
       ctx%pset%frac_obrun = value
    CASE ('uvb')
       ctx%pset%uvb = value
    CASE ('mwr')
       ctx%pset%mwr = value
    CASE ('dust1_index')
       ctx%pset%dust1_index = value
    CASE ('sf_start')
       ctx%pset%sf_start = value
    CASE ('sf_trunc')
       ctx%pset%sf_trunc = value
    CASE ('sf_slope')
       ctx%pset%sf_slope = value
    CASE ('duste_gamma')
       ctx%pset%duste_gamma = value
    CASE ('duste_umin')
       ctx%pset%duste_umin = value
    CASE ('duste_qpah')
       ctx%pset%duste_qpah = value
    CASE ('sigma_smooth')
       ctx%pset%sigma_smooth = value
    CASE ('min_wave_smooth')
       ctx%pset%min_wave_smooth = value
    CASE ('max_wave_smooth')
       ctx%pset%max_wave_smooth = value
    CASE ('gas_logu')
       ctx%pset%gas_logu = value
    CASE ('gas_logz')
       ctx%pset%gas_logz = value
    CASE ('igm_factor')
       ctx%pset%igm_factor = value
    CASE ('fagn')
       ctx%pset%fagn = value
    CASE ('agn_tau')
       ctx%pset%agn_tau = value
    CASE DEFAULT
       status = FSPS_ERR_UNKNOWN_FLOAT_PARAM
    END SELECT
  END SUBROUTINE fsps_context_set_param_float

  SUBROUTINE fsps_context_set_param_str(ctx, key, value, status)
    TYPE(fsps_context_t), INTENT(INOUT) :: ctx
    CHARACTER(LEN=*), INTENT(IN) :: key
    CHARACTER(LEN=*), INTENT(IN) :: value
    INTEGER, INTENT(OUT) :: status

    status = 0
    SELECT CASE (TRIM(key))
    CASE ('imf_filename')
       ctx%pset%imf_filename = TRIM(value)
    CASE ('sfh_filename')
       ctx%pset%sfh_filename = TRIM(value)
    CASE DEFAULT
       status = FSPS_ERR_UNKNOWN_STRING_PARAM
    END SELECT
  END SUBROUTINE fsps_context_set_param_str


   SUBROUTINE fsps_context_get_paths(ctx, sps_home_out, data_home_out, output_home_out)
      TYPE(fsps_context_t), INTENT(IN) :: ctx
      CHARACTER(LEN=*), INTENT(OUT) :: sps_home_out
      CHARACTER(LEN=*), INTENT(OUT) :: data_home_out
      CHARACTER(LEN=*), INTENT(OUT) :: output_home_out

      sps_home_out = TRIM(ctx%sps_home)
      data_home_out = TRIM(ctx%data_home)
      output_home_out = TRIM(ctx%output_home)
   END SUBROUTINE fsps_context_get_paths

   SUBROUTINE fsps_context_compute_ssp(ctx, mass_ssp, lbol_ssp, spec_ssp)
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(SP), DIMENSION(:), INTENT(out) :: mass_ssp, lbol_ssp
      REAL(SP), DIMENSION(:,:), INTENT(out) :: spec_ssp

        CALL fsps_context_ensure_setup(ctx)
      CALL fsps_context_prepare_pset(ctx)
      CALL SSP_GEN(ctx, ctx%pset, mass_ssp, lbol_ssp, spec_ssp)
   END SUBROUTINE fsps_context_compute_ssp

   SUBROUTINE fsps_context_compute_csp(ctx, write_compsp, nzin, outfile, mass_ssp, lbol_ssp, spec_ssp, ocompsp)
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      INTEGER, INTENT(IN) :: write_compsp, nzin
      CHARACTER(LEN=*), INTENT(IN) :: outfile
      REAL(SP), DIMENSION(:,:), INTENT(IN) :: mass_ssp, lbol_ssp
      REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: spec_ssp
      TYPE(COMPSPOUT), DIMENSION(:), INTENT(INOUT) :: ocompsp

        CALL fsps_context_ensure_setup(ctx)
      CALL fsps_context_prepare_pset(ctx)
      CALL COMPSP(ctx, write_compsp, nzin, outfile, mass_ssp, lbol_ssp, spec_ssp, ctx%pset, ocompsp)
   END SUBROUTINE fsps_context_compute_csp

   SUBROUTINE fsps_context_prepare_pset(ctx)
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx

      IF (.NOT. ALLOCATED(ctx%pset%mag_compute)) THEN
           ALLOCATE(ctx%pset%mag_compute(ctx%state%nbands))
          ctx%pset%mag_compute = 1
      END IF
      IF (.NOT. ALLOCATED(ctx%pset%ssp_gen_age)) THEN
           ALLOCATE(ctx%pset%ssp_gen_age(ctx%state%nt))
          ctx%pset%ssp_gen_age = 1
      END IF
   END SUBROUTINE fsps_context_prepare_pset

END MODULE FSPS_CONTEXT
