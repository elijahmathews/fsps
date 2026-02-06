MODULE SPS_UTILS

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

  INTERFACE
     SUBROUTINE COMPSP(ctx, write_compsp, nzin, outfile, mass_ssp, &
          lbol_ssp, spec_ssp, pset, ocompsp)
       USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
       USE fsps_types, ONLY: PARAMS, COMPSPOUT
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       INTEGER, INTENT(in) :: write_compsp,nzin
      REAL(WP), INTENT(in), DIMENSION(:,:) :: lbol_ssp,mass_ssp
      REAL(WP), INTENT(in), DIMENSION(:,:,:) :: spec_ssp
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)   :: pset
       TYPE(COMPSPOUT), INTENT(inout), DIMENSION(:) :: ocompsp
     END SUBROUTINE COMPSP
  END INTERFACE

  INTERFACE
     SUBROUTINE CSP_GEN(ctx, mass_ssp, lbol_ssp, spec_ssp, pset, tage, nzin,&
                        mass_csp, lbol_csp, spec_csp, mdust_csp, emlin_ssp, emlin_csp)
       USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
      USE fsps_constants, ONLY: NEMLINE
       USE fsps_types, ONLY: PARAMS
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(WP), DIMENSION(:,:), INTENT(in) :: mass_ssp, lbol_ssp
      REAL(WP), DIMENSION(:,:,:), INTENT(in) :: spec_ssp
       TYPE(PARAMS), intent(in) :: pset
      REAL(WP), INTENT(in)  :: tage
       INTEGER, INTENT(IN) :: nzin
      REAL(WP), INTENT(out) :: mass_csp, lbol_csp, mdust_csp
      REAL(WP), INTENT(out), DIMENSION(:) :: spec_csp
      REAL(WP), DIMENSION(:,:,:), intent(in) :: emlin_ssp
      REAL(WP), DIMENSION(NEMLINE), intent(out) :: emlin_csp
     END SUBROUTINE CSP_GEN
  END INTERFACE

  INTERFACE
       SUBROUTINE WRITE_ISOCHRONE(ctx, outfile, pset)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          TYPE(PARAMS), INTENT(in) :: pset
          CHARACTER(100), INTENT(in)  :: outfile
     END SUBROUTINE WRITE_ISOCHRONE
  END INTERFACE

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

   SUBROUTINE SPS_TAKEDOWN(ctx)
      USE fsps_cache, ONLY: fsps_cache_release_setup
      USE fsps_context_types, ONLY: fsps_context_t, fsps_context_state_destroy
      IMPLICIT NONE
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx

      IF (ASSOCIATED(ctx%setup_cache)) THEN
         CALL fsps_cache_release_setup(ctx%setup_cache)
         NULLIFY(ctx%setup_cache)
      ENDIF

      CALL fsps_context_state_destroy(ctx%state)
   END SUBROUTINE SPS_TAKEDOWN

END MODULE SPS_UTILS
