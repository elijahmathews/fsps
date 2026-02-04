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
     SUBROUTINE ADD_AGB_DUST(ctx, weight, tspec, mact, logt, logl, logg, &
          zz,tco,lmdot)
       USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(WP), DIMENSION(:), INTENT(inout) :: tspec
      REAL(WP), INTENT(in)  :: weight,mact,logt,logl,logg,zz,tco,lmdot
     END SUBROUTINE ADD_AGB_DUST
  END INTERFACE


   INTERFACE
       SUBROUTINE ADD_DUST(ctx, pset, csp1, csp2, specdust, mdust, ncsp1, ncsp2, nebdust)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_constants, ONLY: NEMLINE
           USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          REAL(WP), INTENT(out) :: mdust
           REAL(WP), DIMENSION(:), INTENT(in) :: csp1,csp2
          TYPE(PARAMS), INTENT(in) :: pset
           REAL(WP), DIMENSION(:), INTENT(out) :: specdust
          REAL(WP), DIMENSION(NEMLINE), INTENT(in) :: ncsp1,ncsp2
          REAL(WP), DIMENSION(NEMLINE), INTENT(out) :: nebdust
       END SUBROUTINE ADD_DUST
   END INTERFACE


  INTERFACE
     FUNCTION AGN_DUST(ctx, lam, spec, pset, lbol_csp)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          REAL(WP), DIMENSION(:), INTENT(in) :: lam,spec
          REAL(WP), INTENT(in) :: lbol_csp
          TYPE(PARAMS), INTENT(in) :: pset
          REAL(WP), DIMENSION(SIZE(lam)) :: agn_dust
     END FUNCTION AGN_DUST
  END INTERFACE

  INTERFACE
     FUNCTION AIRTOVAC(lam)
          USE fsps_precision, ONLY: WP
       REAL(WP), DIMENSION(:), INTENT(in) :: lam
       REAL(WP), DIMENSION(SIZE(lam)) :: airtovac
     END FUNCTION AIRTOVAC
  END INTERFACE

 INTERFACE
   FUNCTION ATTN_CURVE(ctx, lambda, dtype, pset)
      USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
      USE fsps_types, ONLY: PARAMS
      TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      INTEGER, INTENT(in) :: dtype
      REAL(WP), INTENT(in), DIMENSION(:) :: lambda
      TYPE(PARAMS), INTENT(in) :: pset
      REAL(WP), DIMENSION(SIZE(lambda)) :: attn_curve
     END FUNCTION ATTN_CURVE
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
     SUBROUTINE COMPSP_GRID(pset,nti,specout)
      USE fsps_precision, ONLY: WP
       USE fsps_types, ONLY: PARAMS
       TYPE(PARAMS), INTENT(in) :: pset
       INTEGER, INTENT(in) :: nti
      REAL(WP), DIMENSION(:), INTENT(inout) :: specout
     END SUBROUTINE COMPSP_GRID
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
     SUBROUTINE SSP_GEN(ctx, pset, mass_ssp, lbol_ssp, spec_ssp)
       USE fsps_context_types, ONLY: fsps_context_t
      USE fsps_precision, ONLY: WP
       USE fsps_types, ONLY: PARAMS
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       TYPE(PARAMS), INTENT(in) :: pset
      REAL(WP), INTENT(inout), DIMENSION(:) :: mass_ssp, lbol_ssp
      REAL(WP), INTENT(inout), DIMENSION(:,:) :: spec_ssp
     END SUBROUTINE SSP_GEN
  END INTERFACE

  INTERFACE
     SUBROUTINE GETZMET(smass,pos)
      USE fsps_precision, ONLY: WP
       USE fsps_types, ONLY: PARAMS
      REAL(WP), INTENT(in) :: smass
       TYPE(PARAMS), INTENT(inout) :: pos
     END SUBROUTINE GETZMET
  END INTERFACE

  INTERFACE
       SUBROUTINE GETINDX(ctx, lambda, spec, indices)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(WP), INTENT(in), DIMENSION(:) :: spec,lambda
      REAL(WP), INTENT(inout), DIMENSION(:) :: indices
     END SUBROUTINE GETINDX
  END INTERFACE

  INTERFACE
       FUNCTION GET_TUNIV(ctx, z)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(IN) :: ctx
      REAL(WP), INTENT(in) :: z
      REAL(WP) :: get_tuniv
     END FUNCTION GET_TUNIV
  END INTERFACE
 
  INTERFACE
       FUNCTION GET_LUMDIST(ctx, z)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(IN) :: ctx
      REAL(WP), INTENT(in) :: z
      REAL(WP) :: get_lumdist
     END FUNCTION GET_LUMDIST
  END INTERFACE

   INTERFACE
       SUBROUTINE PZ_CONVOL(ctx, yield, zave, spec_pz, lbol_pz, mass_pz)
          USE fsps_context_types, ONLY: fsps_context_t
           USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(IN) :: ctx
          REAL(WP), INTENT(in) :: yield
          REAL(WP), INTENT(out) :: zave
          REAL(WP), INTENT(out), DIMENSION(:,:) :: spec_pz
          REAL(WP), INTENT(out), DIMENSION(:) :: lbol_pz, mass_pz
       END SUBROUTINE PZ_CONVOL
   END INTERFACE
  
  INTERFACE
       SUBROUTINE GETMAGS(ctx, zred, spec, mags, mag_compute)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(WP), INTENT(in) :: zred
      REAL(WP), INTENT(inout), DIMENSION(:) :: spec
      REAL(WP), DIMENSION(:) :: mags
       INTEGER, DIMENSION(:), INTENT(in), OPTIONAL  :: mag_compute
     END SUBROUTINE GETMAGS
  END INTERFACE
  
  INTERFACE
       SUBROUTINE GETSPEC(ctx, pset, mact, logt, lbol, logg, phase, ffco, lmdot, wght, spec)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          REAL(WP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco,wght,lmdot
          TYPE(PARAMS), INTENT(in) :: pset
      REAL(WP), INTENT(inout), DIMENSION(:) :: spec 
     END SUBROUTINE GETSPEC
  END INTERFACE

  INTERFACE
     FUNCTION IGM_ABSORB(lam,spec,zz,factor)
      USE fsps_precision, ONLY: WP
      REAL(WP), DIMENSION(:), INTENT(in) :: lam,spec
      REAL(WP), INTENT(in) :: zz,factor
      REAL(WP), DIMENSION(SIZE(lam)) :: igm_absorb
     END FUNCTION IGM_ABSORB
  END INTERFACE

  INTERFACE
     FUNCTION INTIND(lam,func,lo,hi)
      USE fsps_precision, ONLY: WP
      REAL(WP), INTENT(in), DIMENSION(:) :: lam,func
      REAL(WP), INTENT(in) :: lo,hi
      REAL(WP) :: intind
     END FUNCTION INTIND
  END INTERFACE

  INTERFACE
       FUNCTION INTSFWGHT(ctx, sspind, logt, sfh)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: SFHPARAMS
          TYPE(fsps_context_t), INTENT(IN) :: ctx
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: sspind
      REAL(WP), DIMENSION(2), INTENT(in) :: logt
      REAL(WP) :: intsfwght
     END FUNCTION INTSFWGHT
  END INTERFACE


  INTERFACE
       SUBROUTINE SBF(ctx, pset, outfile)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)    :: pset
     END SUBROUTINE SBF
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP_TABULAR_SFH(ctx, pset, nzin)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          TYPE(PARAMS), INTENT(in) :: pset
          INTEGER, INTENT(in) :: nzin
     END SUBROUTINE SETUP_TABULAR_SFH
  END INTERFACE

  INTERFACE
     SUBROUTINE SFHINFO(ctx, pset, age, mfrac, sfr, frac_linear)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: PARAMS
          TYPE(fsps_context_t), INTENT(IN) :: ctx
          TYPE(PARAMS), INTENT(in) :: pset
          REAL(WP), INTENT(in) :: age
          REAL(WP), INTENT(out) :: mfrac, sfr, frac_linear
     END SUBROUTINE SFHINFO
  END INTERFACE

  INTERFACE
       FUNCTION SFHLIMIT(ctx, tlim, sfh)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: SFHPARAMS
          TYPE(fsps_context_t), INTENT(IN) :: ctx
       TYPE(SFHPARAMS), INTENT(in) :: sfh
      REAL(WP), INTENT(in) :: tlim
      REAL(WP) :: sfhlimit
     END FUNCTION SFHLIMIT
  END INTERFACE

  INTERFACE
       FUNCTION SFH_WEIGHT(ctx, sfh, imin, imax)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          USE fsps_types, ONLY: SFHPARAMS
          TYPE(fsps_context_t), INTENT(IN) :: ctx
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: imin, imax
          REAL(WP), DIMENSION(ctx%state%ntfull) :: sfh_weight
     END FUNCTION SFH_WEIGHT
  END INTERFACE

  INTERFACE
       SUBROUTINE SMOOTHSPEC(ctx, lambda, spec, sigma, minl, maxl, ires)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
      REAL(WP), INTENT(inout), DIMENSION(:) :: spec
      REAL(WP), INTENT(in), DIMENSION(:) :: lambda
      REAL(WP), INTENT(in), DIMENSION(:), OPTIONAL :: ires
      REAL(WP), INTENT(in) :: sigma,minl,maxl
     END SUBROUTINE SMOOTHSPEC
  END INTERFACE

  INTERFACE
     FUNCTION VACTOAIR(lam)
          USE fsps_precision, ONLY: WP
       REAL(WP), DIMENSION(:), INTENT(in) :: lam
       REAL(WP), DIMENSION(SIZE(lam)) :: vactoair
     END FUNCTION VACTOAIR
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

  INTERFACE
       SUBROUTINE ZTINTERP(ctx, zpos, spec, lbol, mass, tpos, zpow)
          USE fsps_context_types, ONLY: fsps_context_t
          USE fsps_precision, ONLY: WP
          TYPE(fsps_context_t), INTENT(IN) :: ctx
      REAL(WP),INTENT(in) :: zpos
      REAL(WP),INTENT(in), OPTIONAL :: tpos,zpow
      REAL(WP),INTENT(inout),DIMENSION(:) :: mass, lbol
      REAL(WP),INTENT(inout),DIMENSION(:,:) :: spec
     END SUBROUTINE ZTINTERP
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
