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
       USE sps_vars
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), DIMENSION(nspec), INTENT(out) :: tspec
       REAL(SP), INTENT(in)  :: weight,mact,logt,logl,logg,zz,tco,lmdot
     END SUBROUTINE ADD_AGB_DUST
  END INTERFACE

  INTERFACE
   SUBROUTINE ADD_BS(ctx, s_bs, t, mini, mact, logl, logt, logg, phase, &
          wght,hb_wght,nmass)
     USE fsps_context_types, ONLY: fsps_context_t
     USE sps_vars
     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,logg,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in) :: hb_wght,s_bs
       INTEGER, INTENT(in)  :: t
       INTEGER, INTENT(inout), DIMENSION(nt)  :: nmass
     END SUBROUTINE ADD_BS
  END INTERFACE

   INTERFACE
       SUBROUTINE ADD_DUST(ctx, pset, csp1, csp2, specdust, mdust, ncsp1, ncsp2, nebdust)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          REAL(SP), INTENT(out) :: mdust
          REAL(SP), DIMENSION(nspec), INTENT(in) :: csp1,csp2
          TYPE(PARAMS), INTENT(in) :: pset
          REAL(SP), DIMENSION(nspec), INTENT(out) :: specdust
          REAL(SP), DIMENSION(nemline), INTENT(in) :: ncsp1,ncsp2
          REAL(SP), DIMENSION(nemline), INTENT(out) :: nebdust
       END SUBROUTINE ADD_DUST
   END INTERFACE

  INTERFACE
       SUBROUTINE ADD_NEBULAR(ctx, pset, sspi, sspo, nebemline)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
       REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
       REAL(SP), INTENT(inout), DIMENSION(nemline,ntfull), OPTIONAL :: nebemline
     END SUBROUTINE ADD_NEBULAR
  END INTERFACE

  INTERFACE
       SUBROUTINE ADD_XRB(ctx, pset, sspi, sspo)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in), DIMENSION(nspec,ntfull)    :: sspi
       REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: sspo
     END SUBROUTINE ADD_XRB
  END INTERFACE

  INTERFACE
       SUBROUTINE ADD_REMNANTS(ctx, mass, maxmass)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(inout) :: mass
       REAL(SP), INTENT(in) :: maxmass
     END SUBROUTINE ADD_REMNANTS
  END INTERFACE

  INTERFACE
     FUNCTION AGN_DUST(lam,spec,pset,lbol_csp)
       USE sps_vars
       REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
       REAL(SP), INTENT(in) :: lbol_csp
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec) :: agn_dust
     END FUNCTION AGN_DUST
  END INTERFACE

  INTERFACE
     FUNCTION AIRTOVAC(lam)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: lam
       REAL(SP), DIMENSION(SIZE(lam)) :: airtovac
     END FUNCTION AIRTOVAC
  END INTERFACE

 INTERFACE
     FUNCTION ATTN_CURVE(lambda,dtype,pset)
       USE sps_vars
       INTEGER, INTENT(in) :: dtype
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), DIMENSION(nspec) :: attn_curve
     END FUNCTION ATTN_CURVE
  END INTERFACE

  INTERFACE
     SUBROUTINE COMPSP(ctx, write_compsp, nzin, outfile, mass_ssp, &
          lbol_ssp, spec_ssp, pset, ocompsp)
       USE fsps_context_types, ONLY: fsps_context_t
       USE sps_vars
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       INTEGER, INTENT(in) :: write_compsp,nzin
       REAL(SP), INTENT(in), DIMENSION(ntfull,nzin) :: lbol_ssp,mass_ssp
       REAL(SP), INTENT(in), DIMENSION(nspec,ntfull,nzin) :: spec_ssp
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)   :: pset
       TYPE(COMPSPOUT), INTENT(inout), DIMENSION(ntfull) :: ocompsp
     END SUBROUTINE COMPSP
  END INTERFACE

  INTERFACE
     SUBROUTINE COMPSP_GRID(pset,nti,specout)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       INTEGER, INTENT(in) :: nti
       REAL(SP), DIMENSION(nspec), INTENT(inout) :: specout
     END SUBROUTINE COMPSP_GRID
  END INTERFACE

  INTERFACE
     SUBROUTINE CSP_GEN(ctx, mass_ssp, lbol_ssp, spec_ssp, pset, tage, nzin,&
                        mass_csp, lbol_csp, spec_csp, mdust_csp, emlin_ssp, emlin_csp)
       USE fsps_context_types, ONLY: fsps_context_t
       USE sps_vars
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), DIMENSION(ntfull, nzin), INTENT(in) :: mass_ssp, lbol_ssp
       REAL(SP), DIMENSION(nspec, ntfull, nzin), INTENT(in) :: spec_ssp
       TYPE(PARAMS), intent(in) :: pset
       REAL(SP), INTENT(in)  :: tage
       INTEGER, INTENT(IN) :: nzin
       REAL(SP), INTENT(out) :: mass_csp, lbol_csp, mdust_csp
       REAL(SP), INTENT(out), DIMENSION(nspec) :: spec_csp
       REAL(SP), DIMENSION(nemline, ntfull, nzin), intent(in) :: emlin_ssp
       REAL(SP), DIMENSION(nemline), intent(out) :: emlin_csp
     END SUBROUTINE CSP_GEN
  END INTERFACE

  INTERFACE
     FUNCTION FUNCINT(func,a,b)
       USE sps_vars
       REAL(SP), INTENT(IN) :: a,b
       REAL(SP) :: funcint
       INTERFACE
          FUNCTION func(x)
            USE sps_vars
            REAL(SP), DIMENSION(:), INTENT(IN) :: x
            REAL(SP), DIMENSION(SIZE(x)) :: func
          END FUNCTION func
       END INTERFACE
     END FUNCTION FUNCINT
  END INTERFACE

  INTERFACE
     SUBROUTINE GETZMET(smass,pos)
       USE sps_vars
       REAL(SP), INTENT(in) :: smass
       TYPE(PARAMS), INTENT(inout) :: pos
     END SUBROUTINE GETZMET
  END INTERFACE

  INTERFACE
       SUBROUTINE GETINDX(ctx, lambda, spec, indices)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(in), DIMENSION(nspec) :: spec,lambda
       REAL(SP), INTENT(inout), DIMENSION(nindx) :: indices
     END SUBROUTINE GETINDX
  END INTERFACE

  INTERFACE
     FUNCTION GET_TUNIV(z)
       USE sps_vars
       REAL(SP), INTENT(in) :: z
       REAL(SP) :: get_tuniv
     END FUNCTION GET_TUNIV
  END INTERFACE
 
  INTERFACE
     FUNCTION GET_LUMDIST(z)
       USE sps_vars
       REAL(SP), INTENT(in) :: z
       REAL(SP) :: get_lumdist
     END FUNCTION GET_LUMDIST
  END INTERFACE
  
  INTERFACE
       SUBROUTINE GETMAGS(ctx, zred, spec, mags, mag_compute)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(in) :: zred
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), DIMENSION(nbands) :: mags
       INTEGER, DIMENSION(nbands), INTENT(in), OPTIONAL  :: mag_compute
     END SUBROUTINE GETMAGS
  END INTERFACE
  
  INTERFACE
       SUBROUTINE GETSPEC(ctx, pset, mact, logt, lbol, logg, phase, ffco, lmdot, wght, spec)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          REAL(SP), INTENT(in) :: mact,logt,lbol,logg,phase,ffco,wght,lmdot
          TYPE(PARAMS), INTENT(in) :: pset
          REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec 
     END SUBROUTINE GETSPEC
  END INTERFACE

  INTERFACE
     FUNCTION IGM_ABSORB(lam,spec,zz,factor)
       USE sps_vars
       REAL(SP), DIMENSION(nspec), INTENT(in) :: lam,spec
       REAL(SP), INTENT(in) :: zz,factor
       REAL(SP), DIMENSION(nspec) :: igm_absorb
     END FUNCTION IGM_ABSORB
  END INTERFACE

  INTERFACE
     FUNCTION INTIND(lam,func,lo,hi)
       USE sps_vars
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lam,func
       REAL(SP), INTENT(in) :: lo,hi
       REAL(SP) :: intind
     END FUNCTION INTIND
  END INTERFACE

  INTERFACE
     FUNCTION INTSFWGHT(sspind, logt, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: sspind
       REAL(SP), DIMENSION(2), INTENT(in) :: logt
       REAL(SP) :: intsfwght
     END FUNCTION INTSFWGHT
  END INTERFACE

  INTERFACE
     FUNCTION IMF(mass)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: mass
       REAL(SP), DIMENSION(size(mass)) :: imf
     END FUNCTION IMF
  END INTERFACE 

  INTERFACE
       SUBROUTINE IMF_WEIGHT(ctx, mini, wght, nmass)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), INTENT(in), DIMENSION(nm)    :: mini
       INTEGER, INTENT(in) :: nmass
     END SUBROUTINE IMF_WEIGHT
  END INTERFACE 

  INTERFACE
     FUNCTION LINTERP(xin,yin,xout)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP), INTENT(in)  :: xout
       REAL(SP) :: linterp
     END FUNCTION LINTERP
  END INTERFACE

  INTERFACE
     FUNCTION LINTERPARR(xin,yin,xout)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP), INTENT(in), DIMENSION(:) :: xout
       REAL(SP), DIMENSION(SIZE(xout)) :: linterparr
     END FUNCTION LINTERPARR
  END INTERFACE

  INTERFACE
     FUNCTION LOCATE(xx,x)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(IN) :: xx
       REAL(SP), INTENT(IN) :: x
       INTEGER :: locate
     END FUNCTION LOCATE
  END INTERFACE

  INTERFACE
     SUBROUTINE MOD_GB(ctx, zz, t, age, delt, dell, pagb, redgb, agb, &
          nn, logl, logt, phase, wght)
       USE fsps_context_types, ONLY: fsps_context_t
       USE sps_vars
       TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       INTEGER,  INTENT(in) :: t, nn,zz
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: logl,logt
       REAL(SP), INTENT(in), DIMENSION(nt,nm)    :: phase
       REAL(SP), INTENT(inout), DIMENSION(nm)    :: wght
       REAL(SP), INTENT(in) :: delt, dell, pagb,redgb, agb
       REAL(SP), INTENT(in), DIMENSION(nt) :: age
     END SUBROUTINE MOD_GB
  END INTERFACE

  INTERFACE
   SUBROUTINE MOD_HB(ctx, f_bhb, t, mini, mact, logl, logt, logg, phase, &
      wght, hb_wght, nmass, hbtime)
     USE fsps_context_types, ONLY: fsps_context_t
     USE sps_vars
     TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(inout), DIMENSION(nt,nm) :: mini,mact,&
            logl,logt,logg,phase
       REAL(SP), INTENT(inout), DIMENSION(nm) :: wght
       REAL(SP), DIMENSION(nm) :: tphase=0.0
       INTEGER, INTENT(inout), DIMENSION(nt) :: nmass
       REAL(SP), INTENT(inout) :: hb_wght
       INTEGER, INTENT(in) :: t
       REAL(SP), INTENT(in) :: f_bhb, hbtime
     END SUBROUTINE MOD_HB
  END INTERFACE

  INTERFACE
       SUBROUTINE SBF(ctx, pset, outfile)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       CHARACTER(100), INTENT(in) :: outfile
       TYPE(PARAMS), INTENT(in)    :: pset
     END SUBROUTINE SBF
  END INTERFACE

  INTERFACE
     SUBROUTINE SETUP_TABULAR_SFH(pset, nzin)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       INTEGER, INTENT(in) :: nzin
     END SUBROUTINE SETUP_TABULAR_SFH
  END INTERFACE

  INTERFACE
     SUBROUTINE SFHINFO(pset, age, mfrac, sfr, frac_linear)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(in) :: age
       REAL(SP), INTENT(out) :: mfrac, sfr, frac_linear
     END SUBROUTINE SFHINFO
  END INTERFACE

  INTERFACE
     FUNCTION SFHLIMIT(tlim, sfh)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       REAL(SP), INTENT(in) :: tlim
       REAL(SP) :: sfhlimit
     END FUNCTION SFHLIMIT
  END INTERFACE

  INTERFACE
     FUNCTION SFH_WEIGHT(sfh, imin, imax)
       USE sps_vars
       TYPE(SFHPARAMS), INTENT(in) :: sfh
       INTEGER, INTENT(in) :: imin, imax
       REAL(SP), DIMENSION(ntfull) :: sfh_weight
     END FUNCTION SFH_WEIGHT
  END INTERFACE

  INTERFACE
     FUNCTION TSUM(xin,yin)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: xin,yin
       REAL(SP) :: tsum
     END FUNCTION TSUM
  END INTERFACE

  INTERFACE
       SUBROUTINE SMOOTHSPEC(ctx, lambda, spec, sigma, minl, maxl, ires)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
       REAL(SP), INTENT(inout), DIMENSION(nspec) :: spec
       REAL(SP), INTENT(in), DIMENSION(nspec) :: lambda
       REAL(SP), INTENT(in), DIMENSION(nspec), OPTIONAL :: ires
       REAL(SP), INTENT(in) :: sigma,minl,maxl
     END SUBROUTINE SMOOTHSPEC
  END INTERFACE

  INTERFACE
     FUNCTION VACTOAIR(lam)
       USE sps_vars
       REAL(SP), DIMENSION(:), INTENT(in) :: lam
       REAL(SP), DIMENSION(SIZE(lam)) :: vactoair
     END FUNCTION VACTOAIR
  END INTERFACE

  INTERFACE
       SUBROUTINE WRITE_ISOCHRONE(ctx, outfile, pset)
          USE fsps_context_types, ONLY: fsps_context_t
          USE sps_vars
          TYPE(fsps_context_t), INTENT(INOUT) :: ctx
          TYPE(PARAMS), INTENT(in) :: pset
          CHARACTER(100), INTENT(in)  :: outfile
     END SUBROUTINE WRITE_ISOCHRONE
  END INTERFACE

  INTERFACE
     SUBROUTINE ZTINTERP(zpos,spec,lbol,mass,tpos,zpow)
       USE sps_vars
       REAL(SP),INTENT(in) :: zpos
       REAL(SP),INTENT(in), OPTIONAL :: tpos,zpow
       REAL(SP),INTENT(inout),DIMENSION(:) :: mass, lbol
       REAL(SP),INTENT(inout),DIMENSION(:,:) :: spec
     END SUBROUTINE ZTINTERP
  END INTERFACE

CONTAINS

  LOGICAL FUNCTION fsps_data_exists(path)
    USE sps_vars
    CHARACTER(LEN=*), INTENT(IN) :: path
    LOGICAL :: ok

    ok = .FALSE.
    INQUIRE(FILE=TRIM(path)//'/allfilters.dat', EXIST=ok)
    IF (.NOT. ok) INQUIRE(FILE=TRIM(path)//'/FILTER_LIST', EXIST=ok)
    fsps_data_exists = ok
  END FUNCTION fsps_data_exists

  SUBROUTINE fsps_resolve_paths()
    USE sps_vars
    CHARACTER(250) :: env, candidate
    LOGICAL :: system_prefix

    DATA_HOME = ''
    OUTPUT_HOME = ''

    CALL getenv('SPS_HOME', SPS_HOME)
    IF (LEN_TRIM(SPS_HOME) > 0) THEN
       candidate = TRIM(SPS_HOME)
       IF (.NOT. fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = ''
    ENDIF

    CALL getenv('FSPS_DATA_HOME', env)
    IF (LEN_TRIM(env) > 0) THEN
       candidate = TRIM(env)
       IF (fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = candidate
    ENDIF

    IF (LEN_TRIM(SPS_HOME) == 0) THEN
       CALL getenv('XDG_DATA_HOME', env)
       IF (LEN_TRIM(env) > 0) THEN
          candidate = TRIM(env)//'/fsps'
          IF (fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = candidate
       ENDIF
    ENDIF

    IF (LEN_TRIM(SPS_HOME) == 0) THEN
       CALL getenv('HOME', env)
       IF (LEN_TRIM(env) > 0) THEN
          candidate = TRIM(env)//'/.local/share/fsps'
          IF (fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = candidate
       ENDIF
    ENDIF

    IF (LEN_TRIM(SPS_HOME) == 0) THEN
       candidate = '/usr/share/fsps'
       IF (fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = candidate
    ENDIF

    IF (LEN_TRIM(SPS_HOME) == 0) THEN
       candidate = '/usr/local/share/fsps'
       IF (fsps_data_exists(TRIM(candidate)//'/data')) SPS_HOME = candidate
    ENDIF

    IF (LEN_TRIM(SPS_HOME) == 0) THEN
       WRITE(*,*) 'SPS_SETUP ERROR: FSPS data path not found. Set FSPS_DATA_HOME or SPS_HOME.'
       STOP
    ENDIF

    DATA_HOME = TRIM(SPS_HOME)//'/data'

    CALL getenv('FSPS_OUTPUT_HOME', env)
    IF (LEN_TRIM(env) > 0) THEN
       OUTPUT_HOME = TRIM(env)
    ELSE
       system_prefix = .FALSE.
       IF (LEN_TRIM(SPS_HOME) >= 5) THEN
          IF (SPS_HOME(1:5) == '/usr/') system_prefix = .TRUE.
       ENDIF
       IF (LEN_TRIM(SPS_HOME) >= 10) THEN
          IF (SPS_HOME(1:10) == '/usr/local') system_prefix = .TRUE.
       ENDIF
       IF (system_prefix) THEN
          CALL getenv('HOME', env)
          IF (LEN_TRIM(env) > 0) THEN
             OUTPUT_HOME = TRIM(env)//'/.local/share/fsps'
          ELSE
             OUTPUT_HOME = '.'
          ENDIF
       ELSE
          OUTPUT_HOME = TRIM(SPS_HOME)
       ENDIF
    ENDIF
  END SUBROUTINE fsps_resolve_paths

   SUBROUTINE SPS_TAKEDOWN(ctx)
      USE sps_vars
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
