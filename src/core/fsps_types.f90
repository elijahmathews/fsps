MODULE FSPS_TYPES

  USE fsps_constants, ONLY: SP, NPZPHOT

  IMPLICIT NONE
  SAVE

  !structure for the set of parameters necessary to generate a model
  TYPE PARAMS
     REAL(SP) :: pagb=1.0,dell=0.,delt=0.,fbhb=0.,sbss=0.,tau=1.0,&
          const=0.,tage=0.,fburst=0.,tburst=11.0,dust1=0.,dust2=0.,&
          logzsol=0.,zred=0.,pmetals=0.02,imf1=1.3,imf2=2.3,imf3=2.3,&
          vdmc=0.08,dust_clumps=-99.,frac_nodust=0.,dust_index=-0.7,&
          dust_tesc=7.0,frac_obrun=0.,uvb=1.0,mwr=3.1,redgb=1.0,agb=1.0,&
          dust1_index=-1.0,mdave=0.5,sf_start=0.,sf_trunc=0.,sf_slope=0.,&
          duste_gamma=0.01,duste_umin=1.0,duste_qpah=3.5,fcstar=1.0,&
          masscut=150.0,sigma_smooth=0.,agb_dust=1.0,min_wave_smooth=1E3,&
          max_wave_smooth=1E4,gas_logu=-2.0,gas_logz=0.,igm_factor=1.0,&
          fagn=0.0,agn_tau=10.0,frac_xrb=1.0,dust3=0.
     INTEGER :: zmet=1,sfh=0,wgp1=1,wgp2=1,wgp3=1,evtype=-1
     INTEGER, ALLOCATABLE :: mag_compute(:)
     INTEGER, ALLOCATABLE :: ssp_gen_age(:)
     CHARACTER(50) :: imf_filename='', sfh_filename=''
  END TYPE PARAMS

  !structure for the output of the compsp routine
  TYPE COMPSPOUT
     REAL(SP) :: age=0.,mass_csp=0.,lbol_csp=0.,sfr=0.,mdust=0.,mformed=0.
     REAL(SP), ALLOCATABLE  :: mags(:)
     REAL(SP), ALLOCATABLE   :: spec(:)
     REAL(SP), ALLOCATABLE   :: indx(:)
     REAL(SP), ALLOCATABLE :: emlines(:)
  END TYPE COMPSPOUT

  ! A structure to hold SFH params converted to intrinsic units
  TYPE SFHPARAMS
     REAL(SP) :: tau=1.0,tage=0.,tburst=0.,sf_trunc=0.,sf_slope=0.,&
          tq=0.,t0=0.,tb=0.
     INTEGER :: type=0,use_simha_limits=0
  END TYPE SFHPARAMS

  TYPE TLSF
     REAL(SP), ALLOCATABLE :: lsf(:)
     REAL(SP) :: minlam=0.,maxlam=0.
  END TYPE TLSF

  !structure for observational data
  TYPE OBSDAT
     REAL(SP)                    :: zred=0.,logsmass=0.
     REAL(SP), ALLOCATABLE :: mags(:),magerr(:)
     REAL(SP), ALLOCATABLE  :: spec(:),specerr(:)
  END TYPE OBSDAT

  TYPE TPZPHOT
     REAL(SP), DIMENSION(NPZPHOT) :: zz=0.,pz=0.
  END TYPE TPZPHOT

END MODULE FSPS_TYPES
