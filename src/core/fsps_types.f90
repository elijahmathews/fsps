MODULE FSPS_TYPES

  IMPLICIT NONE
  SAVE

  !note that "SP" actually means double precision; this is a hack
  !to turn the nr routines into DP
  INTEGER, PARAMETER :: SP = KIND(1.d0)

  !controls the level of output
  !0 = minimal output to screen.
  !1 = lots of output to screen.  useful for debugging.
  INTEGER, PARAMETER :: verbose=0

  !turn-on time for BHB and SBS phases, time is in log(yrs)
  REAL(SP), PARAMETER :: bhb_sbs_time=9.5

  !the factor by which we increase the time array
  !this should no longer need to be set to anything other than 1
  INTEGER, PARAMETER :: time_res_incr=1

  !Use Aringer et al. (2009) Carbon star library if set
  !otherwise use Lancon & Wood (2002) empirical spectra
  INTEGER, PARAMETER :: cstar_aringer=1

  !flag indicating the type of normalization used in the BaSeL library
  !pdva = normalized to Padova isochrones
  !wlbc = normalized to Teff-color relations
  !NB: currently only the wlbc option is included in the public release
  CHARACTER(4), PARAMETER :: basel_str = 'wlbc'

  !---------Dimensions of various arrays----------!

  !max dimension of array for each isochrone
  INTEGER, PARAMETER :: nm=2000
  !max number of lines to read in
  INTEGER, PARAMETER ::  nlines=1000000
  !max number of lines in tabulated SFH, LSF
  INTEGER, PARAMETER :: ntabmax=20000
  !dimensions of BaSeL library
  INTEGER, PARAMETER :: ndim_logt=68, ndim_logg=19
  !number of O-rich, C-rich AGB spectra (and Aringer C-rich spec)
  INTEGER, PARAMETER :: n_agb_o=9, n_agb_c=5, n_agb_car=9
  !number of post-AGB spectra
  INTEGER, PARAMETER :: ndim_pagb=14
  !number of WR spectra
  INTEGER, PARAMETER :: ndim_wr=12
  !dimensions of WMBasic grid
  INTEGER, PARAMETER :: ndim_wmb_logt=11,ndim_wmb_logg=3
  !parameters for circumstellar dust models
  INTEGER, PARAMETER :: ntau_dagb=50, nteff_dagb=6
  !number of emission lines and continuum emission points
  INTEGER, PARAMETER :: nemline=166, nlam_nebcont=1963
  !number of metallicity, age, and ionization parameter points
  INTEGER, PARAMETER :: nebnz=11, nebnage=10, nebnip=7
  !number of optical depths for AGN dust models
  INTEGER, PARAMETER :: nagndust=9
  !number of spectral points in the input library
  INTEGER, PARAMETER :: nagndust_spec=125

  !------------IMF-related Constants--------------!

  !Chabrier 2003 IMF parameters
  REAL(SP), PARAMETER :: chab_mc=0.08, chab_sigma2=0.69*0.69,&
       chab_ind=1.3
  !van Dokkum 2008 IMF parameters
  REAL(SP), PARAMETER :: vd_sigma2=0.69*0.69, vd_ah=0.0443,&
       vd_ind=1.3, vd_al=0.14, vd_nc=25.

  !-------------Physical Constants---------------!
  !-------in cgs units where applicable----------!

  !constant such that g = C MT^4/L
  REAL(SP), PARAMETER :: gsig4pi = 1/4.13E10
  !pi
  REAL(SP), PARAMETER :: mypi    = 3.14159265
  !hc/k (Ang*K)
  REAL(SP), PARAMETER :: hck     = 1.43878E8
  !speed of light (Ang/s)
  REAL(SP), PARAMETER :: clight  = 2.9979E18
  !hc^2/sigma_SB
  REAL(SP), PARAMETER :: hc2sig  = 0.105021
  !Solar mass in grams
  REAL(SP), PARAMETER :: msun    = 1.989E33
  !Solar radius in cm
  REAL(SP), PARAMETER :: rsun    = 6.955E10
  !Solar luminosity in erg/s
  REAL(SP), PARAMETER :: lsun    = 3.839E33
  !Newton's constant
  REAL(SP), PARAMETER :: newton  = 6.67428E-8
  !cm in a pc
  REAL(SP), PARAMETER :: pc2cm   = 3.08568E18
  !seconds per year
  REAL(SP), PARAMETER :: yr2sc   = 3.15569E7
  !Planck's constant
  REAL(SP), PARAMETER :: hplank  = 6.6261E-27
   !constant to convert mags into propert units (see spec_mags.f90)
  REAL(SP), PARAMETER :: mag2cgs = LOG10(lsun/4.0/mypi/(pc2cm*pc2cm)/100.0)

  !define large and small numbers.  numbers whose abs values
  !are less than tiny_number are treated as equal to 0.0
  REAL(SP), PARAMETER :: huge_number = 10**(70.d0)
  REAL(SP), PARAMETER :: tiny_number = 10**(-70.d0)
  REAL(SP), PARAMETER :: tiny30      = 10**(-30.0)

  !-----the following structures are not used in the public code-----!
  !--they are included here because some users of FSPS utilize them--!

  !structure for using P(z) in chi2
  INTEGER, PARAMETER :: npzphot   = 200

  !------------Define TYPE structures-------------!

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
     REAL(SP), DIMENSION(npzphot) :: zz=0.,pz=0.
  END TYPE TPZPHOT

END MODULE FSPS_TYPES
