SUBROUTINE SPS_SETUP(zin, input_isoc_type, input_spec_type)

  !read in isochrones and spectral libraries for all metallicities.
  !read in band-pass info and the spectrum for Vega.
  !Arrays are stored in a common block defined in sps_vars.f90

  !If zin=-1 then all metallicities are read in, otherwise only the
  !metallicity corresponding to zin in the look-up table zlegend.dat
  !is read.  Specifying only the metallicity of interest results
  !in a much faster setup.

  USE sps_vars
  USE sps_utils, ONLY: linterp, linterparr, locate, tsum, get_tuniv, get_lumdist, airtovac
  IMPLICIT NONE
  INTEGER, INTENT(in) :: zin
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: input_isoc_type
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: input_spec_type

  INTEGER :: stat=1,n,i,j,m,jj,k,i1,i2,stat2=1
  INTEGER, PARAMETER :: ntlam=1221,nspec_agb=6146,nspec_aringer=9032
  INTEGER, PARAMETER :: nlamwr=1963,nspec_pagb=9281
  INTEGER, PARAMETER :: nzwmb=12, nspec_wmb=5508
  INTEGER :: n_isoc,z,zmin,zmax,nlam
  CHARACTER(1) :: char,sqpah
  CHARACTER(6) :: zstype
  CHARACTER(5) :: zstype5
  REAL(SP) :: dumr1,d1,d2,logage,x,a,zero=0.0,d,one=1.0,dz,dlam
  
  ! Local allocatable arrays for those that depend on dynamic sizes
  CHARACTER(5), ALLOCATABLE :: zlegend_str(:)
  REAL(SP), ALLOCATABLE :: tspec(:)

  CHARACTER(5), DIMENSION(nz_xrb) :: zz_str_xrb=''
  
  REAL(SP), DIMENSION(ntlam) :: tvega_lam=0.,tvega_spec=0.
  REAL(SP), DIMENSION(ntlam) :: tsun_lam=0.,tsun_spec=0.
  REAL(SP), DIMENSION(nlamwr) :: tlamwr=0.,tspecwr=0.
  REAL(SP), DIMENSION(nlamwr,ndim_wr,5) :: twrc=0.,twrn=0.
  REAL(SP), DIMENSION(5) :: twrzmet=0.
  REAL(SP), DIMENSION(50000) :: readlamb=0.,readband=0.
  REAL(SP), DIMENSION(25) :: wglam=0.
  REAL(SP), DIMENSION(25,18,2,6) :: wgtmp=0.
  REAL(SP), DIMENSION(10000) :: lambda_dagb=0.,fluxin_dagb=0.
  REAL(SP), DIMENSION(14) :: lami=0.
  INTEGER,  DIMENSION(14) :: ind
  REAL(SP), DIMENSION(nlam_nebcont) :: readlambneb=0.,readcontneb=0.
  REAL(SP), DIMENSION(22,n_agb_o)   :: tagb_logt_o
  REAL(SP), DIMENSION(22)           :: tagb_logz_o
  REAL(SP), DIMENSION(nspec_pagb) :: pagb_lam=0.0
  REAL(SP), DIMENSION(nspec_pagb,ndim_pagb,2) :: pagb_specinit=0.
  REAL(SP), DIMENSION(nspec_agb)  :: agb_lam=0.0
  REAL(SP), DIMENSION(nspec_aringer)  :: aringer_lam=0.0
  REAL(SP), DIMENSION(nspec_agb,n_agb_o) :: agb_specinit_o=0.
  REAL(SP), DIMENSION(nspec_agb,n_agb_c) :: agb_specinit_c=0.
  REAL(SP), DIMENSION(nspec_aringer,n_agb_car) :: aringer_specinit=0.
  REAL(SP), DIMENSION(nagndust_spec)           :: agndust_lam=0.
  REAL(SP), DIMENSION(nagndust_spec,nagndust)  :: agndust_specinit=0.
  
  REAL(KIND(1.0)), allocatable :: speclibinit(:,:,:,:)
  REAL(KIND(1.0)), ALLOCATABLE :: spec_buffer(:,:,:)
  REAL(SP), ALLOCATABLE :: wmbsi(:,:,:,:)
  REAL(SP), DIMENSION(nzwmb)     :: zwmb=0.
  REAL(SP), DIMENSION(nspec_wmb) :: wmb_lam=0.
  REAL(SP), DIMENSION(nspec_wmb,ndim_wmb_logt,ndim_wmb_logg) :: wmb_specinit=0.
  REAL(SP), DIMENSION(ntabmax)   :: lsflam=0.,lsfsig=0.
  REAL(SP), DIMENSION(30) :: g03lam=0., g03smc=0.
  REAL(SP), DIMENSION(nspec_xrb) :: tspec_xrb

  CHARACTER(256) :: isoc_dir
  CHARACTER(256) :: file_path
  INTEGER :: ios
  CHARACTER(1024) :: line_buffer

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  IF (verbose.EQ.1) THEN
     WRITE(*,*)
     WRITE(*,*) '    Setting up SPS...'
  ENDIF

  !----------------------------------------------------------------!
  !--------------Confirm that variables are properly set-----------!
  !----------------------------------------------------------------!

  CALL getenv('SPS_HOME',SPS_HOME)
  IF (LEN_TRIM(SPS_HOME).EQ.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: spsdir environment variable not set!'
     STOP
  ENDIF

  ! 1. Determine Isochrone Type
  IF (PRESENT(input_isoc_type)) THEN
     isoc_type = input_isoc_type
  ELSE
     isoc_type = 'mist' ! Default
  ENDIF

  ! Map isoc_type to directory and set zsol
  SELECT CASE (TRIM(isoc_type))
  CASE ('mist')
     isoc_dir = 'MIST'
     zsol = 0.0142
  CASE ('pdva')
     isoc_dir = 'Padova/Padova2007'
     zsol = 0.019
  CASE ('prsc')
     isoc_dir = 'PARSEC'
     zsol = 0.01524
  CASE ('bsti')
     isoc_dir = 'BaSTI'
     zsol = 0.020
  CASE ('gnva')
     isoc_dir = 'Geneva'
     zsol = 0.020
  CASE ('bpss')
     isoc_dir = 'BPASS'
     zsol = 0.020
  CASE DEFAULT
     WRITE(*,*) 'SPS_SETUP ERROR: Unknown isoc_type: ', TRIM(isoc_type)
     STOP
  END SELECT

  ! 2. Determine Spectral Library Type
  IF (PRESENT(input_spec_type)) THEN
     spec_type = input_spec_type
  ELSEIF (TRIM(isoc_type) == 'bpss') THEN
     spec_type = 'bpass'
  ELSE
     spec_type = 'miles' ! Default
  ENDIF

  ! Set zsol_spec
  IF (TRIM(spec_type) == 'miles') THEN
     zsol_spec = 0.019
  ELSEIF (TRIM(spec_type) == 'basel') THEN
     zsol_spec = 0.020
  ELSEIF (INDEX(TRIM(spec_type), 'c3k') > 0) THEN
     zsol_spec = 0.0134
  ELSEIF (TRIM(spec_type) == 'bpass') THEN
     zsol_spec = 0.020
  ENDIF

  ! 3. Determine nbands (Count filters in allfilters.dat)
  nbands = 0
  IF (TRIM(alt_filter_file).EQ.'') THEN
     file_path = TRIM(SPS_HOME)//'/data/allfilters.dat'
  ELSE
     file_path = TRIM(SPS_HOME)//'/data/'//TRIM(alt_filter_file)
  ENDIF
  
  OPEN(99, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Filter file cannot be opened: ', TRIM(file_path)
     STOP
  ENDIF
  DO
     READ(99, '(A)', IOSTAT=ios) line_buffer
     ! Check for valid line content before checking EOF
     IF (LEN_TRIM(line_buffer) > 0) THEN
        IF (line_buffer(1:1) == '#') nbands = nbands + 1
     ENDIF
     IF (ios /= 0) EXIT
  ENDDO
  CLOSE(99)
  
  IF (nbands == 0) THEN
      WRITE(*,*) 'SPS_SETUP ERROR: No filters found in ', TRIM(file_path)
      STOP
  ENDIF

  ! 4. Determine nz (Count lines in zlegend.dat)
  file_path = TRIM(SPS_HOME)//'/ISOCHRONES/'//TRIM(isoc_dir)//'/zlegend.dat'
  OPEN(90, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zlegend.dat cannot be opened: ', TRIM(file_path)
     STOP
  ENDIF
  
  nz = 0
  DO
     READ(90, *, IOSTAT=ios)
     IF (ios /= 0) EXIT
     nz = nz + 1
  ENDDO
  CLOSE(90)

  IF (zin.GT.nz) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zin GT nz', zin,nz
     STOP
  ENDIF

  ! 5. Determine nt (Count blocks in first isochrone file)
  ! We need to read the first metallicity to find the file
  OPEN(90, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  ! Read the first entry to construct filename
  IF (TRIM(isoc_type) == 'mist') THEN
      READ(90, '(A)') line_buffer ! Read as string for MIST
      ! line_buffer holds zlegend_str(1)
      file_path = TRIM(SPS_HOME)//'/ISOCHRONES/MIST/isoc_z'//TRIM(ADJUSTL(line_buffer))//'.dat'
  ELSE
      READ(90, *) dumr1 ! Read as float
      WRITE(zstype, '(F6.4)') dumr1
      file_path = TRIM(SPS_HOME)//'/ISOCHRONES/'//TRIM(isoc_dir)//'/isoc_z'//zstype//'.dat'
  ENDIF
  CLOSE(90)

  nt = 0
  OPEN(97, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  IF (stat.NE.0) THEN
      WRITE(*,*) 'SPS_SETUP ERROR: Cannot open isochrone file to determine nt: ', TRIM(file_path)
      STOP
  ENDIF
  
  DO
     READ(97, '(A)', IOSTAT=ios) char
     IF (ios /= 0) EXIT
     IF (char == '#') nt = nt + 1
  ENDDO
  CLOSE(97)

  ! 6. Determine nspec (Count lines in .lambda file)
  IF (TRIM(isoc_type) == 'bpss') THEN
      file_path = TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass.lambda'
  ELSE
      IF (spec_type.EQ.'basel') THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel.lambda'
      ELSE IF (spec_type.EQ.'miles') THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/MILES/miles.lambda'
      ELSE IF (INDEX(spec_type, 'c3k') > 0) THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/C3K/'//TRIM(spec_type)//'.lambda'
      ENDIF
  ENDIF

  OPEN(91, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: wavelength grid cannot be opened: ', TRIM(file_path)
     STOP
  ENDIF
  nspec = 0
  DO
     READ(91, *, IOSTAT=ios)
     IF (ios /= 0) EXIT
     nspec = nspec + 1
  ENDDO
  CLOSE(91)

  PRINT *, "DEBUG: Detected nspec =", nspec
  PRINT *, "DEBUG: spec_type =", TRIM(spec_type)

  ! 7. Determine nzinit (Count lines in spectral zlegend.dat)
  IF (TRIM(isoc_type) == 'bpss') THEN
      nzinit = 1
  ELSE
      IF (spec_type.EQ.'basel') THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/zlegend.dat'
      ELSE IF (spec_type.EQ.'miles') THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/MILES/zlegend.dat'
      ELSE IF (INDEX(spec_type, 'c3k') > 0) THEN
         file_path = TRIM(SPS_HOME)//'/SPECTRA/C3K/zlegend.dat'
      ENDIF
      
      OPEN(93, FILE=file_path, STATUS='OLD', IOSTAT=stat, ACTION='READ')
      IF (stat.NE.0) THEN
         WRITE(*,*) 'SPS_SETUP ERROR: spectral zlegend cannot be opened: ', TRIM(file_path)
         STOP
      ENDIF
      nzinit = 0
      DO
         READ(93, *, IOSTAT=ios)
         IF (ios /= 0) EXIT
         nzinit = nzinit + 1
      ENDDO
      CLOSE(93)
  ENDIF

  ! 8. Calculate ntfull
  ntfull = time_res_incr * nt

  !----------------------------------------------------------------!
  !-------------------ALLOCATE ARRAYS------------------------------!
  !----------------------------------------------------------------!
  
  ! Local arrays
  IF (ALLOCATED(zlegend_str)) DEALLOCATE(zlegend_str)
  ALLOCATE(zlegend_str(nz), stat=stat)
  IF (ALLOCATED(tspec)) DEALLOCATE(tspec)
  ALLOCATE(tspec(nspec), stat=stat)
  tspec = 0.0
  
  ! sps_vars arrays
  IF (ALLOCATED(mact_isoc)) DEALLOCATE(mact_isoc, logl_isoc, logt_isoc, &
       logg_isoc, ffco_isoc, phase_isoc, mini_isoc, lmdot_isoc)
  ALLOCATE(mact_isoc(nz,nt,nm), logl_isoc(nz,nt,nm), logt_isoc(nz,nt,nm), &
           logg_isoc(nz,nt,nm), ffco_isoc(nz,nt,nm), phase_isoc(nz,nt,nm), &
           mini_isoc(nz,nt,nm), lmdot_isoc(nz,nt,nm), stat=stat)
  IF (stat /= 0) STOP 'Allocation failed for isochrone arrays'

  mact_isoc = 0.0; logl_isoc = 0.0; logt_isoc = 0.0
  logg_isoc = 0.0; ffco_isoc = 0.0; phase_isoc = 0.0
  mini_isoc = 0.0; lmdot_isoc = 0.0

  IF (ALLOCATED(nmass_isoc)) DEALLOCATE(nmass_isoc, timestep_isoc)
  ALLOCATE(nmass_isoc(nz,nt), timestep_isoc(nz,nt), stat=stat)
  IF (ALLOCATED(zlegend)) DEALLOCATE(zlegend, zlegendinit)
  ALLOCATE(zlegend(nz), zlegendinit(nzinit), stat=stat)
  nmass_isoc = 0; timestep_isoc = 0.0
  
  IF (ALLOCATED(spec_ssp_zz)) DEALLOCATE(spec_ssp_zz, mass_ssp_zz, lbol_ssp_zz, time_full)
  ALLOCATE(spec_ssp_zz(nspec,ntfull,nz), mass_ssp_zz(ntfull,nz), &
           lbol_ssp_zz(ntfull,nz), time_full(ntfull), stat=stat)
  spec_ssp_zz = 0.0; mass_ssp_zz = 0.0 
  lbol_ssp_zz = 0.0; time_full = 0.0

  IF (ALLOCATED(weight_ssp)) DEALLOCATE(weight_ssp)
  ALLOCATE(weight_ssp(ntfull,nz), stat=stat)
  weight_ssp = 0.0

  IF (ALLOCATED(spec_young)) DEALLOCATE(spec_young, spec_old)
  ALLOCATE(spec_young(nspec), spec_old(nspec), stat=stat)
  spec_young = 0.0; spec_old = 0.0
  
  IF (ALLOCATED(bpass_spec_ssp)) DEALLOCATE(bpass_spec_ssp, bpass_mass_ssp)
  ALLOCATE(bpass_spec_ssp(nspec,nt,nz), bpass_mass_ssp(nt,nz), stat=stat)
  bpass_spec_ssp = 0.0; bpass_mass_ssp = 0.0

  IF (ALLOCATED(spec_xrb)) DEALLOCATE(spec_xrb)
  ALLOCATE(spec_xrb(nspec,nt_xrb,nz_xrb), stat=stat)
  spec_xrb = 0.0
  
  IF (ALLOCATED(bands)) DEALLOCATE(bands)
  ALLOCATE(bands(nspec,nbands), stat=stat)
  bands = 0.0

  IF (ALLOCATED(magsun)) DEALLOCATE(magsun, magvega, filter_leff)
  ALLOCATE(magsun(nbands), magvega(nbands), filter_leff(nbands), stat=stat)
  magsun = 0.0; magvega = 0.0; filter_leff = 0.0

  IF (ALLOCATED(vega_spec)) DEALLOCATE(vega_spec, sun_spec, spec_lambda, spec_nu, spec_res)
  ALLOCATE(vega_spec(nspec), sun_spec(nspec), spec_lambda(nspec), &
           spec_nu(nspec), spec_res(nspec), stat=stat)
  vega_spec = 0.0; sun_spec = 0.0; spec_lambda = 0.0
  spec_nu = 0.0; spec_res = 0.0
           
  IF (ALLOCATED(speclib)) DEALLOCATE(speclib)
  ALLOCATE(speclib(nspec,nz,ndim_logt,ndim_logg), stat=stat)
  speclib = 0.0

  IF (ALLOCATED(speclibinit)) DEALLOCATE(speclibinit)
  ALLOCATE(speclibinit(nspec,nzinit,ndim_logt,ndim_logg), stat=stat)
  speclibinit = 0.0

  IF (ALLOCATED(wmbsi)) DEALLOCATE(wmbsi)
  ALLOCATE(wmbsi(nspec,nzwmb,ndim_wmb_logt,ndim_wmb_logg), stat=stat)
  wmbsi = 0.0

  IF (ALLOCATED(wmb_spec)) DEALLOCATE(wmb_spec)
  ALLOCATE(wmb_spec(nspec,nz,ndim_wmb_logt,ndim_wmb_logg), stat=stat)
  wmb_spec = 0.0
  
  IF (ALLOCATED(agb_spec_o)) DEALLOCATE(agb_spec_o, agb_logt_o)
  ALLOCATE(agb_spec_o(nspec,n_agb_o), agb_logt_o(nz,n_agb_o), stat=stat)
  agb_spec_o = 0.0; agb_logt_o = 0.0

  IF (ALLOCATED(agb_spec_c)) DEALLOCATE(agb_spec_c, agb_logt_c)
  ALLOCATE(agb_spec_c(nspec,n_agb_c), agb_logt_c(n_agb_c), stat=stat)
  agb_spec_c = 0.0; agb_logt_c = 0.0

  IF (ALLOCATED(agb_logt_car)) DEALLOCATE(agb_logt_car, agb_spec_car)
  ALLOCATE(agb_logt_car(n_agb_car), agb_spec_car(nspec,n_agb_car), stat=stat)
  agb_logt_car = 0.0; agb_spec_car = 0.0
  
  IF (ALLOCATED(pagb_spec)) DEALLOCATE(pagb_spec)
  ALLOCATE(pagb_spec(nspec,ndim_pagb,2), stat=stat)
  pagb_spec = 0.0

  IF (ALLOCATED(wrn_spec)) DEALLOCATE(wrn_spec, wrc_spec)
  ALLOCATE(wrn_spec(nspec,ndim_wr,nz), wrc_spec(nspec,ndim_wr,nz), stat=stat)
  wrn_spec = 0.0; wrc_spec = 0.0
  
  IF (ALLOCATED(dustem2_dustem)) DEALLOCATE(dustem2_dustem)
  ALLOCATE(dustem2_dustem(nspec,nqpah_dustem,numin_dustem*2), stat=stat)
  dustem2_dustem = 0.0

  IF (ALLOCATED(flux_dagb)) DEALLOCATE(flux_dagb)
  ALLOCATE(flux_dagb(nspec,2,nteff_dagb,ntau_dagb), stat=stat)
  flux_dagb = 0.0
  
  IF (ALLOCATED(nebem_cont)) DEALLOCATE(nebem_cont, xnebem_cont)
  ALLOCATE(nebem_cont(nspec,nebnz,nebnage,nebnip), &
           xnebem_cont(nspec,nebnz,nebnage,nebnip), stat=stat)
  nebem_cont = 0.0; xnebem_cont = 0.0   

  IF (ALLOCATED(neb_res_min)) DEALLOCATE(neb_res_min, gaussnebarr)
  ALLOCATE(neb_res_min(nspec), gaussnebarr(nspec,nemline), stat=stat)
  neb_res_min = 0.0; gaussnebarr = 0.0
  
  IF (ALLOCATED(agndust_spec)) DEALLOCATE(agndust_spec)
  ALLOCATE(agndust_spec(nspec,nagndust), stat=stat)
  agndust_spec = 0.0
  
  IF (ALLOCATED(mwdindex)) DEALLOCATE(mwdindex, wgdust, g03smcextn)
  ALLOCATE(mwdindex(nspec), wgdust(nspec,18,6,2), g03smcextn(nspec), stat=stat)
  mwdindex = 0; wgdust = 0.0; g03smcextn = 0.0
  
  ! Type allocations
  IF (ALLOCATED(lsfinfo%lsf)) DEALLOCATE(lsfinfo%lsf)
  ALLOCATE(lsfinfo%lsf(nspec), stat=stat)
  lsfinfo%lsf = 0.0
  
  ! Check allocations
  IF (stat /= 0) THEN
      WRITE(*,*) 'SPS_SETUP ERROR: Allocation failed.'
      STOP
  ENDIF

  ! Initialize
  wmbsi = 0.
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
  n_isoc        = 0
  m             = 1

  IF (basel_str.NE.'pdva'.AND.basel_str.NE.'wlbc') THEN
     WRITE(*,*) 'SPS_SETUP ERROR: basel_str var set to invalid type: ',basel_str
     STOP
  ENDIF

  !----------------------------------------------------------------!
  !----------------Read in metallicity values----------------------!
  !----------------------------------------------------------------!

  ! Re-open zlegend file
  file_path = TRIM(SPS_HOME)//'/ISOCHRONES/'//TRIM(isoc_dir)//'/zlegend.dat'
  OPEN(90,FILE=file_path,STATUS='OLD',iostat=stat,ACTION='READ')
  
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: zlegend.dat cannot be opened'
     STOP
  ENDIF

  IF (TRIM(isoc_type).EQ.'mist') THEN
     DO z=1,nz
        READ(90,'(A5)') zlegend_str(z)
        zstype5 = zlegend_str(z)
        READ(zstype5(2:5),'(F4.2)') zlegend(z)
        IF (zstype5(1:1).EQ.'m') THEN
           zlegend(z) = 10**(-1*zlegend(z)) * zsol
        ELSE
           zlegend(z) = 10**(1*zlegend(z)) * zsol
        ENDIF
     ENDDO
  ELSE
     DO z=1,nz
        READ(90,'(F6.4)') zlegend(z)
     ENDDO
  ENDIF

  CLOSE(90)

  ! --- DEBUG: Verify Metallicity Parsing ---
  PRINT *, "DEBUG: zsol =", zsol
  PRINT *, "DEBUG: zlegend(1) =", zlegend(1)
  PRINT *, "DEBUG: zlegend_str(1) =", zlegend_str(1)
  PRINT *, "DEBUG: zlegendinit(1) =", zlegendinit(1)
  ! -----------------------------------------

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

  IF (TRIM(isoc_type).EQ.'bpss') THEN

     IF (time_res_incr.NE.1) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: cannot have time_res_incr>1 w/ BPASS models'
        STOP
     ENDIF

     OPEN(91,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: wavelength grid cannot be opened'
        STOP
     ENDIF
     DO i=1,nspec
        READ(91,*) spec_lambda(i)
     ENDDO
     CLOSE(91)

     OPEN(92,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass.mass',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: wavelength grid cannot be opened'
        STOP
     ENDIF
     DO i=1,nt
        READ(92,*) time_full(i),bpass_mass_ssp(i,:)
     ENDDO
     CLOSE(92)

     OPEN(93,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass_v2.2_salpeter100'&
          //'.ssp.bin',FORM='UNFORMATTED',&
          STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
          recl=nspec*nt*nz*8)
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: failed to find bpass ssp models'
        STOP
     ENDIF

     READ(93,rec=1) bpass_spec_ssp
     CLOSE(93)

  ENDIF

  !----------------------------------------------------------------!
  !-----------------Read in spectral libraries---------------------!
  !----------------------------------------------------------------!

  IF (TRIM(isoc_type).NE.'bpss') THEN

  !read in wavelength array and spectral metallicity grid
  IF (TRIM(spec_type).EQ.'basel') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/zlegend.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel.res',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (TRIM(spec_type).EQ.'miles') THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/miles.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/zlegend.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/miles.res',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE IF (INDEX(TRIM(spec_type), 'c3k') > 0) THEN
     OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/'//TRIM(spec_type)//'.lambda',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/zlegend.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
     OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/'//TRIM(spec_type)//'.res',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: wavelength grid cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec
     READ(91,*) spec_lambda(i)
  ENDDO
  CLOSE(91)

  !read in spectral resolution
  DO i=1,nspec
     READ(94,*) spec_res(i)
  ENDDO
  CLOSE(94)

  !read in primary logg and logt arrays
  !NB: these are the same for all spectral libraries
  OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_logt.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logt
     READ(91,*) speclib_logt(i)
  ENDDO
  CLOSE(91)
  OPEN(91,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_logg.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,ndim_logg
     READ(91,*) speclib_logg(i)
  ENDDO
  CLOSE(91)

  !read in each metallicity
  DO z=1,nzinit

     READ(93,*) zlegendinit(z)
     WRITE(zstype,'(F6.4)') zlegendinit(z)

     !read in the spectral library
     IF (TRIM(spec_type).EQ.'basel') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'//basel_str//&
             '_z'//zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ELSE IF (TRIM(spec_type).EQ.'miles') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
        ! --- DEBUG PRINT ---
        IF (z == 1) THEN
           PRINT *, "DEBUG: Opening MILES binary. Z=", zstype
           PRINT *, "DEBUG: nspec =", nspec
           PRINT *, "DEBUG: ndim_logt =", ndim_logt
           PRINT *, "DEBUG: ndim_logg =", ndim_logg
           PRINT *, "DEBUG: RECL =", nspec*ndim_logg*ndim_logt*4
        ENDIF
        ! -------------------
     ELSE IF (INDEX(TRIM(spec_type), 'c3k') > 0) THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/C3K/'//spec_type//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='OLD',iostat=stat,ACTION='READ',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
     ENDIF
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: Library cannot be opened'
        STOP
     ENDIF
     
     ! --- FIX: Read into contiguous buffer first ---
     ! Allocate buffer for exactly ONE metallicity (contiguous memory)
     IF (ALLOCATED(spec_buffer)) DEALLOCATE(spec_buffer)
     ALLOCATE(spec_buffer(nspec, ndim_logt, ndim_logg))
     
     ! Read into the buffer (This is safe because it matches the file struct)
     READ(92, rec=1, IOSTAT=stat2) spec_buffer
     
     IF (stat2 /= 0) THEN
         WRITE(*,*) 'SPS_SETUP ERROR: Read failed for Z index', z, 'IOSTAT=', stat2
         STOP
     ENDIF
     
     ! Copy buffer into the fragmented main array
     speclibinit(:,z,:,:) = spec_buffer
     
     DEALLOCATE(spec_buffer)
     ! --- DEBUG PRINT ---
     IF (stat2 /= 0) THEN
         PRINT *, "DEBUG: READ FAILED with IOSTAT =", stat2
     ELSE IF (z == 1) THEN
         PRINT *, "DEBUG: Read Successful."
         PRINT *, "DEBUG: speclibinit(1,1,1,1) =", speclibinit(1,z,1,1)
         PRINT *, "DEBUG: speclibinit(nspec,1,1,1) =", speclibinit(nspec,z,1,1)
         PRINT *, "DEBUG: Sum of spectrum =", SUM(speclibinit(:,z,1,1))
     ENDIF
     ! -------------------
     CLOSE(92)

  ENDDO

  CLOSE(93)

  !interpolate the input spectral library to the isochrone grid
  !notice that we're interpolating at fixed Z/Zsol even in cases
  !where the isochrones and spectra might have different Zsol. This might
  !in fact be the best thing to do.  Either way, its not ideal.
  DO z=1,nz

     i1 = MIN(MAX(locate(LOG10(zlegendinit/zsol_spec),&
          LOG10(zlegend(z)/zsol)),1),nzinit-1)
     dz = (LOG10(zlegend(z)/zsol)-LOG10(zlegendinit(i1)/zsol_spec)) / &
          (LOG10(zlegendinit(i1+1)/zsol_spec)-LOG10(zlegendinit(i1)/zsol_spec))
     dz = MIN(MAX(dz,0.0),1.0) !no extrapolation!

     speclib(:,z,:,:) = (1-dz)*LOG10(speclibinit(:,i1,:,:)+tiny_number) + &
          dz*LOG10(speclibinit(:,i1+1,:,:)+tiny_number)
     speclib(:,z,:,:) = 10**speclib(:,z,:,:)

  ENDDO

  DEALLOCATE(speclibinit)

  !--------------Read WMBasic Grid from JJ Eldridge----------------;

  !read in Teff array
  OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/WMBASIC.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/Hot_spectra/'//&
          'WMBASIC.teff cannot be opened'
     STOP
  ENDIF
  DO i=1,ndim_wmb_logt
     READ(93,*) wmb_logt(i)
  ENDDO
  CLOSE(93)

  !logg for WMB grid
  wmb_logg = (/3.5,4.0,4.5/)

  OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/WMBASIC_zlegend.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')

  DO z=1,nzwmb

     READ(93,*) zwmb(z)
     WRITE(zstype,'(F6.4)') zwmb(z)

     OPEN(95,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/WMBASIC_z'//&
          zstype//'.spec',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: /Hot_spectra/'//&
          'WMBASIC_z'//zstype//'.spec '//'cannot be opened'
        STOP
     ENDIF
     DO i=1,nspec_wmb
        READ(95,*) wmb_lam(i),wmb_specinit(i,:,1),wmb_specinit(i,:,2),&
             wmb_specinit(i,:,3)
     ENDDO
     CLOSE(95)

     !interpolate to the main spectral grid
     !NB: should be smoothing the models first to the resolution of the
     !input spectral grid
     DO i=1,ndim_wmb_logt
        DO j=1,ndim_wmb_logg
           wmbsi(:,z,i,j) = MAX(linterparr(wmb_lam,wmb_specinit(:,i,j),&
                spec_lambda),tiny_number)
        ENDDO
     ENDDO

  ENDDO

  CLOSE(93)

  !Now interpolate the input spectral library to the isochrone grid
  !notice that we're interpolating at fixed Z/Zsol even in cases
  !where the isochrones and spectra might have different Zsol. This might
  !in fact be the best thing to do.  Either way, its not ideal.
  DO z=1,nz

     i1 = MIN(MAX(locate(LOG10(zwmb/zsol_spec),&
          LOG10(zlegend(z)/zsol)),1),nzwmb-1)
     dz = (LOG10(zlegend(z)/zsol)-LOG10(zwmb(i1)/zsol_spec)) / &
          (LOG10(zwmb(i1+1)/zsol_spec)-LOG10(zwmb(i1)/zsol_spec))
     dz = MIN(MAX(dz,0.0),1.0) !no extrapolation!

     wmb_spec(:,z,:,:) = (1-dz)*LOG10(wmbsi(:,i1,:,:)+tiny_number) + &
          dz*LOG10(wmbsi(:,i1+1,:,:)+tiny_number)
     wmb_spec(:,z,:,:) = 10**wmb_spec(:,z,:,:)

  ENDDO

  !-----------Read in TP-AGB Library from Lancon & Wood------------;

  !read in AGB Teff array for O-rich spectra
  OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Orich.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/AGB_spectra/'//&
          'Orich.teff cannot be opened'
     STOP
  ENDIF
  !burn the header
  READ(93,*) char
  READ(93,*) dumr1, tagb_logz_o
  DO i=1,n_agb_o
     READ(93,*) dumr1, tagb_logt_o(:,i)
  ENDDO
  CLOSE(93)

  !now interpolate the master Teff array to the particular Z array
  DO i=1,nz
     i1 = MIN(MAX(locate(tagb_logz_o,LOG10(zlegend(i)/zsol_spec)),1),22-1)
     dz = (LOG10(zlegend(i)/zsol_spec)-tagb_logz_o(i1)) / &
          (tagb_logz_o(i1+1)-tagb_logz_o(i1))
     agb_logt_o(i,:) = (1-dz)*tagb_logt_o(i1,:)+dz*tagb_logt_o(i1+1,:)

  ENDDO
  agb_logt_o = LOG10(agb_logt_o)

  !read in AGB Teff array for C-rich spectra
  OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Crich.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/AGB_spectra/'//&
          'Crich.teff cannot be opened'
     STOP
  ENDIF
  !burn the header
  READ(94,*) char
  DO i=1,n_agb_c
     READ(94,*) dumr1, agb_logt_c(i)
  ENDDO
  CLOSE(94)
  agb_logt_c = LOG10(agb_logt_c)

  !read in TP-AGB O-rich spectra
  OPEN(95,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Orich.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /AGB_spectra/'//&
          'Orich.spec '//'cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec_agb
     READ(95,*) agb_lam(i),agb_specinit_o(i,:)
  ENDDO
  CLOSE(95)
  !interpolate to the main spectral grid
  DO i=1,n_agb_o
     agb_spec_o(:,i) = MAX(linterparr(agb_lam,agb_specinit_o(:,i),&
          spec_lambda),tiny_number)
  ENDDO

  !read in TP-AGB C-rich spectra
  OPEN(96,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Crich.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /AGB_spectra/'//&
          'Crich.spec '//'cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec_agb
     READ(96,*) agb_lam(i),agb_specinit_c(i,:)
  ENDDO
  CLOSE(96)
  !interpolate to the main spectral grid
  DO i=1,n_agb_c
     agb_spec_c(:,i) = MAX(linterparr(agb_lam,agb_specinit_c(:,i),&
          spec_lambda),tiny_number)
  ENDDO

  ! --- DEBUG: AGB Arrays ---
  PRINT *, "DEBUG: AGB_SPEC_O (1,1) (UV) =", agb_spec_o(1,1)
  PRINT *, "DEBUG: AGB_SPEC_O Sum =", SUM(agb_spec_o)
  PRINT *, "DEBUG: AGB_SPEC_C Sum =", SUM(agb_spec_c)
  ! -------------------------

  !---------Read in Aringer carbon star library---------!

  !read in Aringer C-rich Teff grid
  OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/AGB_spectra/Crich_Aringer.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/AGB_spectra/'//&
          'Crich_Aringer.teff cannot be opened'
     STOP
  ENDIF
  !burn the header
  DO i=1,n_agb_car
     READ(94,*) agb_logt_car(i)
  ENDDO
  CLOSE(94)
  agb_logt_car = LOG10(agb_logt_car)

  !read in Aringer C-rich spectra
  OPEN(96,FILE=TRIM(SPS_HOME)//&
       '/SPECTRA/AGB_spectra/Crich_Aringer.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /AGB_spectra/'//&
          'Crich_Aringer.spec '//'cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec_aringer
     READ(96,*) aringer_lam(i),aringer_specinit(i,:)
  ENDDO
  CLOSE(96)
  !interpolate to the main spectral grid
  DO i=1,n_agb_car
     agb_spec_car(:,i) = MAX(linterparr(aringer_lam,aringer_specinit(:,i),&
          spec_lambda),tiny_number)
  ENDDO

  !------------read in post-AGB spectra from Rauch 2003------------;

  !read in post-AGB Teff array
  OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/ipagb.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/ipagb.teff cannot be opened'
     STOP
  ENDIF
  DO i=1,ndim_pagb
     READ(94,*) pagb_logt(i)
  ENDDO
  CLOSE(94)
  pagb_logt = LOG10(pagb_logt)

  !read in solar metallicity post-AGB spectra
  OPEN(97,FILE=TRIM(SPS_HOME)//'&
       /SPECTRA/Hot_spectra/ipagb_solar.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/Hot_spectra/'//&
          'ipagb.spec_solar cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec_pagb
     READ(97,*) pagb_lam(i),pagb_specinit(i,:,2)
  ENDDO
  CLOSE(97)

  !read in halo metallicity post-AGB spectra
  OPEN(97,FILE=TRIM(SPS_HOME)//&
       '/SPECTRA/Hot_spectra/ipagb_halo.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: /SPECTRA/Hot_spectra/'//&
          'ipagb.spec_halo cannot be opened'
     STOP
  ENDIF
  DO i=1,nspec_pagb
     READ(97,*) pagb_lam(i),pagb_specinit(i,:,1)
  ENDDO
  CLOSE(97)

  !interpolate to the main spectral array
  DO j=1,2
     DO i=1,ndim_pagb
        pagb_spec(:,i,j) = MAX(linterparr(pagb_lam,pagb_specinit(:,i,j),&
             spec_lambda),tiny_number)
     ENDDO
  ENDDO

  !--------------read in WR spectra from Smith et al.--------------;

  !read in WR-N Teff array
  OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WN.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/CMFGEN_WN.teff cannot be opened'
     STOP
  ENDIF
  DO i=1,ndim_wr
     READ(94,*) wrn_logt(i)
  ENDDO
  CLOSE(94)

  !read in WR-C Teff array
  OPEN(94,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WC.teff',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/CMFGEN_WC.teff cannot be opened'
     STOP
  ENDIF
  DO i=1,ndim_wr
     READ(94,*) wrc_logt(i)
  ENDDO
  CLOSE(94)

  !read in WR-N spectra
  OPEN(97,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WN_Zall'//&
       '.spec',STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/CMFGEN_WN_*.spec '//&
          'cannot be opened'
     STOP
  ENDIF
  READ(97,*) tlamwr
  DO j=1,5
     DO i=1,ndim_wr
        READ(97,*) d1,twrzmet(j)
        READ(97,*) twrn(:,i,j)
     ENDDO
  ENDDO
  CLOSE(97)
  twrzmet = LOG10(twrzmet/zsol_spec)

  !interpolate to the main array
  DO j=1,nz
     i1 = MIN(MAX(locate(twrzmet,LOG10(zlegend(j)/zsol_spec)),1),SIZE(twrzmet)-1)
     dz = (LOG10(zlegend(j)/zsol_spec)-twrzmet(i1))/(twrzmet(i1+1)-twrzmet(i1))
     dz = MIN(MAX(dz,0.0),1.)
     DO i=1,ndim_wr
        tspecwr = (1-dz)*LOG10(twrn(:,i,i1)+tiny_number) + &
             dz*LOG10(twrn(:,i,i1+1)+tiny_number)
        wrn_spec(:,i,j) = 10**linterparr(LOG10(tlamwr),tspecwr,&
             LOG10(spec_lambda))-tiny_number
     ENDDO
  ENDDO

  !read in WR-C spectra
  OPEN(97,FILE=TRIM(SPS_HOME)//'/SPECTRA/Hot_spectra/CMFGEN_WC_Zall'//&
       '.spec',STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: Hot_spectra/CMFGEN_WC_*.spec '//&
          'cannot be opened'
     STOP
  ENDIF
  READ(97,*) tlamwr
  DO j=1,5
     DO i=1,ndim_wr
        READ(97,*) d1,twrzmet(j)
        READ(97,*) twrc(:,i,j)
     ENDDO
  ENDDO
  CLOSE(97)
  twrzmet = LOG10(twrzmet/zsol_spec)

  !interpolate to the main array
  DO j=1,nz
     i1 = MIN(MAX(locate(twrzmet,LOG10(zlegend(j)/zsol_spec)),1),SIZE(twrzmet)-1)
     dz = (LOG10(zlegend(j)/zsol_spec)-twrzmet(i1))/(twrzmet(i1+1)-twrzmet(i1))
     dz = MIN(MAX(dz,0.0),1.)
     DO i=1,ndim_wr
        tspecwr = (1-dz)*LOG10(twrc(:,i,i1)+tiny_number) + &
             dz*LOG10(twrc(:,i,i1+1)+tiny_number)
        wrc_spec(:,i,j) = 10**linterparr(LOG10(tlamwr),tspecwr,&
             LOG10(spec_lambda))-tiny_number
     ENDDO
  ENDDO

  !----------------------------------------------------------------!
  !--------------------Read in isochrones--------------------------!
  !----------------------------------------------------------------!

  !read in all metallicities
  DO z=zmin,zmax

     n_isoc = 0
     WRITE(zstype,'(F6.4)') zlegend(z)

     !open Padova isochrones
     IF (TRIM(isoc_type).EQ.'pdva') OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/Padova/Padova2007/isoc_z'//&
          zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     !open PARSEC isochrones
     IF (TRIM(isoc_type).EQ.'prsc') OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/PARSEC/isoc_z'//&
          zstype//'.dat',STATUS='OLD', IOSTAT=stat,ACTION='READ')
     !open MIST isochrones
     IF (TRIM(isoc_type).EQ.'mist') OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/MIST/isoc_z'//zlegend_str(z)//'.dat',STATUS='OLD',&
          IOSTAT=stat,ACTION='READ')
     !open BaSTI isochrones
     IF (TRIM(isoc_type).EQ.'bsti') OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/BaSTI/isoc_z'//zstype//'.dat',STATUS='OLD',&
          IOSTAT=stat,ACTION='READ')
     !open Geneva isochrones
     IF (TRIM(isoc_type).EQ.'gnva') OPEN(97,FILE=TRIM(SPS_HOME)//&
          '/ISOCHRONES/Geneva/isoc_z'//zstype//'.dat',STATUS='OLD',&
          IOSTAT=stat,ACTION='READ')

     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: isochrone files cannot be opened'
        STOP
     END IF

     DO i=1,nlines

        READ(97,*,IOSTAT=stat) char
        IF (stat.NE.0) GOTO 20

        IF (char.EQ.'#') THEN
           m = 1
        ELSE

           IF (m.EQ.1) n_isoc = n_isoc+1
           BACKSPACE(97)
           IF (m.GT.nm) THEN
              WRITE(*,*) 'SPS_SETUP ERROR: number of mass points GT nm'
              STOP
           ENDIF
           IF (TRIM(isoc_type).EQ.'mist') THEN
              READ(97,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                   mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                   logt_isoc(z,n_isoc,m),logg_isoc(z,n_isoc,m),&
                   ffco_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m),&
                   lmdot_isoc(z,n_isoc,m)
           ELSE
              READ(97,*,IOSTAT=stat) logage,mini_isoc(z,n_isoc,m),&
                   mact_isoc(z,n_isoc,m),logl_isoc(z,n_isoc,m),&
                   logt_isoc(z,n_isoc,m),logg_isoc(z,n_isoc,m),&
                   ffco_isoc(z,n_isoc,m),phase_isoc(z,n_isoc,m)
              lmdot_isoc(z,n_isoc,m)=-99.
           ENDIF
           IF (stat.NE.0) GOTO 20
           IF (m.EQ.1) timestep_isoc(z,n_isoc) = logage
           nmass_isoc(z,n_isoc) = nmass_isoc(z,n_isoc)+1
           m = m+1

        ENDIF

     ENDDO

     WRITE(*,*) 'SPS_SETUP ERROR: didnt finish reading in the isochrones!'
     STOP

20   CONTINUE
     CLOSE(97)

     ! For dynamically determined nt, we should verify, but here nt is set by reading the first file.
     ! If subsequent files have different number of ages, that's an error.
     IF (n_isoc.NE.nt) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: number of isochrones NE nt',n_isoc,nt
        STOP
     ENDIF

  ENDDO

  !this is necessary because Geneva does not extend below 1.0 Msun
  !see imf_weight.f90 for details
  IF (TRIM(isoc_type).EQ.'gnva') THEN
     imf_lower_bound = MINVAL(mini_isoc(zmin,1,1:nmass_isoc(zmin,1)))*0.99
  ELSE
     imf_lower_bound = imf_lower_limit
  ENDIF

  ENDIF

  !----------------------------------------------------------------!
  !--------Read in dust emission spectra from Draine & Li----------!
  !----------------------------------------------------------------!

  DO k=1,nqpah_dustem
     WRITE(sqpah,'(I1)') k-1
     IF (k-1.EQ.10) THEN
        OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dustem/'//TRIM(str_dustem)//&
             '_MW3.1_100.dat',STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE
        OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dustem/'//TRIM(str_dustem)//&
             '_MW3.1_'//sqpah//'0.dat',STATUS='OLD',iostat=stat,ACTION='READ')
     ENDIF
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: error opening dust emission file'
        STOP
     ENDIF
     DO i=1,2  !burn the header
        READ(99,*)
     ENDDO
     DO i=1,ndim_dustem
        READ(99,*,IOSTAT=stat) lambda_dustem(i),dustem_dustem(i,:)
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR: error during dust emission read'
           STOP
        ENDIF
     ENDDO
     CLOSE(99)
     lambda_dustem = lambda_dustem*1E4  !convert to Ang

     !now interpolate the dust spectra onto the master wavelength array
     DO j=1,numin_dustem*2
        !the dust models only extend to 1um
        jj = locate(spec_lambda/1E4,one)
        dustem2_dustem(jj:,k,j) = linterparr(lambda_dustem,&
             dustem_dustem(:,j),spec_lambda(jj:))
     ENDDO

  ENDDO

  !----------------------------------------------------------------!
  !-------------Read in circumstellar AGB dust models--------------!
  !----------------------------------------------------------------!

  !O-rich spectra
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dusty/Orich_dusty.spec',&
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

  DO i=1,nteff_dagb
     DO j=1,ntau_dagb
        READ(99,*,IOSTAT=stat) teff_dagb(1,i), tau1_dagb(1,j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
        jj = locate(spec_lambda,lambda_dagb(1))
        flux_dagb(jj:,1,i,j) = linterparr(lambda_dagb(1:nlam),&
             fluxin_dagb(1:nlam),spec_lambda(jj:))
     ENDDO
  ENDDO
  CLOSE(99)

  ! --- DEBUG: Dusty AGB Arrays ---
  PRINT *, "DEBUG: FLUX_DAGB Sum =", SUM(flux_dagb)
  ! -------------------------------

  !C-rich spectra
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/dusty/Crich_dusty.spec',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: error opening dusty models'
     STOP
  ENDIF

  !number of wavelength points in the AGB grid
  READ(99,*) nlam
  !read in the wavelength grid
  READ(99,*,IOSTAT=stat) lambda_dagb(1:nlam)

  DO i=1,nteff_dagb
     DO j=1,ntau_dagb
        READ(99,*,IOSTAT=stat) teff_dagb(2,i), tau1_dagb(2,j)
        READ(99,*,IOSTAT=stat) fluxin_dagb(1:nlam)
        IF (stat.NE.0) THEN
           WRITE(*,*) 'SPS_SETUP ERROR: error reading dusty models'
           STOP
        ENDIF
        !interpolate the dust spectra onto the master wavelength array
        jj = locate(spec_lambda,lambda_dagb(1))
        flux_dagb(jj:,2,i,j) = linterparr(lambda_dagb(1:nlam),&
             fluxin_dagb(1:nlam),spec_lambda(jj:))
     ENDDO
  ENDDO
  CLOSE(99)

  !----------------------------------------------------------------!
  !--------------------Set up AGN dust model-----------------------!
  !----------------------------------------------------------------!

  !models from Nenkova et al. 2008

  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/Nenkova08_y010_torusg_n10_q2.0.dat',&
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

  DO i=1,nagndust_spec
     READ(99,*) agndust_lam(i),agndust_specinit(i,:)
  ENDDO

  i1 = locate(spec_lambda,agndust_lam(1))
  i2 = locate(spec_lambda,agndust_lam(nagndust_spec))
  DO i=1,nagndust
     agndust_spec(i1:i2,i) = 10**linterparr(LOG10(agndust_lam),&
          LOG10(agndust_specinit(:,i)+tiny30),LOG10(spec_lambda(i1:i2)))-tiny30
  ENDDO


  !----------------------------------------------------------------!
  !----------------Set up nebular emission arrays------------------!
  !----------------------------------------------------------------!

  IF (TRIM(isoc_type).EQ.'mist'.OR.TRIM(isoc_type).EQ.'pdva'.OR.&
     TRIM(isoc_type).EQ.'prsc'.OR.TRIM(isoc_type).EQ.'bpss') THEN

     !read in nebular continuum arrays.  Units are Lsun/Hz/Q
     IF (cloudy_dust.EQ.1) THEN
        OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WD_'//TRIM(isoc_type)//'.cont',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE
        OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_ND_'//TRIM(isoc_type)//'.cont',&
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
     DO i=1,nebnz
        DO j=1,nebnage
           DO k=1,nebnip
              READ(99,*,iostat=stat) nebem_logz(i),nebem_age(j),nebem_logu(k)
              READ(99,*,iostat=stat) readcontneb
              !interpolate onto the main wavelength grid
              !some values in the table are 0.0, set a floor of 1E-95
              nebem_cont(:,i,j,k) = linterparr(readlambneb,&
                   LOG10(readcontneb+10**(-95.d0)),spec_lambda)
           ENDDO
        ENDDO
     ENDDO
     CLOSE(99)

     !read in nebular emission line luminosities.  Units are Lsun/Q
     IF (cloudy_dust.EQ.1) THEN
        OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WD_'//TRIM(isoc_type)//'.lines',&
             STATUS='OLD',iostat=stat,ACTION='READ')
     ELSE
        OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_ND_'//TRIM(isoc_type)//'.lines',&
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
     DO i=1,nebnz
        DO j=1,nebnage
           DO k=1,nebnip
              READ(99,*,iostat=stat) nebem_logz(i),nebem_age(j),nebem_logu(k)
              READ(99,*,iostat=stat) nebem_line(:,i,j,k)
           ENDDO
        ENDDO
     ENDDO
     CLOSE(99)

     !convert the nebem_age array to log(age), and log the emission arrays
     nebem_age  = LOG10(nebem_age)
     nebem_line = LOG10(nebem_line + 10**(-95.d0))

     !define the minimum resolution of the emission lines
     !based on the resolution of the spectral library
     DO i=1,nemline
        j = MIN(MAX(locate(spec_lambda,nebem_line_pos(i)),1),nspec-1)
        neb_res_min(i) = spec_lambda(j+1)-spec_lambda(j)
     ENDDO

     !set up a "master" array of normalized Gaussians
     !this makes the code much faster
     IF (setup_nebular_gaussians.EQ.1) THEN
        DO i=1,nemline
           IF (smooth_velocity.EQ.1) THEN
              !smoothing variable is km/s
              dlam = nebem_line_pos(i)*nebular_smooth_init/clight*1E13
           ELSE
              !smoothing variable is A
              dlam = nebular_smooth_init
           ENDIF
           !broaden the line to at least the resolution element
           !of the spectrum (x2).
           dlam = MAX(dlam,neb_res_min(i)*2)
           gaussnebarr(:,i) = 1/SQRT(2*mypi)/dlam*&
                EXP(-(spec_lambda-nebem_line_pos(i))**2/2/dlam**2)  / &
                clight*nebem_line_pos(i)**2
        ENDDO
     ENDIF

  ENDIF

  !----------------------------------------------------------------!
  !------------------Set up X-ray nebular --------------------!
  !----------------------------------------------------------------!

  IF (TRIM(isoc_type).EQ.'bpss') THEN
      !read in nebular continuum arrays.  Units are Lsun/Hz/Q
      IF (cloudy_dust.EQ.1) THEN
         OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WX_WD_'//TRIM(isoc_type)//'.cont',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ELSE
         OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WX_ND_'//TRIM(isoc_type)//'.cont',&
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
      DO i=1,nebnz
         DO j=1,nebnage
            DO k=1,nebnip
               READ(99,*,iostat=stat) nebem_logz(i),nebem_age(j),nebem_logu(k)
               READ(99,*,iostat=stat) readcontneb
               !interpolate onto the main wavelength grid
               !some values in the table are 0.0, set a floor of 1E-95
               xnebem_cont(:,i,j,k) = linterparr(readlambneb,&
                     LOG10(readcontneb+10**(-95.d0)),spec_lambda)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(99)

      !read in nebular emission line luminosities.  Units are Lsun/Q
      IF (cloudy_dust.EQ.1) THEN
         OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WX_WD_'//TRIM(isoc_type)//'.lines',&
               STATUS='OLD',iostat=stat,ACTION='READ')
      ELSE
         OPEN(99,FILE=TRIM(SPS_HOME)//'/nebular/ZAU_WX_ND_'//TRIM(isoc_type)//'.lines',&
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
      DO i=1,nebnz
         DO j=1,nebnage
            DO k=1,nebnip
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

  OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/xrb/xsp.lambda',&
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
     OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/xrb/xsp_feh'//zz_str_xrb(j)&
          //'.spec',STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: xsp_feh'//zz_str_xrb(j)//&
             '.spec cannot be opened'
        STOP
     ENDIF
     DO i=1,nt_xrb
        READ(98,*) tspec_xrb
        !interpolate to the main wavelength array
        spec_xrb(:,i,j) = MAX(linterparr(lam_xrb,tspec_xrb,spec_lambda),tiny_number)
     ENDDO
     CLOSE(98)
  ENDDO

  !convert to Lsun/Hz/Msun
  spec_xrb = spec_xrb * lsun

  !----------------------------------------------------------------!
  !-------------------Set up magnitude info------------------------!
  !----------------------------------------------------------------!

  !read in Vega-like star (lambda, Flambda)
  !(this is actually a Kurucz (1992) model for Vega)
  OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/A0V_KURUCZ_92.SED',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: SPECTRA/A0V_KURUCZ_92.SED cannot be opened'
     STOP
  ENDIF
  !burn the header
  READ(98,*)
  DO i=1,ntlam
     READ(98,*) tvega_lam(i), tvega_spec(i)
  ENDDO
  CLOSE(98)

  !interpolate the Vega spectrum onto the wavelength grid
  !and convert to fnu
  jj = locate(vega_spec,tvega_lam(ntlam))
  vega_spec(:jj) = 10**linterparr(LOG10(tvega_lam),&
          LOG10(tvega_spec+tiny_number),LOG10(spec_lambda(:jj)))
  vega_spec = vega_spec*spec_lambda**2
  vega_spec(jj+1:) = tiny_number

  !read in Solar spectrum; units are fnu, flux is appropriate for
  !deriving absolute magnitudes.  spectrum from STScI, extrapolated
  !beyond 2.5um with a blackbody.
  OPEN(98,FILE=TRIM(SPS_HOME)//'/SPECTRA/SUN_STScI.SED',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  IF (stat.NE.0) THEN
     WRITE(*,*) 'SPS_SETUP ERROR: SPECTRA/SUN_STScI.SED cannot be opened'
     STOP
  END IF
  DO i=1,ntlam
     READ(98,*) tsun_lam(i), tsun_spec(i)
  ENDDO
  CLOSE(98)

  !interpolate the Solar spectrum onto the wavelength grid
  jj = locate(sun_spec,tsun_lam(ntlam))
  sun_spec(:jj) = 10**linterparr(LOG10(tsun_lam),&
          LOG10(tsun_spec+tiny_number),LOG10(spec_lambda(:jj)))
  sun_spec(jj+1:) = tiny_number

  !read in and set up band-pass filters
  IF (TRIM(alt_filter_file).EQ.'') THEN
     OPEN(99,FILE=TRIM(SPS_HOME)//'/data/allfilters.dat',&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ELSE
     OPEN(99,FILE=TRIM(SPS_HOME)//'/data/'//TRIM(alt_filter_file),&
          STATUS='OLD',iostat=stat,ACTION='READ')
  ENDIF
  READ(99,*)


  !loop over all the transmission filters
  DO i=1,nbands

     jj=0
     readlamb = 0.0
     readband = 0.0
     DO j=1,50000
        READ(99,*,iostat=stat) d1,d2
        IF (stat.NE.0) GOTO 909
        IF (jj.EQ.0) THEN
           jj = jj+1
           readlamb(jj) = d1
           !force the transmission to be GE 0
           readband(jj) = MAX(d2,0.0)
        ELSE
           !only read unique lambda points
           IF (readlamb(jj).NE.d1) THEN
              jj = jj+1
              readlamb(jj) = d1
              !force the transmission to be GE 0
              readband(jj) = MAX(d2,0.0)
           ENDIF
        ENDIF
     ENDDO

909  CONTINUE

     IF (j.GE.50000) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: did not finish reading in filter ',i
        STOP
     ENDIF
     IF (jj.EQ.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: error during filter definition read-in',i
        STOP
     ENDIF

     !interpolate the filter onto the master wavelength array
     i1 = MAX(locate(spec_lambda,readlamb(1)),1)
     i2 = locate(spec_lambda,readlamb(jj))
     IF (i1.NE.i2) bands(i1:i2,i) = &
          linterparr(readlamb(1:jj),readband(1:jj),spec_lambda(i1:i2))

     !normalize
     dumr1 = TSUM(spec_lambda,bands(:,i)/spec_lambda)
     !in this case the band is entirely outside the wavelength array
     IF (dumr1.LE.tiny_number) dumr1=1.0
     bands(:,i) = bands(:,i) / dumr1
     bands(:,i) = MAX(bands(:,i),0.0)  !force no negative values

     !compute absolute magnitude of the Sun
     magsun(i) = TSUM(spec_lambda,sun_spec*bands(:,i)/spec_lambda)
     IF (magsun(i).LT.2*tiny_number) THEN
        magsun(i) = 99.0
     ELSE
        magsun(i) = -2.5*LOG10(magsun(i)) - 48.60
     ENDIF

     !compute mags of Vega
     magvega(i) = TSUM(spec_lambda,vega_spec*bands(:,i)/spec_lambda)
     IF (magvega(i).LE.tiny_number) THEN
        magvega(i) = 99.0
     ELSE
        magvega(i) = -2.5 * LOG10(magvega(i)) - 48.60
     ENDIF

     !put Sun magnitudes in the Vega system if keyword is set
     IF (compute_vega_mags.EQ.1.AND.magsun(i).NE.99.0) &
          magsun(i) = (magsun(i)-magsun(1)) - &
          (magvega(i)-magvega(1)) + magsun(1)

  ENDDO
  CLOSE(99)

  !only execute this loop for the standard filter list
  IF (TRIM(alt_filter_file).EQ.'') THEN
     !normalize the IRAC, PACS, SPIRE, and IRAS photometry to nu*fnu=const
     !Note: this turns out to be irrelevant and is the result of rather
     !confusing documentation on the IRAC website
     lami = (/3.550,4.493,5.731,7.872,70.0,100.0,160.0,250.0,350.0,500.0,&
          12.0,25.0,60.0,100.0/)*1E4
     ind=(/53,54,55,56,95,96,97,98,99,100,101,102,103,104/)
     DO j=1,14
        IF (ind(j).GT.nbands) THEN
           ! WRITE(*,*) 'SPS_SETUP ERROR: trying to index a filter that does not exist!'
           ! We just skip if it does not exist now that nbands is dynamic
           CYCLE
        ENDIF
        d = TSUM(spec_lambda,(spec_lambda/lami(j))**(-1.0)*bands(:,ind(j))/&
             spec_lambda)
        bands(:,ind(j)) = bands(:,ind(j)) / MAX(d,tiny_number)
     ENDDO

     !normalize the MIPS photometry to a BB (beta=2)
     !this part is *not* irrelevant.
     lami(1:3) = (/23.68,71.42,155.9/)*1E4
     ind(1:3)  = (/90,91,92/)
     DO j=1,3
        IF (ind(j).GT.nbands) THEN
           ! WRITE(*,*) 'SPS_SETUP ERROR: trying to index a filter that does not exist!'
           CYCLE
        ENDIF
        d = TSUM(spec_lambda,(spec_lambda/lami(j))**(-2.0)*bands(:,ind(j))/&
             spec_lambda)
        bands(:,ind(j)) = bands(:,ind(j)) / MAX(d,tiny_number)
     ENDDO
  ENDIF

  !compute the effective wavelength of each filter
  !NB: These are sometimes referred to as "pivot" wavelengths
  ! in the literature.  See Bessell & Murphy 2012 A.2.1 for details'
  DO i=1,nbands
     filter_leff(i) = TSUM(spec_lambda,spec_lambda*bands(:,i)) / &
          TSUM(spec_lambda,bands(:,i)/spec_lambda)
     filter_leff(i) = SQRT(filter_leff(i))
  ENDDO

  !----------------------------------------------------------------!
  !---------------Set up extinction curve indices------------------!
  !----------------------------------------------------------------!

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
  !----------Set up Witt & Gordon 2000 attenuation curves----------!
  !----------------------------------------------------------------!

  !read in the WG00 dust model attenuation curves
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/alldirty_h.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(99,*)
  READ(99,*)
  READ(99,*)
  DO i=1,18
     DO j=1,25
        READ(99,*) wglam(j), d1,wgtmp(j,i,1,:)
     ENDDO
  ENDDO
  CLOSE(99)
  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/alldirty_c.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  READ(99,*)
  READ(99,*)
  READ(99,*)
  DO i=1,18
     DO j=1,25
        READ(99,*) wglam(j), d1,wgtmp(j,i,2,:)
     ENDDO
  ENDDO
  CLOSE(99)

  !interpolate the WG00 models onto the spectral grid
  DO k=1,2
     DO i=1,18
        DO j=1,6
           DO n=1,nspec
              IF (spec_lambda(n).GT.wglam(25)) THEN
                 wgdust(n,i,j,k)=0.0
              ELSE IF (spec_lambda(n).LT.wglam(1)) THEN
                 wgdust(n,i,j,k) = wgtmp(1,i,k,j)
              ELSE
                 wgdust(n,i,j,k) = linterp(wglam,wgtmp(:,i,k,j),&
                      spec_lambda(n))
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !set up Gordon et al. (2003) SMC bar extinction curve

  OPEN(99,FILE=TRIM(SPS_HOME)//'/dust/Gordon03_table4.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,30
     !the data are in reverse wavelength order
     READ(99,*) g03lam(30-i+1),d1,g03smc(30-i+1)
  ENDDO
  CLOSE(99)
  g03lam = g03lam*1E4 !convert from um to A

  DO n=1,nspec
     IF (spec_lambda(n).GT.g03lam(30)) THEN
        g03smcextn(n)=0.0
     ELSE IF (spec_lambda(n).LT.g03lam(1)) THEN
        g03smcextn(n) = g03smc(1)
     ELSE
        g03smcextn(n) = linterp(g03lam,g03smc,spec_lambda(n))
     ENDIF
     !write(34,*) spec_lambda(n),g03smcextn(n)
  ENDDO

  !----------------------------------------------------------------!
  !--------------set up the redshift-age-DL relations--------------!
  !----------------------------------------------------------------!

  DO i=1,500
     a = (i-1)/499.*(1-1/1001.)+1/1001.
     cosmospl(i,1) = 1/a-1  !redshift
     cosmospl(i,2) = get_tuniv(cosmospl(i,1))   ! Tuniv in Gyr
     cosmospl(i,3) = get_lumdist(cosmospl(i,1)) ! Lum Dist in pc
  ENDDO

  !set Tuniv
  tuniv = get_tuniv(zero)

  !----------------------------------------------------------------!
  !-----------------read in index definitions----------------------!
  !----------------------------------------------------------------!

  OPEN(99,FILE=TRIM(SPS_HOME)//'/data/allindices.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
  DO i=1,4  !burn the header
     READ(99,*)
  ENDDO
  DO i=1,nindx
     READ(99,*,IOSTAT=stat) indexdefined(:,i)
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: error during index defintion read'
        STOP
     ENDIF
     !convert the Lick indices from air to vacuum wavelengths
     IF (i.LE.25) THEN
        indexdefined(1:6,i) = airtovac(indexdefined(1:6,i))
     ENDIF
  ENDDO
  CLOSE(99)

  !----------------------------------------------------------------!
  !-----------------set up expanded time array---------------------!
  !----------------------------------------------------------------!

  IF (TRIM(isoc_type).NE.'bpss') THEN

     DO i=1,ntfull
        IF (MOD(i-1,time_res_incr).EQ.0) THEN
           time_full(i) = timestep_isoc(zmin,(i-1)/time_res_incr+1)
        ELSE
           IF ((i-1)/time_res_incr+2.LT.nt) THEN
              d1 = (timestep_isoc(zmin,(i-1)/time_res_incr+2)-&
                   timestep_isoc(zmin,(i-1)/time_res_incr+1))/time_res_incr
           ENDIF
           time_full(i) = timestep_isoc(zmin,(i-1)/time_res_incr+1)+d1
           time_full(i) = time_full(i-1)+d1
        ENDIF
     ENDDO

  ENDIF

  !----------------------------------------------------------------!
  !------------------------set up the LSF--------------------------!
  !----------------------------------------------------------------!

  IF (smooth_lsf.EQ.1) THEN

     OPEN(99,FILE=TRIM(SPS_HOME)//'/data/lsf.dat',&
       STATUS='OLD',iostat=stat,ACTION='READ')
     IF (stat.NE.0) THEN
        WRITE(*,*) 'SPS_SETUP ERROR: lsf.dat cannot be opened'
        STOP
     ENDIF

     DO i=1,ntabmax
        READ(99,*,iostat=stat) lsflam(i),lsfsig(i)
        IF (stat.NE.0) GOTO 910
        IF (i.EQ.1) lsfinfo%minlam=lsflam(i)
     ENDDO

     WRITE(*,*) 'SPS_SETUP ERROR: read to end of lsf.dat file'
     STOP

910  CONTINUE

     lsfinfo%maxlam = lsflam(i-1)

     !interpolate onto the main wavelength array
     DO n=1,nspec
        IF (spec_lambda(n).GE.lsfinfo%minlam.AND.&
             spec_lambda(n).LE.lsfinfo%maxlam) THEN
           lsfinfo%lsf(n) = linterp(lsflam(1:i-1),lsfsig(1:i-1),&
                spec_lambda(n))
        ENDIF
     ENDDO

  ENDIF


  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  whlam5000 = locate(spec_lambda,5000.d0)
  whlylim   = locate(spec_lambda,912.d0)
  !define the frequency array
  spec_nu   = clight / spec_lambda

  ! --- DEBUG: Wavelength Grid ---
  whlam5000 = locate(spec_lambda,5000.d0)
  PRINT *, "DEBUG: whlam5000 (Index of 5000A) =", whlam5000
  PRINT *, "DEBUG: spec_lambda(1) =", spec_lambda(1)
  PRINT *, "DEBUG: spec_lambda(50) =", spec_lambda(50)
  PRINT *, "DEBUG: spec_lambda(whlam5000) =", spec_lambda(whlam5000)
  ! ------------------------------

  !set flag indicating that sps_setup has been run, initializing
  !important common block vars/arrays
  check_sps_setup = 1

  ! --- DEBUG: Verify Final Spectral Library ---
  PRINT *, "DEBUG: End of Setup."
  PRINT *, "DEBUG: speclib sum =", SUM(speclib)
  PRINT *, "DEBUG: speclib(1,1,1,1) =", speclib(1,1,1,1)
  PRINT *, "DEBUG: speclib(nspec,1,1,1) =", speclib(nspec,1,1,1)
  ! --------------------------------------------

  IF (verbose.EQ.1) THEN
     WRITE(*,*) '      ...done'
     WRITE(*,*)
  ENDIF

END SUBROUTINE SPS_SETUP
