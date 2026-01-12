PROGRAM AUTOSPS

  USE sps_vars
  USE sps_utils
  
  IMPLICIT NONE

  INTEGER :: z
  REAL(SP), ALLOCATABLE, DIMENSION(:,:)  :: spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:)        :: mass_ssp,lbol_ssp
  CHARACTER(100) :: file1='',aux
  CHARACTER(3)  :: str
  TYPE(PARAMS)  :: pset
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: ocompsp
  INTEGER :: alloc_stat

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Initialize FSPS with Padova Isochrones (Required for autosps)
  ! We call this early to populate zlegend and nz dimensions
  WRITE(6,*) 'Initializing SPS with Padova isochrones...'
  CALL SPS_SETUP(-1, input_isoc_type='pdva')

  ! Verify that the initialization succeeded with the correct library
  IF (TRIM(isoc_type).NE.'pdva') THEN
     WRITE(*,*) 'ERROR: autosps only works with the "Padova" isochrones'
     STOP
  ENDIF

  ! Allocate arrays dependent on runtime dimensions
  ALLOCATE(spec_ssp(ntfull,nspec), stat=alloc_stat)
  ALLOCATE(mass_ssp(ntfull), lbol_ssp(ntfull), stat=alloc_stat)
  ALLOCATE(ocompsp(ntfull), stat=alloc_stat)
  IF (alloc_stat /= 0) STOP 'Allocation failed in AUTOSPS'

  ! Allocate PARAMS components
  IF (.NOT. ALLOCATED(pset%mag_compute)) ALLOCATE(pset%mag_compute(nbands))
  IF (.NOT. ALLOCATED(pset%ssp_gen_age)) ALLOCATE(pset%ssp_gen_age(nt))
  pset%mag_compute = 1
  pset%ssp_gen_age = 1
  
  !set SFH
  WRITE(6,*)
  WRITE(6,*)  'Specify SFH [0-2, def:0]'
  WRITE(6,*)  '(0=SSP, 1=CSP, 2=tabulated)'
  READ(5,'(A)')  aux
  IF (len(trim(aux)).EQ.0) THEN
     pset%sfh = 0
  ELSE
     READ(aux,'(I1)') pset%sfh
  ENDIF
  IF (pset%sfh.EQ.0) WRITE(6,'(" ---> Computing an SSP")') 
  IF (pset%sfh.EQ.1) WRITE(6,'(" ---> Computing a CSP")') 
  IF (pset%sfh.EQ.2) WRITE(6,'(" ---> Computing a tabulated SFH")') 

  IF (pset%sfh.EQ.1) THEN
     WRITE(6,*)
     WRITE(6,*) 'input parameters for CSP: tau, const, age, fburst, tburst'
     WRITE(6,*) ' - tau in Gyr' 
     WRITE(6,*) ' - const as fraction of mass formed in constant component'
     WRITE(6,*) ' - age of the system in Gyr.  i.e. results span the time 0<t<age'
     WRITE(6,*) ' - fburst as fraction of mass formed in an instantaneous burst'
     WRITE(6,*) ' - tburst as time of burst, with tburst<age'
     READ(5,*) pset%tau,pset%const,tuniv,pset%fburst,pset%tburst
     pset%tage = 0.0
     IF (pset%const.LT.0.OR.pset%const.GT.1.0) THEN
        WRITE(6,*) 'ERROR, const out of bounds: ',pset%const
        STOP
     ENDIF
     IF (pset%fburst.LT.0.OR.pset%fburst.GT.1.0) THEN
        WRITE(6,*) 'ERROR, fburst out of bounds: ',pset%fburst
        STOP
     ENDIF
     WRITE(6,'(" ---> (tau const age fburst tburst)=(",5(F6.2),")")') &
          pset%tau,pset%const,tuniv,pset%fburst,pset%tburst

  ENDIF

  IF (pset%sfh.NE.2) THEN
     !set metallicity
     WRITE(6,*)
     WRITE(6,*)  'enter metallicity [1-22; def:20]:'
     READ(5,'(A)')  aux
     IF (len(trim(aux)).EQ.0) THEN
        pset%zmet = 20
     ELSE
        READ(aux,'(I2)') pset%zmet
     ENDIF
     IF (pset%zmet.LT.1.OR.pset%zmet.GT.nz) THEN
        WRITE(*,*) 'ERROR: Z out of bounds: ',pset%zmet
        STOP
     ENDIF
     WRITE(6,'(" ---> Using metallicity",1x,I2," corresponding to log(Z/Zsol)=",1x,F5.2)') &
          pset%zmet,LOG10(zlegend(pset%zmet)/zsol)
  ENDIF

  !set dust
  WRITE(6,*)
  WRITE(6,*)  'Include default dust model? [yes/no, def:no]'
  WRITE(6,*)  '(default: tau1=1.0, tau2=0.3, MW extinction)'
  READ(5,'(A)')  aux
  IF (len(trim(aux)).NE.0) THEN
     READ(aux,'(A3)') str
     IF (str(1:1).EQ.'y') THEN
        dust_type  = 1
        pset%dust1 = 1.0
        pset%dust2 = 0.3
     ENDIF
  ENDIF
  WRITE(6,'(" ---> tau1=",1x,F5.2,", tau2=",1x,F5.2)') pset%dust1,pset%dust2


  !set filename
  WRITE(6,*)
  WRITE(6,*)  'Enter filename [def: "CSP.out"]'
  READ(5,'(A)')  aux
  IF (len(trim(aux)).EQ.0) THEN
     file1 = 'CSP.out'
  ELSE
     READ(aux,'(A)') file1
  ENDIF
  WRITE(6,'(" ---> Output filename:",1x,A100)') file1


  WRITE(6,'(" ---> Running model.......")')
  
  IF (pset%sfh.EQ.2) THEN
     DO z=1,nz
        pset%zmet=z
        CALL SSP_GEN(pset,mass_ssp_zz(:,z),lbol_ssp_zz(:,z),spec_ssp_zz(:,:,z))
     ENDDO
     CALL COMPSP(3,nz,file1,mass_ssp_zz,lbol_ssp_zz,spec_ssp_zz,pset,ocompsp)
  ELSE
     CALL SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
     CALL COMPSP(3,1,file1,mass_ssp,lbol_ssp,spec_ssp,pset,ocompsp)
  ENDIF


END PROGRAM AUTOSPS 
