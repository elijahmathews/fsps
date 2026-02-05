PROGRAM AUTOSPS

   USE fsps_precision, ONLY: WP
   USE fsps_types, ONLY: PARAMS, COMPSPOUT
   USE sps_utils
   USE fsps_ssp, ONLY: generate_ssp_grid
   USE fsps_context, ONLY: fsps_context_create
   USE fsps_context_types, ONLY: fsps_context_t
  
  IMPLICIT NONE

   INTEGER :: z
   TYPE(fsps_context_t) :: ctx

   REAL(WP), ALLOCATABLE :: spec_ssp(:,:,:)
   REAL(WP), ALLOCATABLE :: mass_ssp(:,:),lbol_ssp(:,:)
  TYPE(COMPSPOUT), ALLOCATABLE :: ocompsp(:)

  CHARACTER(100) :: file1='',aux
  CHARACTER(3)  :: str
  TYPE(PARAMS)  :: pset
   REAL(WP) :: tuniv

  ! Variables for library selection
  CHARACTER(10) :: iso_in, spec_in

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  WRITE(6,*) '=================================================='
  WRITE(6,*) '             FSPS INTERACTIVE MODE                '
  WRITE(6,*) '=================================================='

  ! --- Ask for libraries at runtime ---
  WRITE(6,*) 'Choose Isochrones [def: mist]:'
  WRITE(6,*) '(mist, pdva, parsec, basti, bpass, etc.)'
  READ(5,'(A)') aux
  IF (LEN_TRIM(aux).EQ.0) THEN
     iso_in = 'mist'
  ELSE
     READ(aux,'(A)') iso_in
  ENDIF
  
  WRITE(6,*) 'Choose Spectral Library [def: miles]:'
  WRITE(6,*) '(miles, basel, bpass, c3k_afe+0.0, etc.)'
  READ(5,'(A)') aux
  IF (LEN_TRIM(aux).EQ.0) THEN
     spec_in = 'miles'
  ELSE
     READ(aux,'(A)') spec_in
  ENDIF

  ! --- Initialize Environment Immediately ---
   WRITE(6,*) 'Initializing libraries...'
   CALL fsps_context_create(ctx)
  ! We load ALL metallicities (-1) so we can query nz and zlegend
   CALL SPS_SETUP(ctx, -1, TRIM(iso_in), TRIM(spec_in))

  ! --- Allocate Memory ---
  IF (.NOT. ALLOCATED(spec_ssp)) THEN
     ALLOCATE(spec_ssp(ctx%state%nspec, ctx%state%ntfull, ctx%state%nz))
     ALLOCATE(mass_ssp(ctx%state%ntfull, ctx%state%nz))
     ALLOCATE(lbol_ssp(ctx%state%ntfull, ctx%state%nz))
     ALLOCATE(ocompsp(ctx%state%ntfull))
  END IF

  !set IMF
  WRITE(6,*)
  WRITE(6,*)  'enter IMF [0-5; def:0]:'
  WRITE(6,*) ' (0=Salpeter, 1=Chabrier 2003, 2=Kroupa 2001, '//&
       '3=van Dokkum 2008, 4=Dave 2008, 5=tabulated)'
  READ(5,'(A)')  aux
  IF (LEN(TRIM(aux)).EQ.0) THEN
     ctx%imf_type_val = 0
  ELSE
     READ(aux,'(I1)') ctx%imf_type_val
  ENDIF
  IF (ctx%imf_type_val.LT.0.OR.ctx%imf_type_val.GT.5) THEN
     WRITE(*,*) 'ERROR: imf out of bounds: ',ctx%imf_type_val
     STOP
  ENDIF
  WRITE(6,'(" ---> Using IMF",1x,I1)') ctx%imf_type_val

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
     ! Dynamic prompt using the loaded nz
     WRITE(6,'("enter metallicity index [1-",I0,"; def:",I0,"]:")') ctx%state%nz, ctx%state%nz
     READ(5,'(A)')  aux
     IF (len(trim(aux)).EQ.0) THEN
        pset%zmet = ctx%state%nz  ! Default to the last index (usually safest/solar-ish)
     ELSE
        READ(aux,'(I2)') pset%zmet
     ENDIF
     IF (pset%zmet.LT.1.OR.pset%zmet.GT.ctx%state%nz) THEN
        WRITE(*,*) 'ERROR: Z out of bounds: ',pset%zmet
        STOP
     ENDIF
     WRITE(6,'(" ---> Using metallicity",1x,I2," corresponding to log(Z/Zsol)=",1x,F5.2)') &
          pset%zmet,LOG10(ctx%state%zlegend(pset%zmet)/ctx%state%zsol)
  ENDIF

  !set dust
  WRITE(6,*)
  WRITE(6,*)  'Include default dust model? [yes/no, def:no]'
  WRITE(6,*)  '(default: tau1=1.0, tau2=0.3, MW extinction)'
  READ(5,'(A)')  aux
  IF (len(trim(aux)).NE.0) THEN
     READ(aux,'(A3)') str
     IF (str(1:1).EQ.'y') THEN
      ctx%dust_type_val  = 1
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
     ! We already called SPS_SETUP(-1) at the top, so variables are ready.
     DO z=1,ctx%state%nz
        pset%zmet=z
      CALL generate_ssp_grid(ctx, pset, mass_ssp(:,z), lbol_ssp(:,z), spec_ssp(:,:,z))
     ENDDO
     CALL COMPSP(ctx, 3, ctx%state%nz, file1, mass_ssp, lbol_ssp, spec_ssp, pset, ocompsp)
  ELSE
     ! We already called SPS_SETUP(-1), so speclib is populated.
   CALL generate_ssp_grid(ctx, pset, mass_ssp(:,pset%zmet), lbol_ssp(:,pset%zmet), spec_ssp(:,:,pset%zmet))
   CALL COMPSP(ctx, 3, 1, file1, mass_ssp(:,pset%zmet:pset%zmet), lbol_ssp(:,pset%zmet:pset%zmet), &
        spec_ssp(:,:,pset%zmet:pset%zmet), pset, ocompsp)
  ENDIF

  ! Clean up
  IF (ALLOCATED(spec_ssp)) DEALLOCATE(spec_ssp)
  IF (ALLOCATED(mass_ssp)) DEALLOCATE(mass_ssp)
  IF (ALLOCATED(lbol_ssp)) DEALLOCATE(lbol_ssp)
  IF (ALLOCATED(ocompsp))  DEALLOCATE(ocompsp)

END PROGRAM AUTOSPS
