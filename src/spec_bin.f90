PROGRAM SPEC_BIN

  !routine to convert ascii spectral files to binary
  !must be run twice for each value of isoc_type var

  USE sps_vars
  USE sps_utils, ONLY: sps_setup ! To access setup if needed, or rely on vars
  IMPLICIT NONE
  INTEGER  :: z,dumi1,i,j
  REAL(SP) :: dumr1,d2,d3
  CHARACTER(6) :: zstype
  INTEGER :: alloc_stat

  !----------------------------------------------------------------!

  CALL GETENV('SPS_HOME',SPS_HOME)
   
   ! Initialize variables and dimensions using default libraries.
   ! This sets nspec, nzinit, and allocates speclib arrays.
   CALL SPS_SETUP(-1)
   
   ! Close the zlegend file if it was left open (safety check)
   CLOSE(90, IOSTAT=dumi1) 

   DO z=1,nzinit
     
     WRITE(zstype,'(F6.4)') zlegendinit(z)

     IF (TRIM(spec_type).EQ.'basel') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'&
             //basel_str//'_z'//zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/BaSeL3.1/basel_'&
             //basel_str//'_z'//zstype//'.spectra.bin',&
             FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)

     ELSE IF (TRIM(spec_type).EQ.'miles') THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/MILES/imiles_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)

     ELSE IF (INDEX(TRIM(spec_type), 'c3k') > 0) THEN
        OPEN(92,FILE=TRIM(SPS_HOME)//'/SPECTRA/CKC14/'//TRIM(spec_type)//'_z'&
             //zstype//'.spectra',FORM='FORMATTED',&
             STATUS='OLD',ACTION='READ')
        OPEN(93,FILE=TRIM(SPS_HOME)//'/SPECTRA/CKC14/'//TRIM(spec_type)//'_z'&
             //zstype//'.spectra.bin',FORM='UNFORMATTED',&
             STATUS='REPLACE',access='direct',&
             recl=nspec*ndim_logg*ndim_logt*4)
  
    ENDIF

     DO i=1,ndim_logg
        DO j=1,ndim_logt
           READ(92,*) dumi1,d2,d3,dumr1
           READ(92,*) speclib(:,z,j,i)
        ENDDO
     ENDDO

     WRITE(93,rec=1) speclib(:,z,:,:)
     CLOSE(92)
     CLOSE(93)
     
  ENDDO
  
END PROGRAM SPEC_BIN
