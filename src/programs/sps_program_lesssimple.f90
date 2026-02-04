PROGRAM LESSSIMPLE

  !set up modules
      USE fsps_constants, ONLY: SP
      USE fsps_types, ONLY: PARAMS, COMPSPOUT
        USE sps_utils
        USE fsps_context, ONLY: fsps_context_create
        USE fsps_context_types, ONLY: fsps_context_t
  
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90
        INTEGER :: i
        TYPE(fsps_context_t) :: ctx
  !define variable for SSP spectrum
      REAL(SP), ALLOCATABLE :: spec_pz(:,:)
      REAL(SP), ALLOCATABLE :: spec_pz_zz(:,:,:)
  !define variables for Mass and Lbol info
      REAL(SP), ALLOCATABLE :: mass_pz(:),lbol_pz(:)
      REAL(SP), ALLOCATABLE :: mass_pz_zz(:,:), lbol_pz_zz(:,:)
  CHARACTER(100) :: file2=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), ALLOCATABLE :: ocompsp(:)
  REAL(SP) :: zave

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Now we're going to show you how to use full  
  ! metallicity-dependent info                   
  
        CALL fsps_context_create(ctx)

      ctx%imf_type_val = 0    ! Salpeter IMF
  pset%sfh = 0    ! compute SSP

  !here we have to read in all the librarries
        CALL SPS_SETUP(ctx, -1)

  ! Allocate memory now that SPS_SETUP has defined ntfull/nspec
  IF (.NOT. ALLOCATED(spec_pz)) THEN
        ALLOCATE(spec_pz(ctx%state%nspec, ctx%state%ntfull))
        ALLOCATE(spec_pz_zz(ctx%state%nspec, ctx%state%ntfull, 1))
        ALLOCATE(mass_pz(ctx%state%ntfull))
        ALLOCATE(lbol_pz(ctx%state%ntfull))
        ALLOCATE(mass_pz_zz(ctx%state%ntfull, 1))
        ALLOCATE(lbol_pz_zz(ctx%state%ntfull, 1))
        ALLOCATE(ocompsp(ctx%state%ntfull))
  END IF

  !compute all SSPs (i.e. at all Zs)
  !nz and the various *ssp_zz arrays are stored 
  !in the common block set up in sps_vars.f90
  DO i=1,ctx%state%nz
     pset%zmet = i
        CALL SSP_GEN(ctx, pset, ctx%state%mass_ssp_zz(:,i),&
              ctx%state%lbol_ssp_zz(:,i),ctx%state%spec_ssp_zz(:,:,i))
  ENDDO

  !define the yield for a closed box distribution
  pset%pmetals = 0.02
  !compute SSP convolved with a closed box  
        CALL PZ_CONVOL(ctx, pset%pmetals, zave, spec_pz, lbol_pz, mass_pz)
  spec_pz_zz(:,:,1) = spec_pz
  mass_pz_zz(:,1) = mass_pz
  lbol_pz_zz(:,1) = lbol_pz
  file2    = 'SSP_pz.out'
  !now compute magnitudes for this SSP
        CALL COMPSP(ctx, 1, 1, file2, mass_pz_zz, lbol_pz_zz, spec_pz_zz, pset, ocompsp)

  !run compsp for a tabulated sfh with a metallicity history
  !NB: one must have setup all the SSPs, as was done in the DO-loop above
  pset%sfh = 2
  file2    = 'CSP_tabsfh.out'
  CALL COMPSP(ctx, 1, ctx%state%nz, file2, ctx%state%mass_ssp_zz, ctx%state%lbol_ssp_zz, &
        ctx%state%spec_ssp_zz, pset, ocompsp)

  ! Clean up memory before exiting
      IF (ALLOCATED(spec_pz)) DEALLOCATE(spec_pz)
      IF (ALLOCATED(spec_pz_zz)) DEALLOCATE(spec_pz_zz)
  IF (ALLOCATED(mass_pz)) DEALLOCATE(mass_pz)
  IF (ALLOCATED(lbol_pz)) DEALLOCATE(lbol_pz)
      IF (ALLOCATED(mass_pz_zz)) DEALLOCATE(mass_pz_zz)
      IF (ALLOCATED(lbol_pz_zz)) DEALLOCATE(lbol_pz_zz)
  IF (ALLOCATED(ocompsp))  DEALLOCATE(ocompsp)

END PROGRAM LESSSIMPLE
