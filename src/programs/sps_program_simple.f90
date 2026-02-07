 PROGRAM SIMPLE

  !set up modules
  USE fsps_precision, ONLY: WP
  USE fsps_types, ONLY: PARAMS, COMPSPOUT
  USE sps_utils
  USE fsps_csp, ONLY: compute_csp_scenario
  USE fsps_io, ONLY: write_csp_output_files
  USE fsps_ssp, ONLY: generate_ssp_grid
  USE fsps_context, ONLY: fsps_context_create
  use fsps_sfh, only: compute_sfh_statistics
  USE fsps_context_types, ONLY: fsps_context_t
  IMPLICIT NONE

  !NB: the various structure types are defined in sps_vars.f90
  !    variables not explicitly defined here are defined in sps_vars.f90

  TYPE(fsps_context_t) :: ctx
  !define variable for SSP spectrum
  REAL(WP), ALLOCATABLE :: spec_ssp(:,:,:)
  !define variables for Mass and Lbol info
  REAL(WP), ALLOCATABLE :: mass_ssp(:,:), lbol_ssp(:,:)
  CHARACTER(100) :: file1=''
  !structure containing all necessary parameters
  TYPE(PARAMS) :: pset
  !define structure for CSP spectrum
  TYPE(COMPSPOUT), ALLOCATABLE :: results(:)
  REAL(WP) :: ssfr_stats(3), ave_age

  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  
  ! Lets compute an SSP, solar metallicity, with a Chabrier IMF
  ! with no dust, and the 'default' assumptions regarding the 
  ! locations of the isochrones

  CALL fsps_context_create(ctx)

  ctx%imf_type_val  = 0             !define the IMF (1=Chabrier 2003)
                            !see sps_vars.f90 for details of this var
  pset%zmet = 10            !define the metallicity (see the manual)
                            !20 = solar metallacity

  CALL SPS_SETUP(ctx, pset%zmet) !read in the isochrones and spectral libraries

  ! Allocate memory now that SPS_SETUP has defined ntfull/nspec
      IF (.NOT. ALLOCATED(spec_ssp)) THEN
        ALLOCATE(spec_ssp(ctx%state%nspec, ctx%state%ntfull, 1))
        ALLOCATE(mass_ssp(ctx%state%ntfull, 1))
        ALLOCATE(lbol_ssp(ctx%state%ntfull, 1))
          ! results is allocated by compute_csp_scenario
  END IF

  !define the parameter set.  These are the default values, specified 
  !in sps_vars.f90, but are explicitly included here for transparency
  pset%sfh   = 0     !set SFH to "SSP"
  pset%const = 1.
  pset%zred  = 0.0   !redshift  
  pset%dust1 = 0.0   !dust parameter 1
  pset%dust2 = 0.0   !dust parameter 2

  pset%dell  = 0.0   !shift in log(L) for TP-AGB stars
  pset%delt  = 0.0   !shift in log(Teff) for TP-AGB stars
  pset%fbhb  = 0.0   !fraction of blue HB stars
  pset%sbss  = 0.0   !specific frequency of BS stars

  ctx%add_neb_emission_val=1

  !compute the SSP
  CALL generate_ssp_grid(ctx, pset, mass_ssp(:,1), lbol_ssp(:,1), spec_ssp(:,:,1))
  !compute mags and write out mags and spec for SSP
  file1 = 'SSP_BPASS.out'
  CALL compute_csp_scenario(ctx, pset, 1, spec_ssp, mass_ssp, lbol_ssp, results)
  CALL write_csp_output_files(ctx, pset, results, file1, 3)
  DEALLOCATE(results)



  ! Now lets compute a 1 Gyr tau model SFH with a Salpeter IMF
  ! with a simple dust model, at a particular time 
  ! (rather than outputing all the time info)

  ctx%imf_type_val  = 0                !define the IMF (0=Salpeter)
                               !see sps_vars.f90 for details of this var

  !NB: you only need to re-run SPS_SETUP if you have changed the metallicity
  !    or, even better (but slower), you can call SPS_SETUP(-1) and this will
  !    set up all the metallicities at once.
  CALL SPS_SETUP(ctx, pset%zmet)    !read in the isochrones and spectral libraries

  !define the parameter set. 
  pset%sfh   = 1     !set SFH to "CSP"; sfh=1 means a normal tau model
  pset%tau   = 2.0   !tau units are Gyr
  pset%dust1 = 1.0   !dust parameter 1
  pset%dust2 = 0.3   !dust parameter 2
  !pset%tage  = 12.5  !age at which we want the mags

  !compute the CSP
  CALL generate_ssp_grid(ctx, pset, mass_ssp(:,1), lbol_ssp(:,1), spec_ssp(:,:,1))
  !compute mags, and write out mags and spec for CSP
  file1 = 'CSP.out'
  CALL compute_csp_scenario(ctx, pset, 1, spec_ssp, mass_ssp, lbol_ssp, results)
  CALL write_csp_output_files(ctx, pset, results, file1, 3)

  !compute basic SFH statistics for the last entry in the ocompsp array
  !results are returned in the variable ssfr_stats (array) and ave_age
  CALL compute_sfh_statistics(ctx, pset, results(1), ssfr_stats, ave_age)

  ! Clean up memory before exiting
  IF (ALLOCATED(spec_ssp)) DEALLOCATE(spec_ssp)
  IF (ALLOCATED(mass_ssp)) DEALLOCATE(mass_ssp)
  IF (ALLOCATED(lbol_ssp)) DEALLOCATE(lbol_ssp)
  IF (ALLOCATED(results))  DEALLOCATE(results)


END PROGRAM SIMPLE
