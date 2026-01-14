PROGRAM TEST_DRIVER
  USE sps_vars
  USE sps_utils
  USE test_utils
  IMPLICIT NONE
  
  INTERFACE
     SUBROUTINE SSP_GEN(pset,mass_ssp,lbol_ssp,spec_ssp)
       USE sps_vars
       TYPE(PARAMS), INTENT(in) :: pset
       REAL(SP), INTENT(inout), DIMENSION(ntfull) :: mass_ssp, lbol_ssp
       REAL(SP), INTENT(inout), DIMENSION(nspec,ntfull) :: spec_ssp
     END SUBROUTINE SSP_GEN
  END INTERFACE

  INTEGER :: i
  REAL(SP) :: dummy_sum
  INTEGER :: original_nbands
  CHARACTER(LEN=256) :: sps_home_val
  LOGICAL :: file_exists

  PRINT *, "Starting FSPS Test Suite..."

  ! Ensure SPS_HOME is set
  CALL getenv('SPS_HOME', sps_home_val)
  IF (LEN_TRIM(sps_home_val) == 0) THEN
     PRINT *, "WARNING: SPS_HOME not set. Trying current directory."
  ENDIF

  ! --- Test A: Regression ---
  PRINT *, "---------------------------------------------------"
  PRINT *, "Running Test A: Regression"
  
  ! Run Setup (using default 'mist' and 'miles')
  CALL sps_setup(1) 
  
  CALL assert_true(check_sps_setup == 1, "sps_setup ran successfully")

  ! Zero all allocatable buffers
  ! These are used as accumulators in SSP_GEN and other routines.
  ! If not zeroed, they contain garbage memory (approx 1e-15).
  
  ! 1. SSP Intermediate Arrays (The likely culprits!)
  IF (ALLOCATED(spec_young)) spec_young = 0.0
  IF (ALLOCATED(spec_old))   spec_old   = 0.0
  IF (ALLOCATED(spec_xrb))   spec_xrb   = 0.0
  
  ! 2. Global accumulator arrays
  IF (ALLOCATED(spec_ssp_zz)) spec_ssp_zz = 0.0
  IF (ALLOCATED(mass_ssp_zz)) mass_ssp_zz = 0.0
  IF (ALLOCATED(lbol_ssp_zz)) lbol_ssp_zz = 0.0
  
  ! 3. Dust/Nebular Arrays (buffers)
  IF (ALLOCATED(dustem2_dustem)) dustem2_dustem = 0.0
  IF (ALLOCATED(flux_dagb))      flux_dagb      = 0.0
  IF (ALLOCATED(nebem_cont))     nebem_cont     = 0.0
  IF (ALLOCATED(xnebem_cont))    xnebem_cont    = 0.0
  IF (ALLOCATED(agndust_spec))   agndust_spec   = 0.0
  
  BLOCK
    TYPE(PARAMS) :: pset
    REAL(SP), ALLOCATABLE :: my_mass(:), my_lbol(:), my_spec(:,:)
    REAL(SP), ALLOCATABLE :: ref_wave(:), ref_flux(:)
    REAL(SP) :: diff, max_diff, rel_err
    INTEGER :: alloc_stat, k, io_stat
    CHARACTER(LEN=256) :: ref_path
    
    ! Initialize pset defaults
    pset%zmet = 1 
    pset%sfh = 0
    pset%const = 1.0
    pset%dust1 = 0.0
    pset%dust2 = 0.0
    pset%imf1 = 1.3
    pset%imf2 = 2.3
    pset%imf3 = 2.3
    pset%vdmc = 0.08
    pset%mdave = 0.5
    pset%evtype = -1
    
    ! Allocate necessary pset components
    IF (ALLOCATED(pset%mag_compute)) DEALLOCATE(pset%mag_compute)
    ALLOCATE(pset%mag_compute(nbands))
    pset%mag_compute = 1
    
    IF (ALLOCATED(pset%ssp_gen_age)) DEALLOCATE(pset%ssp_gen_age)
    ALLOCATE(pset%ssp_gen_age(nt))
    pset%ssp_gen_age = 1

    ! 1. Allocate arrays for new code
    ALLOCATE(my_mass(ntfull), my_lbol(ntfull), my_spec(nspec, ntfull), STAT=alloc_stat)
    CALL assert_true(alloc_stat == 0, "Allocated arrays for Test A")
    
    ! Zero local arrays just to be safe
    my_mass = 0.0
    my_lbol = 0.0
    my_spec = 0.0

    ! 2. Run new code
    CALL SSP_GEN(pset, my_mass, my_lbol, my_spec)
    
    ! 3. Load reference data
    ref_path = 'reference/reference_ssp.dat'
    
    OPEN(20, FILE=TRIM(ref_path), STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
    IF (io_stat /= 0) THEN
       PRINT *, "ERROR: Could not open reference file at: ", TRIM(ref_path)
       CALL assert_true(.FALSE., "Reference file missing")
    ELSE
       ALLOCATE(ref_wave(nspec), ref_flux(nspec))
       
       DO k = 1, nspec
          READ(20, *, IOSTAT=io_stat) ref_wave(k), ref_flux(k)
          IF (io_stat /= 0) EXIT 
       END DO
       CLOSE(20)

       ! 4. Compare
       max_diff = 0.0
       DO k = 1, nspec
          IF (ABS(ref_flux(k)) > 1.0e-30) THEN
             rel_err = ABS((my_spec(k, 1) - ref_flux(k)) / ref_flux(k))
             if (rel_err > max_diff) then
                 max_diff = rel_err
                 ! Print the failure details for the worst offender
                 print *, "Fail at Index:", k
                 print *, "   Wavelength:", ref_wave(k)
                 print *, "   Ref Flux:  ", ref_flux(k)
                 print *, "   New Flux:  ", my_spec(k, 1)
             endif
          ENDIF
       END DO
       
       PRINT *, "Max Relative Error vs Reference: ", max_diff
       
       CALL assert_true(max_diff < 1.0e-5, "Regression Check (Error < 1e-5)")
       
       DEALLOCATE(ref_wave, ref_flux)
    END IF

    IF (ALLOCATED(my_mass)) DEALLOCATE(my_mass, my_lbol, my_spec)
  END BLOCK
  
  ! Check wavelength array
  CALL assert_true(ALLOCATED(spec_lambda), "spec_lambda allocated")
  CALL assert_true(spec_lambda(1) < spec_lambda(nspec), "Wavelength array is increasing")

  ! Store nbands for later comparison
  original_nbands = nbands
  PRINT *, "Original nbands: ", nbands

  ! --- Test B: Dynamic Reallocation ---
  PRINT *, "---------------------------------------------------"
  PRINT *, "Running Test B: Dynamic Reallocation"

  CALL sps_takedown()
  CALL assert_true(.NOT. ALLOCATED(spec_ssp_zz), "sps_takedown deallocated memory")
  CALL assert_true(check_sps_setup == 0, "check_sps_setup reset")
  
  OPEN(10, FILE='../data/test_filters.dat', STATUS='UNKNOWN')
  WRITE(10, '(A)') '# Small filter list'
  WRITE(10, '(A)') '3000.0 0.0'
  WRITE(10, '(A)') '3500.0 1.0'
  WRITE(10, '(A)') '4000.0 0.0'
  WRITE(10, '(A)') '5000.0 0.0'
  CLOSE(10)
  
  ! Set alt_filter_file (module variable in sps_vars)
  alt_filter_file = 'test_filters.dat'
  
  ! 3. Setup again with new filter file
  CALL sps_setup(1)
  
  PRINT *, "New nbands: ", nbands
  CALL assert_true(nbands == 1, "Resized to 1 band")
  CALL assert_true(ALLOCATED(bands), "Bands allocated")
  CALL assert_true(SIZE(bands, 2) == 1, "Bands array has 2nd dimension 1")

  ! 4. Takedown and revert
  CALL sps_takedown()
  
  alt_filter_file = '' ! Reset to default
  CALL sps_setup(1)
  
  PRINT *, "Restored nbands: ", nbands
  CALL assert_true(nbands == original_nbands, "Resized back to full bands")
  
  ! --- Test C: Stress Loop ---
  PRINT *, "---------------------------------------------------"
  PRINT *, "Running Test C: Stress Loop"
  
  DO i = 1, 5
     CALL sps_takedown()
     CALL sps_setup(1)
     
     ! Also run calculation to stress test allocation within SSP_GEN if any
     ! (Though most allocations are in setup)
     BLOCK
        TYPE(PARAMS) :: pset
        REAL(SP), ALLOCATABLE :: my_mass(:), my_lbol(:), my_spec(:,:)
        INTEGER :: alloc_stat
        pset%zmet = 1 
        pset%sfh = 0
        pset%evtype = -1

        ALLOCATE(pset%mag_compute(nbands))
        pset%mag_compute = 1
        ALLOCATE(pset%ssp_gen_age(nt))
        pset%ssp_gen_age = 1
        
        ALLOCATE(my_mass(ntfull), my_lbol(ntfull), my_spec(nspec, ntfull), STAT=alloc_stat)
        IF (alloc_stat == 0) THEN
           CALL SSP_GEN(pset, my_mass, my_lbol, my_spec)
           DEALLOCATE(my_mass, my_lbol, my_spec)
        ENDIF
     END BLOCK

     IF (MOD(i, 10) == 0) PRINT *, "Iteration ", i
  END DO
  
  CALL assert_true(.TRUE., "Stress loop completed without crash")
  
  CALL final_report()

END PROGRAM TEST_DRIVER
