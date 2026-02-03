PROGRAM TEST_RUNNER

  ! FSPS Regression Test Runner
  ! Reads reference data from disk, runs current FSPS code with identical
  ! parameters, and compares outputs with a relative tolerance.

  USE, INTRINSIC :: IEEE_ARITHMETIC
      USE fsps_types, ONLY: SP, PARAMS, COMPSPOUT, nemline
      USE sps_utils
   USE fsps_context_types, ONLY: fsps_context_t
         USE fsps_context, ONLY: fsps_context_create, fsps_context_set_pset
  IMPLICIT NONE

  ! Exit codes
  INTEGER, PARAMETER :: EXIT_SUCCESS = 0
  INTEGER, PARAMETER :: EXIT_FAILURE = 1
  
  ! Default tolerance (can be overridden by env var FSPS_TEST_RTOL)
  REAL(SP), PARAMETER :: DEFAULT_RTOL = 1.0E-5

  ! Test arrays (allocatable)
  ! Reference data (read from disk)
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: ref_spec_ssp
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ref_mass_ssp, ref_lbol_ssp
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: ref_ocompsp

  ! New data (computed on the fly)
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: new_spec_ssp
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: new_spec_ssp_ctx
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: new_spec_ssp_cmp
   REAL(SP), ALLOCATABLE, DIMENSION(:,:,:) :: new_spec_ssp3
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: new_mass_ssp2, new_lbol_ssp2
  REAL(SP), ALLOCATABLE, DIMENSION(:)   :: new_mass_ssp, new_lbol_ssp
  TYPE(COMPSPOUT), ALLOCATABLE, DIMENSION(:) :: new_ocompsp

  ! Control variables
   TYPE(fsps_context_t) :: ctx
   TYPE(PARAMS) :: pset
   INTEGER :: i, unit_in, status, arg_count
  CHARACTER(LEN=255) :: filename_in, env_buffer, arg_val
  CHARACTER(LEN=100) :: csp_dummy_file
  CHARACTER(LEN=20) :: isoc_arg, spec_arg, dust_arg
  LOGICAL :: isoc_set, spec_set, dust_set
  
  ! Dimensions read from file
   INTEGER :: file_nspec, file_ntfull, file_nbands
   INTEGER :: nspec_ctx, ntfull_ctx, nbands_ctx, nt_ctx, nindx_ctx

  ! Comparison stats
   REAL(SP) :: rtol
   LOGICAL :: test_passed
   INTEGER :: nfail
   LOGICAL :: verbose_output
   INTEGER :: max_fail_print, fail_printed
   LOGICAL :: fail_suppression_noted


   ! Configuration
   unit_in = 40
   csp_dummy_file = 'dummy_csp.out'

   CALL ENSURE_OUTPUT_DIR('OUTPUTS/')

   test_passed = .TRUE.
   nfail = 0
   verbose_output = .FALSE.
   max_fail_print = 20
   fail_printed = 0
   fail_suppression_noted = .FALSE.
   isoc_set = .FALSE.
   spec_set = .FALSE.
   dust_set = .FALSE.
   isoc_arg = 'mist'
   spec_arg = 'miles'
   dust_arg = 'DL07'
   filename_in = ''

  WRITE(*,*) '========================================='
  WRITE(*,*) 'FSPS TEST RUNNER'
  WRITE(*,*) '========================================='

  ! Parse command line arguments
  arg_count = COMMAND_ARGUMENT_COUNT()
  i = 1
  DO WHILE (i <= arg_count)
     CALL GET_COMMAND_ARGUMENT(i, arg_val, STATUS=status)
     IF (status /= 0) EXIT
     
     IF (TRIM(arg_val) == '--isoc') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, isoc_arg, STATUS=status)
           isoc_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --isoc requires an argument'
           STOP EXIT_FAILURE
        END IF
     ELSE IF (TRIM(arg_val) == '--spec') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, spec_arg, STATUS=status)
           spec_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --spec requires an argument'
           STOP EXIT_FAILURE
        END IF
     ELSE IF (TRIM(arg_val) == '--dust') THEN
        i = i + 1
        IF (i <= arg_count) THEN
           CALL GET_COMMAND_ARGUMENT(i, dust_arg, STATUS=status)
           dust_set = .TRUE.
        ELSE
           WRITE(*,*) 'ERROR: --dust requires an argument'
           STOP EXIT_FAILURE
        END IF
     ELSE
        ! Assume it's the filename if it doesn't start with --
        ! Or if we haven't found one yet.
        IF (LEN_TRIM(filename_in) == 0) THEN
           filename_in = TRIM(arg_val)
        END IF
     END IF
     i = i + 1
  END DO

  IF (LEN_TRIM(filename_in) == 0) THEN
     WRITE(*,*) 'ERROR: Must provide reference filename as argument.'
     WRITE(*,*) 'Usage: ./test_runner [--isoc type] [--spec type] [--dust type] tests/data/sps_ref_XXX.bin'
     STOP EXIT_FAILURE
  END IF

  ! Read tolerance from environment variable
  CALL GET_ENVIRONMENT_VARIABLE("FSPS_TEST_RTOL", VALUE=env_buffer, STATUS=status)
  IF (status == 0) THEN
     READ(env_buffer, *) rtol
     WRITE(*,*) 'Using RTOL from environment: ', rtol
  ELSE
     rtol = DEFAULT_RTOL
     WRITE(*,*) 'Using default RTOL: ', rtol
  END IF

  ! Optional verbose output
  CALL GET_ENVIRONMENT_VARIABLE("FSPS_TEST_VERBOSE", VALUE=env_buffer, STATUS=status)
  IF (status == 0) THEN
     IF (LEN_TRIM(env_buffer) > 0) THEN
        SELECT CASE (env_buffer(1:1))
        CASE ('1','t','T','y','Y')
           verbose_output = .TRUE.
        END SELECT
     END IF
  END IF

  ! Optional maximum printed failures
  CALL GET_ENVIRONMENT_VARIABLE("FSPS_TEST_MAXFAIL", VALUE=env_buffer, STATUS=status)
  IF (status == 0) THEN
     READ(env_buffer, *, IOSTAT=status) max_fail_print
     IF (status /= 0) max_fail_print = 20
  END IF

  ! Initialize FSPS and check dimensions
  ! Note: We must initialize FSPS before allocating, but we must read the 
  ! file header before we know if dimensions match.
  
   CALL fsps_context_create(ctx)
   ctx%imf_type_val = 1
   pset%zmet = 10
  
  WRITE(*,*) 'Initializing FSPS...'
  ! Always provide defaults that match the reference generator unless overridden
   CALL SPS_SETUP(ctx, pset%zmet, isoc_type_in=TRIM(isoc_arg), spec_type_in=TRIM(spec_arg), dust_type_in=TRIM(dust_arg))
   IF (verbose_output) CALL DUMP_STATE('AFTER SPS_SETUP', ctx, pset)

  nspec_ctx = ctx%state%nspec
  ntfull_ctx = ctx%state%ntfull
  nbands_ctx = ctx%state%nbands
  nt_ctx = ctx%state%nt
  nindx_ctx = ctx%state%nindx

   ! Global dimensions are synchronized via fsps_context_apply_globals

  ! Allocate pset allocatable components
  IF (ALLOCATED(pset%mag_compute)) DEALLOCATE(pset%mag_compute)
   ALLOCATE(pset%mag_compute(nbands_ctx))
  pset%mag_compute = 1
  
  IF (ALLOCATED(pset%ssp_gen_age)) DEALLOCATE(pset%ssp_gen_age)
   ALLOCATE(pset%ssp_gen_age(nt_ctx))
  pset%ssp_gen_age = 1

  ! Open the reference file
  OPEN(UNIT=unit_in, FILE=TRIM(filename_in), STATUS='OLD', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
  IF (status /= 0) THEN
     WRITE(*,*) 'ERROR: Could not open file: ', TRIM(filename_in)
     STOP EXIT_FAILURE
  END IF

  ! Read Header
  READ(unit_in) file_nspec
  READ(unit_in) file_ntfull
  READ(unit_in) file_nbands

   WRITE(*,*) 'Reference Dimensions: nspec=', file_nspec, ' nt=', file_ntfull
   WRITE(*,*) 'Compiled Dimensions:  nspec=', nspec_ctx,  ' nt=', ntfull_ctx

  ! Strict dimension check
   IF (file_nspec /= nspec_ctx .OR. file_ntfull /= ntfull_ctx) THEN
     WRITE(*,*) 'FATAL: Binary dimensions do not match compiled FSPS dimensions.'
     WRITE(*,*) 'Ensure you are running the test with the same flags/arguments used to generate the data.'
     STOP EXIT_FAILURE
  END IF

  ! Allocate and read reference data
   ALLOCATE(ref_spec_ssp(nspec_ctx,ntfull_ctx))
   ALLOCATE(ref_mass_ssp(ntfull_ctx))
   ALLOCATE(ref_lbol_ssp(ntfull_ctx))
   ALLOCATE(ref_ocompsp(ntfull_ctx))

   ALLOCATE(new_spec_ssp(ntfull_ctx,nspec_ctx))
   ALLOCATE(new_spec_ssp_ctx(nspec_ctx,ntfull_ctx))
   ALLOCATE(new_spec_ssp_cmp(nspec_ctx,ntfull_ctx))
   ALLOCATE(new_mass_ssp(ntfull_ctx))
   ALLOCATE(new_lbol_ssp(ntfull_ctx))
   ALLOCATE(new_mass_ssp2(ntfull_ctx,1))
   ALLOCATE(new_lbol_ssp2(ntfull_ctx,1))
   ALLOCATE(new_spec_ssp3(nspec_ctx,ntfull_ctx,1))
   ALLOCATE(new_ocompsp(ntfull_ctx))

  WRITE(*,*) 'Reading SSP reference data...'
  READ(unit_in) ref_mass_ssp
  READ(unit_in) ref_lbol_ssp
  READ(unit_in) ref_spec_ssp

  WRITE(*,*) 'Reading CSP reference data...'
   DO i = 1, ntfull_ctx
     ! Allocate components of derived type before reading
   ALLOCATE(ref_ocompsp(i)%mags(nbands_ctx))
   ALLOCATE(ref_ocompsp(i)%spec(nspec_ctx))
   ALLOCATE(ref_ocompsp(i)%indx(nindx_ctx))
     ALLOCATE(ref_ocompsp(i)%emlines(nemline))
     
     READ(unit_in) ref_ocompsp(i)%age
     READ(unit_in) ref_ocompsp(i)%mass_csp
     READ(unit_in) ref_ocompsp(i)%lbol_csp
     READ(unit_in) ref_ocompsp(i)%sfr
     READ(unit_in) ref_ocompsp(i)%mdust
     READ(unit_in) ref_ocompsp(i)%mformed
     READ(unit_in) ref_ocompsp(i)%mags      
     READ(unit_in) ref_ocompsp(i)%spec      
     READ(unit_in) ref_ocompsp(i)%indx      
     READ(unit_in) ref_ocompsp(i)%emlines   
  END DO
  CLOSE(unit_in)

  ! Generate new data
  WRITE(*,*) 'Generating new SSP data...'
  pset%sfh   = 0     
  pset%const = 0.0
  pset%zred  = 0.0
  pset%dust1 = 0.0
  pset%dust2 = 0.0
   ctx%add_neb_emission_val = 1
   CALL fsps_context_set_pset(ctx, pset)
   IF (verbose_output) CALL DUMP_STATE('BEFORE SSP_GEN (SSP)', ctx, pset)
   CALL SSP_GEN(ctx, pset, new_mass_ssp, new_lbol_ssp, new_spec_ssp_ctx)
   IF (verbose_output) CALL DUMP_SSP_SUMMARY('AFTER SSP_GEN (SSP)', ctx, new_mass_ssp, new_lbol_ssp)
   new_spec_ssp_cmp = new_spec_ssp_ctx

  WRITE(*,*) 'Generating new CSP data...'
  pset%sfh   = 1    
  pset%tau   = 2.0   
  pset%dust1 = 1.0   
  pset%dust2 = 0.3
   CALL fsps_context_set_pset(ctx, pset)
   IF (verbose_output) CALL DUMP_STATE('BEFORE SSP_GEN (CSP)', ctx, pset)
   CALL SSP_GEN(ctx, pset, new_mass_ssp, new_lbol_ssp, new_spec_ssp_ctx)
   IF (verbose_output) CALL DUMP_SSP_SUMMARY('AFTER SSP_GEN (CSP)', ctx, new_mass_ssp, new_lbol_ssp)

  DO i = 1, ntfull_ctx
     IF (.NOT. ALLOCATED(new_ocompsp(i)%mags)) ALLOCATE(new_ocompsp(i)%mags(nbands_ctx))
     IF (.NOT. ALLOCATED(new_ocompsp(i)%spec)) ALLOCATE(new_ocompsp(i)%spec(nspec_ctx))
     IF (.NOT. ALLOCATED(new_ocompsp(i)%indx)) ALLOCATE(new_ocompsp(i)%indx(nindx_ctx))
     IF (.NOT. ALLOCATED(new_ocompsp(i)%emlines)) ALLOCATE(new_ocompsp(i)%emlines(nemline))
  END DO

   new_mass_ssp2(:,1) = new_mass_ssp
   new_lbol_ssp2(:,1) = new_lbol_ssp
    ! Match COMPSP input layout: (nspec, ntfull, nzin)
    new_spec_ssp3(:,:,1) = new_spec_ssp_ctx

   IF (verbose_output) CALL DUMP_STATE('BEFORE COMPSP', ctx, pset)
   CALL COMPSP(ctx, 3, 1, csp_dummy_file, new_mass_ssp2, new_lbol_ssp2, new_spec_ssp3, pset, new_ocompsp)
   IF (verbose_output) CALL DUMP_CSP_SUMMARY('AFTER COMPSP', new_ocompsp)


  ! Compare results
  WRITE(*,*) 'Verifying results (RTOL = ', rtol, ')...'

  ! Helper internal subroutine to check arrays
   CALL CHECK_ARRAY_2D("SSP Spectra", ref_spec_ssp, new_spec_ssp_cmp, nspec_ctx, ntfull_ctx)
   CALL CHECK_ARRAY_1D("SSP Mass", ref_mass_ssp, new_mass_ssp, ntfull_ctx)
   CALL CHECK_ARRAY_1D("SSP Lbol", ref_lbol_ssp, new_lbol_ssp, ntfull_ctx)

  ! Check CSP structure components manually
   DO i = 1, ntfull_ctx
     ! Check Scalars
     CALL CHECK_VAL("CSP Lbol", i, ref_ocompsp(i)%lbol_csp, new_ocompsp(i)%lbol_csp)
     CALL CHECK_VAL("CSP Mass", i, ref_ocompsp(i)%mass_csp, new_ocompsp(i)%mass_csp)
     CALL CHECK_VAL("CSP SFR", i, ref_ocompsp(i)%sfr, new_ocompsp(i)%sfr)
     CALL CHECK_VAL("CSP Dust Mass", i, ref_ocompsp(i)%mdust, new_ocompsp(i)%mdust)
     CALL CHECK_VAL("CSP Mass Formed", i, ref_ocompsp(i)%mformed, new_ocompsp(i)%mformed)

     ! Check Arrays for EVERY time step
   CALL CHECK_MAGS_1D("CSP Mags (flux)", ref_ocompsp(i)%mags, new_ocompsp(i)%mags, nbands_ctx)
   CALL CHECK_ARRAY_1D("CSP Indx", ref_ocompsp(i)%indx, new_ocompsp(i)%indx, nindx_ctx)
     CALL CHECK_ARRAY_1D("CSP Emlines", ref_ocompsp(i)%emlines, new_ocompsp(i)%emlines, nemline)
     
     ! Check Spectrum
   CALL CHECK_ARRAY_1D("CSP Spec", ref_ocompsp(i)%spec, new_ocompsp(i)%spec, nspec_ctx)
  END DO

  ! Report results
  WRITE(*,*) '--------------------------------------------------'  

   CALL SPS_TAKEDOWN(ctx)

   WRITE(*,*) 'Total failures:', nfail
   IF (test_passed) THEN
     WRITE(*,*) 'TEST RESULT: PASS'
     STOP EXIT_SUCCESS
  ELSE
     WRITE(*,*) 'TEST RESULT: FAIL'
     STOP EXIT_FAILURE
  END IF


CONTAINS

  SUBROUTINE CHECK_ARRAY_2D(label, ref, new, d1, d2)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: d1, d2
    REAL(SP), DIMENSION(d1,d2), INTENT(IN) :: ref, new
   REAL(SP) :: delta, threshold, max_delta, max_rel, ref_val
   INTEGER :: j, k, mj, mk

   max_delta = 0.0
   max_rel = 0.0
   mj = 1
   mk = 1

    DO k = 1, d2
       DO j = 1, d1
          IF (IEEE_IS_NAN(ref(j,k)) .OR. IEEE_IS_NAN(new(j,k))) THEN
             WRITE(*,*) 'FAIL: ', label, ' contains NaN at index (',j,',',k,')'
             test_passed = .FALSE.
             RETURN 
          END IF

          ref_val = ref(j,k)
          delta = ABS(ref_val - new(j,k))
          if (ABS(ref_val) > 0.0_SP) then
             max_rel = MAX(max_rel, delta / ABS(ref_val))
          end if
          IF (delta > max_delta) THEN
             max_delta = delta
             mj = j
             mk = k
          END IF
          ! If ref is close to zero, use absolute tolerance, else relative
          threshold = MAX(ABS(ref_val) * rtol, 1.0E-30) 
          IF (delta > threshold) THEN
             test_passed = .FALSE.
             nfail = nfail + 1
             IF (verbose_output .OR. fail_printed < max_fail_print) THEN
                WRITE(*,*) 'FAIL: ', label, ' mismatch at index (',j,',',k,')'
                WRITE(*,*) '  Ref:', ref_val, ' New:', new(j,k), ' Diff:', delta
                fail_printed = fail_printed + 1
             ELSEIF (.NOT. fail_suppression_noted) THEN
                WRITE(*,*) 'NOTE: Further failure details suppressed (set FSPS_TEST_VERBOSE=1 to expand).' 
                fail_suppression_noted = .TRUE.
             END IF
          END IF
       END DO
    END DO
    WRITE(*,*) 'SUMMARY: ', label, ' max abs diff=', max_delta, ' at (', mj, ',', mk, ')', &
         ' max rel diff=', max_rel
  END SUBROUTINE CHECK_ARRAY_2D

  SUBROUTINE CHECK_ARRAY_1D(label, ref, new, d1)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: d1
    REAL(SP), DIMENSION(d1), INTENT(IN) :: ref, new
   REAL(SP) :: delta, threshold, max_delta, max_rel, ref_val
   INTEGER :: j, mj

   max_delta = 0.0
   max_rel = 0.0
   mj = 1

    DO j = 1, d1
       IF (IEEE_IS_NAN(ref(j)) .OR. IEEE_IS_NAN(new(j))) THEN
          WRITE(*,*) 'FAIL: ', label, ' contains NaN at index (',j,')'
          test_passed = .FALSE.
          RETURN
       END IF
       
       ref_val = ref(j)
       delta = ABS(ref_val - new(j))
       if (ABS(ref_val) > 0.0_SP) then
          max_rel = MAX(max_rel, delta / ABS(ref_val))
       end if
       IF (delta > max_delta) THEN
          max_delta = delta
          mj = j
       END IF
       threshold = MAX(ABS(ref_val) * rtol, 1.0E-30)
       IF (delta > threshold) THEN
          test_passed = .FALSE.
          nfail = nfail + 1
          IF (verbose_output .OR. fail_printed < max_fail_print) THEN
             WRITE(*,*) 'FAIL: ', label, ' mismatch at index (',j,')'
             WRITE(*,*) '  Ref:', ref_val, ' New:', new(j), ' Diff:', delta
             fail_printed = fail_printed + 1
          ELSEIF (.NOT. fail_suppression_noted) THEN
             WRITE(*,*) 'NOTE: Further failure details suppressed (set FSPS_TEST_VERBOSE=1 to expand).' 
             fail_suppression_noted = .TRUE.
          END IF
       END IF
    END DO
    WRITE(*,*) 'SUMMARY: ', label, ' max abs diff=', max_delta, ' at (', mj, ')', &
         ' max rel diff=', max_rel
  END SUBROUTINE CHECK_ARRAY_1D

  SUBROUTINE CHECK_MAGS_1D(label, ref_mag, new_mag, d1)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: d1
    REAL(SP), DIMENSION(d1), INTENT(IN) :: ref_mag, new_mag
    REAL(SP) :: delta, threshold, max_delta, max_rel, ref_val
    REAL(SP) :: flux_ref, flux_new, exp_ref, exp_new
    INTEGER :: j, mj

    max_delta = 0.0_SP
    max_rel = 0.0_SP
    mj = 1

    DO j = 1, d1
       IF (IEEE_IS_NAN(ref_mag(j)) .OR. IEEE_IS_NAN(new_mag(j))) THEN
          WRITE(*,*) 'FAIL: ', label, ' contains NaN at index (',j,')'
          test_passed = .FALSE.
          RETURN
       END IF

       ! Convert magnitudes to linear flux units (relative scale)
       exp_ref = -0.4_SP * ref_mag(j) * LOG(10.0_SP)
       exp_new = -0.4_SP * new_mag(j) * LOG(10.0_SP)
       IF (exp_ref < -700.0_SP) THEN
          flux_ref = 0.0_SP
       ELSE IF (exp_ref > 700.0_SP) THEN
          flux_ref = HUGE(1.0_SP)
       ELSE
          flux_ref = EXP(exp_ref)
       END IF
       IF (exp_new < -700.0_SP) THEN
          flux_new = 0.0_SP
       ELSE IF (exp_new > 700.0_SP) THEN
          flux_new = HUGE(1.0_SP)
       ELSE
          flux_new = EXP(exp_new)
       END IF

       ref_val = flux_ref
       delta = ABS(flux_ref - flux_new)
       IF (ABS(ref_val) > 0.0_SP) THEN
          max_rel = MAX(max_rel, delta / ABS(ref_val))
       END IF
       IF (delta > max_delta) THEN
          max_delta = delta
          mj = j
       END IF
       threshold = MAX(ABS(ref_val) * rtol, 1.0E-30_SP)
       IF (delta > threshold) THEN
          test_passed = .FALSE.
          nfail = nfail + 1
          IF (verbose_output .OR. fail_printed < max_fail_print) THEN
             WRITE(*,*) 'FAIL: ', label, ' mismatch at index (',j,')'
             WRITE(*,*) '  Ref mag:', ref_mag(j), ' New mag:', new_mag(j)
             WRITE(*,*) '  Ref flux:', flux_ref, ' New flux:', flux_new, ' Diff:', delta
             fail_printed = fail_printed + 1
          ELSEIF (.NOT. fail_suppression_noted) THEN
             WRITE(*,*) 'NOTE: Further failure details suppressed (set FSPS_TEST_VERBOSE=1 to expand).'
             fail_suppression_noted = .TRUE.
          END IF
       END IF
    END DO
    WRITE(*,*) 'SUMMARY: ', label, ' max abs diff=', max_delta, ' at (', mj, ')', &
         ' max rel diff=', max_rel
  END SUBROUTINE CHECK_MAGS_1D

  SUBROUTINE CHECK_VAL(label, idx, r, n)
    CHARACTER(*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: idx
    REAL(SP), INTENT(IN) :: r, n
    REAL(SP) :: delta, threshold

    IF (IEEE_IS_NAN(r) .OR. IEEE_IS_NAN(n)) THEN
       WRITE(*,*) 'FAIL: ', label, ' contains NaN at step ', idx
       test_passed = .FALSE.
       RETURN
    END IF

    delta = ABS(r - n)
    threshold = MAX(ABS(r) * rtol, 1.0E-30)

      IF (delta > threshold) THEN
       test_passed = .FALSE.
          nfail = nfail + 1
       IF (verbose_output .OR. fail_printed < max_fail_print) THEN
          WRITE(*,*) 'FAIL: ', label, ' mismatch at step ', idx
          WRITE(*,*) '  Ref:', r, ' New:', n, ' Diff:', delta
          fail_printed = fail_printed + 1
       ELSEIF (.NOT. fail_suppression_noted) THEN
          WRITE(*,*) 'NOTE: Further failure details suppressed (set FSPS_TEST_VERBOSE=1 to expand).' 
          fail_suppression_noted = .TRUE.
       END IF
    END IF
  END SUBROUTINE CHECK_VAL

   SUBROUTINE DUMP_STATE(label, ctx, pset)
      CHARACTER(*), INTENT(IN) :: label
      TYPE(fsps_context_t), INTENT(IN) :: ctx
      TYPE(PARAMS), INTENT(IN) :: pset
      WRITE(*,*) '--- STATE:', TRIM(label)
      WRITE(*,*) '  imf_type=', ctx%imf_type_val
      WRITE(*,*) '  imf_lower_limit=', ctx%state%imf_lower_limit, ' imf_upper_limit=', ctx%state%imf_upper_limit
      WRITE(*,*) '  dust_type=', ctx%dust_type_val, ' add_dust_emission=', ctx%add_dust_emission_val
      WRITE(*,*) '  add_neb_emission=', ctx%add_neb_emission_val, ' nebemlineinspec=', ctx%nebemlineinspec_val
      WRITE(*,*) '  interpolation_type=', ctx%interpolation_type_val, ' tiny_logt=', ctx%tiny_logt_val
      WRITE(*,*) '  pset: sfh=', pset%sfh, ' tau=', pset%tau, ' const=', pset%const, ' fburst=', pset%fburst
      WRITE(*,*) '  pset: sf_start=', pset%sf_start, ' sf_trunc=', pset%sf_trunc, ' tburst=', pset%tburst
      WRITE(*,*) '  pset: dust1=', pset%dust1, ' dust2=', pset%dust2, ' zred=', pset%zred
      WRITE(*,*) '  dims: ntfull=', ctx%state%ntfull, ' nspec=', ctx%state%nspec, ' nbands=', ctx%state%nbands
   END SUBROUTINE DUMP_STATE

   SUBROUTINE DUMP_SSP_SUMMARY(label, ctx, mass_ssp, lbol_ssp)
      CHARACTER(*), INTENT(IN) :: label
      TYPE(fsps_context_t), INTENT(IN) :: ctx
      REAL(SP), DIMENSION(:), INTENT(IN) :: mass_ssp, lbol_ssp
      INTEGER :: n
      n = SIZE(mass_ssp)
      WRITE(*,*) '--- SSP SUMMARY:', TRIM(label)
      WRITE(*,*) '  mass_ssp(1)=', mass_ssp(1), ' mass_ssp(n)=', mass_ssp(n)
      WRITE(*,*) '  lbol_ssp(1)=', lbol_ssp(1), ' lbol_ssp(n)=', lbol_ssp(n)
      WRITE(*,*) '  time_full(1)=', ctx%state%time_full(1), ' time_full(n)=', ctx%state%time_full(n)
   END SUBROUTINE DUMP_SSP_SUMMARY

   SUBROUTINE DUMP_CSP_SUMMARY(label, ocompsp)
      CHARACTER(*), INTENT(IN) :: label
      TYPE(COMPSPOUT), DIMENSION(:), INTENT(IN) :: ocompsp
      INTEGER :: n
      n = SIZE(ocompsp)
      WRITE(*,*) '--- CSP SUMMARY:', TRIM(label)
      WRITE(*,*) '  age(1)=', ocompsp(1)%age, ' age(n)=', ocompsp(n)%age
      WRITE(*,*) '  mass_csp(1)=', ocompsp(1)%mass_csp, ' mass_csp(n)=', ocompsp(n)%mass_csp
      WRITE(*,*) '  lbol_csp(1)=', ocompsp(1)%lbol_csp, ' lbol_csp(n)=', ocompsp(n)%lbol_csp
   END SUBROUTINE DUMP_CSP_SUMMARY

   SUBROUTINE ENSURE_OUTPUT_DIR(path)
      CHARACTER(*), INTENT(IN) :: path
      INTEGER :: last_slash, ierr
      LOGICAL :: exists
      CHARACTER(:), ALLOCATABLE :: dir_path

      ! If the path provided is the directory (ends in /)
      IF (path(LEN_TRIM(path):LEN_TRIM(path)) == '/') THEN
         dir_path = path(1:LEN_TRIM(path)-1)
      ELSE
         ! Original logic for "path/to/file" strings
         last_slash = SCAN(path, '/', BACK=.TRUE.)
         IF (last_slash > 0) THEN
            dir_path = path(1:last_slash-1)
         ELSE
            RETURN ! No directory part found
         END IF
      END IF
         
      INQUIRE(FILE=dir_path, EXIST=exists)
      
      IF (.NOT. exists) THEN
         WRITE(*,*) 'Pre-test check: Creating missing directory: ', dir_path
         CALL EXECUTE_COMMAND_LINE('mkdir -p ' // dir_path, EXITSTAT=ierr)
         IF (ierr /= 0) THEN
            WRITE(*,*) 'FATAL: Could not create output directory.'
            STOP 1
         END IF
      END IF
   END SUBROUTINE ENSURE_OUTPUT_DIR

END PROGRAM TEST_RUNNER
