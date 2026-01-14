MODULE TEST_UTILS
  IMPLICIT NONE
  INTEGER, SAVE :: num_failures = 0

CONTAINS

  SUBROUTINE assert_true(condition, message)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: message
    
    IF (.NOT. condition) THEN
       PRINT *, "FAILED: ", TRIM(message)
       num_failures = num_failures + 1
    ELSE
       PRINT *, "PASSED: ", TRIM(message)
    ENDIF
  END SUBROUTINE assert_true

  SUBROUTINE final_report()
    IF (num_failures > 0) THEN
       PRINT *, "Total Failures: ", num_failures
       CALL EXIT(1)
    ELSE
       PRINT *, "All tests passed!"
       CALL EXIT(0)
    ENDIF
  END SUBROUTINE final_report

END MODULE TEST_UTILS
