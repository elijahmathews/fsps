PROGRAM TEST_CACHE
  USE fsps_context_types, ONLY: fsps_context_t
  USE fsps_context, ONLY: fsps_context_create, fsps_context_destroy
  USE sps_utils, ONLY: SPS_SETUP
  IMPLICIT NONE

  TYPE(fsps_context_t) :: ctx1, ctx2, ctx3
  LOGICAL :: same_cache, diff_cache

  WRITE(*,*) '========================================='
  WRITE(*,*) 'FSPS CACHE REUSE TEST'
  WRITE(*,*) '========================================='

  CALL fsps_context_create(ctx1)
  CALL fsps_context_create(ctx2)
  CALL fsps_context_create(ctx3)

  CALL SPS_SETUP(ctx1, -1, 'mist', 'miles', 'DL07')
  CALL SPS_SETUP(ctx2, -1, 'mist', 'miles', 'DL07')
  CALL SPS_SETUP(ctx3, -1, 'mist', 'basel', 'DL07')

  same_cache = ASSOCIATED(ctx1%setup_cache, ctx2%setup_cache)
  diff_cache = .NOT.ASSOCIATED(ctx1%setup_cache, ctx3%setup_cache)

  IF (.NOT.(same_cache .AND. diff_cache)) THEN
     WRITE(*,*) 'Cache reuse test FAIL'
     CALL fsps_context_destroy(ctx1)
     CALL fsps_context_destroy(ctx2)
     CALL fsps_context_destroy(ctx3)
     STOP 1
  ENDIF

  CALL fsps_context_destroy(ctx1)
  CALL fsps_context_destroy(ctx2)
  CALL fsps_context_destroy(ctx3)

  WRITE(*,*) 'Cache reuse test PASS'
END PROGRAM TEST_CACHE
