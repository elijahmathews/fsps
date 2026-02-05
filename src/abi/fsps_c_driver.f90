MODULE FSPS_C_DRIVER
    USE ISO_C_BINDING
      USE fsps_precision, ONLY: WP
      USE fsps_constants, ONLY: NEMLINE
      USE fsps_types, ONLY: PARAMS, COMPSPOUT
      USE sps_utils
      USE fsps_spectral_library, ONLY: get_stellar_spectrum
      USE fsps_smoothing, ONLY: apply_smoothing
      USE fsps_cosmology, ONLY: vacuum_to_air
      USE fsps_interpolation, ONLY: find_interval
     USE fsps_context, ONLY: fsps_context_t, fsps_context_create, fsps_context_setup, &
        fsps_context_destroy, fsps_context_set_param_int, fsps_context_set_param_float, &
        fsps_context_set_param_str, fsps_context_compute_ssp, fsps_context_get_paths, &
        fsps_context_prepare_pset, fsps_context_ensure_setup
  IMPLICIT NONE

  ! 1. GLOBAL STATE POINTERS
  TYPE(PARAMS), POINTER :: global_pset => NULL()
  TYPE(COMPSPOUT), POINTER :: global_ocompsp(:) => NULL()
   INTEGER, ALLOCATABLE :: has_ssp(:)
   INTEGER, ALLOCATABLE :: has_ssp_age(:,:)
   ! Context pool for handle-based API
   TYPE(fsps_context_t), ALLOCATABLE :: ctx_pool(:)
   LOGICAL, ALLOCATABLE :: ctx_inuse(:)
   ! Driver error state
   INTEGER :: fsps_last_status = 0
   CHARACTER(LEN=256) :: fsps_last_error = ''
   INTEGER :: fsps_debug = 0
   INTEGER :: fsps_lock_state = 0
   TYPE(fsps_context_t), SAVE :: fsps_default_ctx
   LOGICAL :: fsps_default_ctx_ready = .FALSE.
   INTEGER, PARAMETER :: fsps_driver_version_major = 1
   INTEGER, PARAMETER :: fsps_driver_version_minor = 0
   INTEGER, PARAMETER :: fsps_driver_version_patch = 0

CONTAINS

   SUBROUTINE fsps_ensure_default_ctx()
      IF (.NOT. fsps_default_ctx_ready) THEN
         CALL fsps_context_create(fsps_default_ctx)
         fsps_default_ctx_ready = .TRUE.
      END IF
   END SUBROUTINE fsps_ensure_default_ctx

   SUBROUTINE fsps_context_alloc_slot(slot)
      INTEGER, INTENT(OUT) :: slot
      INTEGER :: i, n
      TYPE(fsps_context_t), ALLOCATABLE :: new_pool(:)
      LOGICAL, ALLOCATABLE :: new_inuse(:)

      IF (.NOT. ALLOCATED(ctx_pool)) THEN
         ALLOCATE(ctx_pool(1))
         ALLOCATE(ctx_inuse(1))
         ctx_inuse = .FALSE.
      END IF

      slot = 0
      DO i = 1, SIZE(ctx_pool)
         IF (.NOT. ctx_inuse(i)) THEN
            slot = i
            EXIT
         END IF
      END DO

      IF (slot == 0) THEN
         n = SIZE(ctx_pool)
         ALLOCATE(new_pool(n+1))
         ALLOCATE(new_inuse(n+1))
         new_pool(1:n) = ctx_pool
         new_inuse(1:n) = ctx_inuse
         new_inuse(n+1) = .FALSE.
         CALL MOVE_ALLOC(new_pool, ctx_pool)
         CALL MOVE_ALLOC(new_inuse, ctx_inuse)
         slot = n + 1
      END IF
   END SUBROUTINE fsps_context_alloc_slot

   SUBROUTINE fsps_ensure_legacy_state()
      INTEGER :: i
      INTEGER :: n_bands, n_t, n_tfull, n_spec, n_indx, n_z

      CALL fsps_ensure_default_ctx()
      n_bands = fsps_default_ctx%state%nbands
      n_t = fsps_default_ctx%state%nt
      n_tfull = fsps_default_ctx%state%ntfull
      n_spec = fsps_default_ctx%state%nspec
      n_indx = fsps_default_ctx%state%nindx
      n_z = fsps_default_ctx%state%nz

      IF (.NOT. ASSOCIATED(global_pset)) THEN
         ALLOCATE(global_pset)
      END IF

      IF (.NOT. ALLOCATED(global_pset%mag_compute)) THEN
         ALLOCATE(global_pset%mag_compute(n_bands))
         global_pset%mag_compute = 1
      ELSE IF (SIZE(global_pset%mag_compute) /= n_bands) THEN
         DEALLOCATE(global_pset%mag_compute)
         ALLOCATE(global_pset%mag_compute(n_bands))
         global_pset%mag_compute = 1
      END IF

      IF (.NOT. ALLOCATED(global_pset%ssp_gen_age)) THEN
         ALLOCATE(global_pset%ssp_gen_age(n_t))
         global_pset%ssp_gen_age = 1
      ELSE IF (SIZE(global_pset%ssp_gen_age) /= n_t) THEN
         DEALLOCATE(global_pset%ssp_gen_age)
         ALLOCATE(global_pset%ssp_gen_age(n_t))
         global_pset%ssp_gen_age = 1
      END IF

      IF (ASSOCIATED(global_ocompsp)) THEN
         IF (SIZE(global_ocompsp) /= n_tfull) THEN
            DO i = 1, SIZE(global_ocompsp)
               IF (ALLOCATED(global_ocompsp(i)%mags))    DEALLOCATE(global_ocompsp(i)%mags)
               IF (ALLOCATED(global_ocompsp(i)%spec))    DEALLOCATE(global_ocompsp(i)%spec)
               IF (ALLOCATED(global_ocompsp(i)%indx))    DEALLOCATE(global_ocompsp(i)%indx)
               IF (ALLOCATED(global_ocompsp(i)%emlines)) DEALLOCATE(global_ocompsp(i)%emlines)
            END DO
            DEALLOCATE(global_ocompsp)
         END IF
      END IF

      IF (.NOT. ASSOCIATED(global_ocompsp)) THEN
         ALLOCATE(global_ocompsp(n_tfull))
         DO i = 1, n_tfull
            ALLOCATE(global_ocompsp(i)%mags(n_bands))
            ALLOCATE(global_ocompsp(i)%spec(n_spec))
            ALLOCATE(global_ocompsp(i)%indx(n_indx))
            ALLOCATE(global_ocompsp(i)%emlines(NEMLINE))
         END DO
      END IF

      IF (.NOT. ALLOCATED(has_ssp)) ALLOCATE(has_ssp(n_z))
      IF (.NOT. ALLOCATED(has_ssp_age)) ALLOCATE(has_ssp_age(n_z, n_t))
      has_ssp = 0
      has_ssp_age = 0
   END SUBROUTINE fsps_ensure_legacy_state

   SUBROUTINE fsps_copy_pset_from_ctx(ctx)
      TYPE(fsps_context_t), INTENT(IN) :: ctx
      global_pset = ctx%pset
   END SUBROUTINE fsps_copy_pset_from_ctx

   ! Record a driver error or warning.
   SUBROUTINE fsps_set_error(status, message)
      INTEGER, INTENT(IN) :: status
      CHARACTER(LEN=*), INTENT(IN) :: message
      fsps_last_status = status
      fsps_last_error = message
      IF (fsps_debug /= 0) THEN
          WRITE(*,*) TRIM(message)
      END IF
   END SUBROUTINE fsps_set_error

   ! Clear error state.
   SUBROUTINE fsps_clear_error() BIND(C, name="fsps_clear_error")
      fsps_last_status = 0
      fsps_last_error = ''
   END SUBROUTINE fsps_clear_error

   ! Retrieve the last error message into a C buffer.
   SUBROUTINE fsps_get_last_error(status, c_msg, c_len) &
          BIND(C, name="fsps_get_last_error")
      INTEGER(C_INT), INTENT(OUT) :: status
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_msg
      INTEGER(C_INT), VALUE :: c_len

      status = fsps_last_status
      CALL f_to_c_string(fsps_last_error, c_msg, c_len)
   END SUBROUTINE fsps_get_last_error

   ! Enable or disable driver debug prints.
   SUBROUTINE fsps_set_debug(flag) BIND(C, name="fsps_set_debug")
      INTEGER(C_INT), VALUE :: flag
      fsps_debug = flag
   END SUBROUTINE fsps_set_debug

   SUBROUTINE fsps_lock(status) BIND(C, name="fsps_lock")
      INTEGER(C_INT), INTENT(OUT) :: status
      IF (fsps_lock_state /= 0) THEN
          status = 1
          CALL fsps_set_error(401, "[FSPS-C] Error: fsps_lock already held")
          RETURN
      END IF
      fsps_lock_state = 1
      status = 0
   END SUBROUTINE fsps_lock

   SUBROUTINE fsps_unlock(status) BIND(C, name="fsps_unlock")
      INTEGER(C_INT), INTENT(OUT) :: status
      IF (fsps_lock_state == 0) THEN
          status = 1
          CALL fsps_set_error(402, "[FSPS-C] Error: fsps_unlock without lock")
          RETURN
      END IF
      fsps_lock_state = 0
      status = 0
   END SUBROUTINE fsps_unlock

   SUBROUTINE fsps_get_driver_version(major, minor, patch) &
          BIND(C, name="fsps_get_driver_version")
      INTEGER(C_INT), INTENT(OUT) :: major, minor, patch
      major = fsps_driver_version_major
      minor = fsps_driver_version_minor
      patch = fsps_driver_version_patch
   END SUBROUTINE fsps_get_driver_version

     ! -------------------------------------------------------------------------
     ! CONTEXT-BASED C API
     ! -------------------------------------------------------------------------
     SUBROUTINE fsps_context_create_handle(handle, status) &
            BIND(C, name="fsps_context_create")
       INTEGER(C_INT), INTENT(OUT) :: handle
       INTEGER(C_INT), INTENT(OUT) :: status
       INTEGER :: slot

       CALL fsps_context_alloc_slot(slot)
       CALL fsps_context_create(ctx_pool(slot))
       ctx_inuse(slot) = .TRUE.
       handle = slot
       status = 0
     END SUBROUTINE fsps_context_create_handle

     SUBROUTINE fsps_context_destroy_handle(handle, status) &
            BIND(C, name="fsps_context_destroy")
       INTEGER(C_INT), VALUE :: handle
       INTEGER(C_INT), INTENT(OUT) :: status

       status = 0
       IF (.NOT. ALLOCATED(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (handle < 1 .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

       CALL fsps_context_destroy(ctx_pool(handle))
       ctx_inuse(handle) = .FALSE.
     END SUBROUTINE fsps_context_destroy_handle

     SUBROUTINE fsps_context_setup_handle(zin, c_isoc, c_spec, c_dust, handle, status) &
            BIND(C, name="fsps_context_setup")
       INTEGER(C_INT), VALUE :: zin
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_isoc
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_spec
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_dust
       INTEGER(C_INT), VALUE :: handle
       INTEGER(C_INT), INTENT(OUT) :: status
       CHARACTER(LEN=64) :: isoc_type_in
       CHARACTER(LEN=64) :: spec_type_in
       CHARACTER(LEN=64) :: dust_type_in

       status = 0
       IF (.NOT. ALLOCATED(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (handle < 1 .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

       CALL c_to_f_string(c_isoc, isoc_type_in)
       CALL c_to_f_string(c_spec, spec_type_in)
       CALL c_to_f_string(c_dust, dust_type_in)

       IF (LEN_TRIM(isoc_type_in) == 0 .AND. LEN_TRIM(spec_type_in) == 0 .AND. &
           LEN_TRIM(dust_type_in) == 0) THEN
          CALL fsps_context_setup(ctx_pool(handle), zin)
       ELSE IF (LEN_TRIM(spec_type_in) == 0 .AND. LEN_TRIM(dust_type_in) == 0) THEN
          CALL fsps_context_setup(ctx_pool(handle), zin, TRIM(isoc_type_in))
       ELSE IF (LEN_TRIM(dust_type_in) == 0) THEN
          CALL fsps_context_setup(ctx_pool(handle), zin, TRIM(isoc_type_in), TRIM(spec_type_in))
       ELSE
          CALL fsps_context_setup(ctx_pool(handle), zin, TRIM(isoc_type_in), TRIM(spec_type_in), TRIM(dust_type_in))
       END IF
     END SUBROUTINE fsps_context_setup_handle

     SUBROUTINE fsps_context_set_int_handle(handle, c_key, value, status) &
            BIND(C, name="fsps_context_set_int")
       INTEGER(C_INT), VALUE :: handle
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
       INTEGER(C_INT), VALUE :: value
       INTEGER(C_INT), INTENT(OUT) :: status
       CHARACTER(LEN=64) :: key

       status = 0
       IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

       CALL c_to_f_string(c_key, key)
       CALL fsps_context_set_param_int(ctx_pool(handle), TRIM(key), value, status)
     END SUBROUTINE fsps_context_set_int_handle

     SUBROUTINE fsps_context_set_float_handle(handle, c_key, value, status) &
            BIND(C, name="fsps_context_set_float")
       INTEGER(C_INT), VALUE :: handle
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
       REAL(C_DOUBLE), VALUE :: value
       INTEGER(C_INT), INTENT(OUT) :: status
       CHARACTER(LEN=64) :: key

       status = 0
       IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

       CALL c_to_f_string(c_key, key)
       CALL fsps_context_set_param_float(ctx_pool(handle), TRIM(key), value, status)
     END SUBROUTINE fsps_context_set_float_handle

     SUBROUTINE fsps_context_set_str_handle(handle, c_key, c_val, status) &
            BIND(C, name="fsps_context_set_str")
       INTEGER(C_INT), VALUE :: handle
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_val
       INTEGER(C_INT), INTENT(OUT) :: status
       CHARACTER(LEN=64) :: key
       CHARACTER(LEN=128) :: val

       status = 0
       IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

       CALL c_to_f_string(c_key, key)
       CALL c_to_f_string(c_val, val)
       CALL fsps_context_set_param_str(ctx_pool(handle), TRIM(key), TRIM(val), status)
     END SUBROUTINE fsps_context_set_str_handle

     SUBROUTINE fsps_context_compute_ssp_handle(handle, c_spec, c_mass, c_lbol, status) &
            BIND(C, name="fsps_context_compute_ssp")
       INTEGER(C_INT), VALUE :: handle
       TYPE(C_PTR), VALUE :: c_spec
       TYPE(C_PTR), VALUE :: c_mass
       TYPE(C_PTR), VALUE :: c_lbol
       INTEGER(C_INT), INTENT(OUT) :: status
       REAL(C_DOUBLE), POINTER :: spec_ptr(:,:)
       REAL(C_DOUBLE), POINTER :: mass_ptr(:)
       REAL(C_DOUBLE), POINTER :: lbol_ptr(:)
      INTEGER :: n_spec, n_time

       status = 0
       IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
          status = 1
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          status = 1
          RETURN
       END IF

      n_spec = ctx_pool(handle)%state%nspec
      n_time = ctx_pool(handle)%state%ntfull
      CALL c_f_pointer(c_spec, spec_ptr, [n_spec, n_time])
      CALL c_f_pointer(c_mass, mass_ptr, [n_time])
      CALL c_f_pointer(c_lbol, lbol_ptr, [n_time])
       CALL fsps_context_compute_ssp(ctx_pool(handle), mass_ptr, lbol_ptr, spec_ptr)
     END SUBROUTINE fsps_context_compute_ssp_handle

     SUBROUTINE fsps_context_get_paths_handle(handle, c_sps, sps_len, c_data, data_len, c_out, out_len) &
            BIND(C, name="fsps_context_get_paths")
       INTEGER(C_INT), VALUE :: handle
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_sps
       INTEGER(C_INT), VALUE :: sps_len
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_data
       INTEGER(C_INT), VALUE :: data_len
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_out
       INTEGER(C_INT), VALUE :: out_len
       INTEGER(C_INT) :: status
       CHARACTER(LEN=250) :: sps_path, data_path, out_path

       status = 0
       IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
          CALL f_to_c_string('', c_sps, sps_len)
          CALL f_to_c_string('', c_data, data_len)
          CALL f_to_c_string('', c_out, out_len)
          RETURN
       END IF
       IF (.NOT. ctx_inuse(handle)) THEN
          CALL f_to_c_string('', c_sps, sps_len)
          CALL f_to_c_string('', c_data, data_len)
          CALL f_to_c_string('', c_out, out_len)
          RETURN
       END IF

       CALL fsps_context_get_paths(ctx_pool(handle), sps_path, data_path, out_path)
       CALL f_to_c_string(TRIM(sps_path), c_sps, sps_len)
       CALL f_to_c_string(TRIM(data_path), c_data, data_len)
       CALL f_to_c_string(TRIM(out_path), c_out, out_len)
     END SUBROUTINE fsps_context_get_paths_handle

    SUBROUTINE fsps_context_get_dims_handle(handle, n_spec, n_time, status) &
           BIND(C, name="fsps_context_get_dims")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_spec, n_time
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_spec = 0
      n_time = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_spec = ctx_pool(handle)%state%nspec
      n_time = ctx_pool(handle)%state%ntfull
    END SUBROUTINE fsps_context_get_dims_handle

    SUBROUTINE fsps_context_get_nspec_handle(handle, n_spec, status) &
           BIND(C, name="fsps_context_get_nspec")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_spec
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_spec = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_spec = ctx_pool(handle)%state%nspec
    END SUBROUTINE fsps_context_get_nspec_handle

    SUBROUTINE fsps_context_get_ntfull_handle(handle, n_time, status) &
           BIND(C, name="fsps_context_get_ntfull")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_time
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_time = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_time = ctx_pool(handle)%state%ntfull
    END SUBROUTINE fsps_context_get_ntfull_handle

    SUBROUTINE fsps_context_get_nbands_handle(handle, n_bands, status) &
           BIND(C, name="fsps_context_get_nbands")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_bands
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_bands = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_bands = ctx_pool(handle)%state%nbands
    END SUBROUTINE fsps_context_get_nbands_handle

    SUBROUTINE fsps_context_get_nindx_handle(handle, n_indices, status) &
           BIND(C, name="fsps_context_get_nindx")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_indices
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_indices = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_indices = ctx_pool(handle)%state%nindx
    END SUBROUTINE fsps_context_get_nindx_handle

    SUBROUTINE fsps_context_get_nz_handle(handle, n_z, status) &
           BIND(C, name="fsps_context_get_nz")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_z
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_z = 0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      n_z = ctx_pool(handle)%state%nz
    END SUBROUTINE fsps_context_get_nz_handle

    SUBROUTINE fsps_context_get_nemline_handle(handle, n_line, status) &
           BIND(C, name="fsps_context_get_nemline")
      INTEGER(C_INT), VALUE :: handle
      INTEGER(C_INT), INTENT(OUT) :: n_line
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      n_line = NEMLINE
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         n_line = 0
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         n_line = 0
         RETURN
      END IF
    END SUBROUTINE fsps_context_get_nemline_handle

    SUBROUTINE fsps_context_get_zsol_handle(handle, z_sol, status) &
           BIND(C, name="fsps_context_get_zsol")
      INTEGER(C_INT), VALUE :: handle
      REAL(C_DOUBLE), INTENT(OUT) :: z_sol
      INTEGER(C_INT), INTENT(OUT) :: status

      status = 0
      z_sol = 0.0
      IF (handle < 1 .OR. .NOT. ALLOCATED(ctx_pool) .OR. handle > SIZE(ctx_pool)) THEN
         status = 1
         RETURN
      END IF
      IF (.NOT. ctx_inuse(handle)) THEN
         status = 1
         RETURN
      END IF

      z_sol = ctx_pool(handle)%state%zsol
    END SUBROUTINE fsps_context_get_zsol_handle

#ifdef FSPS_ENABLE_LEGACY

  ! -------------------------------------------------------------------------
  ! INITIALIZATION
  ! -------------------------------------------------------------------------
   ! Initialize with default libraries.
   SUBROUTINE fsps_initialize(zin) BIND(C, name="fsps_initialize")
    INTEGER(C_INT), VALUE :: zin

      ! Call standard FSPS setup
         CALL fsps_ensure_default_ctx()
         CALL SPS_SETUP(fsps_default_ctx, zin, 'mist', 'miles', 'DL07')
      CALL fsps_initialize_state(zin)
  END SUBROUTINE fsps_initialize

   ! Initialize with explicit library selections.
   SUBROUTINE fsps_initialize_full(zin, compute_vega_mags0, vactoair_flag0, &
                                  c_isoc, c_spec, c_dust) &
       BIND(C, name="fsps_initialize_full")
    INTEGER(C_INT), VALUE :: zin
    INTEGER(C_INT), VALUE :: compute_vega_mags0
    INTEGER(C_INT), VALUE :: vactoair_flag0
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_isoc
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_spec
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_dust
    CHARACTER(LEN=64) :: isoc_type_in
    CHARACTER(LEN=64) :: spec_type_in
    CHARACTER(LEN=64) :: dust_type_in

      CALL fsps_ensure_default_ctx()
      fsps_default_ctx%compute_vega_mags_val = compute_vega_mags0
      fsps_default_ctx%vactoair_flag_val = vactoair_flag0

    CALL c_to_f_string(c_isoc, isoc_type_in)
    CALL c_to_f_string(c_spec, spec_type_in)
    CALL c_to_f_string(c_dust, dust_type_in)

    IF (LEN_TRIM(isoc_type_in) == 0 .AND. LEN_TRIM(spec_type_in) == 0 .AND. &
        LEN_TRIM(dust_type_in) == 0) THEN
       CALL SPS_SETUP(fsps_default_ctx, zin)
    ELSE IF (LEN_TRIM(spec_type_in) == 0 .AND. LEN_TRIM(dust_type_in) == 0) THEN
       CALL SPS_SETUP(fsps_default_ctx, zin, TRIM(isoc_type_in))
    ELSE IF (LEN_TRIM(dust_type_in) == 0) THEN
       CALL SPS_SETUP(fsps_default_ctx, zin, TRIM(isoc_type_in), TRIM(spec_type_in))
    ELSE
       CALL SPS_SETUP(fsps_default_ctx, zin, TRIM(isoc_type_in), TRIM(spec_type_in), TRIM(dust_type_in))
    END IF

    CALL fsps_initialize_state(zin)
  END SUBROUTINE fsps_initialize_full

  ! -------------------------------------------------------------------------
  ! PARAMETER CONTROL
  ! -------------------------------------------------------------------------
   ! Set integer parameters by name.
   SUBROUTINE fsps_set_int(c_key, val) BIND(C, name="fsps_set_int")
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
    INTEGER(C_INT), VALUE :: val
    
    CHARACTER(LEN=64) :: key
    CALL c_to_f_string(c_key, key)
      CALL fsps_ensure_default_ctx()

    SELECT CASE (TRIM(key))
    ! Globals
    CASE ('imf_type')
       fsps_default_ctx%imf_type_val = val
    CASE ('tpagb_norm_type')
       fsps_default_ctx%tpagb_norm_type_val = val
    CASE ('pzcon')
       fsps_default_ctx%pzcon_val = val
    CASE ('interpolation_type')
       fsps_default_ctx%interpolation_type_val = val
    CASE ('add_agb_dust_model')
       fsps_default_ctx%add_agb_dust_model_val = val
    CASE ('add_stellar_remnants')
       fsps_default_ctx%add_stellar_remnants_val = val
    CASE ('add_agn_dust')
       fsps_default_ctx%add_agn_dust_val = val
    CASE ('use_wr_spectra')
       fsps_default_ctx%use_wr_spectra_val = val
    CASE ('add_xrb_emission')
       fsps_default_ctx%add_xrb_emission_val = val
    CASE ('smooth_lsf')
       fsps_default_ctx%smooth_lsf_val = val
    CASE ('smoothspec_fast')
       fsps_default_ctx%smoothspec_fast_val = val
    CASE ('dust_type')
       fsps_default_ctx%dust_type_val = val
    CASE ('add_dust_emission')
       fsps_default_ctx%add_dust_emission_val = val
    CASE ('add_neb_emission')
       fsps_default_ctx%add_neb_emission_val = val
    CASE ('add_neb_continuum')
       fsps_default_ctx%add_neb_continuum_val = val
    CASE ('cloudy_dust')
       fsps_default_ctx%cloudy_dust_val = val
    CASE ('add_igm_absorption')
       fsps_default_ctx%add_igm_absorption_val = val
    CASE ('nebemlineinspec')
       fsps_default_ctx%nebemlineinspec_val = val
    CASE ('smooth_velocity')
       fsps_default_ctx%smooth_velocity_val = val
    CASE ('redshift_colors')
       fsps_default_ctx%redshift_colors_val = val
    CASE ('compute_light_ages')
       fsps_default_ctx%compute_light_ages_val = val
    CASE ('compute_vega_mags')
       fsps_default_ctx%compute_vega_mags_val = val
    CASE ('vactoair_flag')
       fsps_default_ctx%vactoair_flag_val = val
    CASE ('use_isoc_mdot')
       fsps_default_ctx%use_isoc_mdot_val = val
    CASE ('setup_nebular_gaussians')
       fsps_default_ctx%setup_nebular_gaussians_val = val

    ! PARAMS Members
    CASE ('evtype')
       global_pset%evtype = val
    CASE ('sfh')
       global_pset%sfh = val
    CASE ('zmet')
       global_pset%zmet = val
    CASE ('wgp1')
       global_pset%wgp1 = val
    CASE ('wgp2')
       global_pset%wgp2 = val
    CASE ('wgp3')
       global_pset%wgp3 = val

    CASE DEFAULT
       CALL fsps_set_error(101, "[FSPS-C] Warning: Unknown integer parameter: "//TRIM(key))
    END SELECT
  END SUBROUTINE fsps_set_int

   ! Set floating-point parameters by name.
   SUBROUTINE fsps_set_float(c_key, val) BIND(C, name="fsps_set_float")
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
    REAL(C_DOUBLE), VALUE :: val
    
    CHARACTER(LEN=64) :: key
    CALL c_to_f_string(c_key, key)

    SELECT CASE (TRIM(key))
    ! Globals
    CASE ('om0')
      om0 = REAL(val, WP)
    CASE ('ol0')
      ol0 = REAL(val, WP)
    CASE ('H0')
      H0 = REAL(val, WP)
    CASE ('tiny_logt')
      tiny_logt = REAL(val, WP)
    CASE ('imf_upper_limit')
      imf_upper_limit = REAL(val, WP)
    CASE ('imf_lower_limit')
      imf_lower_limit = REAL(val, WP)
    CASE ('logt_wmb_hot')
      logt_wmb_hot = REAL(val, WP)
    CASE ('nebular_smooth_init')
      nebular_smooth_init = REAL(val, WP)

    ! PARAMS Members - SSP
    CASE ('imf1')
      global_pset%imf1 = REAL(val, WP)
    CASE ('imf2')
      global_pset%imf2 = REAL(val, WP)
    CASE ('imf3')
      global_pset%imf3 = REAL(val, WP)
    CASE ('vdmc')
      global_pset%vdmc = REAL(val, WP)
    CASE ('mdave')
      global_pset%mdave = REAL(val, WP)
    CASE ('dell')
      global_pset%dell = REAL(val, WP)
    CASE ('delt')
      global_pset%delt = REAL(val, WP)
    CASE ('sbss')
      global_pset%sbss = REAL(val, WP)
    CASE ('fbhb')
      global_pset%fbhb = REAL(val, WP)
    CASE ('pagb')
      global_pset%pagb = REAL(val, WP)
    CASE ('agb_dust')
      global_pset%agb_dust = REAL(val, WP)
    CASE ('redgb')
      global_pset%redgb = REAL(val, WP)
    CASE ('agb')
      global_pset%agb = REAL(val, WP)
    CASE ('masscut')
      global_pset%masscut = REAL(val, WP)
    CASE ('fcstar')
      global_pset%fcstar = REAL(val, WP)
    CASE ('frac_xrb')
      global_pset%frac_xrb = REAL(val, WP)

    ! PARAMS Members - CSP
    CASE ('logzsol')
      global_pset%logzsol = REAL(val, WP)
    CASE ('tau')
      global_pset%tau = REAL(val, WP)
    CASE ('const')
      global_pset%const = REAL(val, WP)
    CASE ('tage')
      global_pset%tage = REAL(val, WP)
    CASE ('fburst')
      global_pset%fburst = REAL(val, WP)
    CASE ('tburst')
      global_pset%tburst = REAL(val, WP)
    CASE ('dust1')
      global_pset%dust1 = REAL(val, WP)
    CASE ('dust2')
      global_pset%dust2 = REAL(val, WP)
    CASE ('dust3')
      global_pset%dust3 = REAL(val, WP)
    CASE ('zred')
      global_pset%zred = REAL(val, WP)
    CASE ('pmetals')
      global_pset%pmetals = REAL(val, WP)
    CASE ('dust_clumps')
      global_pset%dust_clumps = REAL(val, WP)
    CASE ('frac_nodust')
      global_pset%frac_nodust = REAL(val, WP)
    CASE ('dust_index')
      global_pset%dust_index = REAL(val, WP)
    CASE ('dust_tesc')
      global_pset%dust_tesc = REAL(val, WP)
    CASE ('frac_obrun')
      global_pset%frac_obrun = REAL(val, WP)
    CASE ('uvb')
      global_pset%uvb = REAL(val, WP)
    CASE ('mwr')
      global_pset%mwr = REAL(val, WP)
    CASE ('dust1_index')
      global_pset%dust1_index = REAL(val, WP)
    CASE ('sf_start')
      global_pset%sf_start = REAL(val, WP)
    CASE ('sf_trunc')
      global_pset%sf_trunc = REAL(val, WP)
    CASE ('sf_slope')
      global_pset%sf_slope = REAL(val, WP)
    CASE ('duste_gamma')
      global_pset%duste_gamma = REAL(val, WP)
    CASE ('duste_umin')
      global_pset%duste_umin = REAL(val, WP)
    CASE ('duste_qpah')
      global_pset%duste_qpah = REAL(val, WP)
    CASE ('sigma_smooth')
      global_pset%sigma_smooth = REAL(val, WP)
    CASE ('min_wave_smooth')
      global_pset%min_wave_smooth = REAL(val, WP)
    CASE ('max_wave_smooth')
      global_pset%max_wave_smooth = REAL(val, WP)
    CASE ('gas_logu')
      global_pset%gas_logu = REAL(val, WP)
    CASE ('gas_logz')
      global_pset%gas_logz = REAL(val, WP)
    CASE ('igm_factor')
      global_pset%igm_factor = REAL(val, WP)
    CASE ('fagn')
      global_pset%fagn = REAL(val, WP)
    CASE ('agn_tau')
      global_pset%agn_tau = REAL(val, WP)

    CASE DEFAULT
       CALL fsps_set_error(102, "[FSPS-C] Warning: Unknown float parameter: "//TRIM(key))
    END SELECT
  END SUBROUTINE fsps_set_float

   ! Validate common parameter constraints.
   SUBROUTINE fsps_validate_params(status) BIND(C, name="fsps_validate_params")
    INTEGER(C_INT), INTENT(OUT) :: status
    REAL(WP) :: sumcb

    status = 0
       CALL fsps_ensure_default_ctx()
       IF (global_pset%zmet < 1 .OR. global_pset%zmet > fsps_default_ctx%state%nz) THEN
       status = 1
       CALL fsps_set_error(201, "[FSPS-C] Warning: zmet out of range")
    END IF
       IF (fsps_default_ctx%dust_type_val < 0 .OR. fsps_default_ctx%dust_type_val > 6) THEN
       status = 2
       CALL fsps_set_error(202, "[FSPS-C] Warning: dust_type out of range")
    END IF
       IF (fsps_default_ctx%imf_type_val < 0 .OR. fsps_default_ctx%imf_type_val > 5) THEN
       status = 3
       CALL fsps_set_error(203, "[FSPS-C] Warning: imf_type out of range")
    END IF
    IF (global_pset%tage > 0.0 .AND. global_pset%sf_start > global_pset%tage) THEN
       status = 4
       CALL fsps_set_error(204, "[FSPS-C] Warning: sf_start > tage")
    END IF
    sumcb = global_pset%const + global_pset%fburst
    IF (sumcb > 1.0) THEN
       status = 5
       CALL fsps_set_error(205, "[FSPS-C] Warning: const + fburst > 1")
    END IF
  END SUBROUTINE fsps_validate_params

   ! Set string parameters by name.
   SUBROUTINE fsps_set_str(c_key, c_val) BIND(C, name="fsps_set_str")
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_key
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_val
    CHARACTER(LEN=64) :: key
    CHARACTER(LEN=256) :: val

    CALL c_to_f_string(c_key, key)
    CALL c_to_f_string(c_val, val)

    SELECT CASE (TRIM(key))
    CASE ('imf_filename')
       global_pset%imf_filename = TRIM(val)
    CASE ('sfh_filename')
       global_pset%sfh_filename = TRIM(val)
    CASE DEFAULT
       CALL fsps_set_error(103, "[FSPS-C] Warning: Unknown string parameter: "//TRIM(key))
    END SELECT
  END SUBROUTINE fsps_set_str

  ! -------------------------------------------------------------------------
  ! COMPUTATION
  ! -------------------------------------------------------------------------
   ! Compute SSP or CSP depending on `sfh`.
   SUBROUTINE fsps_compute(c_spec) BIND(C, name="fsps_compute")
    TYPE(C_PTR), VALUE :: c_spec
    REAL(WP), POINTER :: f_spec(:,:) 
    
    ! SSP Workspace
    REAL(WP), ALLOCATABLE, TARGET :: ssp_mass(:), ssp_lbol(:)
    REAL(WP), ALLOCATABLE, TARGET :: ssp_spec(:,:)
      REAL(WP), ALLOCATABLE :: ssp_mass_zz(:,:), ssp_lbol_zz(:,:), ssp_spec_zz(:,:,:)
    
    INTEGER :: i
      INTEGER :: n_spec, n_time
    CHARACTER(LEN=128) :: junk_file = 'fsps.out'

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(301, "[FSPS-C] Error: fsps_compute called before initialize!")
       RETURN
    END IF

   CALL fsps_ensure_default_ctx()
   n_spec = fsps_default_ctx%state%nspec
   n_time = fsps_default_ctx%state%ntfull

   CALL C_F_POINTER(c_spec, f_spec, [n_spec, n_time])

    ! 1. Calculate the SSP for the current global_pset%zmet
    ! We allocate workspace because we might need to feed this into COMPSP
   ALLOCATE(ssp_mass(n_time))
   ALLOCATE(ssp_lbol(n_time))
   ALLOCATE(ssp_spec(n_spec, n_time))
    
   CALL SSP_GEN(fsps_default_ctx, global_pset, ssp_mass, ssp_lbol, ssp_spec)

    IF (global_pset%sfh .EQ. 0) THEN
       ! --- SSP Mode ---
       ! Copy directly to output buffer
       f_spec = ssp_spec
       
       ! Populate global_ocompsp so that get_mags works
         DO i=1, n_time
           global_ocompsp(i)%spec = ssp_spec(:,i)
           global_ocompsp(i)%mags = 0.0 ! Will be calc'd by get_mags
       END DO
    ELSE
       ! --- CSP Mode ---
       ! Use COMPSP to handle the complex CSP generation.
       ! This avoids guessing the changing signature of CSP_GEN.
       ! args: (ztype, n_z, outfile, mass_in, lbol_in, spec_in, pset, ocompsp_out)
       ! ztype=0 (Single Z), n_z=1
            ALLOCATE(ssp_mass_zz(n_time, 1))
            ALLOCATE(ssp_lbol_zz(n_time, 1))
            ALLOCATE(ssp_spec_zz(n_spec, n_time, 1))
            ssp_mass_zz(:,1) = ssp_mass
            ssp_lbol_zz(:,1) = ssp_lbol
            ssp_spec_zz(:,:,1) = ssp_spec
            CALL COMPSP(fsps_default_ctx, 0, 1, junk_file, ssp_mass_zz, ssp_lbol_zz, ssp_spec_zz, &
               global_pset, global_ocompsp)
       
       ! Copy result to output buffer
       DO i=1, n_time
           f_spec(:,i) = global_ocompsp(i)%spec
       END DO
    END IF

    DEALLOCATE(ssp_mass)
    DEALLOCATE(ssp_lbol)
    DEALLOCATE(ssp_spec)
    IF (ALLOCATED(ssp_mass_zz)) DEALLOCATE(ssp_mass_zz)
    IF (ALLOCATED(ssp_lbol_zz)) DEALLOCATE(ssp_lbol_zz)
    IF (ALLOCATED(ssp_spec_zz)) DEALLOCATE(ssp_spec_zz)

  END SUBROUTINE fsps_compute

   ! Compute CSP with metallicity interpolation.
   SUBROUTINE fsps_compute_csp(zcontinuous) BIND(C, name="fsps_compute_csp")
    INTEGER(C_INT), VALUE :: zcontinuous

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(302, "[FSPS-C] Error: fsps_compute_csp called before initialize!")
       RETURN
    END IF
    CALL fsps_ensure_default_ctx()
    IF (zcontinuous == 3 .AND. fsps_default_ctx%add_neb_emission_val /= 0) THEN
       CALL fsps_set_error(206, "[FSPS-C] Warning: zcontinuous=3 with nebular emission enabled")
    END IF

    CALL fsps_compute_zdep(zcontinuous)
  END SUBROUTINE fsps_compute_csp

   ! Compute and cache a single SSP at a metallicity index.
   SUBROUTINE fsps_compute_ssp(zin) BIND(C, name="fsps_compute_ssp")
    INTEGER(C_INT), VALUE :: zin
    INTEGER :: zidx
    INTEGER :: old_z

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(303, "[FSPS-C] Error: fsps_compute_ssp called before initialize!")
       RETURN
    END IF

    CALL fsps_ensure_default_ctx()
    zidx = zin
    IF (zidx < 1 .OR. zidx > fsps_default_ctx%state%nz) THEN
       CALL fsps_set_error(304, "[FSPS-C] Error: z index out of bounds in fsps_compute_ssp")
       RETURN
    END IF

    old_z = global_pset%zmet
    global_pset%zmet = zidx
   CALL SSP_GEN(fsps_default_ctx, global_pset, fsps_default_ctx%state%mass_ssp_zz(:,zidx), &
                fsps_default_ctx%state%lbol_ssp_zz(:,zidx), fsps_default_ctx%state%spec_ssp_zz(:,:,zidx))
    has_ssp(zidx) = 1
    has_ssp_age(zidx,:) = global_pset%ssp_gen_age
    global_pset%zmet = old_z
  END SUBROUTINE fsps_compute_ssp

   ! Compute and cache SSPs across the full Z grid.
   SUBROUTINE fsps_compute_ssps() BIND(C, name="fsps_compute_ssps")
    INTEGER :: zidx
      CALL fsps_ensure_default_ctx()
      DO zidx = 1, fsps_default_ctx%state%nz
       CALL fsps_compute_ssp(zidx)
    END DO
  END SUBROUTINE fsps_compute_ssps

   ! Compute CSP with metallicity interpolation mode (0-3).
   SUBROUTINE fsps_compute_zdep(ztype) BIND(C, name="fsps_compute_zdep")
    INTEGER(C_INT), VALUE :: ztype
    REAL(WP), ALLOCATABLE :: mass(:), lbol(:)
    REAL(WP), ALLOCATABLE :: spec(:,:)
      REAL(WP), ALLOCATABLE :: mass_zz(:,:), lbol_zz(:,:), spec_zz(:,:,:)
    REAL(WP) :: zpos
      INTEGER :: zlo, zmet
      INTEGER :: n_spec, n_time
    CHARACTER(LEN=128) :: junk_file = 'fsps.out'

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(305, "[FSPS-C] Error: fsps_compute_zdep called before initialize!")
       RETURN
    END IF

   CALL fsps_ensure_default_ctx()
   n_spec = fsps_default_ctx%state%nspec
   n_time = fsps_default_ctx%state%ntfull
   ALLOCATE(mass(n_time))
   ALLOCATE(lbol(n_time))
   ALLOCATE(spec(n_spec, n_time))
   ALLOCATE(mass_zz(n_time, 1))
   ALLOCATE(lbol_zz(n_time, 1))
   ALLOCATE(spec_zz(n_spec, n_time, 1))

    SELECT CASE (ztype)
    CASE (0)
       zmet = global_pset%zmet
       IF (has_ssp(zmet) == 0) CALL fsps_compute_ssp(zmet)
      CALL COMPSP(fsps_default_ctx, 0, 1, junk_file, &
              fsps_default_ctx%state%mass_ssp_zz(:,zmet:zmet), &
              fsps_default_ctx%state%lbol_ssp_zz(:,zmet:zmet), &
              fsps_default_ctx%state%spec_ssp_zz(:,:,zmet:zmet), global_pset, global_ocompsp)
    CASE (1)
       zpos = global_pset%logzsol
        zlo = MAX(MIN(find_interval(LOG10(fsps_default_ctx%state%zlegend/fsps_default_ctx%state%zsol), zpos), &
            fsps_default_ctx%state%nz-1), 1)
       DO zmet = zlo, zlo+1
          IF (has_ssp(zmet) == 0) CALL fsps_compute_ssp(zmet)
       END DO
      CALL ztinterp(fsps_default_ctx, zpos, spec, lbol, mass)
      mass_zz(:,1) = mass
      lbol_zz(:,1) = lbol
      spec_zz(:,:,1) = spec
      CALL COMPSP(fsps_default_ctx, 0, 1, junk_file, mass_zz, lbol_zz, spec_zz, global_pset, global_ocompsp)
    CASE (2)
       zpos = global_pset%logzsol
       DO zmet = 1, fsps_default_ctx%state%nz
          IF (has_ssp(zmet) == 0) CALL fsps_compute_ssp(zmet)
       END DO
      CALL ztinterp(fsps_default_ctx, zpos, spec, lbol, mass, zpow=global_pset%pmetals)
      mass_zz(:,1) = mass
      lbol_zz(:,1) = lbol
      spec_zz(:,:,1) = spec
      CALL COMPSP(fsps_default_ctx, 0, 1, junk_file, mass_zz, lbol_zz, spec_zz, global_pset, global_ocompsp)
    CASE (3)
       DO zmet = 1, fsps_default_ctx%state%nz
          IF (has_ssp(zmet) == 0) CALL fsps_compute_ssp(zmet)
       END DO
      CALL COMPSP(fsps_default_ctx, 0, fsps_default_ctx%state%nz, junk_file, &
              fsps_default_ctx%state%mass_ssp_zz, fsps_default_ctx%state%lbol_ssp_zz, &
              fsps_default_ctx%state%spec_ssp_zz, global_pset, global_ocompsp)
    CASE DEFAULT
       CALL fsps_set_error(306, "[FSPS-C] Error: Unknown ztype in fsps_compute_zdep")
    END SELECT

    DEALLOCATE(mass)
    DEALLOCATE(lbol)
    DEALLOCATE(spec)
   DEALLOCATE(mass_zz)
   DEALLOCATE(lbol_zz)
   DEALLOCATE(spec_zz)
  END SUBROUTINE fsps_compute_zdep

   ! Interpolate an SSP to a target metallicity and age.
   SUBROUTINE fsps_interp_ssp(zpos, tpos, c_spec, c_mass, c_lbol) &
       BIND(C, name="fsps_interp_ssp")
    REAL(C_DOUBLE), VALUE :: zpos, tpos
    TYPE(C_PTR), VALUE :: c_spec, c_mass, c_lbol
    REAL(WP), POINTER :: f_spec(:,:), f_mass(:), f_lbol(:)
      REAL(WP), ALLOCATABLE :: time(:)
      INTEGER :: zlo, zmet, tlo, n_spec, n_t

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(307, "[FSPS-C] Error: fsps_interp_ssp called before initialize!")
       RETURN
    END IF

       CALL fsps_ensure_default_ctx()
       n_t = fsps_default_ctx%state%nt
       n_spec = fsps_default_ctx%state%nspec
       ALLOCATE(time(n_t))
       zlo = MAX(MIN(find_interval(LOG10(fsps_default_ctx%state%zlegend/fsps_default_ctx%state%zsol), &
          REAL(zpos, WP)), fsps_default_ctx%state%nz-1), 1)
       time = fsps_default_ctx%state%timestep_isoc(zlo,:)
      tlo = MAX(MIN(find_interval(time, REAL(tpos, WP)), n_t-1), 1)

    DO zmet = zlo, zlo+1
       IF (has_ssp_age(zmet,tlo) == 0 .OR. has_ssp_age(zmet,tlo+1) == 0) THEN
          global_pset%ssp_gen_age = 0
          global_pset%ssp_gen_age(tlo:tlo+1) = 1
          CALL fsps_compute_ssp(zmet)
          global_pset%ssp_gen_age = 1
       END IF
    END DO

   CALL C_F_POINTER(c_spec, f_spec, [n_spec, 1])
    CALL C_F_POINTER(c_mass, f_mass, [1])
    CALL C_F_POINTER(c_lbol, f_lbol, [1])
   CALL ztinterp(fsps_default_ctx, REAL(zpos, WP), f_spec, f_lbol, f_mass, tpos=REAL(tpos, WP))
   DEALLOCATE(time)
  END SUBROUTINE fsps_interp_ssp

  ! -------------------------------------------------------------------------
  ! PHOTOMETRY
  ! -------------------------------------------------------------------------
   ! Compute magnitudes for all bands.
   SUBROUTINE fsps_get_mags(zred, c_mags) BIND(C, name="fsps_get_mags")
     REAL(C_DOUBLE), VALUE :: zred
     TYPE(C_PTR), VALUE :: c_mags
     REAL(WP), POINTER :: f_mags(:,:) ! (nbands, ntfull)
     
     INTEGER :: i
        INTEGER :: n_spec, n_bands, n_time
        REAL(WP), ALLOCATABLE :: tspec(:)
        INTEGER, ALLOCATABLE :: all_bands(:)
     
        CALL fsps_ensure_default_ctx()
        n_spec = fsps_default_ctx%state%nspec
        n_bands = fsps_default_ctx%state%nbands
        n_time = fsps_default_ctx%state%ntfull
        CALL C_F_POINTER(c_mags, f_mags, [n_bands, n_time])

        ALLOCATE(tspec(n_spec))
        ALLOCATE(all_bands(n_bands))
     
     all_bands = 1
     
        DO i = 1, n_time
        tspec = global_ocompsp(i)%spec
      CALL GETMAGS(fsps_default_ctx, REAL(zred, WP), tspec, f_mags(:,i), all_bands)
     END DO

        DEALLOCATE(tspec)
        DEALLOCATE(all_bands)
     
  END SUBROUTINE fsps_get_mags

   ! Compute magnitudes for a band mask.
   SUBROUTINE fsps_get_mags_mask(zred, c_mags, c_mc) BIND(C, name="fsps_get_mags_mask")
     REAL(C_DOUBLE), VALUE :: zred
     TYPE(C_PTR), VALUE :: c_mags
     TYPE(C_PTR), VALUE :: c_mc
     REAL(WP), POINTER :: f_mags(:,:) ! (nbands, ntfull)
     INTEGER(C_INT), POINTER :: f_mc(:)

     INTEGER :: i
       INTEGER :: n_spec, n_bands, n_time
       REAL(WP), ALLOCATABLE :: tspec(:)

       CALL fsps_ensure_default_ctx()
       n_spec = fsps_default_ctx%state%nspec
       n_bands = fsps_default_ctx%state%nbands
       n_time = fsps_default_ctx%state%ntfull
       CALL C_F_POINTER(c_mags, f_mags, [n_bands, n_time])
       CALL C_F_POINTER(c_mc, f_mc, [n_bands])

       ALLOCATE(tspec(n_spec))

       DO i = 1, n_time
        tspec = global_ocompsp(i)%spec
      CALL GETMAGS(fsps_default_ctx, REAL(zred, WP), tspec, f_mags(:,i), f_mc)
     END DO

       DEALLOCATE(tspec)

  END SUBROUTINE fsps_get_mags_mask

   ! Return spectra from the compsp output buffer.
   SUBROUTINE fsps_get_spec(c_spec) BIND(C, name="fsps_get_spec")
     TYPE(C_PTR), VALUE :: c_spec
     REAL(WP), POINTER :: f_spec(:,:)
        INTEGER :: i
        INTEGER :: n_spec, n_time

        CALL fsps_ensure_default_ctx()
        n_spec = fsps_default_ctx%state%nspec
        n_time = fsps_default_ctx%state%ntfull
        CALL C_F_POINTER(c_spec, f_spec, [n_spec, n_time])

        DO i = 1, n_time
        f_spec(:,i) = global_ocompsp(i)%spec
     END DO
  END SUBROUTINE fsps_get_spec

   ! Return spectra converted to Lsun/Angstrom.
   SUBROUTINE fsps_get_spec_peraa(c_spec) BIND(C, name="fsps_get_spec_peraa")
     TYPE(C_PTR), VALUE :: c_spec
     REAL(WP), POINTER :: f_spec(:,:)
     REAL(WP) :: lam
        REAL(WP) :: lamarr(1)
        INTEGER :: i, j
        INTEGER :: n_spec, n_time

        CALL fsps_ensure_default_ctx()
        n_spec = fsps_default_ctx%state%nspec
        n_time = fsps_default_ctx%state%ntfull
        CALL C_F_POINTER(c_spec, f_spec, [n_spec, n_time])

        DO j = 1, n_time
           DO i = 1, n_spec
              IF (fsps_default_ctx%vactoair_flag_val == 1) THEN
                 lamarr = vacuum_to_air(fsps_default_ctx%state%spec_lambda(i:i))
              lam = lamarr(1)
           ELSE
                 lam = fsps_default_ctx%state%spec_lambda(i)
           END IF
           f_spec(i,j) = global_ocompsp(j)%spec(i) * (3.0e18_WP / (lam*lam))
        END DO
     END DO
  END SUBROUTINE fsps_get_spec_peraa

   ! Return CSP statistics and emission lines.
   SUBROUTINE fsps_get_stats(c_age, c_mass, c_lbol, c_sfr, c_mdust, c_mformed, &
                            c_emlines) BIND(C, name="fsps_get_stats")
     TYPE(C_PTR), VALUE :: c_age, c_mass, c_lbol, c_sfr, c_mdust, c_mformed
     TYPE(C_PTR), VALUE :: c_emlines
     REAL(WP), POINTER :: f_age(:), f_mass(:), f_lbol(:), f_sfr(:), f_mdust(:), f_mformed(:)
     REAL(WP), POINTER :: f_emlines(:,:)
     INTEGER :: i
     INTEGER :: n_time

     CALL fsps_ensure_default_ctx()
     n_time = fsps_default_ctx%state%ntfull
     CALL C_F_POINTER(c_age, f_age, [n_time])
     CALL C_F_POINTER(c_mass, f_mass, [n_time])
     CALL C_F_POINTER(c_lbol, f_lbol, [n_time])
     CALL C_F_POINTER(c_sfr, f_sfr, [n_time])
     CALL C_F_POINTER(c_mdust, f_mdust, [n_time])
     CALL C_F_POINTER(c_mformed, f_mformed, [n_time])
     CALL C_F_POINTER(c_emlines, f_emlines, [NEMLINE, n_time])

     DO i = 1, n_time
        f_age(i) = global_ocompsp(i)%age
        f_mass(i) = global_ocompsp(i)%mass_csp
        f_lbol(i) = global_ocompsp(i)%lbol_csp
        f_sfr(i) = global_ocompsp(i)%sfr
        f_mdust(i) = global_ocompsp(i)%mdust
        f_mformed(i) = global_ocompsp(i)%mformed
        f_emlines(:,i) = global_ocompsp(i)%emlines
     END DO
  END SUBROUTINE fsps_get_stats

  SUBROUTINE fsps_get_indices(c_spec, c_indices) BIND(C, name="fsps_get_indices")
     TYPE(C_PTR), VALUE :: c_spec
     TYPE(C_PTR), VALUE :: c_indices
     REAL(WP), POINTER :: f_spec(:)
     REAL(WP), POINTER :: f_indices(:)
     REAL(WP), ALLOCATABLE :: lamarr(:)
     INTEGER :: n_spec, n_indx

     CALL fsps_ensure_default_ctx()
     n_spec = fsps_default_ctx%state%nspec
     n_indx = fsps_default_ctx%state%nindx
     ALLOCATE(lamarr(n_spec))

     CALL C_F_POINTER(c_spec, f_spec, [n_spec])
     CALL C_F_POINTER(c_indices, f_indices, [n_indx])

     IF (fsps_default_ctx%vactoair_flag_val == 1) THEN
        lamarr = vacuum_to_air(fsps_default_ctx%state%spec_lambda)
     ELSE
        lamarr = fsps_default_ctx%state%spec_lambda
     END IF

   CALL GETINDX(fsps_default_ctx, lamarr, f_spec, f_indices)
    DEALLOCATE(lamarr)
  END SUBROUTINE fsps_get_indices

   ! Return a stellar spectrum for stellar parameters.
   SUBROUTINE fsps_stellar_spectrum(mact, logt, lbol, logg, phase, ffco, lmdot, &
                                   wght, c_spec) BIND(C, name="fsps_stellar_spectrum")
    REAL(C_DOUBLE), VALUE :: mact, logt, lbol, logg, phase, ffco, lmdot, wght
    TYPE(C_PTR), VALUE :: c_spec
    REAL(WP), POINTER :: f_spec(:)

    IF (.NOT. ASSOCIATED(global_pset)) THEN
       CALL fsps_set_error(308, "[FSPS-C] Error: fsps_stellar_spectrum called before initialize!")
       RETURN
    END IF

   CALL fsps_ensure_default_ctx()
   CALL C_F_POINTER(c_spec, f_spec, [fsps_default_ctx%state%nspec])
   CALL get_stellar_spectrum(fsps_default_ctx, global_pset, REAL(mact, WP), REAL(logt, WP), REAL(lbol, WP), &
             REAL(logg, WP), REAL(phase, WP), REAL(ffco, WP), REAL(lmdot, WP), f_spec)
  END SUBROUTINE fsps_stellar_spectrum

  ! -------------------------------------------------------------------------
  ! TEARDOWN
  ! -------------------------------------------------------------------------
  SUBROUTINE fsps_finalize() BIND(C, name="fsps_finalize")
     INTEGER :: i
     
     IF (ASSOCIATED(global_ocompsp)) THEN
        DO i = 1, SIZE(global_ocompsp)
           IF (ALLOCATED(global_ocompsp(i)%mags))    DEALLOCATE(global_ocompsp(i)%mags)
           IF (ALLOCATED(global_ocompsp(i)%spec))    DEALLOCATE(global_ocompsp(i)%spec)
           IF (ALLOCATED(global_ocompsp(i)%indx))    DEALLOCATE(global_ocompsp(i)%indx)
           IF (ALLOCATED(global_ocompsp(i)%emlines)) DEALLOCATE(global_ocompsp(i)%emlines)
        END DO
        DEALLOCATE(global_ocompsp)
        global_ocompsp => NULL()
     END IF

     IF (ASSOCIATED(global_pset)) THEN
        IF (ALLOCATED(global_pset%mag_compute)) DEALLOCATE(global_pset%mag_compute)
        IF (ALLOCATED(global_pset%ssp_gen_age)) DEALLOCATE(global_pset%ssp_gen_age)
        DEALLOCATE(global_pset)
        global_pset => NULL()
     END IF

     IF (ALLOCATED(has_ssp)) DEALLOCATE(has_ssp)
     IF (ALLOCATED(has_ssp_age)) DEALLOCATE(has_ssp_age)

     IF (fsps_default_ctx_ready) THEN
        CALL SPS_TAKEDOWN(fsps_default_ctx)
     END IF
  END SUBROUTINE fsps_finalize

  ! -------------------------------------------------------------------------
  ! UTILS
  ! -------------------------------------------------------------------------
  SUBROUTINE fsps_get_dims(n_spec, n_time) BIND(C, name="fsps_get_dims")
    INTEGER(C_INT), INTENT(OUT) :: n_spec, n_time
      CALL fsps_ensure_default_ctx()
      n_spec = fsps_default_ctx%state%nspec
      n_time = fsps_default_ctx%state%ntfull
  END SUBROUTINE fsps_get_dims
  
  SUBROUTINE fsps_get_nbands(n_bands) BIND(C, name="fsps_get_nbands")
    INTEGER(C_INT), INTENT(OUT) :: n_bands
      CALL fsps_ensure_default_ctx()
      n_bands = fsps_default_ctx%state%nbands
  END SUBROUTINE fsps_get_nbands

   SUBROUTINE fsps_get_nindx(n_indices) BIND(C, name="fsps_get_nindx")
      INTEGER(C_INT), INTENT(OUT) :: n_indices
      CALL fsps_ensure_default_ctx()
      n_indices = fsps_default_ctx%state%nindx
   END SUBROUTINE fsps_get_nindx

   ! Dimension getters
   SUBROUTINE fsps_get_nspec(n_spec) BIND(C, name="fsps_get_nspec")
      INTEGER(C_INT), INTENT(OUT) :: n_spec
      CALL fsps_ensure_default_ctx()
      n_spec = fsps_default_ctx%state%nspec
   END SUBROUTINE fsps_get_nspec

   SUBROUTINE fsps_get_ntfull(n_time) BIND(C, name="fsps_get_ntfull")
      INTEGER(C_INT), INTENT(OUT) :: n_time
      CALL fsps_ensure_default_ctx()
      n_time = fsps_default_ctx%state%ntfull
   END SUBROUTINE fsps_get_ntfull

   SUBROUTINE fsps_get_nt(n_time) BIND(C, name="fsps_get_nt")
      INTEGER(C_INT), INTENT(OUT) :: n_time
      CALL fsps_ensure_default_ctx()
      n_time = fsps_default_ctx%state%nt
   END SUBROUTINE fsps_get_nt

   SUBROUTINE fsps_get_nm(n_mass) BIND(C, name="fsps_get_nm")
      INTEGER(C_INT), INTENT(OUT) :: n_mass
      n_mass = NM
   END SUBROUTINE fsps_get_nm

   SUBROUTINE fsps_get_ntabmax(n_tabmax) BIND(C, name="fsps_get_ntabmax")
      INTEGER(C_INT), INTENT(OUT) :: n_tabmax
      n_tabmax = NTABMAX
   END SUBROUTINE fsps_get_ntabmax

   SUBROUTINE fsps_get_nz(n_z) BIND(C, name="fsps_get_nz")
      INTEGER(C_INT), INTENT(OUT) :: n_z
      CALL fsps_ensure_default_ctx()
      n_z = fsps_default_ctx%state%nz
   END SUBROUTINE fsps_get_nz

   SUBROUTINE fsps_get_nemline(n_line) BIND(C, name="fsps_get_nemline")
      INTEGER(C_INT), INTENT(OUT) :: n_line
      n_line = NEMLINE
   END SUBROUTINE fsps_get_nemline

   ! Isochrone metadata
   SUBROUTINE fsps_get_isochrone_dimensions(n_age, n_mass) &
          BIND(C, name="fsps_get_isochrone_dimensions")
      INTEGER(C_INT), INTENT(OUT) :: n_age, n_mass
      CALL fsps_ensure_default_ctx()
      n_age = fsps_default_ctx%state%nt
      n_mass = NM
   END SUBROUTINE fsps_get_isochrone_dimensions

   SUBROUTINE fsps_get_nmass_isochrone(z_idx, t_idx, n_mass) &
          BIND(C, name="fsps_get_nmass_isochrone")
      INTEGER(C_INT), VALUE :: z_idx, t_idx
      INTEGER(C_INT), INTENT(OUT) :: n_mass
         CALL fsps_ensure_default_ctx()
         IF (z_idx < 1 .OR. z_idx > fsps_default_ctx%state%nz .OR. t_idx < 1 .OR. t_idx > fsps_default_ctx%state%nt) THEN
          n_mass = -1
          CALL fsps_set_error(207, "[FSPS-C] Warning: get_nmass_isochrone index out of range")
          RETURN
       END IF
         n_mass = fsps_default_ctx%state%nmass_isoc(z_idx, t_idx)
   END SUBROUTINE fsps_get_nmass_isochrone

   SUBROUTINE fsps_get_zsol(z_sol) BIND(C, name="fsps_get_zsol")
      REAL(C_DOUBLE), INTENT(OUT) :: z_sol
      CALL fsps_ensure_default_ctx()
      z_sol = fsps_default_ctx%state%zsol
   END SUBROUTINE fsps_get_zsol

   SUBROUTINE fsps_get_zlegend(c_zlegend) BIND(C, name="fsps_get_zlegend")
      TYPE(C_PTR), VALUE :: c_zlegend
      REAL(WP), POINTER :: f_zlegend(:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_zlegend, f_zlegend, [fsps_default_ctx%state%nz])
      f_zlegend = fsps_default_ctx%state%zlegend
   END SUBROUTINE fsps_get_zlegend

   SUBROUTINE fsps_get_timefull(c_timefull) BIND(C, name="fsps_get_timefull")
      TYPE(C_PTR), VALUE :: c_timefull
      REAL(WP), POINTER :: f_timefull(:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_timefull, f_timefull, [fsps_default_ctx%state%ntfull])
      f_timefull = fsps_default_ctx%state%time_full
   END SUBROUTINE fsps_get_timefull

   SUBROUTINE fsps_get_lambda(c_lambda) BIND(C, name="fsps_get_lambda")
      TYPE(C_PTR), VALUE :: c_lambda
      REAL(WP), POINTER :: f_lambda(:)
        CALL fsps_ensure_default_ctx()
        CALL C_F_POINTER(c_lambda, f_lambda, [fsps_default_ctx%state%nspec])
        IF (fsps_default_ctx%vactoair_flag_val == 1) THEN
           f_lambda = vacuum_to_air(fsps_default_ctx%state%spec_lambda)
      ELSE
           f_lambda = fsps_default_ctx%state%spec_lambda
      END IF
   END SUBROUTINE fsps_get_lambda

   SUBROUTINE fsps_get_emlambda(c_emlambda) BIND(C, name="fsps_get_emlambda")
      TYPE(C_PTR), VALUE :: c_emlambda
      REAL(WP), POINTER :: f_emlambda(:)
      CALL C_F_POINTER(c_emlambda, f_emlambda, [NEMLINE])
        CALL fsps_ensure_default_ctx()
        IF (fsps_default_ctx%vactoair_flag_val == 1) THEN
           f_emlambda = vacuum_to_air(fsps_default_ctx%state%nebem_line_pos)
      ELSE
           f_emlambda = fsps_default_ctx%state%nebem_line_pos
      END IF
   END SUBROUTINE fsps_get_emlambda

   SUBROUTINE fsps_get_res(c_res) BIND(C, name="fsps_get_res")
      TYPE(C_PTR), VALUE :: c_res
      REAL(WP), POINTER :: f_res(:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_res, f_res, [fsps_default_ctx%state%nspec])
      f_res = fsps_default_ctx%state%spec_res
   END SUBROUTINE fsps_get_res

   SUBROUTINE fsps_get_filter_data(c_wave_eff, c_mag_vega, c_mag_sun) &
          BIND(C, name="fsps_get_filter_data")
      TYPE(C_PTR), VALUE :: c_wave_eff, c_mag_vega, c_mag_sun
      REAL(WP), POINTER :: f_wave_eff(:), f_mag_vega(:), f_mag_sun(:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_wave_eff, f_wave_eff, [fsps_default_ctx%state%nbands])
      CALL C_F_POINTER(c_mag_vega, f_mag_vega, [fsps_default_ctx%state%nbands])
      CALL C_F_POINTER(c_mag_sun, f_mag_sun, [fsps_default_ctx%state%nbands])
      f_wave_eff = fsps_default_ctx%state%filter_leff
      f_mag_vega = fsps_default_ctx%state%magvega - fsps_default_ctx%state%magvega(1)
      f_mag_sun = fsps_default_ctx%state%magsun
   END SUBROUTINE fsps_get_filter_data

   SUBROUTINE fsps_get_ssp_weights(c_wghts) BIND(C, name="fsps_get_ssp_weights")
      TYPE(C_PTR), VALUE :: c_wghts
      REAL(WP), POINTER :: f_wghts(:,:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_wghts, f_wghts, [fsps_default_ctx%state%ntfull, fsps_default_ctx%state%nz])
      f_wghts = fsps_default_ctx%state%weight_ssp
   END SUBROUTINE fsps_get_ssp_weights

   SUBROUTINE fsps_get_csp_components(c_young, c_old) BIND(C, name="fsps_get_csp_components")
      TYPE(C_PTR), VALUE :: c_young, c_old
      REAL(WP), POINTER :: f_young(:), f_old(:)
      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_young, f_young, [fsps_default_ctx%state%nspec])
      CALL C_F_POINTER(c_old, f_old, [fsps_default_ctx%state%nspec])
      f_young = fsps_default_ctx%state%spec_young
      f_old = fsps_default_ctx%state%spec_old
   END SUBROUTINE fsps_get_csp_components

   SUBROUTINE fsps_get_ssp_spec(c_spec, c_mass, c_lbol) BIND(C, name="fsps_get_ssp_spec")
      TYPE(C_PTR), VALUE :: c_spec, c_mass, c_lbol
      REAL(WP), POINTER :: f_spec(:,:,:)
      REAL(WP), POINTER :: f_mass(:,:), f_lbol(:,:)
      INTEGER :: zidx

      CALL fsps_ensure_default_ctx()
      DO zidx = 1, fsps_default_ctx%state%nz
          IF (has_ssp(zidx) == 0) CALL fsps_compute_ssp(zidx)
      END DO

      CALL C_F_POINTER(c_spec, f_spec, [fsps_default_ctx%state%nspec, fsps_default_ctx%state%ntfull, fsps_default_ctx%state%nz])
      CALL C_F_POINTER(c_mass, f_mass, [fsps_default_ctx%state%ntfull, fsps_default_ctx%state%nz])
      CALL C_F_POINTER(c_lbol, f_lbol, [fsps_default_ctx%state%ntfull, fsps_default_ctx%state%nz])
      f_spec = fsps_default_ctx%state%spec_ssp_zz
      f_mass = fsps_default_ctx%state%mass_ssp_zz
      f_lbol = fsps_default_ctx%state%lbol_ssp_zz
   END SUBROUTINE fsps_get_ssp_spec

   ! Set a tabular SFH.
   SUBROUTINE fsps_set_sfh_tab(ntab, c_age, c_sfr, c_met) BIND(C, name="fsps_set_sfh_tab")
      INTEGER(C_INT), VALUE :: ntab
      TYPE(C_PTR), VALUE :: c_age, c_sfr, c_met
      REAL(WP), POINTER :: f_age(:), f_sfr(:), f_met(:)

      CALL C_F_POINTER(c_age, f_age, [ntab])
      CALL C_F_POINTER(c_sfr, f_sfr, [ntab])
      CALL C_F_POINTER(c_met, f_met, [ntab])

      CALL fsps_ensure_default_ctx()
      fsps_default_ctx%state%ntabsfh = ntab
      fsps_default_ctx%state%sfh_tab(1,1:fsps_default_ctx%state%ntabsfh) = f_age
      fsps_default_ctx%state%sfh_tab(2,1:fsps_default_ctx%state%ntabsfh) = f_sfr
      fsps_default_ctx%state%sfh_tab(3,1:fsps_default_ctx%state%ntabsfh) = f_met
   END SUBROUTINE fsps_set_sfh_tab

   ! Set band computation mask for mags.
   SUBROUTINE fsps_set_mag_compute(n_bands, c_mask) BIND(C, name="fsps_set_mag_compute")
       INTEGER(C_INT), VALUE :: n_bands
       TYPE(C_PTR), VALUE :: c_mask
       INTEGER(C_INT), POINTER :: f_mask(:)

       CALL C_F_POINTER(c_mask, f_mask, [n_bands])
       IF (.NOT. ASSOCIATED(global_pset)) THEN
          CALL fsps_set_error(309, "[FSPS-C] Error: fsps_set_mag_compute called before initialize!")
          RETURN
       END IF
       CALL fsps_ensure_default_ctx()
       IF (n_bands < 1 .OR. n_bands > fsps_default_ctx%state%nbands) THEN
          CALL fsps_set_error(208, "[FSPS-C] Warning: fsps_set_mag_compute size out of range")
          RETURN
       END IF
       IF (.NOT. ALLOCATED(global_pset%mag_compute)) THEN
          ALLOCATE(global_pset%mag_compute(fsps_default_ctx%state%nbands))
          global_pset%mag_compute = 0
       END IF
       global_pset%mag_compute(1:n_bands) = f_mask(1:n_bands)
     END SUBROUTINE fsps_set_mag_compute

   ! Set SSP age computation mask.
   SUBROUTINE fsps_set_ssp_gen_age(n_age, c_mask) BIND(C, name="fsps_set_ssp_gen_age")
       INTEGER(C_INT), VALUE :: n_age
       TYPE(C_PTR), VALUE :: c_mask
       INTEGER(C_INT), POINTER :: f_mask(:)

       CALL C_F_POINTER(c_mask, f_mask, [n_age])
       IF (.NOT. ASSOCIATED(global_pset)) THEN
          CALL fsps_set_error(310, "[FSPS-C] Error: fsps_set_ssp_gen_age called before initialize!")
          RETURN
       END IF
       CALL fsps_ensure_default_ctx()
       IF (n_age < 1 .OR. n_age > fsps_default_ctx%state%nt) THEN
          CALL fsps_set_error(209, "[FSPS-C] Warning: fsps_set_ssp_gen_age size out of range")
          RETURN
       END IF
       IF (.NOT. ALLOCATED(global_pset%ssp_gen_age)) THEN
          ALLOCATE(global_pset%ssp_gen_age(fsps_default_ctx%state%nt))
          global_pset%ssp_gen_age = 0
       END IF
       global_pset%ssp_gen_age(1:n_age) = f_mask(1:n_age)
     END SUBROUTINE fsps_set_ssp_gen_age

   ! Provide a wavelength-dependent LSF (sigma in km/s).
   SUBROUTINE fsps_set_ssp_lsf(nsv, c_sigma, wlo, whi) BIND(C, name="fsps_set_ssp_lsf")
      INTEGER(C_INT), VALUE :: nsv
      TYPE(C_PTR), VALUE :: c_sigma
      REAL(C_DOUBLE), VALUE :: wlo, whi
      REAL(WP), POINTER :: f_sigma(:)

      CALL C_F_POINTER(c_sigma, f_sigma, [nsv])

      CALL fsps_ensure_default_ctx()
      fsps_default_ctx%state%lsfinfo%minlam = REAL(wlo, WP)
      fsps_default_ctx%state%lsfinfo%maxlam = REAL(whi, WP)
      IF (ALLOCATED(fsps_default_ctx%state%lsfinfo%lsf)) DEALLOCATE(fsps_default_ctx%state%lsfinfo%lsf)
      ALLOCATE(fsps_default_ctx%state%lsfinfo%lsf(nsv))
      fsps_default_ctx%state%lsfinfo%lsf = f_sigma
   END SUBROUTINE fsps_set_ssp_lsf

   ! Smooth a spectrum using a Gaussian kernel.
   SUBROUTINE fsps_smooth_spectrum(c_wave, c_spec, sigma_broad, minw, maxw) &
          BIND(C, name="fsps_smooth_spectrum")
      TYPE(C_PTR), VALUE :: c_wave, c_spec
      REAL(C_DOUBLE), VALUE :: sigma_broad, minw, maxw
      REAL(WP), POINTER :: f_wave(:), f_spec(:)

      CALL fsps_ensure_default_ctx()
      CALL C_F_POINTER(c_wave, f_wave, [fsps_default_ctx%state%nspec])
      CALL C_F_POINTER(c_spec, f_spec, [fsps_default_ctx%state%nspec])

      CALL apply_smoothing(fsps_default_ctx, f_wave, f_spec, REAL(sigma_broad, WP), &
         REAL(minw, WP), REAL(maxw, WP))
   END SUBROUTINE fsps_smooth_spectrum

   ! Write isochrone data to a .cmd file.
   SUBROUTINE fsps_write_isochrone(c_outfile) BIND(C, name="fsps_write_isochrone")
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_outfile
      CHARACTER(LEN=100) :: outfile
      CALL c_to_f_string(c_outfile, outfile)
      CALL WRITE_ISOCHRONE(fsps_default_ctx, TRIM(outfile), global_pset)
   END SUBROUTINE fsps_write_isochrone

   SUBROUTINE fsps_get_setup_vars(cvms, vta_flag) BIND(C, name="fsps_get_setup_vars")
      INTEGER(C_INT), INTENT(OUT) :: cvms, vta_flag
      CALL fsps_ensure_default_ctx()
      cvms = fsps_default_ctx%compute_vega_mags_val
      vta_flag = fsps_default_ctx%vactoair_flag_val
   END SUBROUTINE fsps_get_setup_vars

   ! Return the currently-selected library names.
   SUBROUTINE fsps_get_libraries(c_isoc, c_isoc_len, c_spec, c_spec_len, &
                                                c_dust, c_dust_len) &
          BIND(C, name="fsps_get_libraries")
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_isoc
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_spec
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_dust
      INTEGER(C_INT), VALUE :: c_isoc_len, c_spec_len, c_dust_len

      CALL fsps_ensure_default_ctx()
      CALL f_to_c_string(fsps_default_ctx%state%isoc_type, c_isoc, c_isoc_len)
      CALL f_to_c_string(fsps_default_ctx%state%spec_type, c_spec, c_spec_len)
      CALL f_to_c_string(fsps_default_ctx%state%str_dustem, c_dust, c_dust_len)
   END SUBROUTINE fsps_get_libraries

#endif

   ! Convert a null-terminated C string into a Fortran fixed-length string.
   SUBROUTINE c_to_f_string(c_ptr, f_str)
    CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: c_ptr
    CHARACTER(LEN=*), INTENT(OUT) :: f_str
    INTEGER :: i
    
    f_str = ''
    i = 1
    DO WHILE (c_ptr(i) /= C_NULL_CHAR .AND. i <= LEN(f_str))
       f_str(i:i) = c_ptr(i)
       i = i + 1
    END DO
  END SUBROUTINE c_to_f_string

   ! Copy a Fortran string into a C buffer with null termination.
   SUBROUTINE f_to_c_string(f_str, c_ptr, c_len)
      CHARACTER(LEN=*), INTENT(IN) :: f_str
      CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) :: c_ptr
      INTEGER(C_INT), VALUE :: c_len
      INTEGER :: i, ncopy

      ncopy = MIN(LEN_TRIM(f_str), c_len - 1)
      DO i = 1, ncopy
          c_ptr(i) = f_str(i:i)
      END DO
      IF (ncopy + 1 <= c_len) c_ptr(ncopy+1) = C_NULL_CHAR
      IF (ncopy + 2 <= c_len) c_ptr(ncopy+2:c_len) = C_NULL_CHAR
   END SUBROUTINE f_to_c_string

END MODULE FSPS_C_DRIVER