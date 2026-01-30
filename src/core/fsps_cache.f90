MODULE FSPS_CACHE
  USE sps_vars, ONLY: SP, ndim_logt, ndim_logg, ndim_wmb_logt, ndim_wmb_logg, &
       n_agb_o, n_agb_c, n_agb_car, ndim_pagb, ndim_wr, ntau_dagb, nteff_dagb, &
       nemline, nebnz, nebnage, nebnip, nagndust, nm
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: fsps_setup_cache_t, fsps_cache_get_setup, fsps_cache_release_setup

  TYPE :: fsps_setup_cache_t
     CHARACTER(LEN=512) :: key = ''
     INTEGER :: refcount = 0
     INTEGER :: nz = 0
     INTEGER :: nt = 0
     INTEGER :: nspec = 0
     INTEGER :: nzinit = 0
     INTEGER :: nbands = 0
     INTEGER :: nindx = 0
     INTEGER :: ntfull = 0
     INTEGER :: nspec_xrb = 0
     INTEGER :: nt_xrb = 0
     INTEGER :: nz_xrb = 0
     CHARACTER(30) :: alt_filter_file = ''
     CHARACTER(LEN=64) :: isoc_type = ''
     CHARACTER(LEN=64) :: spec_type = ''
     CHARACTER(LEN=64) :: dust_type = ''
     INTEGER :: setup_nebular_gaussians = 0
     INTEGER :: smooth_velocity = 0
     INTEGER :: add_neb_emission = 0
     INTEGER :: add_neb_continuum = 0
     INTEGER :: add_dust_emission = 0
     INTEGER :: add_agn_dust = 0
     INTEGER :: add_xrb_emission = 0
     INTEGER :: add_agb_dust_model = 0
     INTEGER :: use_wr_spectra = 0

   REAL(SP), POINTER :: indexdefined(:,:) => NULL()
   REAL(SP), POINTER :: wgdust(:,:,:,:) => NULL()
   REAL(SP), POINTER :: g03smcextn(:) => NULL()
   REAL(SP), POINTER :: bands(:,:) => NULL()
   REAL(SP), POINTER :: magsun(:) => NULL()
   REAL(SP), POINTER :: magvega(:) => NULL()
   REAL(SP), POINTER :: filter_leff(:) => NULL()
   REAL(SP), POINTER :: vega_spec(:) => NULL()
   REAL(SP), POINTER :: sun_spec(:) => NULL()
   REAL(SP), POINTER :: spec_lambda(:) => NULL()
   REAL(SP), POINTER :: spec_nu(:) => NULL()
   REAL(SP), POINTER :: spec_res(:) => NULL()
   REAL(KIND(1.0)), POINTER :: speclib(:,:,:,:) => NULL()
   REAL(KIND(1.0)), POINTER :: wmb_spec(:,:,:,:) => NULL()
   REAL(SP), POINTER :: agb_spec_o(:,:) => NULL()
   REAL(SP), POINTER :: agb_logt_o(:,:) => NULL()
   REAL(SP), POINTER :: agb_spec_c(:,:) => NULL()
   REAL(SP), POINTER :: agb_logt_c(:) => NULL()
   REAL(SP), POINTER :: agb_spec_car(:,:) => NULL()
   REAL(SP), POINTER :: pagb_spec(:,:,:) => NULL()
   REAL(SP), POINTER :: wrn_spec(:,:,:) => NULL()
   REAL(SP), POINTER :: wrc_spec(:,:,:) => NULL()
   REAL(SP), POINTER :: qpaharr(:) => NULL()
   REAL(SP), POINTER :: uminarr(:) => NULL()
   REAL(SP), POINTER :: lambda_dustem(:) => NULL()
   REAL(SP), POINTER :: dustem_dustem(:,:) => NULL()
   REAL(SP), POINTER :: dustem2_dustem(:,:,:) => NULL()
   REAL(SP), POINTER :: flux_dagb(:,:,:,:) => NULL()
   REAL(SP), POINTER :: nebem_cont(:,:,:,:) => NULL()
   REAL(SP), POINTER :: xnebem_cont(:,:,:,:) => NULL()
   REAL(SP), POINTER :: neb_res_min(:) => NULL()
   REAL(SP), POINTER :: gaussnebarr(:,:) => NULL()
   REAL(SP), POINTER :: agndust_spec(:,:) => NULL()
   REAL(SP), POINTER :: mact_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: logl_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: logt_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: logg_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: ffco_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: phase_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: mini_isoc(:,:,:) => NULL()
   REAL(SP), POINTER :: lmdot_isoc(:,:,:) => NULL()
   INTEGER, POINTER :: nmass_isoc(:,:) => NULL()
   REAL(SP), POINTER :: timestep_isoc(:,:) => NULL()
   REAL(SP), POINTER :: zlegend(:) => NULL()
   REAL(SP), POINTER :: zlegendinit(:) => NULL()
   REAL(SP), POINTER :: bpass_spec_ssp(:,:,:) => NULL()
   REAL(SP), POINTER :: bpass_mass_ssp(:,:) => NULL()
   REAL(SP), POINTER :: lam_xrb(:) => NULL()
   REAL(SP), POINTER :: spec_xrb(:,:,:) => NULL()
   REAL(SP), POINTER :: ages_xrb(:) => NULL()
   REAL(SP), POINTER :: zmet_xrb(:) => NULL()
   REAL(SP), POINTER :: time_full(:) => NULL()
  END TYPE fsps_setup_cache_t

  TYPE(fsps_setup_cache_t), TARGET, ALLOCATABLE :: setup_cache(:)

CONTAINS

  SUBROUTINE fsps_cache_get_setup(key, entry, is_new)
    CHARACTER(LEN=*), INTENT(IN) :: key
    TYPE(fsps_setup_cache_t), POINTER :: entry
    LOGICAL, INTENT(OUT) :: is_new

    INTEGER :: i, empty_slot
    TYPE(fsps_setup_cache_t), ALLOCATABLE :: tmp(:)

    is_new = .FALSE.
    empty_slot = 0
    NULLIFY(entry)

    IF (ALLOCATED(setup_cache)) THEN
       DO i=1,SIZE(setup_cache)
          IF (LEN_TRIM(setup_cache(i)%key) == 0) THEN
             IF (empty_slot == 0) empty_slot = i
          ELSE IF (TRIM(setup_cache(i)%key) == TRIM(key)) THEN
             setup_cache(i)%refcount = setup_cache(i)%refcount + 1
             entry => setup_cache(i)
             RETURN
          ENDIF
       ENDDO
    ENDIF

    is_new = .TRUE.

    IF (ALLOCATED(setup_cache)) THEN
       IF (empty_slot > 0) THEN
          setup_cache(empty_slot)%key = TRIM(key)
          setup_cache(empty_slot)%refcount = 1
          entry => setup_cache(empty_slot)
          RETURN
       ENDIF
       ALLOCATE(tmp(SIZE(setup_cache)+1))
       tmp(1:SIZE(setup_cache)) = setup_cache
       CALL MOVE_ALLOC(tmp, setup_cache)
    ELSE
       ALLOCATE(setup_cache(1))
    ENDIF

    setup_cache(SIZE(setup_cache))%key = TRIM(key)
    setup_cache(SIZE(setup_cache))%refcount = 1
    entry => setup_cache(SIZE(setup_cache))
  END SUBROUTINE fsps_cache_get_setup

  SUBROUTINE fsps_cache_release_setup(entry)
    TYPE(fsps_setup_cache_t), POINTER :: entry

    IF (.NOT.ASSOCIATED(entry)) RETURN

    entry%refcount = entry%refcount - 1
    IF (entry%refcount <= 0) THEN
       CALL fsps_cache_setup_clear(entry)
    ENDIF

    NULLIFY(entry)
  END SUBROUTINE fsps_cache_release_setup

  SUBROUTINE fsps_cache_setup_clear(entry)
    TYPE(fsps_setup_cache_t), INTENT(INOUT) :: entry

   IF (ASSOCIATED(entry%indexdefined)) DEALLOCATE(entry%indexdefined)
   IF (ASSOCIATED(entry%wgdust)) DEALLOCATE(entry%wgdust)
   IF (ASSOCIATED(entry%g03smcextn)) DEALLOCATE(entry%g03smcextn)
   IF (ASSOCIATED(entry%bands)) DEALLOCATE(entry%bands)
   IF (ASSOCIATED(entry%magsun)) DEALLOCATE(entry%magsun)
   IF (ASSOCIATED(entry%magvega)) DEALLOCATE(entry%magvega)
   IF (ASSOCIATED(entry%filter_leff)) DEALLOCATE(entry%filter_leff)
   IF (ASSOCIATED(entry%vega_spec)) DEALLOCATE(entry%vega_spec)
   IF (ASSOCIATED(entry%sun_spec)) DEALLOCATE(entry%sun_spec)
   IF (ASSOCIATED(entry%spec_lambda)) DEALLOCATE(entry%spec_lambda)
   IF (ASSOCIATED(entry%spec_nu)) DEALLOCATE(entry%spec_nu)
   IF (ASSOCIATED(entry%spec_res)) DEALLOCATE(entry%spec_res)
   IF (ASSOCIATED(entry%speclib)) DEALLOCATE(entry%speclib)
   IF (ASSOCIATED(entry%wmb_spec)) DEALLOCATE(entry%wmb_spec)
   IF (ASSOCIATED(entry%agb_spec_o)) DEALLOCATE(entry%agb_spec_o)
   IF (ASSOCIATED(entry%agb_logt_o)) DEALLOCATE(entry%agb_logt_o)
   IF (ASSOCIATED(entry%agb_spec_c)) DEALLOCATE(entry%agb_spec_c)
   IF (ASSOCIATED(entry%agb_logt_c)) DEALLOCATE(entry%agb_logt_c)
   IF (ASSOCIATED(entry%agb_spec_car)) DEALLOCATE(entry%agb_spec_car)
   IF (ASSOCIATED(entry%pagb_spec)) DEALLOCATE(entry%pagb_spec)
   IF (ASSOCIATED(entry%wrn_spec)) DEALLOCATE(entry%wrn_spec)
   IF (ASSOCIATED(entry%wrc_spec)) DEALLOCATE(entry%wrc_spec)
   IF (ASSOCIATED(entry%qpaharr)) DEALLOCATE(entry%qpaharr)
   IF (ASSOCIATED(entry%uminarr)) DEALLOCATE(entry%uminarr)
   IF (ASSOCIATED(entry%lambda_dustem)) DEALLOCATE(entry%lambda_dustem)
   IF (ASSOCIATED(entry%dustem_dustem)) DEALLOCATE(entry%dustem_dustem)
   IF (ASSOCIATED(entry%dustem2_dustem)) DEALLOCATE(entry%dustem2_dustem)
   IF (ASSOCIATED(entry%flux_dagb)) DEALLOCATE(entry%flux_dagb)
   IF (ASSOCIATED(entry%nebem_cont)) DEALLOCATE(entry%nebem_cont)
   IF (ASSOCIATED(entry%xnebem_cont)) DEALLOCATE(entry%xnebem_cont)
   IF (ASSOCIATED(entry%neb_res_min)) DEALLOCATE(entry%neb_res_min)
   IF (ASSOCIATED(entry%gaussnebarr)) DEALLOCATE(entry%gaussnebarr)
   IF (ASSOCIATED(entry%agndust_spec)) DEALLOCATE(entry%agndust_spec)
   IF (ASSOCIATED(entry%mact_isoc)) DEALLOCATE(entry%mact_isoc)
   IF (ASSOCIATED(entry%logl_isoc)) DEALLOCATE(entry%logl_isoc)
   IF (ASSOCIATED(entry%logt_isoc)) DEALLOCATE(entry%logt_isoc)
   IF (ASSOCIATED(entry%logg_isoc)) DEALLOCATE(entry%logg_isoc)
   IF (ASSOCIATED(entry%ffco_isoc)) DEALLOCATE(entry%ffco_isoc)
   IF (ASSOCIATED(entry%phase_isoc)) DEALLOCATE(entry%phase_isoc)
   IF (ASSOCIATED(entry%mini_isoc)) DEALLOCATE(entry%mini_isoc)
   IF (ASSOCIATED(entry%lmdot_isoc)) DEALLOCATE(entry%lmdot_isoc)
   IF (ASSOCIATED(entry%nmass_isoc)) DEALLOCATE(entry%nmass_isoc)
   IF (ASSOCIATED(entry%timestep_isoc)) DEALLOCATE(entry%timestep_isoc)
   IF (ASSOCIATED(entry%zlegend)) DEALLOCATE(entry%zlegend)
   IF (ASSOCIATED(entry%zlegendinit)) DEALLOCATE(entry%zlegendinit)
   IF (ASSOCIATED(entry%bpass_spec_ssp)) DEALLOCATE(entry%bpass_spec_ssp)
   IF (ASSOCIATED(entry%bpass_mass_ssp)) DEALLOCATE(entry%bpass_mass_ssp)
   IF (ASSOCIATED(entry%lam_xrb)) DEALLOCATE(entry%lam_xrb)
   IF (ASSOCIATED(entry%spec_xrb)) DEALLOCATE(entry%spec_xrb)
   IF (ASSOCIATED(entry%ages_xrb)) DEALLOCATE(entry%ages_xrb)
   IF (ASSOCIATED(entry%zmet_xrb)) DEALLOCATE(entry%zmet_xrb)
   IF (ASSOCIATED(entry%time_full)) DEALLOCATE(entry%time_full)

    entry%key = ''
    entry%refcount = 0
    entry%nz = 0
    entry%nt = 0
    entry%nspec = 0
    entry%nzinit = 0
    entry%nbands = 0
    entry%nindx = 0
    entry%ntfull = 0
    entry%nspec_xrb = 0
    entry%nt_xrb = 0
    entry%nz_xrb = 0
  END SUBROUTINE fsps_cache_setup_clear

END MODULE FSPS_CACHE
