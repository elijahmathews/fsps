PROGRAM LEGACY_BATCH
  ! Legacy FSPS (upstream Fortran) evaluated on the parameter vectors of demo_batch_sed.py, one
  ! vector at a time, as python-fsps does with zcontinuous = 1: the SSPs of all metallicities are
  ! computed once, then per vector ZTINTERP to its logzsol and COMPSP on that single-metallicity
  ! set (nebular emission, SSP sum with the given SSP weights, dust, dust emission, IGM, AGN,
  ! magnitudes at its redshift), and GETINDX on the result.
  !
  ! The SSP weights enter through the use_given_weights patch of csp_gen (weights.sed): they
  ! replace the weights csp_gen would compute from the SFH parameters; nothing else changes.
  !
  ! Usage: legacy_batch settings.txt vectors.bin out.bin
  !   settings.txt  key=value lines (sps_vars globals, PARAMS fields), applied before SPS_SETUP
  !   vectors.bin   int32 n, int32 ntfull, float64 params(12, n) (order of FSPS_ARRAYS in
  !                 demo_batch_sed.py), float64 weights(ntfull, n)
  !   out.bin       int32 n, nspec, nbands, nindx, nemline, ntfull; float64 time_full,
  !                 spec_lambda, mags(nbands, n), spec(nspec, n), indx(nindx, n),
  !                 emlines(nemline, n), mass(n), lbol(n), mdust(n), then timings (s): setup,
  !                 SSPs of all metallicities, per-vector ztinterp, compsp, getindx (totals)

  USE sps_vars
  USE sps_utils
  IMPLICIT NONE

  INTEGER, PARAMETER :: NPAR = 12
  REAL(SP), DIMENSION(nspec,ntfull,1) :: spec1
  REAL(SP), DIMENSION(ntfull,1) :: mass1, lbol1
  REAL(SP), ALLOCATABLE :: par(:,:), wts(:,:), mags(:,:), spec(:,:), indx(:,:), eml(:,:)
  REAL(SP), ALLOCATABLE :: mass(:), lbol(:), mdust(:)
  TYPE(PARAMS) :: pset, base
  TYPE(COMPSPOUT), DIMENSION(ntfull) :: ocompsp
  CHARACTER(256) :: fset, fvec, fout, line
  CHARACTER(100) :: dummy_out = ''
  INTEGER(4) :: n, nt_in
  INTEGER :: u, m, z, ios
  INTEGER(8) :: c0, c1, rate
  REAL(8) :: t_setup, t_ssp, t_zt = 0, t_csp = 0, t_idx = 0

  CALL GET_COMMAND_ARGUMENT(1, fset)
  CALL GET_COMMAND_ARGUMENT(2, fvec)
  CALL GET_COMMAND_ARGUMENT(3, fout)
  CALL SYSTEM_CLOCK(count_rate=rate)

  ! Settings (globals and the PARAMS fields shared by every vector)
  OPEN(NEWUNIT=u, FILE=TRIM(fset), STATUS='OLD')
  DO
     READ(u, '(A)', IOSTAT=ios) line
     IF (ios /= 0) EXIT
     IF (LEN_TRIM(line) == 0 .OR. line(1:1) == '#') CYCLE
     CALL apply(TRIM(line))
  END DO
  CLOSE(u)

  OPEN(NEWUNIT=u, FILE=TRIM(fvec), STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
  READ(u) n, nt_in
  IF (nt_in /= ntfull) THEN
     WRITE(*,*) 'legacy_batch: ntfull of the input (', nt_in, ') differs from legacy (', ntfull, ')'
     STOP 1
  END IF
  ALLOCATE(par(NPAR, n), wts(ntfull, n))
  READ(u) par
  READ(u) wts
  CLOSE(u)
  ALLOCATE(mags(nbands, n), spec(nspec, n), indx(nindx, n), eml(nemline, n), mass(n), lbol(n), mdust(n))

  ! Setup and the SSPs of every metallicity (once, as python-fsps with zcontinuous = 1)
  CALL SYSTEM_CLOCK(c0)
  CALL SPS_SETUP(-1)
  CALL SYSTEM_CLOCK(c1); t_setup = REAL(c1 - c0, 8) / rate
  CALL SYSTEM_CLOCK(c0)
  DO z = 1, nz
     pset%zmet = z
     CALL SSP_GEN(pset, mass_ssp_zz(:,z), lbol_ssp_zz(:,z), spec_ssp_zz(:,:,z))
  END DO
  CALL SYSTEM_CLOCK(c1); t_ssp = REAL(c1 - c0, 8) / rate

  ! Per vector: an SSP-weight CSP observed at its own redshift (sfh = 0 with tage > 0 so that
  ! COMPSP returns one output; the weights replace the SFH's in csp_gen)
  base = pset
  base%sfh = 0
  base%tage = 10**(time_full(ntfull) - 9)
  use_given_weights = 1
  DO m = 1, n
     pset = base
     pset%dust2       = par(2, m)
     pset%dust1       = par(3, m)
     pset%dust_index  = par(4, m)
     pset%fagn        = par(5, m)
     pset%agn_tau     = par(6, m)
     pset%gas_logz    = par(7, m)
     pset%gas_logu    = par(8, m)
     pset%duste_umin  = par(9, m)
     pset%duste_gamma = par(10, m)
     pset%duste_qpah  = par(11, m)
     pset%zred        = par(12, m)
     pset%logzsol     = par(1, m)
     given_weights = wts(:, m)

     CALL SYSTEM_CLOCK(c0)
     CALL ZTINTERP(pset%logzsol, spec1(:,:,1), lbol1(:,1), mass1(:,1))
     CALL SYSTEM_CLOCK(c1); t_zt = t_zt + REAL(c1 - c0, 8) / rate
     c0 = c1
     CALL COMPSP(0, 1, dummy_out, mass1, lbol1, spec1, pset, ocompsp)
     CALL SYSTEM_CLOCK(c1); t_csp = t_csp + REAL(c1 - c0, 8) / rate
     c0 = c1
     CALL GETINDX(spec_lambda, ocompsp(1)%spec, ocompsp(1)%indx)
     CALL SYSTEM_CLOCK(c1); t_idx = t_idx + REAL(c1 - c0, 8) / rate

     mags(:, m) = ocompsp(1)%mags
     spec(:, m) = ocompsp(1)%spec
     indx(:, m) = ocompsp(1)%indx
     eml(:, m) = ocompsp(1)%emlines
     mass(m) = ocompsp(1)%mass_csp
     lbol(m) = ocompsp(1)%lbol_csp
     mdust(m) = ocompsp(1)%mdust
  END DO

  OPEN(NEWUNIT=u, FILE=TRIM(fout), STATUS='REPLACE', FORM='UNFORMATTED', ACCESS='STREAM')
  WRITE(u) n, INT(nspec, 4), INT(nbands, 4), INT(nindx, 4), INT(nemline, 4), INT(ntfull, 4)
  WRITE(u) REAL(time_full, 8), REAL(spec_lambda, 8)
  WRITE(u) REAL(mags, 8), REAL(spec, 8), REAL(indx, 8), REAL(eml, 8)
  WRITE(u) REAL(mass, 8), REAL(lbol, 8), REAL(mdust, 8)
  WRITE(u) t_setup, t_ssp, t_zt, t_csp, t_idx
  CLOSE(u)
  WRITE(*,'(A,I0,A,F10.3,A)') 'legacy: ', n, ' vectors, ', 1e3 * (t_zt + t_csp + t_idx) / n, ' ms per vector'

CONTAINS

  !> One key=value: an sps_vars global or a PARAMS field (applied to pset).
  SUBROUTINE apply(kv)
    CHARACTER(*), INTENT(in) :: kv
    CHARACTER(64) :: key
    REAL(SP) :: r
    INTEGER :: k, eq

    eq = INDEX(kv, '=')
    key = kv(1:eq-1)
    READ(kv(eq+1:), *) r
    k = NINT(r)
    SELECT CASE (TRIM(key))
    CASE ('imf_type');             imf_type = k
    CASE ('dust_type');            dust_type = k
    CASE ('add_dust_emission');    add_dust_emission = k
    CASE ('add_neb_emission');     add_neb_emission = k
    CASE ('add_neb_continuum');    add_neb_continuum = k
    CASE ('nebemlineinspec');      nebemlineinspec = k
    CASE ('add_agn_dust');         add_agn_dust = k
    CASE ('add_igm_absorption');   add_igm_absorption = k
    CASE ('add_stellar_remnants'); add_stellar_remnants = k
    CASE ('add_agb_dust_model');   add_agb_dust_model = k
    CASE ('use_wr_spectra');       use_wr_spectra = k
    CASE ('tpagb_norm_type');      tpagb_norm_type = k
    CASE ('compute_vega_mags');    compute_vega_mags = k
    CASE ('vactoair_flag');        vactoair_flag = k
    CASE ('smooth_velocity');      smooth_velocity = k
    CASE ('compute_indices')       ! the fork's switch; GETINDX is called below
    CASE ('H0');                   H0 = r
    CASE ('om0');                  om0 = r
    CASE ('ol0');                  ol0 = r
    CASE ('dust_tesc');            pset%dust_tesc = r
    CASE ('dust1_index');          pset%dust1_index = r
    CASE ('frac_nodust');          pset%frac_nodust = r
    CASE ('frac_obrun');           pset%frac_obrun = r
    CASE ('igm_factor');           pset%igm_factor = r
    CASE ('sigma_smooth');         pset%sigma_smooth = r
    CASE ('min_wave_smooth');      pset%min_wave_smooth = r
    CASE ('max_wave_smooth');      pset%max_wave_smooth = r
    CASE DEFAULT
       WRITE(*,*) 'legacy_batch: unknown setting ', TRIM(key)
       STOP 1
    END SELECT
  END SUBROUTINE apply

END PROGRAM LEGACY_BATCH
