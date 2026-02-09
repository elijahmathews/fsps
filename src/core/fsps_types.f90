module fsps_types
    !> @brief
    !> Core data structures used across FSPS.
    !>
    !> @details
    !> Defines parameter sets and output containers used throughout the library.

    use fsps_precision, only: WP
    use fsps_constants, only: NPZPHOT

    implicit none
    save
    private

    public :: PARAMS
    public :: COMPSPOUT
    public :: SFHPARAMS
    public :: TLSF
    public :: OBSDAT
    public :: TPZPHOT

    ! ---------------------------------------------------------------------
    ! Derived types
    ! ---------------------------------------------------------------------

    !> @brief Parameter set used to generate a model.
    type PARAMS
        real(WP) :: pagb = 1.0, dell = 0.0, delt = 0.0, fbhb = 0.0, sbss = 0.0, tau = 1.0, &
                    const = 0.0, tage = 0.0, fburst = 0.0, tburst = 11.0, dust1 = 0.0, dust2 = 0.0, &
                    logzsol = 0.0, zred = 0.0, pmetals = 0.02, imf1 = 1.3, imf2 = 2.3, imf3 = 2.3, &
                    vdmc = 0.08, dust_clumps = -99.0, frac_nodust = 0.0, dust_index = -0.7, &
                    dust_tesc = 7.0, frac_obrun = 0.0, uvb = 1.0, mwr = 3.1, redgb = 1.0, agb = 1.0, &
                    dust1_index = -1.0, mdave = 0.5, sf_start = 0.0, sf_trunc = 0.0, sf_slope = 0.0, &
                    duste_gamma = 0.01, duste_umin = 1.0, duste_qpah = 3.5, fcstar = 1.0, &
                    masscut = 150.0, sigma_smooth = 0.0, agb_dust = 1.0, min_wave_smooth = 1.0e3, &
                    max_wave_smooth = 1.0e4, gas_logu = -2.0, gas_logz = 0.0, igm_factor = 1.0, &
                    fagn = 0.0, agn_tau = 10.0, frac_xrb = 1.0, dust3 = 0.0
        integer :: zmet = 1, sfh = 0, wgp1 = 1, wgp2 = 1, wgp3 = 1, evtype = -1
        integer, allocatable :: mag_compute(:)
        integer, allocatable :: ssp_gen_age(:)
        character(50) :: imf_filename = '', sfh_filename = ''
    end type PARAMS

    !> @brief Output container for CSP results.
    type COMPSPOUT
        real(WP) :: age = 0.0, mass_csp = 0.0, lbol_csp = 0.0, sfr = 0.0, mdust = 0.0, mformed = 0.0
        real(WP), allocatable :: mags(:)
        real(WP), allocatable :: spec(:)
        real(WP), allocatable :: indx(:)
        real(WP), allocatable :: emlines(:)
    end type COMPSPOUT

    !> @brief SFH parameters converted to intrinsic units.
    type SFHPARAMS
        real(WP) :: tau = 1.0, tage = 0.0, tburst = 0.0, sf_trunc = 0.0, sf_slope = 0.0, &
                    tq = 0.0, t0 = 0.0, tb = 0.0
        integer :: type = 0, use_simha_limits = 0
    end type SFHPARAMS

    !> @brief Line-spread function container.
    type TLSF
        real(WP), allocatable :: lsf(:)
        real(WP) :: minlam = 0.0, maxlam = 0.0
    end type TLSF

    !> @brief Observational data container.
    type OBSDAT
        real(WP) :: zred = 0.0, logsmass = 0.0
        real(WP), allocatable :: mags(:), magerr(:)
        real(WP), allocatable :: spec(:), specerr(:)
    end type OBSDAT

    !> @brief Photometric redshift grid container.
    type TPZPHOT
        real(WP), dimension(NPZPHOT) :: zz = 0.0, pz = 0.0
    end type TPZPHOT

end module fsps_types
