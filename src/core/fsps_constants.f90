module fsps_constants
    !> @brief
    !> Defines global constants, precision parameters, and array dimensions for FSPS.
    !>
    !> @details
    !> This module consolidates all physical constants, array limits, and configuration
    !> parameters used throughout the FSPS codebase. It serves as the root dependency
    !> for precision definitions (SP) and is used by almost all other modules.

    ! Use the modern ISO standard for precision definitions
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    public

    ! ------------------------------------------------------------------------
    ! PRECISION & CONFIGURATION
    ! ------------------------------------------------------------------------

    !> @brief Precision definition (64-bit Double Precision).
    !> @details We map the legacy FSPS name "SP" to the standard real64.
    integer, parameter :: SP = real64

    !> @brief Verbosity Flags
    integer, parameter :: VERBOSE_MINIMAL = 0
    integer, parameter :: VERBOSE_DEBUG   = 1

    !> @brief Controls the level of output to screen.
    integer, parameter :: VERBOSE = VERBOSE_MINIMAL

    ! ------------------------------------------------------------------------
    ! PHYSICAL CONSTANTS (FUNDAMENTAL)
    ! ------------------------------------------------------------------------
    ! Defined first so they can be used to calculate derived constants below.

    !> @brief π (machine precision).
    real(sp), parameter :: PI = acos(-1.0_sp)

    !> @brief Planck constant [erg s] (2022 CODATA; exact).
    real(sp), parameter :: H_PLANCK = 6.62607015e-27_sp

    !> @brief Boltzmann constant [erg K⁻¹] (2022 CODATA; exact).
    real(sp), parameter :: K_BOLTZMANN = 1.380649e-16_sp

    !> @brief Speed of light in vacuum [Å s⁻¹] (2022 CODATA; exact).
    real(sp), parameter :: C_LIGHT = 2.99792458e18_sp

    !> @brief Newtonian constant of gravitation [cm³ g⁻¹ s⁻²] (2022 CODATA; estimated).
    real(sp), parameter :: G_NEWTON = 6.6743e-8_sp

    !> @brief Nominal solar radius [cm] (IAU 2015 Resolution B3; exact).
    real(sp), parameter :: R_SOL = 6.957e10_sp

    !> @brief Nominal solar luminosity [erg s⁻¹] (IAU 2015 Resolution B3; exact).
    real(sp), parameter :: L_SOL = 3.828e33_sp

    !> @brief Nominal solar mass parameter [cm³ s⁻²] (IAU 2015 Resolution B3; exact).
    real(sp), parameter :: M_SOL_PAR = 1.3271244e26_sp

    !> @brief Astronomical unit [cm] (IAU 2012 Resolution B2; exact).
    real(sp), parameter :: A_U = 1.495978707e13_sp

    ! ------------------------------------------------------------------------
    ! PHYSICAL CONSTANTS (DERIVED)
    ! ------------------------------------------------------------------------

    !> @brief Parsec [cm].
    real(sp), parameter :: PARSEC = (648000.0_sp / PI) * A_U

    !> @brief Solar mass [g].
    real(sp), parameter :: M_SOL = M_SOL_PAR / G_NEWTON

    !> @brief Stefan-Boltzmann constant [erg cm⁻² s⁻¹ K⁻⁴].
    real(sp), parameter :: SIGMA_SB = (2.0_sp * PI**5 * K_BOLTZMANN**4) / &
                                      (15.0_sp * (C_LIGHT * 1.0e-8_sp)**2 * H_PLANCK**3)

    !> @brief Coefficient to calculate Surface Gravity (g) from L, M, and Teff.
    !> Derivation: From g = GM/R² and L = 4πR²σT⁴, we substitute R² to get:
    !>             g = (4π * G * σ_SB) * (M * T⁴ / L)
    !> This constant contains the (4π * G * σ_SB) term, adjusted for Solar units.
    real(sp), parameter :: GRAVITY_L_M_T_COEFF = 4.0_sp * PI * SIGMA_SB * G_NEWTON * (M_SOL / L_SOL)

    !> @brief Reference log-flux for Absolute Bolometric Magnitude [log10(erg s⁻¹ cm⁻²)].
    !> Represents the log10 flux of 1.0 Solar Luminosity at a distance of 10 parsecs.
    !> Used to convert theoretical Luminosity to Absolute Magnitude (M_bol).
    real(sp), parameter :: ABS_MAG_ZEROPOINT_LOG = log10(L_SOL / (4.0_sp * PI * (10.0_sp * PARSEC)**2))

    !> @brief Julian year [s].
    real(sp), parameter :: YEAR_TO_SECOND = 3.15576e7_sp

    ! ------------------------------------------------------------------------
    ! FSPS SPECIFIC PARAMETERS
    ! ------------------------------------------------------------------------

    !> @brief Turn-on time for BHB and SBS phases, time is in log(yrs).
    real(sp), parameter :: BHB_SBS_TIME = 9.5_sp

    !> @brief The factor by which we increase the time array.
    integer, parameter :: TIME_RES_INCR = 1

    !> @brief Carbon star library selection.
    integer, parameter :: CSTAR_ARINGER = 1

    !> @brief BaSeL library normalization: 'wlbc' = Teff-color relations.
    character(4), parameter :: BASEL_STR = 'wlbc'

    ! ------------------------------------------------------------------------
    ! ARRAY DIMENSIONS
    ! ------------------------------------------------------------------------

    !> @brief Max dimension of array for each isochrone.
    integer, parameter :: NM            = 2000

    !> @brief Max number of lines to read in.
    integer, parameter :: NLINES        = 1000000

    !> @brief Max number of lines in tabulated SFH, LSF.
    integer, parameter :: NTABMAX       = 20000

    !> @brief Dimensions of BaSeL library (Log T and Log g).
    integer, parameter :: NDIM_LOGT     = 68
    integer, parameter :: NDIM_LOGG     = 19

    !> @brief Number of O-rich, C-rich AGB spectra, and Aringer C-rich spectra.
    integer, parameter :: N_AGB_O       = 9
    integer, parameter :: N_AGB_C       = 5
    integer, parameter :: N_AGB_CAR     = 9

    !> @brief Number of post-AGB spectra.
    integer, parameter :: NDIM_PAGB     = 14

    !> @brief Number of WR spectra.
    integer, parameter :: NDIM_WR       = 12

    !> @brief Dimensions of WMBasic grid.
    integer, parameter :: NDIM_WMB_LOGT = 11
    integer, parameter :: NDIM_WMB_LOGG = 3

    !> @brief Parameters for circumstellar dust models (tau and Teff points).
    integer, parameter :: NTAU_DAGB     = 50
    integer, parameter :: NTEFF_DAGB    = 6

    !> @brief Number of emission lines and continuum emission points.
    integer, parameter :: NEMLINE       = 166
    integer, parameter :: NLAM_NEBCONT  = 1963
    integer, parameter :: NEBNZ         = 11

    !> @brief Number of metallicity, age, and ionization parameter points (nebular).
    integer, parameter :: NEBNAGE       = 10
    integer, parameter :: NEBNIP        = 7
    integer, parameter :: NAGNDUST      = 9

    !> @brief Number of spectral points in the input library (AGN).
    integer, parameter :: NAGNDUST_SPEC = 125

    !> @brief Structure for using P(z) in chi2 (legacy/unused).
    integer, parameter :: NPZPHOT = 200

    ! ------------------------------------------------------------------------
    ! NUMERICAL LIMITS
    ! ------------------------------------------------------------------------

    real(sp), parameter :: SAFE_FLOOR   = tiny(1.0_sp)

end module fsps_constants