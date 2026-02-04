module fsps_dust
    !> @brief
    !> Handles dust attenuation, absorption, and re-emission physics.
    !>
    !> @details
    !> This module consolidates all dust-related physics into a single interface.
    !> It provides routines to:
    !> 1. Compute attenuation curves (optical depth vs wavelength) for various
    !>    standard models (Calzetti, Cardelli/CCM89, Witt & Gordon, etc.).
    !> 2. Apply attenuation to stellar and nebular spectra, accounting for
    !>    differential attenuation between birth clouds and the diffuse ISM.
    !> 3. Calculate infrared dust re-emission based on energy balance principles
    !>    (absorbing UV/optical flux and re-radiating in the IR).
    !> 4. Model specialized dust environments, including circumstellar shells
    !>    around AGB stars and dusty tori surrounding AGN.

    use fsps_constants, only: SP, C_LIGHT, NEMLINE, GRAVITY_L_M_T_COEFF, &
                              M_SOL, G_NEWTON, R_SOL, YEAR_TO_SECOND, &
                              SAFE_FLOOR, PI
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: interpolate_linear, find_interval
    use sps_utils, only: smoothspec
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    implicit none
    private

    ! Public Interface
    public :: apply_dust_attenuation_and_emission
    public :: compute_attenuation_curve
    public :: apply_agb_dust_screen
    public :: apply_agn_dust_emission
    public :: interpolate_draine_li_dust_model
    public :: calculate_dust_self_absorption
    public :: compute_circumstellar_optical_depth

    ! ------------------------------------------------------------------------
    ! CONSTANTS: General Dust Parameters
    ! ------------------------------------------------------------------------
    real(SP), parameter :: V_BAND_ANGSTROMS  = 5500.0_sp
    real(SP), parameter :: UV_BUMP_CENTER    = 2175.0_sp
    real(SP), parameter :: CALZETTI_BREAK    = 6300.0_sp

    ! ------------------------------------------------------------------------
    ! CONSTANTS: Cardelli, Clayton, & Mathis (1989) Extinction Curve Parameters
    ! ------------------------------------------------------------------------

    ! Region Boundaries (in inverse microns, x = 1/lambda)
    real(SP), parameter :: CCM_X_IR_MIN  = 0.3_sp
    real(SP), parameter :: CCM_X_OPT_MIN = 1.1_sp
    real(SP), parameter :: CCM_X_NUV_MIN = 3.3_sp
    real(SP), parameter :: CCM_X_MUV_MIN = 5.9_sp
    real(SP), parameter :: CCM_X_FUV_MIN = 8.0_sp

    ! FUV Cutoff (Clamp values for x > 12.0, i.e., lambda < 833 A)
    real(SP), parameter :: CCM_X_CUTOFF = 12.0_sp
    
    ! Infrared Parameters (0.3 <= x < 1.1)
    real(SP), parameter :: CCM_IR_A_SCALE = 0.574_sp
    real(SP), parameter :: CCM_IR_B_SCALE = -0.527_sp
    real(SP), parameter :: CCM_IR_EXP     = 1.61_sp

    ! Optical Parameters (1.1 <= x < 3.3)
    ! Polynomial coefficients for y = x - 1.82 (Powers 0 through 7)
    real(SP), parameter :: CCM_OPT_Y_SHIFT = 1.82_sp
    real(SP), parameter :: CCM_OPT_A_COEFFS(0:7) = [ &
        1.0_sp,      0.17699_sp, -0.50447_sp, -0.02427_sp, &
        0.72085_sp,  0.01979_sp, -0.77530_sp,  0.32999_sp ]
    real(SP), parameter :: CCM_OPT_B_COEFFS(0:7) = [ &
        0.0_sp,      1.41338_sp,  2.28305_sp,  1.07233_sp, &
       -5.38434_sp, -0.62251_sp,  5.30260_sp, -2.09002_sp ]

    ! Near-UV Parameters (3.3 <= x < 5.9)
    ! Base Linear Terms: C1 + C2*x
    real(SP), parameter :: CCM_NUV_A_BASE(2) = [ 1.752_sp, -0.316_sp] 
    real(SP), parameter :: CCM_NUV_B_BASE(2) = [-3.09_sp,   1.825_sp]
    ! Drude Bump Terms: Scale / ((x-Pos)**2 + Width)
    real(SP), parameter :: CCM_NUV_A_BUMP(3) = [-0.104_sp, 4.67_sp, 0.341_sp]
    real(SP), parameter :: CCM_NUV_B_BUMP(3) = [ 1.206_sp, 4.62_sp, 0.263_sp]

    ! Mid-UV Parameters (5.9 <= x < 8.0) - Additions to NUV
    ! F(x) = C2*(x-5.9)**2 + C3*(x-5.9)**3
    real(SP), parameter :: CCM_MUV_A_POLY(2) = [-0.04473_sp, -0.009779_sp]
    real(SP), parameter :: CCM_MUV_B_POLY(2) = [ 0.2130_sp,   0.1207_sp]

    ! Far-UV Parameters (x >= 8.0)
    ! Polynomials in (x - 8.0)
    real(SP), parameter :: CCM_FUV_A_COEFFS(0:3) = [-1.073_sp, -0.628_sp, 0.137_sp, -0.070_sp]
    real(SP), parameter :: CCM_FUV_B_COEFFS(0:3) = [ 13.67_sp,  4.257_sp, -0.42_sp,  0.374_sp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Calzetti et al. (2000) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(SP), parameter :: CALZ_LAM_UV_MIN = 1200.0_sp  ! 0.12 microns
    real(SP), parameter :: CALZ_LAM_BREAK  = 6300.0_sp  ! 0.63 microns
    real(SP), parameter :: CALZ_LAM_IR_MAX = 22000.0_sp ! 2.20 microns

    ! General Parameters
    real(SP), parameter :: CALZ_R_V   = 4.05_sp
    real(SP), parameter :: CALZ_SCALE = 2.659_sp      ! Scaling factor k'

    ! UV/Optical Polynomial Coefficients (0.12 <= lambda < 0.63 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(SP), parameter :: CALZ_UV_COEFFS(0:3) = [-2.156_sp, 1.509_sp, -0.198_sp, 0.011_sp]

    ! Optical/NIR Linear Coefficients (0.63 <= lambda <= 2.2 um)
    ! Linear in (1/lambda_microns): c0 + c1*x
    real(SP), parameter :: CALZ_OPT_COEFFS(0:1) = [-1.857_sp, 1.040_sp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Kriek & Conroy (2013) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------
    
    ! Drude Profile Parameters
    real(SP), parameter :: KC13_BUMP_WIDTH = 350.0_sp  ! Delta lambda (Angstroms)
    
    ! Bump Amplitude Relationship: E_b = 0.85 - 1.9 * delta
    real(SP), parameter :: KC13_AMPL_INTERCEPT = 0.85_sp
    real(SP), parameter :: KC13_AMPL_SLOPE     = 1.9_sp

    ! The base model is tied to Calzetti's specific R_V
    real(SP), parameter :: KC13_R_V_BASE = 4.05_sp

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Reddy et al. (2015) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(SP), parameter :: REDDY_LAM_UV_MIN = 1500.0_sp
    real(SP), parameter :: REDDY_LAM_BREAK  = 6000.0_sp
    real(SP), parameter :: REDDY_LAM_IR_MAX = 28500.0_sp

    ! Region Boundary (Inverse Microns)
    real(SP), parameter :: REDDY_X_UV_MAX   = 1.0e4_sp / REDDY_LAM_UV_MIN

    ! General Parameters
    real(SP), parameter :: REDDY_R_V        = 2.505_sp
    real(SP), parameter :: REDDY_OFFSET     = 2.505_sp ! Base offset added to both curves

    ! UV Polynomial Coefficients (0.15 <= lambda < 0.60 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(SP), parameter :: REDDY_UV_COEFFS(0:3) = [-5.726_sp, 4.004_sp, -0.525_sp, 0.029_sp]
    
    ! Blueward Extrapolation Value (< 0.15 um)
    real(SP), parameter :: REDDY_UV_EXTRAP_VAL = 10.36_sp

    ! Optical/NIR Polynomial Coefficients (0.60 <= lambda < 2.85 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(SP), parameter :: REDDY_OPT_COEFFS(0:3) = [-2.672_sp, -0.010_sp, 1.532_sp, -0.412_sp]
    
    ! Continuity Correction for Optical Range
    real(SP), parameter :: REDDY_OPT_CORRECTION = -0.036221981_sp

    ! --------------------------------------------------------------------------
    ! CONSTANTS: AGB Circumstellar Dust Parameters
    ! --------------------------------------------------------------------------
    
    ! Dust-to-Gas Ratios (delta)
    real(SP), parameter :: AGB_DELTA_C_RICH = 0.0025_sp
    real(SP), parameter :: AGB_DELTA_O_RICH = 0.01_sp

    ! Extinction Coefficients (kappa)
    real(SP), parameter :: AGB_KAPPA_C_RICH = 3200.0_sp ! AmC + SiC
    real(SP), parameter :: AGB_KAPPA_O_RICH = 3000.0_sp ! Silicates

    ! Inner Radius Factors (cm * L^-0.5)
    real(SP), parameter :: AGB_RIN_FACTOR_C = 1.92E12_sp ! Td = 1100 K
    real(SP), parameter :: AGB_RIN_FACTOR_O = 4.74E12_sp ! Td = 700 K

    ! Villaume et al. (2015) Period Relation Coefficients
    ! logP = A + B*logR + C*logM
    real(SP), parameter :: AGB_PER_INTERCEPT = -2.07_sp
    real(SP), parameter :: AGB_PER_SLOPE_R   = 1.94_sp
    real(SP), parameter :: AGB_PER_SLOPE_M   = -0.9_sp

    ! Expansion Velocity Parameters
    real(SP), parameter :: AGB_VEXP_INTERCEPT = -13.5_sp
    real(SP), parameter :: AGB_VEXP_SLOPE     = 0.056_sp
    real(SP), parameter :: AGB_VEXP_MIN       = 3.0_sp  ! km/s
    real(SP), parameter :: AGB_VEXP_MAX       = 15.0_sp ! km/s

    ! Vassiliadis & Wood (1993) Mass Loss Parameters
    real(SP), parameter :: VW93_MDOT_LIMIT_ISO = 1.0e-4_sp
    real(SP), parameter :: VW93_PER_THRESH     = 500.0_sp ! Days
    real(SP), parameter :: VW93_MASS_THRESH    = 2.5_sp   ! Solar Masses
    real(SP), parameter :: VW93_BASE_INTERCEPT = -11.4_sp
    real(SP), parameter :: VW93_BASE_SLOPE     = 0.0123_sp
    real(SP), parameter :: VW93_HIGH_SLOPE     = 0.0125_sp
    real(SP), parameter :: VW93_SUPERWIND_NORM = 1.93e3_sp

    ! Dust-to-Gas Ratio Scaling
    real(SP), parameter :: AGB_DTG_VEL_NORM    = 225.0_sp
    real(SP), parameter :: AGB_DTG_LUM_NORM    = 1.0e4_sp
    real(SP), parameter :: AGB_DTG_LUM_EXP     = -0.6_sp

    ! Smoothing Constants
    real(SP), parameter :: AGB_SMOOTH_SIGMA      = 1.0e4_sp
    real(SP), parameter :: AGB_SMOOTH_LAMBDA_MIN = 3.0e4_sp
    real(SP), parameter :: AGB_SMOOTH_LAMBDA_MAX = 1.0e8_sp

contains

    !> @brief
    !> Main driver for applying dust physics (attenuation and emission) to the spectrum.
    !>
    !> @details
    !> 1. Calculates the attenuation curves for birth clouds (young stars) and diffuse ISM.
    !> 2. Attenuates stellar and nebular spectra based on population age (young vs old).
    !> 3. If enabled, calculates the IR dust emission (Draine & Li 2007 models) based on 
    !>    energy balance (absorbed UV/optical flux = emitted IR flux).
    !> 4. Handles self-absorption of IR dust emission iteratively.
    !>
    !> @param[in]    ctx              FSPS context.
    !> @param[in]    settings         FSPS parameter structure.
    !> @param[in]    spec_young       Spectrum of young stars (birth cloud + diffuse).
    !> @param[in]    spec_old         Spectrum of old stars (diffuse only).
    !> @param[in]    neb_flux_young   Nebular emission line fluxes for young population.
    !> @param[in]    neb_flux_old     Nebular emission line fluxes for old population.
    !> @param[out]   spec_total_out   Final combined spectrum (L_sol/Hz).
    !> @param[out]   dust_mass        Total dust mass (M_sol).
    !> @param[out]   neb_flux_out     Final attenuated nebular line fluxes.
    subroutine apply_dust_attenuation_and_emission(ctx, settings, &
                                                   spec_young, spec_old, &
                                                   neb_flux_young, neb_flux_old, &
                                                   spec_total_out, dust_mass, neb_flux_out)
        
        type(fsps_context_t), intent(in)       :: ctx
        type(params), intent(in)               :: settings
        real(SP), dimension(:), intent(in)     :: spec_young, spec_old
        real(SP), dimension(:), intent(in)     :: neb_flux_young, neb_flux_old
        real(SP), dimension(:), intent(out)    :: spec_total_out
        real(SP), intent(out)                  :: dust_mass
        real(SP), dimension(:), intent(out)    :: neb_flux_out

        ! Local Variables
        real(SP), dimension(size(spec_young)) :: attenuation_curve_diffuse
        real(SP), dimension(size(spec_young)) :: transmission_diffuse
        real(SP), dimension(size(spec_young)) :: transmission_birth_cloud
        real(SP), dimension(size(spec_young)) :: spec_attenuated_sum
        real(SP), dimension(size(spec_young)) :: frequencies
        real(SP), dimension(size(neb_flux_young)) :: transmission_diffuse_neb
        
        real(SP) :: lum_bol_intrinsic, lum_bol_attenuated, lum_absorbed_total
        real(SP), dimension(size(spec_young)) :: dust_emission_shape, dust_emission_final
        real(SP) :: emission_norm_factor

        ! 0. Input Validation
        ! -------------------
        if (settings%uvb < 0.0_sp) return ! Should trigger error handling upstream
        if (settings%wgp1 < 1 .or. settings%wgp2 < 1) return 

        ! 1. Calculate Attenuation Curves & Transmissivities
        ! --------------------------------------------------
        
        ! A. Diffuse ISM (affects all stars)
        attenuation_curve_diffuse = compute_attenuation_curve(ctx%state%spec_lambda, &
                                                              ctx%dust_type_val, settings, ctx)
        
        if (ctx%dust_type_val == 3) then
            ! Witt & Gordon models are self-normalized
            transmission_diffuse = exp(-attenuation_curve_diffuse)
        else
            transmission_diffuse = exp(-settings%dust2 * attenuation_curve_diffuse)
        end if

        ! B. Birth Clouds (affects young stars only)
        ! Standard power-law attenuation centered at 5500A
        transmission_birth_cloud = exp(-settings%dust1 * &
                                   (ctx%state%spec_lambda / V_BAND_ANGSTROMS)**settings%dust1_index)


        ! 2. Apply Attenuation to Stellar Spectra
        ! ---------------------------------------
        ! Young Stars: Part obscured by birth cloud (1-frac_obrun), part runaways (frac_obrun).
        !              ALL young stars see diffuse dust (in standard model).
        ! Old Stars:   Only see diffuse dust (optionally multiplied by dust3 factor).
        
        spec_attenuated_sum = &
            (spec_young * transmission_birth_cloud * (1.0_sp - settings%frac_obrun) + &
             spec_young * settings%frac_obrun) + &
            (spec_old * exp(-settings%dust3 * attenuation_curve_diffuse))

        ! Apply final diffuse screen (allowing for 'frac_nodust' holes in the ISM)
        spec_total_out = spec_attenuated_sum * transmission_diffuse * (1.0_sp - settings%frac_nodust) + &
                         spec_attenuated_sum * settings%frac_nodust


        ! 3. Apply Attenuation to Nebular Lines
        ! -------------------------------------
        ! Note: We must interpolate the diffuse transmission to the line wavelengths
        transmission_diffuse_neb = interpolate_linear(ctx%state%spec_lambda, &
                                                      transmission_diffuse, &
                                                      ctx%state%nebem_line_pos)
        
        ! Handle extrapolation/NaNs safely
        where (ieee_is_nan(transmission_diffuse_neb)) transmission_diffuse_neb = 1.0_sp

        ! Apply birth cloud attenuation to young nebular lines
        neb_flux_out = (neb_flux_young * &
                        exp(-settings%dust1 * (ctx%state%nebem_line_pos / V_BAND_ANGSTROMS)**settings%dust1_index) * &
                        (1.0_sp - settings%frac_obrun) + &
                        neb_flux_young * settings%frac_obrun + &
                        neb_flux_old) 
        
        ! Apply diffuse screen to total nebular flux
        neb_flux_out = neb_flux_out * transmission_diffuse_neb * (1.0_sp - settings%frac_nodust) + &
                       neb_flux_out * settings%frac_nodust


        ! 4. Add Dust Emission (Energy Balance)
        ! -------------------------------------
        if (ctx%add_dust_emission_val == 1 .and. &
            (settings%dust1 > SAFE_FLOOR .or. settings%dust2 > SAFE_FLOOR)) then
            
            frequencies = C_LIGHT / ctx%state%spec_lambda

            ! Calculate Bolometric Luminosities (L_bol)
            ! -----------------------------------------
            ! Intrinsic (Pre-Dust)
            lum_bol_intrinsic = integrate_trapezoid_array(frequencies, spec_young + spec_old)
            if (ctx%nebemlineinspec_val == 0) then
                 lum_bol_intrinsic = lum_bol_intrinsic + sum(neb_flux_young) + sum(neb_flux_old)
            end if

            ! Attenuated (Post-Dust)
            lum_bol_attenuated = integrate_trapezoid_array(frequencies, spec_total_out)
            if (ctx%nebemlineinspec_val == 0) then
                 lum_bol_attenuated = lum_bol_attenuated + sum(neb_flux_out)
            end if
            
            ! Total Energy Absorbed by Dust
            lum_absorbed_total = lum_bol_intrinsic - lum_bol_attenuated

            ! Get Dust Emission Template (Draine & Li 2007)
            ! ---------------------------------------------
            call interpolate_draine_li_dust_model(ctx, settings, dust_emission_shape)
            
            ! Normalize template area
            emission_norm_factor = integrate_trapezoid_array(frequencies, dust_emission_shape)
            
            if (emission_norm_factor <= SAFE_FLOOR) then
                dust_mass = SAFE_FLOOR
                return
            end if

            ! Calculate Self-Absorption & Final Emission
            ! ------------------------------------------
            call calculate_dust_self_absorption(frequencies, dust_emission_shape, transmission_diffuse, &
                                                lum_absorbed_total, dust_emission_final)

            ! Add to total spectrum
            spec_total_out = spec_total_out + dust_emission_final

            ! Estimate Dust Mass (Factor from Draine & Li MW3.1 model)
            ! 3.21e-3 converts Luminosity/Norm to Mass (Solar Units) roughly
            dust_mass = 3.21e-3_sp / (4.0_sp * PI) * (lum_absorbed_total / emission_norm_factor)

        else
            dust_mass = SAFE_FLOOR
        end if

    end subroutine apply_dust_attenuation_and_emission

    !> @brief
    !> Computes the attenuation curve (optical depth shape) for a given dust type.
    !>
    !> @details
    !> Calculates the optical depth \tau_\lambda normalized to the V-band (or similar
    !> normalization depending on the specific paper) for the requested dust model.
    !>
    !> Supported Dust Types:
    !> - 0: Power Law (Charlot & Fall 2000 style)
    !> - 1: Milky Way (Cardelli, Clayton, & Mathis 1989) with variable UV bump
    !> - 2: Calzetti et al. (2000) starburst attenuation curve
    !> - 3: Witt & Gordon (2000) SMC/LMC grids (Interpolated from context tables)
    !> - 4: Kriek & Conroy (2013) (Calzetti + Variable UV Bump)
    !> - 5: SMC Bar (Gordon et al. 2003) (Interpolated from context tables)
    !> - 6: Reddy et al. (2015) MOSDEF survey curve
    !>
    !> @param[in] wavelengths    Vector of wavelengths in angstroms.
    !> @param[in] dust_type_id   Integer ID of the dust model to apply.
    !> @param[in] settings       FSPS parameter structure (containing indexes, UV bump strengths, etc).
    !> @param[in] ctx            FSPS context (containing pre-loaded tables for WG00 and SMC).
    !>
    !> @return attenuation_curve Vector of optical depths (dimension matching wavelengths).
    pure function compute_attenuation_curve(wavelengths, dust_type_id, settings, ctx) result(attenuation_curve)
        real(SP), dimension(:), intent(in) :: wavelengths
        integer, intent(in)                :: dust_type_id
        type(params), intent(in)           :: settings
        type(fsps_context_t), intent(in)   :: ctx
        real(SP), dimension(size(wavelengths)) :: attenuation_curve

        ! Initialize to zero
        attenuation_curve = 0.0_sp

        select case (dust_type_id)
        
        ! --- Power Law Attenuation ---
        case (0)
            attenuation_curve = (wavelengths / V_BAND_ANGSTROMS)**settings%dust_index

        ! --- Cardelli, Clayton, & Mathis (1989) Milky Way Curve ---
        case (1)
            attenuation_curve = get_ccm89_curve(wavelengths, settings%mwr, settings%uvb)

        ! --- Calzetti et al. (2000) ---
        case (2)
            attenuation_curve = get_calzetti_curve(wavelengths)

        ! --- Witt & Gordon (2000) [Table Lookup] ---
        case (3)
            ! Direct table lookup from context. 
            ! Note: The original code implies the table matches the input wavelength grid.
            ! This dependency is retained here.
            attenuation_curve = ctx%state%wgdust(:, settings%wgp1, settings%wgp2, settings%wgp3)

        ! --- Kriek & Conroy (2013) ---
        case (4)
            attenuation_curve = get_kriek_conroy_curve(wavelengths, settings%dust_index)

        ! --- Gordon et al. (2003) SMC [Table Lookup] ---
        case (5)
            ! Direct table lookup
            attenuation_curve = ctx%state%g03smcextn

        ! --- Reddy et al. (2015) ---
        case (6)
            attenuation_curve = get_reddy_curve(wavelengths)

        case default
            ! Return NaN to signal invalid configuration in a pure context.
            ! This ensures the error propagates rather than failing silently with 0.0.
            attenuation_curve = ieee_value(1.0_sp, ieee_quiet_nan)
        end select

    end function compute_attenuation_curve

    !> @brief
    !> Applies a circumstellar dust screen to AGB stars (DUSTY models).
    !>
    !> @details
    !> 1. Estimates the optical depth (tau_1um) of the circumstellar shell based on
    !>    stellar parameters (Mass, T_eff, Luminosity, Pulsation Period).
    !> 2. Interpolates a "transfer function" (flux_out / flux_in) from the 
    !>    pre-computed DUSTY grid (Villaume et al. 2015).
    !> 3. Multiplies the input stellar spectrum by this transfer function.
    !>
    !> @param[inout]    ctx       FSPS context.
    !> @param[in]    weight       User-defined weight/scaling for the dust effect.
    !> @param[inout] spectrum     The stellar spectrum to modify (L_sol/Hz).
    !> @param[in]    mass_act     Actual stellar mass (M_sol).
    !> @param[in]    log_t        Log10 Effective Temperature (K).
    !> @param[in]    log_l        Log10 Luminosity (L_sol).
    !> @param[in]    log_g        Log10 Surface Gravity (cm/s^2).
    !> @param[in]    z_met        Metallicity (Z).
    !> @param[in]    c_o_ratio    Carbon/Oxygen ratio.
    !> @param[in]    log_mdot     Log10 Mass Loss Rate (M_sol/yr) from isochrone (optional).
    subroutine apply_agb_dust_screen(ctx, weight, spectrum, mass_act, log_t, log_l, &
                                     log_g, c_o_ratio, log_mdot)
        
        type(fsps_context_t), intent(inout)   :: ctx
        real(SP), intent(in)                  :: weight
        real(SP), dimension(:), intent(inout) :: spectrum
        real(SP), intent(in)                  :: mass_act, log_t, log_l, log_g
        real(SP), intent(in)                  :: c_o_ratio, log_mdot
        
        ! Local variables
        integer  :: c_rich_flag ! 0 = O-rich, 1 = C-rich
        integer  :: idx_teff, idx_tau
        integer  :: n_teff_grid, n_tau_grid
        real(SP) :: tau_1um, log_g_local
        real(SP) :: w_teff, w_tau ! Interpolation weights
        real(SP), dimension(size(spectrum)) :: dusty_transfer_function

        ! 1. Determine Chemistry (C-rich vs O-rich)
        ! -----------------------------------------
        if (c_o_ratio > 1.0_sp) then
            c_rich_flag = 1
        else
            c_rich_flag = 0
        end if

        ! 2. Handle Log(g) for BaSTI Isochrones
        ! -------------------------------------
        ! BaSTI does not always tabulate log(g), so we compute it physically.
        if (ctx%state%isoc_type == 'bsti') then
            log_g_local = log10(GRAVITY_L_M_T_COEFF * mass_act / (10.0_sp**log_l)) + 4.0_sp * log_t
        else
            log_g_local = log_g
        end if

        ! 3. Compute Circumstellar Optical Depth (tau at 1 micron)
        ! -------------------------------------------------------
        tau_1um = compute_circumstellar_optical_depth(ctx, c_rich_flag, mass_act, log_l, log_g_local, log_mdot)
        
        ! Apply user scaling weight
        tau_1um = tau_1um * weight

        ! Early exit if optical depth is negligible
        if (tau_1um <= SAFE_FLOOR) return


        ! 4. Interpolate DUSTY Model Grid
        ! -------------------------------
        ! Grid dimensions: flux_dagb(:, cstar+1, i_teff, i_tau)
        ! Axes: teff_dagb(cstar+1, :), tau1_dagb(cstar+1, :)

        n_teff_grid = size(ctx%state%teff_dagb, 2)
        n_tau_grid  = size(ctx%state%tau1_dagb, 2)

        ! A. Interpolate in Teff
        idx_teff = find_interval(ctx%state%teff_dagb(c_rich_flag + 1, :), 10.0_sp**log_t)
        idx_teff = max(1, min(idx_teff, n_teff_grid - 1))
        
        w_teff = (10.0_sp**log_t - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff)) / &
                 (ctx%state%teff_dagb(c_rich_flag + 1, idx_teff + 1) - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff))
        w_teff = max(-1.0_sp, min(w_teff, 1.0_sp)) ! Allow slight extrapolation

        ! B. Interpolate in Tau
        idx_tau = find_interval(ctx%state%tau1_dagb(c_rich_flag + 1, :), log10(tau_1um))
        idx_tau = max(1, min(idx_tau, n_tau_grid - 1))

        w_tau = (log10(tau_1um) - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau)) / &
                (ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau + 1) - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau))
        w_tau = max(-1.0_sp, min(w_tau, 1.0_sp))


        ! C. Bilinear Interpolation
        ! f(x, y) ~ (1-x)(1-y)F00 + x(1-y)F10 + (1-x)yF01 + xyF11
        dusty_transfer_function = &
            (1.0_sp - w_teff) * (1.0_sp - w_tau) * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff, idx_tau) + &
            w_teff            * (1.0_sp - w_tau) * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff + 1, idx_tau) + &
            w_teff            * w_tau            * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff + 1, idx_tau + 1) + &
            (1.0_sp - w_teff) * w_tau            * ctx%state%flux_dagb(:, c_rich_flag + 1, idx_teff, idx_tau + 1)


        ! 5. Apply to Spectrum
        ! --------------------
        spectrum = spectrum * dusty_transfer_function

    end subroutine apply_agb_dust_screen

    !> @brief
    !> Adds AGN emission (accretion disk + torus) to the spectrum.
    !>
    !> @details
    !> 1. Interpolates the AGN spectral templates (CLUMPY models) based on the 
    !>    AGN torus optical depth parameter (`agn_tau`).
    !> 2. Attenuates the intrinsic AGN spectrum by the host galaxy's diffuse dust.
    !> 3. Normalizes the AGN emission relative to the stellar population's 
    !>    bolometric luminosity using `fagn`.
    !> 4. Adds the result to the input spectrum in-place.
    !>
    !> @param[in]    ctx              FSPS context (contains AGN template grids).
    !> @param[in]    settings         FSPS parameter structure.
    !> @param[in]    wavelengths      Vector of wavelengths in Angstroms.
    !> @param[in]    log_lbol_stellar Log10 of the bolometric luminosity of the stellar population (L_sol).
    !> @param[inout] spectrum_inout   The spectral energy distribution to be modified (L_sol/Hz or similar).
    subroutine apply_agn_dust_emission(ctx, settings, wavelengths, log_lbol_stellar, spectrum_inout)
        type(fsps_context_t), intent(in)       :: ctx
        type(params), intent(in)               :: settings
        real(SP), dimension(:), intent(in)     :: wavelengths
        real(SP), intent(in)                   :: log_lbol_stellar
        real(SP), dimension(:), intent(inout)  :: spectrum_inout

        ! Local variables
        real(SP), dimension(size(wavelengths)) :: agn_template_interpolated
        real(SP), dimension(size(wavelengths)) :: galaxy_attenuation_curve
        real(SP) :: tau_agn_param, interpolation_weight
        real(SP) :: luminosity_agn_bolometric
        integer  :: idx_tau_grid, n_agn_grid

        ! 0. Early exit if no AGN contribution is specified
        if (settings%fagn <= tiny(0.0_sp)) return
        
        ! 1. Interpolate AGN Template based on Torus Optical Depth (agn_tau)
        ! ------------------------------------------------------------------
        ! The context stores a grid of AGN spectra varying by torus optical depth.
        ! Grid: ctx%state%agndust_spec(:, i_tau)
        ! Axis: ctx%state%agndust_tau(:)
        
        tau_agn_param = settings%agn_tau
        n_agn_grid    = size(ctx%state%agndust_tau)

        ! Use binary search to find the interval
        idx_tau_grid = find_interval(ctx%state%agndust_tau, tau_agn_param)
        
        ! Clamp index to valid range [1, N-1] for interpolation
        idx_tau_grid = max(1, min(idx_tau_grid, n_agn_grid - 1))

        ! Calculate linear interpolation weight
        interpolation_weight = (tau_agn_param - ctx%state%agndust_tau(idx_tau_grid)) / &
                               (ctx%state%agndust_tau(idx_tau_grid + 1) - ctx%state%agndust_tau(idx_tau_grid))
        
        ! Clamp weight to [0, 1] to prevent extrapolation beyond grid bounds
        interpolation_weight = max(0.0_sp, min(interpolation_weight, 1.0_sp))

        ! Interpolate the template
        agn_template_interpolated = (1.0_sp - interpolation_weight) * ctx%state%agndust_spec(:, idx_tau_grid) + &
                                    interpolation_weight * ctx%state%agndust_spec(:, idx_tau_grid + 1)


        ! 2. Attenuate AGN by Host Galaxy Diffuse Dust
        ! --------------------------------------------
        ! Calculate the shape of the galaxy's attenuation curve
        galaxy_attenuation_curve = compute_attenuation_curve(wavelengths, ctx%dust_type_val, settings, ctx)

        ! Apply the optical depth scalar (dust2)
        ! Note: Witt & Gordon models (Type 3) are self-normalized and do not use dust2.
        if (ctx%dust_type_val == 3) then
            agn_template_interpolated = agn_template_interpolated * exp(-galaxy_attenuation_curve)
        else
            agn_template_interpolated = agn_template_interpolated * exp(-settings%dust2 * galaxy_attenuation_curve)
        end if


        ! 3. Normalize and Add to Spectrum
        ! --------------------------------
        ! L_AGN = f_agn * L_bol_stellar
        ! The AGN templates in FSPS are likely pre-normalized, but we scale by the 
        ! bolometric luminosity of the current CSP generation.
        
        luminosity_agn_bolometric = (10.0_sp**log_lbol_stellar) * settings%fagn

        spectrum_inout = spectrum_inout + (luminosity_agn_bolometric * agn_template_interpolated)

    end subroutine apply_agn_dust_emission

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> Interpolates the Draine & Li (2007) dust emission grid based on Q_PAH, U_min, and Gamma.
    subroutine interpolate_draine_li_dust_model(ctx, settings, emission_spectrum)
        type(fsps_context_t), intent(in)   :: ctx
        type(params), intent(in)           :: settings
        real(SP), dimension(:), intent(out) :: emission_spectrum
        
        integer :: idx_q, idx_u
        real(SP) :: w_q, w_u, gamma_frac
        real(SP), dimension(size(emission_spectrum)) :: spec_u_min, spec_u_max
        integer :: n_qpah, n_umin

        ! Grid Dimensions
        n_qpah = ctx%state%nqpah_dustem
        n_umin = ctx%state%numin_dustem

        ! 1. Interpolate Q_PAH (Polycyclic Aromatic Hydrocarbons fraction)
        idx_q = find_interval(ctx%state%qpaharr, settings%duste_qpah)
        idx_q = max(1, min(idx_q, n_qpah - 1))
        
        w_q = (settings%duste_qpah - ctx%state%qpaharr(idx_q)) / &
              (ctx%state%qpaharr(idx_q + 1) - ctx%state%qpaharr(idx_q))
        w_q = max(0.0_sp, min(w_q, 1.0_sp))

        ! 2. Interpolate U_min (Minimum Radiation Field Intensity)
        idx_u = find_interval(ctx%state%uminarr, settings%duste_umin)
        idx_u = max(1, min(idx_u, n_umin - 1))

        w_u = (settings%duste_umin - ctx%state%uminarr(idx_u)) / &
              (ctx%state%uminarr(idx_u + 1) - ctx%state%uminarr(idx_u))
        w_u = max(0.0_sp, min(w_u, 1.0_sp))

        ! 3. Gamma Fraction (Fraction of dust in high-intensity PDRs)
        gamma_frac = max(0.0_sp, min(settings%duste_gamma, 1.0_sp))

        ! 4. Bilinear Interpolation for U_min component (Low Field)
        ! The grid `dustem2_dustem` layout is (lambda, qpah, 2*umin_index - 1) for U_min part
        spec_u_min = &
            (1-w_q)*(1-w_u) * ctx%state%dustem2_dustem(:, idx_q,   2*idx_u - 1) + &
            w_q    *(1-w_u) * ctx%state%dustem2_dustem(:, idx_q+1, 2*idx_u - 1) + &
            w_q    *w_u     * ctx%state%dustem2_dustem(:, idx_q+1, 2*(idx_u+1) - 1) + &
            (1-w_q)*w_u     * ctx%state%dustem2_dustem(:, idx_q,   2*(idx_u+1) - 1)

        ! 5. Bilinear Interpolation for U_max component (High Field PDR)
        ! The grid uses `2*umin_index` for the U_max part
        spec_u_max = &
            (1-w_q)*(1-w_u) * ctx%state%dustem2_dustem(:, idx_q,   2*idx_u) + &
            w_q    *(1-w_u) * ctx%state%dustem2_dustem(:, idx_q+1, 2*idx_u) + &
            w_q    *w_u     * ctx%state%dustem2_dustem(:, idx_q+1, 2*(idx_u+1)) + &
            (1-w_q)*w_u     * ctx%state%dustem2_dustem(:, idx_q,   2*(idx_u+1))

        ! Combine components
        emission_spectrum = (1.0_sp - gamma_frac) * spec_u_min + gamma_frac * spec_u_max
        
        ! Safety floor
        emission_spectrum = max(emission_spectrum, SAFE_FLOOR)

    end subroutine interpolate_draine_li_dust_model


    !> @brief
    !> Calculates the final dust emission spectrum analytically.
    !>
    !> @details
    !> Replaces the iterative self-absorption loop with an exact analytical solution.
    !> Since the dust emission shape is fixed and energy is conserved, the final 
    !> spectrum is simply the attenuated dust shape (S_int * e^-tau) normalized 
    !> such that its total integrated luminosity equals the total stellar energy 
    !> absorbed.
    !>
    !> derivation:
    !> Final Spectrum = (S_int * e^-tau) * (L_absorbed_stellar / Integrate(S_int * e^-tau))
    !>
    subroutine calculate_dust_self_absorption(nu, shape_intrinsic, transmission_ism, &
                                              lum_absorbed_initial, spec_final)
        
        use fsps_integration, only: integrate_trapezoid_array
        
        real(SP), dimension(:), intent(in)  :: nu
        real(SP), dimension(:), intent(in)  :: shape_intrinsic
        real(SP), dimension(:), intent(in)  :: transmission_ism ! e^-tau
        real(SP), intent(in)                :: lum_absorbed_initial
        real(SP), dimension(:), intent(out) :: spec_final

        real(SP), dimension(size(nu)) :: profile_escaped
        real(SP) :: lum_escaped_profile, normalization_factor
        
        ! 1. Calculate the shape of the dust emission that actually escapes the galaxy.
        !    This is the intrinsic dust emission curve attenuated by the dust itself.
        profile_escaped = shape_intrinsic * transmission_ism
        
        ! 2. Integrate this profile to see how much luminosity it currently represents.
        lum_escaped_profile = integrate_trapezoid_array(nu, profile_escaped)
        
        ! 3. Normalize to ensure Energy Conservation.
        !    The total IR energy leaving the galaxy must equal the total UV/Optical 
        !    energy absorbed by the dust (L_absorbed_stellar).
        if (lum_escaped_profile > SAFE_FLOOR) then
            normalization_factor = lum_absorbed_initial / lum_escaped_profile
            spec_final = profile_escaped * normalization_factor
        else
            spec_final = 0.0_sp
        end if

    end subroutine calculate_dust_self_absorption

    !> Implementation of Cardelli, Clayton, & Mathis (1989) extinction curve.
    !> Includes the "hack" for smooth transitions used in the original FSPS.
    pure function get_ccm89_curve(wavelengths, r_v, uv_bump_strength) result(curve)
        real(SP), dimension(:), intent(in) :: wavelengths
        real(SP), intent(in) :: r_v, uv_bump_strength
        real(SP), dimension(size(wavelengths)) :: curve

        ! Array variables for the main calculation
        real(SP), dimension(size(wavelengths)) :: wavenumbers, wavenumber_term, poly_a, poly_b
        real(SP), dimension(size(wavelengths)) :: temp_curve, wavenumbers_clamped

        ! Scalar variables for calculating the Smoothing Hack (at x=3.3)
        real(SP) :: y_anchor, a_scalar, b_scalar
        real(SP) :: opt_val_at_break, nuv_val_at_break, continuity_correction

        ! 1. PRE-CALCULATION: The "Smoothing Hack" Constant
        ! -------------------------------------------------
        ! The Optical and NUV polynomials in the original paper do not meet perfectly 
        ! at x=3.3. FSPS calculates the gap at this specific point and adds a 
        ! correction term to the NUV region to force continuity.
        
        ! A. Calculate Optical Value at x=3.3
        y_anchor = CCM_X_NUV_MIN - CCM_OPT_Y_SHIFT
        
        a_scalar = CCM_OPT_A_COEFFS(0) + y_anchor*(CCM_OPT_A_COEFFS(1) + y_anchor*(CCM_OPT_A_COEFFS(2) + &
                   y_anchor*(CCM_OPT_A_COEFFS(3) + y_anchor*(CCM_OPT_A_COEFFS(4) + y_anchor*(CCM_OPT_A_COEFFS(5) + &
                   y_anchor*(CCM_OPT_A_COEFFS(6) + y_anchor*CCM_OPT_A_COEFFS(7)))))))
        
        b_scalar = CCM_OPT_B_COEFFS(0) + y_anchor*(CCM_OPT_B_COEFFS(1) + y_anchor*(CCM_OPT_B_COEFFS(2) + &
                   y_anchor*(CCM_OPT_B_COEFFS(3) + y_anchor*(CCM_OPT_B_COEFFS(4) + y_anchor*(CCM_OPT_B_COEFFS(5) + &
                   y_anchor*(CCM_OPT_B_COEFFS(6) + y_anchor*CCM_OPT_B_COEFFS(7)))))))

        opt_val_at_break = a_scalar + b_scalar / r_v

        ! B. Calculate NUV Value at x=3.3
        a_scalar = (CCM_NUV_A_BASE(1) + CCM_NUV_A_BASE(2)*CCM_X_NUV_MIN) + &
                   (CCM_NUV_A_BUMP(1) / ((CCM_X_NUV_MIN - CCM_NUV_A_BUMP(2))**2 + CCM_NUV_A_BUMP(3)) * uv_bump_strength)
                   
        b_scalar = (CCM_NUV_B_BASE(1) + CCM_NUV_B_BASE(2)*CCM_X_NUV_MIN) + &
                   (CCM_NUV_B_BUMP(1) / ((CCM_X_NUV_MIN - CCM_NUV_B_BUMP(2))**2 + CCM_NUV_B_BUMP(3)) * uv_bump_strength)
        
        nuv_val_at_break = a_scalar + b_scalar / r_v

        ! C. Define the offset
        continuity_correction = opt_val_at_break - nuv_val_at_break

        ! 2. ARRAY CALCULATION
        ! --------------------
        wavenumbers = get_wavenumber(wavelengths)
        
        poly_a = 0.0_sp
        poly_b = 0.0_sp
        temp_curve = 0.0_sp

        ! Region 1: Infrared (0.3 < x < 1.1)
        where (wavenumbers >= CCM_X_IR_MIN .and. wavenumbers < CCM_X_OPT_MIN)
            temp_curve = (CCM_IR_A_SCALE * wavenumbers**CCM_IR_EXP) + &
                         (CCM_IR_B_SCALE * wavenumbers**CCM_IR_EXP) / r_v
        end where

        ! Region 2: Optical / Near-IR (1.1 <= x < 3.3)
        where (wavenumbers >= CCM_X_OPT_MIN .and. wavenumbers < CCM_X_NUV_MIN)
            wavenumber_term = wavenumbers - CCM_OPT_Y_SHIFT
            
            poly_a = CCM_OPT_A_COEFFS(0) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(1) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(2) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(3) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(4) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(5) + &
                wavenumber_term * (CCM_OPT_A_COEFFS(6) + &
                wavenumber_term * CCM_OPT_A_COEFFS(7)))))))

            poly_b = CCM_OPT_B_COEFFS(0) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(1) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(2) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(3) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(4) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(5) + &
                wavenumber_term * (CCM_OPT_B_COEFFS(6) + &
                wavenumber_term * CCM_OPT_B_COEFFS(7)))))))
            
            temp_curve = poly_a + poly_b / r_v
        end where

        ! --- Regions 3 & 4: UV Base (3.3 <= x < 8.0) ---
        ! Both NUV and Mid-UV share the same base linear term and Drude profile.
        where (wavenumbers >= CCM_X_NUV_MIN .and. wavenumbers < CCM_X_FUV_MIN)
            ! Base Linear Component
            poly_a = CCM_NUV_A_BASE(1) + CCM_NUV_A_BASE(2)*wavenumbers
            poly_b = CCM_NUV_B_BASE(1) + CCM_NUV_B_BASE(2)*wavenumbers
            
            ! Add UV Bump (Drude Profile)
            poly_a = poly_a + CCM_NUV_A_BUMP(1) / &
                     ((wavenumbers - CCM_NUV_A_BUMP(2))**2 + CCM_NUV_A_BUMP(3)) * uv_bump_strength
                     
            poly_b = poly_b + CCM_NUV_B_BUMP(1) / &
                     ((wavenumbers - CCM_NUV_B_BUMP(2))**2 + CCM_NUV_B_BUMP(3)) * uv_bump_strength
        end where

        ! Specific NUV Modifier (3.3 <= x < 5.9): Apply Continuity Correction
        where (wavenumbers >= CCM_X_NUV_MIN .and. wavenumbers < CCM_X_MUV_MIN)
             temp_curve = poly_a + poly_b / r_v + continuity_correction * (CCM_X_NUV_MIN / wavenumbers)**6
        end where

        ! Specific Mid-UV Modifier (5.9 <= x < 8.0): Apply Curvature Polynomials
        where (wavenumbers >= CCM_X_MUV_MIN .and. wavenumbers < CCM_X_FUV_MIN)
            wavenumber_term = wavenumbers - CCM_X_MUV_MIN
            
            poly_a = poly_a + CCM_MUV_A_POLY(1) * wavenumber_term**2 + CCM_MUV_A_POLY(2) * wavenumber_term**3
            poly_b = poly_b + CCM_MUV_B_POLY(1) * wavenumber_term**2 + CCM_MUV_B_POLY(2) * wavenumber_term**3
                
            temp_curve = poly_a + poly_b / r_v
        end where

        ! Region 5: Far-UV (x >= 8.0)
        ! Clamps input x to 12.0 (lambda = 833 A) to prevent divergence
        where (wavenumbers >= CCM_X_FUV_MIN)
            wavenumbers_clamped = min(wavenumbers, CCM_X_CUTOFF)
            
            wavenumber_term = wavenumbers_clamped - CCM_X_FUV_MIN
            
            poly_a = CCM_FUV_A_COEFFS(0) + wavenumber_term*(CCM_FUV_A_COEFFS(1) + &
                wavenumber_term*(CCM_FUV_A_COEFFS(2) + wavenumber_term*CCM_FUV_A_COEFFS(3)))
            
            poly_b = CCM_FUV_B_COEFFS(0) + wavenumber_term*(CCM_FUV_B_COEFFS(1) + &
                wavenumber_term*(CCM_FUV_B_COEFFS(2) + wavenumber_term*CCM_FUV_B_COEFFS(3)))
            
            temp_curve = poly_a + poly_b / r_v
        end where

        curve = temp_curve
    end function get_ccm89_curve


    !> Implementation of Calzetti et al. (2000) starburst attenuation curve.
    pure function get_calzetti_curve(wavelengths) result(curve)
        real(SP), dimension(:), intent(in) :: wavelengths
        real(SP), dimension(size(wavelengths)) :: curve
        real(SP), dimension(size(wavelengths)) :: wavenumbers, extinction_k

        ! Convert to inverse microns (x = 1/lambda_um)
        wavenumbers = get_wavenumber(wavelengths)
        
        extinction_k = 0.0_sp
        
        ! Optical / NIR (0.63um < lambda <= 2.2um)
        where (wavelengths > CALZ_LAM_BREAK .and. wavelengths <= CALZ_LAM_IR_MAX)
            extinction_k = CALZ_SCALE * (CALZ_OPT_COEFFS(0) + CALZ_OPT_COEFFS(1) * wavenumbers) + CALZ_R_V
        end where
        
        ! UV / Optical (0.12um <= lambda <= 0.63um)
        where (wavelengths >= CALZ_LAM_UV_MIN .and. wavelengths <= CALZ_LAM_BREAK)
            ! Use nested multiplication (Horner's method) for clarity and efficiency
            extinction_k = CALZ_R_V + CALZ_SCALE * ( &
                CALZ_UV_COEFFS(0) + wavenumbers * ( &
                    CALZ_UV_COEFFS(1) + wavenumbers * ( &
                        CALZ_UV_COEFFS(2) + wavenumbers * CALZ_UV_COEFFS(3) &
                    ) &
                ) &
            )
        end where

        ! Result is A_lambda / A_V = k_lambda / R_V
        curve = extinction_k / CALZ_R_V
    end function get_calzetti_curve


    !> Implementation of Kriek & Conroy (2013): Calzetti + UV Bump + Tilt.
    pure function get_kriek_conroy_curve(wavelengths, tilt_index) result(curve)
        real(SP), dimension(:), intent(in) :: wavelengths
        real(SP), intent(in) :: tilt_index
        real(SP), dimension(size(wavelengths)) :: curve
        
        real(SP), dimension(size(wavelengths)) :: base_calzetti, drude_profile
        real(SP) :: bump_amplitude ! E_b in paper

        ! 1. Base Calzetti (normalized to E(B-V), i.e., k_lambda scale)
        ! Note: Our helper `get_calzetti_curve` returns A_lambda/A_V.
        ! We must multiply by R_V to get back to k_lambda.
        base_calzetti = get_calzetti_curve(wavelengths) * KC13_R_V_BASE

        ! 2. UV Bump (Drude Profile)
        ! Kriek & Conroy (2013) Eq 3: E_b = 0.85 - 1.9 * delta
        ! (Note: The `uv_bump_strength` parameter is technically not in the standard 
        !  KC13 definition, but FSPS likely passes it to allow modulating the bump 
        !  further. We apply it here to match original logic if intended, 
        !  though the original snippet didn't use `uv_bump_strength` in the Drude calculation.
        !  Based on your provided snippet, `uv_bump_strength` was unused! 
        !  I will stick strictly to your snippet's logic which calculated amplitude solely from tilt.)
        
        bump_amplitude = KC13_AMPL_INTERCEPT - KC13_AMPL_SLOPE * tilt_index
        
        drude_profile = bump_amplitude * (wavelengths * KC13_BUMP_WIDTH)**2 / &
                        ( (wavelengths**2 - UV_BUMP_CENTER**2)**2 + (wavelengths * KC13_BUMP_WIDTH)**2 )

        ! 3. Combine with Tilt
        ! A_lambda = (k_calz + D_bump) / R_V * (lambda / 5500)^delta
        curve = (base_calzetti + drude_profile) / KC13_R_V_BASE * &
                (wavelengths / V_BAND_ANGSTROMS)**tilt_index

    end function get_kriek_conroy_curve


!> Implementation of Reddy et al. (2015) MOSDEF curve.
    pure function get_reddy_curve(wavelengths) result(curve)
        real(SP), dimension(:), intent(in) :: wavelengths
        real(SP), dimension(size(wavelengths)) :: curve
        real(SP), dimension(size(wavelengths)) :: wavenumbers, extinction_k, wavenumbers_clamped
        
        wavenumbers = get_wavenumber(wavelengths)
        extinction_k = 0.0_sp
        
        ! 1. UV Range (Lambda < 6000 A)
        !    For Lambda < 1500, we clamp the wavenumber to the value at 1500,
        !    effectively extrapolating the curve as a constant value blueward.
        where (wavelengths < REDDY_LAM_BREAK)
            wavenumbers_clamped = min(wavenumbers, REDDY_X_UV_MAX)
            
            extinction_k = REDDY_UV_COEFFS(0) + &
                           wavenumbers_clamped * (REDDY_UV_COEFFS(1) + &
                           wavenumbers_clamped * (REDDY_UV_COEFFS(2) + &
                           wavenumbers_clamped * REDDY_UV_COEFFS(3))) + REDDY_OFFSET
        end where

        ! 2. Optical/NIR Range (0.60um <= lambda < 2.85um)
        where (wavelengths >= REDDY_LAM_BREAK .and. wavelengths < REDDY_LAM_IR_MAX)
             extinction_k = REDDY_OPT_COEFFS(0) + &
                            wavenumbers * (REDDY_OPT_COEFFS(1) + &
                            wavenumbers * (REDDY_OPT_COEFFS(2) + &
                            wavenumbers * REDDY_OPT_COEFFS(3))) + &
                            REDDY_OFFSET + REDDY_OPT_CORRECTION
        end where

        ! Convert k_lambda to A_lambda / A_V assuming R_V = 2.505
        curve = extinction_k / REDDY_R_V
    end function get_reddy_curve

    !> Converts wavelength (Angstroms) to wavenumber (inverse microns).
    !> Used frequently for dust curve parameterizations (CCM89, Calzetti, etc.).
    pure function get_wavenumber(wavelengths) result(wavenumbers)
        real(SP), dimension(:), intent(in) :: wavelengths
        real(SP), dimension(size(wavelengths)) :: wavenumbers
        
        ! x = 1 / lambda_microns = 10000 / lambda_angstroms
        wavenumbers = 1.0e4_sp / wavelengths
    end function get_wavenumber

    !> Computes the circumstellar optical depth (tau_1um) from physical parameters.
    !> See Villaume et al. (2015).
    pure function compute_circumstellar_optical_depth(ctx, c_rich_flag, m_act, log_l, log_g, log_mdot_iso) result(tau)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)  :: c_rich_flag
        real(SP), intent(in) :: m_act, log_l, log_g, log_mdot_iso
        real(SP) :: tau
        
        real(SP) :: radius_solar, period_days, velocity_exp, mdot_sol_yr
        real(SP) :: inner_radius_cm, dust_gas_ratio, kappa_eff

        ! 1. Determine Constants based on Chemistry
        ! -----------------------------------------
        if (c_rich_flag == 1) then
            kappa_eff = AGB_KAPPA_C_RICH
        else
            kappa_eff = AGB_KAPPA_O_RICH
        end if

        ! 2. Stellar Parameters
        ! ---------------------
        ! Radius (R_sun) = sqrt(GM / g) / R_sun_cm
        radius_solar = sqrt(m_act * M_SOL * G_NEWTON / (10.0_sp**log_g)) / R_SOL

        ! Fundamental Pulsation Period (Days) - Villaume et al. (2015) relation
        period_days = 10.0_sp**(AGB_PER_INTERCEPT + &
                                AGB_PER_SLOPE_R * log10(radius_solar) + &
                                AGB_PER_SLOPE_M * log10(m_act))

        ! Expansion Velocity (km/s)
        velocity_exp = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        velocity_exp = max(min(velocity_exp, AGB_VEXP_MAX), AGB_VEXP_MIN) 

        ! 3. Mass Loss Rate (M_sun/yr)
        ! ----------------------------
        if (ctx%use_isoc_mdot_val == 1) then
            ! Use Isochrone value (MIST only), capped
            mdot_sol_yr = min(10.0_sp**log_mdot_iso, VW93_MDOT_LIMIT_ISO)
        else
            ! Vassiliadis & Wood (1993) Prescription
            if (period_days < VW93_PER_THRESH) then
                if (m_act < VW93_MASS_THRESH) then
                    mdot_sol_yr = 10.0_sp**(VW93_BASE_INTERCEPT + VW93_BASE_SLOPE * period_days)
                else
                    mdot_sol_yr = 10.0_sp**(VW93_BASE_INTERCEPT + VW93_HIGH_SLOPE * &
                                           (period_days - 100.0_sp * (m_act - VW93_MASS_THRESH)))
                end if
            else
                ! Superwind phase
                mdot_sol_yr = (10.0_sp**log_l) / velocity_exp * VW93_SUPERWIND_NORM * YEAR_TO_SECOND / C_LIGHT
            end if
        end if

        ! 4. Shell Geometry
        ! -----------------
        ! Inner Radius (cm) scaling with Luminosity
        if (c_rich_flag == 1) then
            inner_radius_cm = AGB_RIN_FACTOR_C * (10.0_sp**log_l)**0.5_sp
        else
            inner_radius_cm = AGB_RIN_FACTOR_O * (10.0_sp**log_l)**0.5_sp
        end if

        ! Dust-to-Gas Ratio (delta) scaling with Velocity and Luminosity
        if (c_rich_flag == 1) then
            dust_gas_ratio = AGB_DELTA_C_RICH
        else
            dust_gas_ratio = AGB_DELTA_O_RICH
        end if
        
        dust_gas_ratio = dust_gas_ratio * (velocity_exp**2 / AGB_DTG_VEL_NORM) * &
                         ((10.0_sp**log_l / AGB_DTG_LUM_NORM)**AGB_DTG_LUM_EXP)

        ! 5. Final Optical Depth Calculation
        ! ----------------------------------
        ! tau = kappa * delta * Mdot / (4 * pi * R_in * v_exp)
        ! Note: v_exp is converted from km/s to cm/s (1E5 factor)
        
        tau = kappa_eff * dust_gas_ratio * (mdot_sol_yr * M_SOL / YEAR_TO_SECOND) / &
              inner_radius_cm / (4.0_sp * PI) / (velocity_exp * 1.0e5_sp)

    end function compute_circumstellar_optical_depth

end module fsps_dust