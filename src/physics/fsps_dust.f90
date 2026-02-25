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

    use fsps_precision, only: WP
    use fsps_constants, only: C_LIGHT, NEMLINE, GRAVITY_L_M_T_COEFF, &
                              M_SOL, G_NEWTON, R_SOL, YEAR_TO_SECOND, &
                              SAFE_FLOOR, PI
    use fsps_types, only: params
    use fsps_context_types, only: fsps_context_t
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: interpolate_linear, find_interval
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none
    private

    ! Public Interface
    public :: apply_dust_attenuation_and_emission
    public :: compute_attenuation_curve_point
    public :: apply_agb_dust_screen
    public :: apply_agn_dust_emission
    public :: interpolate_draine_li_dust_model
    public :: calculate_dust_self_absorption
    public :: compute_circumstellar_optical_depth

    ! ------------------------------------------------------------------------
    ! CONSTANTS: General Dust Parameters
    ! ------------------------------------------------------------------------
    real(WP), parameter :: V_BAND_ANGSTROMS  = 5500.0_wp
    real(WP), parameter :: UV_BUMP_CENTER    = 2175.0_wp
    real(WP), parameter :: CALZETTI_BREAK    = 6300.0_wp

    ! ------------------------------------------------------------------------
    ! CONSTANTS: Cardelli, Clayton, & Mathis (1989) Extinction Curve Parameters
    ! ------------------------------------------------------------------------

    ! Region Boundaries (in inverse microns, x = 1/lambda)
    real(WP), parameter :: CCM_X_IR_MIN  = 0.3_wp
    real(WP), parameter :: CCM_X_OPT_MIN = 1.1_wp
    real(WP), parameter :: CCM_X_NUV_MIN = 3.3_wp
    real(WP), parameter :: CCM_X_MUV_MIN = 5.9_wp
    real(WP), parameter :: CCM_X_FUV_MIN = 8.0_wp

    ! FUV Cutoff (Clamp values for x > 12.0, i.e., lambda < 833 A)
    real(WP), parameter :: CCM_X_CUTOFF = 12.0_wp
    
    ! Infrared Parameters (0.3 <= x < 1.1)
    real(WP), parameter :: CCM_IR_A_SCALE = 0.574_wp
    real(WP), parameter :: CCM_IR_B_SCALE = -0.527_wp
    real(WP), parameter :: CCM_IR_EXP     = 1.61_wp

    ! Optical Parameters (1.1 <= x < 3.3)
    ! Polynomial coefficients for y = x - 1.82 (Powers 0 through 7)
    real(WP), parameter :: CCM_OPT_Y_SHIFT = 1.82_wp
    real(WP), parameter :: CCM_OPT_A_COEFFS(0:7) = [ &
        1.0_wp,      0.17699_wp, -0.50447_wp, -0.02427_wp, &
        0.72085_wp,  0.01979_wp, -0.77530_wp,  0.32999_wp ]
    real(WP), parameter :: CCM_OPT_B_COEFFS(0:7) = [ &
        0.0_wp,      1.41338_wp,  2.28305_wp,  1.07233_wp, &
       -5.38434_wp, -0.62251_wp,  5.30260_wp, -2.09002_wp ]

    ! Near-UV Parameters (3.3 <= x < 5.9)
    ! Base Linear Terms: C1 + C2*x
    real(WP), parameter :: CCM_NUV_A_BASE(2) = [ 1.752_wp, -0.316_wp] 
    real(WP), parameter :: CCM_NUV_B_BASE(2) = [-3.09_wp,   1.825_wp]
    ! Drude Bump Terms: Scale / ((x-Pos)**2 + Width)
    real(WP), parameter :: CCM_NUV_A_BUMP(3) = [-0.104_wp, 4.67_wp, 0.341_wp]
    real(WP), parameter :: CCM_NUV_B_BUMP(3) = [ 1.206_wp, 4.62_wp, 0.263_wp]

    ! Mid-UV Parameters (5.9 <= x < 8.0) - Additions to NUV
    ! F(x) = C2*(x-5.9)**2 + C3*(x-5.9)**3
    real(WP), parameter :: CCM_MUV_A_POLY(2) = [-0.04473_wp, -0.009779_wp]
    real(WP), parameter :: CCM_MUV_B_POLY(2) = [ 0.2130_wp,   0.1207_wp]

    ! Far-UV Parameters (x >= 8.0)
    ! Polynomials in (x - 8.0)
    real(WP), parameter :: CCM_FUV_A_COEFFS(0:3) = [-1.073_wp, -0.628_wp, 0.137_wp, -0.070_wp]
    real(WP), parameter :: CCM_FUV_B_COEFFS(0:3) = [ 13.67_wp,  4.257_wp, -0.42_wp,  0.374_wp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Calzetti et al. (2000) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(WP), parameter :: CALZ_LAM_UV_MIN = 1200.0_wp  ! 0.12 microns
    real(WP), parameter :: CALZ_LAM_BREAK  = 6300.0_wp  ! 0.63 microns
    real(WP), parameter :: CALZ_LAM_IR_MAX = 22000.0_wp ! 2.20 microns

    ! General Parameters
    real(WP), parameter :: CALZ_R_V   = 4.05_wp
    real(WP), parameter :: CALZ_SCALE = 2.659_wp      ! Scaling factor k'

    ! UV/Optical Polynomial Coefficients (0.12 <= lambda < 0.63 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(WP), parameter :: CALZ_UV_COEFFS(0:3) = [-2.156_wp, 1.509_wp, -0.198_wp, 0.011_wp]

    ! Optical/NIR Linear Coefficients (0.63 <= lambda <= 2.2 um)
    ! Linear in (1/lambda_microns): c0 + c1*x
    real(WP), parameter :: CALZ_OPT_COEFFS(0:1) = [-1.857_wp, 1.040_wp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Kriek & Conroy (2013) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------
    
    ! Drude Profile Parameters
    real(WP), parameter :: KC13_BUMP_WIDTH = 350.0_wp  ! Delta lambda (Angstroms)
    
    ! Bump Amplitude Relationship: E_b = 0.85 - 1.9 * delta
    real(WP), parameter :: KC13_AMPL_INTERCEPT = 0.85_wp
    real(WP), parameter :: KC13_AMPL_SLOPE     = 1.9_wp

    ! The base model is tied to Calzetti's specific R_V
    real(WP), parameter :: KC13_R_V_BASE = 4.05_wp

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Reddy et al. (2015) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(WP), parameter :: REDDY_LAM_UV_MIN = 1500.0_wp
    real(WP), parameter :: REDDY_LAM_BREAK  = 6000.0_wp
    real(WP), parameter :: REDDY_LAM_IR_MAX = 28500.0_wp

    ! Region Boundary (Inverse Microns)
    real(WP), parameter :: REDDY_X_UV_MAX   = 1.0e4_wp / REDDY_LAM_UV_MIN

    ! General Parameters
    real(WP), parameter :: REDDY_R_V        = 2.505_wp
    real(WP), parameter :: REDDY_OFFSET     = 2.505_wp ! Base offset added to both curves

    ! UV Polynomial Coefficients (0.15 <= lambda < 0.60 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(WP), parameter :: REDDY_UV_COEFFS(0:3) = [-5.726_wp, 4.004_wp, -0.525_wp, 0.029_wp]
    
    ! Blueward Extrapolation Value (< 0.15 um)
    real(WP), parameter :: REDDY_UV_EXTRAP_VAL = 10.36_wp

    ! Optical/NIR Polynomial Coefficients (0.60 <= lambda < 2.85 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(WP), parameter :: REDDY_OPT_COEFFS(0:3) = [-2.672_wp, -0.010_wp, 1.532_wp, -0.412_wp]
    
    ! Continuity Correction for Optical Range
    real(WP), parameter :: REDDY_OPT_CORRECTION = -0.036221981_wp

    ! --------------------------------------------------------------------------
    ! CONSTANTS: AGB Circumstellar Dust Parameters
    ! --------------------------------------------------------------------------
    
    ! Dust-to-Gas Ratios (delta)
    real(WP), parameter :: AGB_DELTA_C_RICH = 0.0025_wp
    real(WP), parameter :: AGB_DELTA_O_RICH = 0.01_wp

    ! Extinction Coefficients (kappa)
    real(WP), parameter :: AGB_KAPPA_C_RICH = 3200.0_wp ! AmC + SiC
    real(WP), parameter :: AGB_KAPPA_O_RICH = 3000.0_wp ! Silicates

    ! Inner Radius Factors (cm * L^-0.5)
    real(WP), parameter :: AGB_RIN_FACTOR_C = 1.92E12_wp ! Td = 1100 K
    real(WP), parameter :: AGB_RIN_FACTOR_O = 4.74E12_wp ! Td = 700 K

    ! Villaume et al. (2015) Period Relation Coefficients
    ! logP = A + B*logR + C*logM
    real(WP), parameter :: AGB_PER_INTERCEPT = -2.07_wp
    real(WP), parameter :: AGB_PER_SLOPE_R   = 1.94_wp
    real(WP), parameter :: AGB_PER_SLOPE_M   = -0.9_wp

    ! Expansion Velocity Parameters
    real(WP), parameter :: AGB_VEXP_INTERCEPT = -13.5_wp
    real(WP), parameter :: AGB_VEXP_SLOPE     = 0.056_wp
    real(WP), parameter :: AGB_VEXP_MIN       = 3.0_wp  ! km/s
    real(WP), parameter :: AGB_VEXP_MAX       = 15.0_wp ! km/s

    ! Vassiliadis & Wood (1993) Mass Loss Parameters
    real(WP), parameter :: VW93_MDOT_LIMIT_ISO = 1.0e-4_wp
    real(WP), parameter :: VW93_PER_THRESH     = 500.0_wp ! Days
    real(WP), parameter :: VW93_MASS_THRESH    = 2.5_wp   ! Solar Masses
    real(WP), parameter :: VW93_BASE_INTERCEPT = -11.4_wp
    real(WP), parameter :: VW93_BASE_SLOPE     = 0.0123_wp
    real(WP), parameter :: VW93_HIGH_SLOPE     = 0.0125_wp
    real(WP), parameter :: VW93_SUPERWIND_NORM = 1.93e3_wp

    ! Dust-to-Gas Ratio Scaling
    real(WP), parameter :: AGB_DTG_VEL_NORM    = 225.0_wp
    real(WP), parameter :: AGB_DTG_LUM_NORM    = 1.0e4_wp
    real(WP), parameter :: AGB_DTG_LUM_EXP     = -0.6_wp

    ! Smoothing Constants
    real(WP), parameter :: AGB_SMOOTH_SIGMA      = 1.0e4_wp
    real(WP), parameter :: AGB_SMOOTH_LAMBDA_MIN = 3.0e4_wp
    real(WP), parameter :: AGB_SMOOTH_LAMBDA_MAX = 1.0e8_wp

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
        real(WP), dimension(:), intent(in)     :: spec_young, spec_old
        real(WP), dimension(:), intent(in)     :: neb_flux_young, neb_flux_old
        real(WP), dimension(:), intent(out)    :: spec_total_out
        real(WP), intent(out)                  :: dust_mass
        real(WP), dimension(:), intent(out)    :: neb_flux_out

        ! Local Variables
        ! We use automatic arrays. If sizes are large, we might need create.
        ! But here we are inside a kernel-like routine (called from csp).
        ! We assume data is present.
        ! Note: Automatic arrays on device stack might be limited.
        ! However, these are nspec sized.
        ! We can use data create if needed, or rely on compiler.
        ! For resident device, explicit data clauses are safer.
        real(WP), dimension(size(spec_young)) :: transmission_diffuse
        real(WP), dimension(size(spec_young)) :: frequencies
        
        real(WP) :: lum_bol_intrinsic, lum_bol_attenuated, lum_absorbed_total
        real(WP), dimension(size(spec_young)) :: dust_emission_shape, dust_emission_final
        real(WP) :: emission_norm_factor
        
        integer :: nspec, i
        real(WP) :: y1, y2
        real(WP) :: curve, trans_birth, trans_old, spec_sum, trans_diffuse_neb, neb_birth
        real(WP) :: frac_obrun, one_minus_obrun, frac_nodust, one_minus_nodust
        real(WP) :: dust1, dust1_index, dust2, dust3
        logical  :: dust_type_is3

        !$acc enter data create(transmission_diffuse, frequencies)
        !$acc enter data create(dust_emission_shape, dust_emission_final)

        nspec = size(spec_young)

        ! 0. Input Validation
        ! -------------------
        if (settings%uvb < 0.0_wp) return ! Should trigger error handling upstream
        if (settings%wgp1 < 1 .or. settings%wgp2 < 1) return 

        frac_obrun      = settings%frac_obrun
        one_minus_obrun = 1.0_wp - frac_obrun
        frac_nodust     = settings%frac_nodust
        one_minus_nodust = 1.0_wp - frac_nodust
        dust1           = settings%dust1
        dust1_index     = settings%dust1_index
        dust2           = settings%dust2
        dust3           = settings%dust3
        dust_type_is3   = (ctx%dust_type_val == 3)

        ! 1. Calculate Attenuation Curves & Transmissivities
        ! --------------------------------------------------
        ! We parallelize the array operations. compute_attenuation_curve needs to be !acc routine seq/vector.
        
        ! A. Diffuse ISM (affects all stars)
        ! B. Birth Clouds (affects young stars only)
        ! 2. Apply Attenuation to Stellar Spectra
        
        !$acc parallel loop present(ctx, spec_young, spec_old, spec_total_out, transmission_diffuse)
        do i = 1, nspec
            ! A. Diffuse Curve (Inline call or routine seq)
            curve = compute_attenuation_curve_point(ctx%state%spec_lambda(i), i, &
                                                    ctx%dust_type_val, settings, ctx)
            
            if (dust_type_is3) then
                transmission_diffuse(i) = exp(-curve)
            else
                transmission_diffuse(i) = exp(-dust2 * curve)
            end if

            ! B. Birth Clouds
            if (dust1 <= SAFE_FLOOR .or. one_minus_obrun <= SAFE_FLOOR) then
                trans_birth = 1.0_wp
            else
                trans_birth = exp(-dust1 * (ctx%state%spec_lambda(i) / V_BAND_ANGSTROMS)**dust1_index)
            end if

            if (abs(dust3) <= SAFE_FLOOR) then
                trans_old = 1.0_wp
            else
                trans_old = exp(-dust3 * curve)
            end if
                                       
            ! 2. Apply Attenuation
            spec_sum = (spec_young(i) * trans_birth * one_minus_obrun + spec_young(i) * frac_obrun) + &
                       (spec_old(i) * trans_old)

            ! Final diffuse screen
            if (one_minus_nodust <= SAFE_FLOOR) then
                spec_total_out(i) = spec_sum
            else
                spec_total_out(i) = spec_sum * (transmission_diffuse(i) * one_minus_nodust + frac_nodust)
            end if
        end do


        ! 3. Apply Attenuation to Nebular Lines
        ! -------------------------------------
        ! Note: We must interpolate the diffuse transmission to the line wavelengths
        ! interpolate_linear is now !acc routine seq
        
           !$acc parallel loop present(ctx, neb_flux_young, neb_flux_old, neb_flux_out, transmission_diffuse)
        do i = 1, size(neb_flux_young)
             if (one_minus_nodust <= SAFE_FLOOR) then
                 trans_diffuse_neb = 1.0_wp
             else
                 trans_diffuse_neb = interpolate_linear(ctx%state%spec_lambda, &
                                                        transmission_diffuse, &
                                                        ctx%state%nebem_line_pos(i))
             end if
                                                              
               if (trans_diffuse_neb /= trans_diffuse_neb) trans_diffuse_neb = 1.0_wp

             if (dust1 <= SAFE_FLOOR .or. one_minus_obrun <= SAFE_FLOOR) then
                 neb_birth = 1.0_wp
             else
                 neb_birth = exp(-dust1 * (ctx%state%nebem_line_pos(i) / V_BAND_ANGSTROMS)**dust1_index)
             end if
             
             neb_flux_out(i) = (neb_flux_young(i) * &
                        neb_birth * &
                        one_minus_obrun + &
                        neb_flux_young(i) * frac_obrun + &
                        neb_flux_old(i))
                        
               if (one_minus_nodust > SAFE_FLOOR) then
                  neb_flux_out(i) = neb_flux_out(i) * (trans_diffuse_neb * one_minus_nodust + frac_nodust)
               end if
        end do


        ! 4. Add Dust Emission (Energy Balance)
        ! -------------------------------------
        if (ctx%add_dust_emission_val == 1 .and. &
            (settings%dust1 > SAFE_FLOOR .or. settings%dust2 > SAFE_FLOOR)) then
            
            !$acc parallel loop present(ctx, frequencies)
            do i = 1, nspec
                frequencies(i) = C_LIGHT / ctx%state%spec_lambda(i)
            end do

            ! Calculate Bolometric Luminosities (L_bol)
            ! Assumes integrate_trapezoid_array is modified to take raw arrays?
            ! No, it takes assumed-shape. This is hard on device if we want to avoid array creation.
            ! But we can compute array expressions? spec_young + spec_old.
            ! This creates temp array.
            ! We should write a specialized kernel or loop for integration.
            
            ! Intrinsic (Pre-Dust)
            lum_bol_intrinsic = 0.0_wp
            !$acc parallel loop reduction(+:lum_bol_intrinsic) present(frequencies, spec_young, spec_old)
            do i = 1, nspec-1
                y1 = spec_young(i) + spec_old(i)
                y2 = spec_young(i+1) + spec_old(i+1)
                lum_bol_intrinsic = lum_bol_intrinsic + 0.5_wp * abs(frequencies(i+1) - frequencies(i)) * (y1 + y2)
            end do
            
            if (ctx%nebemlineinspec_val == 0) then
                 !$acc kernels present(neb_flux_young, neb_flux_old)
                 lum_bol_intrinsic = lum_bol_intrinsic + sum(neb_flux_young) + sum(neb_flux_old)
                 !$acc end kernels
            end if

            ! Attenuated (Post-Dust)
            lum_bol_attenuated = 0.0_wp
            !$acc parallel loop reduction(+:lum_bol_attenuated) present(frequencies, spec_total_out)
            do i = 1, nspec-1
                y1 = spec_total_out(i)
                y2 = spec_total_out(i+1)
                lum_bol_attenuated = lum_bol_attenuated + 0.5_wp * abs(frequencies(i+1) - frequencies(i)) * (y1 + y2)
            end do
            
            if (ctx%nebemlineinspec_val == 0) then
                 !$acc kernels present(neb_flux_out)
                 lum_bol_attenuated = lum_bol_attenuated + sum(neb_flux_out)
                 !$acc end kernels
            end if
            
            ! Total Energy Absorbed by Dust
            lum_absorbed_total = lum_bol_intrinsic - lum_bol_attenuated

            ! Get Dust Emission Template (Draine & Li 2007)
            ! ---------------------------------------------
            call interpolate_draine_li_dust_model(ctx, settings, dust_emission_shape)
            
            ! Normalize template area
            emission_norm_factor = 0.0_wp
            !$acc parallel loop reduction(+:emission_norm_factor) present(frequencies, dust_emission_shape)
            do i = 1, nspec-1
                y1 = dust_emission_shape(i)
                y2 = dust_emission_shape(i+1)
                emission_norm_factor = emission_norm_factor + 0.5_wp * abs(frequencies(i+1) - frequencies(i)) * (y1 + y2)
            end do
            
            if (emission_norm_factor <= SAFE_FLOOR) then
                dust_mass = SAFE_FLOOR
                ! Cleanup
                !$acc exit data delete(transmission_diffuse, frequencies)
                !$acc exit data delete(dust_emission_shape, dust_emission_final)
                return
            end if

            ! Calculate Self-Absorption & Final Emission
            ! ------------------------------------------
            call calculate_dust_self_absorption(frequencies, dust_emission_shape, transmission_diffuse, &
                                                lum_absorbed_total, dust_emission_final)

            ! Add to total spectrum
            !$acc parallel loop present(spec_total_out, dust_emission_final)
            do i = 1, nspec
                spec_total_out(i) = spec_total_out(i) + dust_emission_final(i)
            end do

            ! Estimate Dust Mass (Factor from Draine & Li MW3.1 model)
            ! 3.21e-3 converts Luminosity/Norm to Mass (Solar Units) roughly
            dust_mass = 3.21e-3_wp / (4.0_wp * PI) * (lum_absorbed_total / emission_norm_factor)

        else
            dust_mass = SAFE_FLOOR
        end if
        
        ! Cleanup
        !$acc exit data delete(transmission_diffuse, frequencies)
        !$acc exit data delete(dust_emission_shape, dust_emission_final)

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
    !> @brief Computes attenuation at a single point (Device Compatible).
    !> We refactor the array function into an elemental/scalar one for the parallel loop.
    !> @return attenuation_val
    pure function compute_attenuation_curve_point(wavelength, idx, dust_type_id, settings, ctx) result(attenuation_val)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        integer, intent(in)                :: idx, dust_type_id
        type(params), intent(in)           :: settings
        type(fsps_context_t), intent(in)   :: ctx
        real(WP) :: attenuation_val
        
        attenuation_val = 0.0_wp

        select case (dust_type_id)
        case (0)
            attenuation_val = (wavelength / V_BAND_ANGSTROMS)**settings%dust_index
        case (1)
            attenuation_val = get_ccm89_curve_point(wavelength, settings%mwr, settings%uvb)
        case (2)
            attenuation_val = get_calzetti_curve_point(wavelength)
        case (3)
            ! Table lookup using index
            attenuation_val = ctx%state%wgdust(idx, settings%wgp1, settings%wgp2, settings%wgp3)
        case (4)
            attenuation_val = get_kriek_conroy_curve_point(wavelength, settings%dust_index)
        case (5)
            ! Table lookup using index
            attenuation_val = ctx%state%g03smcextn(idx)
        case (6)
            attenuation_val = get_reddy_curve_point(wavelength)
        case default
            attenuation_val = ieee_value(1.0_wp, ieee_quiet_nan)
        end select
    end function compute_attenuation_curve_point

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
        !$acc routine seq
        
        type(fsps_context_t), intent(inout)   :: ctx
        real(WP), intent(in)                  :: weight
        real(WP), dimension(:), intent(inout) :: spectrum
        real(WP), intent(in)                  :: mass_act, log_t, log_l, log_g
        real(WP), intent(in)                  :: c_o_ratio, log_mdot
        
        ! Local variables
        integer  :: c_rich_flag ! 0 = O-rich, 1 = C-rich
        integer  :: idx_teff, idx_tau
        integer  :: n_teff_grid, n_tau_grid
        real(WP) :: tau_1um, log_g_local
        real(WP) :: teff_linear, log_tau
        real(WP) :: w_teff, w_tau ! Interpolation weights
        real(WP) :: c00, c10, c11, c01, transfer_val
        integer  :: i

        ! 1. Determine Chemistry (C-rich vs O-rich)
        ! -----------------------------------------
        if (c_o_ratio > 1.0_wp) then
            c_rich_flag = 1
        else
            c_rich_flag = 0
        end if

        ! 2. Handle Log(g) for BaSTI Isochrones
        ! -------------------------------------
        ! BaSTI does not always tabulate log(g), so we compute it physically.
        if (ctx%state%isoc_type == 'bsti') then
            log_g_local = log10(GRAVITY_L_M_T_COEFF * mass_act / (10.0_wp**log_l)) + 4.0_wp * log_t
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
        teff_linear = 10.0_wp**log_t

        idx_teff = find_interval(ctx%state%teff_dagb(c_rich_flag + 1, :), teff_linear)
        idx_teff = max(1, min(idx_teff, n_teff_grid - 1))
        
        w_teff = (teff_linear - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff)) / &
                 (ctx%state%teff_dagb(c_rich_flag + 1, idx_teff + 1) - ctx%state%teff_dagb(c_rich_flag + 1, idx_teff))
        w_teff = max(-1.0_wp, min(w_teff, 1.0_wp)) ! Allow slight extrapolation

        ! B. Interpolate in Tau
        log_tau = log10(tau_1um)
        idx_tau = find_interval(ctx%state%tau1_dagb(c_rich_flag + 1, :), log_tau)
        idx_tau = max(1, min(idx_tau, n_tau_grid - 1))

        w_tau = (log_tau - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau)) / &
                (ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau + 1) - ctx%state%tau1_dagb(c_rich_flag + 1, idx_tau))
        w_tau = max(-1.0_wp, min(w_tau, 1.0_wp))


        ! C. Bilinear Interpolation + In-place application
        ! f(x, y) ~ (1-x)(1-y)F00 + x(1-y)F10 + (1-x)yF01 + xyF11
        c00 = (1.0_wp - w_teff) * (1.0_wp - w_tau)
        c10 = w_teff            * (1.0_wp - w_tau)
        c11 = w_teff            * w_tau
        c01 = (1.0_wp - w_teff) * w_tau

        do i = 1, size(spectrum)
            transfer_val = c00 * ctx%state%flux_dagb(i, c_rich_flag + 1, idx_teff,     idx_tau)     + &
                           c10 * ctx%state%flux_dagb(i, c_rich_flag + 1, idx_teff + 1, idx_tau)     + &
                           c11 * ctx%state%flux_dagb(i, c_rich_flag + 1, idx_teff + 1, idx_tau + 1) + &
                           c01 * ctx%state%flux_dagb(i, c_rich_flag + 1, idx_teff,     idx_tau + 1)
            spectrum(i) = spectrum(i) * transfer_val
        end do

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
        real(WP), dimension(:), intent(in)     :: wavelengths
        real(WP), intent(in)                   :: log_lbol_stellar
        real(WP), dimension(:), intent(inout)  :: spectrum_inout

        ! Local variables
        ! Use max size for stack alloc if needed, or assume kernel mode
        real(WP) :: agn_template_interpolated(size(wavelengths))
        real(WP) :: galaxy_attenuation_curve(size(wavelengths))
        real(WP) :: tau_agn_param, interpolation_weight
        real(WP) :: luminosity_agn_bolometric
        integer  :: idx_tau_grid, n_agn_grid
        integer :: i

        ! 0. Early exit if no AGN contribution is specified
        if (settings%fagn <= tiny(0.0_wp)) return
        
        ! 1. Interpolate AGN Template based on Torus Optical Depth (agn_tau)
        tau_agn_param = settings%agn_tau
        n_agn_grid    = size(ctx%state%agndust_tau)

        ! Use binary search to find the interval
        idx_tau_grid = find_interval(ctx%state%agndust_tau, tau_agn_param)
        idx_tau_grid = max(1, min(idx_tau_grid, n_agn_grid - 1))

        ! Calculate linear interpolation weight
        interpolation_weight = (tau_agn_param - ctx%state%agndust_tau(idx_tau_grid)) / &
                               (ctx%state%agndust_tau(idx_tau_grid + 1) - ctx%state%agndust_tau(idx_tau_grid))
        interpolation_weight = max(0.0_wp, min(interpolation_weight, 1.0_wp))

        ! Interpolate the template AND calculate attenuation
        ! Combined loop for performance
        luminosity_agn_bolometric = (10.0_wp**log_lbol_stellar) * settings%fagn
        
        !$acc parallel loop present(ctx, wavelengths, spectrum_inout) private(agn_template_interpolated, galaxy_attenuation_curve)
        do i = 1, size(wavelengths)
            ! Interpolate Template
            agn_template_interpolated(i) = (1.0_wp - interpolation_weight) * ctx%state%agndust_spec(i, idx_tau_grid) + &
                                           interpolation_weight * ctx%state%agndust_spec(i, idx_tau_grid + 1)
            
            ! Calculate Attenuation
            galaxy_attenuation_curve(i) = compute_attenuation_curve_point(wavelengths(i), i, ctx%dust_type_val, settings, ctx)
            
            ! Apply Attenuation
            if (ctx%dust_type_val == 3) then
                agn_template_interpolated(i) = agn_template_interpolated(i) * exp(-galaxy_attenuation_curve(i))
            else
                agn_template_interpolated(i) = agn_template_interpolated(i) * exp(-settings%dust2 * galaxy_attenuation_curve(i))
            end if
            
            ! Add to Spectrum
            spectrum_inout(i) = spectrum_inout(i) + (luminosity_agn_bolometric * agn_template_interpolated(i))
        end do

    end subroutine apply_agn_dust_emission

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> Interpolates the Draine & Li (2007) dust emission grid based on Q_PAH, U_min, and Gamma.
    subroutine interpolate_draine_li_dust_model(ctx, settings, emission_spectrum)
        type(fsps_context_t), intent(in)   :: ctx
        type(params), intent(in)           :: settings
        real(WP), dimension(:), intent(out) :: emission_spectrum
        
        integer :: idx_q, idx_u
        real(WP) :: w_q, w_u, gamma_frac
        real(WP), dimension(size(emission_spectrum)) :: spec_u_min, spec_u_max
        integer :: n_qpah, n_umin

        ! Grid Dimensions
        n_qpah = ctx%state%nqpah_dustem
        n_umin = ctx%state%numin_dustem

        ! 1. Interpolate Q_PAH (Polycyclic Aromatic Hydrocarbons fraction)
        idx_q = find_interval(ctx%state%qpaharr, settings%duste_qpah)
        idx_q = max(1, min(idx_q, n_qpah - 1))
        
        w_q = (settings%duste_qpah - ctx%state%qpaharr(idx_q)) / &
              (ctx%state%qpaharr(idx_q + 1) - ctx%state%qpaharr(idx_q))
        w_q = max(0.0_wp, min(w_q, 1.0_wp))

        ! 2. Interpolate U_min (Minimum Radiation Field Intensity)
        idx_u = find_interval(ctx%state%uminarr, settings%duste_umin)
        idx_u = max(1, min(idx_u, n_umin - 1))

        w_u = (settings%duste_umin - ctx%state%uminarr(idx_u)) / &
              (ctx%state%uminarr(idx_u + 1) - ctx%state%uminarr(idx_u))
        w_u = max(0.0_wp, min(w_u, 1.0_wp))

        ! 3. Gamma Fraction (Fraction of dust in high-intensity PDRs)
        gamma_frac = max(0.0_wp, min(settings%duste_gamma, 1.0_wp))

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
        emission_spectrum = (1.0_wp - gamma_frac) * spec_u_min + gamma_frac * spec_u_max
        
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
        
        real(WP), dimension(:), intent(in)  :: nu
        real(WP), dimension(:), intent(in)  :: shape_intrinsic
        real(WP), dimension(:), intent(in)  :: transmission_ism ! e^-tau
        real(WP), intent(in)                :: lum_absorbed_initial
        real(WP), dimension(:), intent(out) :: spec_final

        real(WP), dimension(size(nu)) :: profile_escaped
        real(WP) :: lum_escaped_profile, normalization_factor
        
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
            spec_final = 0.0_wp
        end if

    end subroutine calculate_dust_self_absorption

    !> Implementation of Cardelli, Clayton, & Mathis (1989) extinction curve.
    !> Includes the "hack" for smooth transitions used in the original FSPS.
    elemental function get_ccm89_curve_point(wavelength, r_v, uv_bump_strength) result(curve)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        real(WP), intent(in) :: r_v, uv_bump_strength
        real(WP) :: curve

        ! Array variables for the main calculation
        real(WP) :: wavenumber, wavenumber_term, poly_a, poly_b
        real(WP) :: temp_curve, wavenumber_clamped

        ! Scalar variables for calculating the Smoothing Hack (at x=3.3)
        real(WP) :: y_anchor, a_scalar, b_scalar
        real(WP) :: opt_val_at_break, nuv_val_at_break, continuity_correction

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

        ! 2. SCALAR CALCULATION
        ! --------------------
        wavenumber = get_wavenumber_point(wavelength)
        
        poly_a = 0.0_wp
        poly_b = 0.0_wp
        temp_curve = 0.0_wp

        ! Region 1: Infrared (0.3 < x < 1.1)
        if (wavenumber >= CCM_X_IR_MIN .and. wavenumber < CCM_X_OPT_MIN) then
            temp_curve = (CCM_IR_A_SCALE * wavenumber**CCM_IR_EXP) + &
                         (CCM_IR_B_SCALE * wavenumber**CCM_IR_EXP) / r_v
        
        ! Region 2: Optical / Near-IR (1.1 <= x < 3.3)
        else if (wavenumber >= CCM_X_OPT_MIN .and. wavenumber < CCM_X_NUV_MIN) then
            wavenumber_term = wavenumber - CCM_OPT_Y_SHIFT
            
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
        end if

        ! --- Regions 3 & 4: UV Base (3.3 <= x < 8.0) ---
        ! Both NUV and Mid-UV share the same base linear term and Drude profile.
        if (wavenumber >= CCM_X_NUV_MIN .and. wavenumber < CCM_X_FUV_MIN) then
            ! Base Linear Component
            poly_a = CCM_NUV_A_BASE(1) + CCM_NUV_A_BASE(2)*wavenumber
            poly_b = CCM_NUV_B_BASE(1) + CCM_NUV_B_BASE(2)*wavenumber
            
            ! Add UV Bump (Drude Profile)
            poly_a = poly_a + CCM_NUV_A_BUMP(1) / &
                     ((wavenumber - CCM_NUV_A_BUMP(2))**2 + CCM_NUV_A_BUMP(3)) * uv_bump_strength
                     
            poly_b = poly_b + CCM_NUV_B_BUMP(1) / &
                     ((wavenumber - CCM_NUV_B_BUMP(2))**2 + CCM_NUV_B_BUMP(3)) * uv_bump_strength
        end if

        ! Specific NUV Modifier (3.3 <= x < 5.9): Apply Continuity Correction
        if (wavenumber >= CCM_X_NUV_MIN .and. wavenumber < CCM_X_MUV_MIN) then
             temp_curve = poly_a + poly_b / r_v + continuity_correction * (CCM_X_NUV_MIN / wavenumber)**6
        end if

        ! Specific Mid-UV Modifier (5.9 <= x < 8.0): Apply Curvature Polynomials
        if (wavenumber >= CCM_X_MUV_MIN .and. wavenumber < CCM_X_FUV_MIN) then
            wavenumber_term = wavenumber - CCM_X_MUV_MIN
            
            poly_a = poly_a + CCM_MUV_A_POLY(1) * wavenumber_term**2 + CCM_MUV_A_POLY(2) * wavenumber_term**3
            poly_b = poly_b + CCM_MUV_B_POLY(1) * wavenumber_term**2 + CCM_MUV_B_POLY(2) * wavenumber_term**3
                
            temp_curve = poly_a + poly_b / r_v
        end if

        ! Region 5: Far-UV (x >= 8.0)
        ! Clamps input x to 12.0 (lambda = 833 A) to prevent divergence
        if (wavenumber >= CCM_X_FUV_MIN) then
            wavenumber_clamped = min(wavenumber, CCM_X_CUTOFF)
            
            wavenumber_term = wavenumber_clamped - CCM_X_FUV_MIN
            
            poly_a = CCM_FUV_A_COEFFS(0) + wavenumber_term*(CCM_FUV_A_COEFFS(1) + &
                wavenumber_term*(CCM_FUV_A_COEFFS(2) + wavenumber_term*CCM_FUV_A_COEFFS(3)))
            
            poly_b = CCM_FUV_B_COEFFS(0) + wavenumber_term*(CCM_FUV_B_COEFFS(1) + &
                wavenumber_term*(CCM_FUV_B_COEFFS(2) + wavenumber_term*CCM_FUV_B_COEFFS(3)))
            
            temp_curve = poly_a + poly_b / r_v
        end if

        curve = temp_curve
    end function get_ccm89_curve_point


    !> Implementation of Calzetti et al. (2000) starburst attenuation curve.
    elemental function get_calzetti_curve_point(wavelength) result(curve)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        real(WP) :: curve
        real(WP) :: wavenumber, extinction_k

        ! Convert to inverse microns (x = 1/lambda_um)
        wavenumber = get_wavenumber_point(wavelength)
        
        extinction_k = 0.0_wp
        
        ! Optical / NIR (0.63um < lambda <= 2.2um)
        if (wavelength > CALZ_LAM_BREAK .and. wavelength <= CALZ_LAM_IR_MAX) then
            extinction_k = CALZ_SCALE * (CALZ_OPT_COEFFS(0) + CALZ_OPT_COEFFS(1) * wavenumber) + CALZ_R_V
        end if
        
        ! UV / Optical (0.12um <= lambda <= 0.63um)
        if (wavelength >= CALZ_LAM_UV_MIN .and. wavelength <= CALZ_LAM_BREAK) then
            ! Use nested multiplication (Horner's method) for clarity and efficiency
            extinction_k = CALZ_R_V + CALZ_SCALE * ( &
                CALZ_UV_COEFFS(0) + wavenumber * ( &
                    CALZ_UV_COEFFS(1) + wavenumber * ( &
                        CALZ_UV_COEFFS(2) + wavenumber * CALZ_UV_COEFFS(3) &
                    ) &
                ) &
            )
        end if

        ! Result is A_lambda / A_V = k_lambda / R_V
        curve = extinction_k / CALZ_R_V
    end function get_calzetti_curve_point


    !> Implementation of Kriek & Conroy (2013): Calzetti + UV Bump + Tilt.
    elemental function get_kriek_conroy_curve_point(wavelength, tilt_index) result(curve)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        real(WP), intent(in) :: tilt_index
        real(WP) :: curve

        real(WP) :: base_calzetti, drude_profile
        real(WP) :: bump_amplitude ! E_b in paper

        ! 1. Base Calzetti (normalized to E(B-V), i.e., k_lambda scale)
        ! Note: Our helper `get_calzetti_curve` returns A_lambda/A_V.
        ! We must multiply by R_V to get back to k_lambda.
        base_calzetti = get_calzetti_curve_point(wavelength) * KC13_R_V_BASE

        ! 2. UV Bump (Drude Profile)
        ! Kriek & Conroy (2013) Eq 3: E_b = 0.85 - 1.9 * delta
        bump_amplitude = KC13_AMPL_INTERCEPT - KC13_AMPL_SLOPE * tilt_index
        
        drude_profile = bump_amplitude * (wavelength * KC13_BUMP_WIDTH)**2 / &
                        ( (wavelength**2 - UV_BUMP_CENTER**2)**2 + (wavelength * KC13_BUMP_WIDTH)**2 )

        ! 3. Combine with Tilt
        ! A_lambda = (k_calz + D_bump) / R_V * (lambda / 5500)^delta
        curve = (base_calzetti + drude_profile) / KC13_R_V_BASE * &
                (wavelength / V_BAND_ANGSTROMS)**tilt_index

    end function get_kriek_conroy_curve_point


!> Implementation of Reddy et al. (2015) MOSDEF curve.
    elemental function get_reddy_curve_point(wavelength) result(curve)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        real(WP) :: curve
        real(WP) :: wavenumber, extinction_k, wavenumber_clamped
        
        wavenumber = get_wavenumber_point(wavelength)
        extinction_k = 0.0_wp
        
        ! 1. UV Range (Lambda < 6000 A)
        if (wavelength < REDDY_LAM_BREAK) then
            wavenumber_clamped = min(wavenumber, REDDY_X_UV_MAX)
            
            extinction_k = REDDY_UV_COEFFS(0) + &
                           wavenumber_clamped * (REDDY_UV_COEFFS(1) + &
                           wavenumber_clamped * (REDDY_UV_COEFFS(2) + &
                           wavenumber_clamped * REDDY_UV_COEFFS(3))) + REDDY_OFFSET
        end if

        ! 2. Optical/NIR Range (0.60um <= lambda < 2.85um)
        if (wavelength >= REDDY_LAM_BREAK .and. wavelength < REDDY_LAM_IR_MAX) then
             extinction_k = REDDY_OPT_COEFFS(0) + &
                            wavenumber * (REDDY_OPT_COEFFS(1) + &
                            wavenumber * (REDDY_OPT_COEFFS(2) + &
                            wavenumber * REDDY_OPT_COEFFS(3))) + &
                            REDDY_OFFSET + REDDY_OPT_CORRECTION
        end if

        ! Convert k_lambda to A_lambda / A_V assuming R_V = 2.505
        curve = extinction_k / REDDY_R_V
    end function get_reddy_curve_point

    !> Converts wavelength (Angstroms) to wavenumber (inverse microns).
    !> Used frequently for dust curve parameterizations (CCM89, Calzetti, etc.).
    elemental function get_wavenumber_point(wavelength) result(wavenumber)
        !$acc routine seq
        real(WP), intent(in) :: wavelength
        real(WP) :: wavenumber
        
        ! x = 1 / lambda_microns = 10000 / lambda_angstroms
        wavenumber = 1.0e4_wp / wavelength
    end function get_wavenumber_point

    !> Computes the circumstellar optical depth (tau_1um) from physical parameters.
    !> See Villaume et al. (2015).
    pure function compute_circumstellar_optical_depth(ctx, c_rich_flag, m_act, log_l, log_g, log_mdot_iso) result(tau)
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in)  :: c_rich_flag
        real(WP), intent(in) :: m_act, log_l, log_g, log_mdot_iso
        real(WP) :: tau

        real(WP) :: radius_solar, period_days, velocity_exp, mdot_sol_yr
        real(WP) :: inner_radius_cm, dust_gas_ratio, kappa_eff

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
        radius_solar = sqrt(m_act * M_SOL * G_NEWTON / (10.0_wp**log_g)) / R_SOL

        ! Fundamental Pulsation Period (Days) - Villaume et al. (2015) relation
        period_days = 10.0_wp**(AGB_PER_INTERCEPT + &
                                AGB_PER_SLOPE_R * log10(radius_solar) + &
                                AGB_PER_SLOPE_M * log10(m_act))

        ! Expansion Velocity (km/s)
        velocity_exp = AGB_VEXP_INTERCEPT + AGB_VEXP_SLOPE * period_days
        velocity_exp = max(min(velocity_exp, AGB_VEXP_MAX), AGB_VEXP_MIN) 

        ! 3. Mass Loss Rate (M_sun/yr)
        ! ----------------------------
        if (ctx%use_isoc_mdot_val == 1) then
            ! Use Isochrone value (MIST only), capped
            mdot_sol_yr = min(10.0_wp**log_mdot_iso, VW93_MDOT_LIMIT_ISO)
        else
            ! Vassiliadis & Wood (1993) Prescription
            if (period_days < VW93_PER_THRESH) then
                if (m_act < VW93_MASS_THRESH) then
                    mdot_sol_yr = 10.0_wp**(VW93_BASE_INTERCEPT + VW93_BASE_SLOPE * period_days)
                else
                    mdot_sol_yr = 10.0_wp**(VW93_BASE_INTERCEPT + VW93_HIGH_SLOPE * &
                                           (period_days - 100.0_wp * (m_act - VW93_MASS_THRESH)))
                end if
            else
                ! Superwind phase
                mdot_sol_yr = (10.0_wp**log_l) / velocity_exp * VW93_SUPERWIND_NORM * YEAR_TO_SECOND / C_LIGHT
            end if
        end if

        ! 4. Shell Geometry
        ! -----------------
        ! Inner Radius (cm) scaling with Luminosity
        if (c_rich_flag == 1) then
            inner_radius_cm = AGB_RIN_FACTOR_C * (10.0_wp**log_l)**0.5_wp
        else
            inner_radius_cm = AGB_RIN_FACTOR_O * (10.0_wp**log_l)**0.5_wp
        end if

        ! Dust-to-Gas Ratio (delta) scaling with Velocity and Luminosity
        if (c_rich_flag == 1) then
            dust_gas_ratio = AGB_DELTA_C_RICH
        else
            dust_gas_ratio = AGB_DELTA_O_RICH
        end if
        
        dust_gas_ratio = dust_gas_ratio * (velocity_exp**2 / AGB_DTG_VEL_NORM) * &
                         ((10.0_wp**log_l / AGB_DTG_LUM_NORM)**AGB_DTG_LUM_EXP)

        ! 5. Final Optical Depth Calculation
        ! ----------------------------------
        ! tau = kappa * delta * Mdot / (4 * pi * R_in * v_exp)
        ! Note: v_exp is converted from km/s to cm/s (1E5 factor)
        
        tau = kappa_eff * dust_gas_ratio * (mdot_sol_yr * M_SOL / YEAR_TO_SECOND) / &
              inner_radius_cm / (4.0_wp * PI) / (velocity_exp * 1.0e5_wp)

    end function compute_circumstellar_optical_depth

end module fsps_dust