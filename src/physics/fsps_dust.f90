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

    use fsps_types, only: sp, params, clight, mypi, tiny_number, nemline
    use fsps_context_types, only: fsps_context_t
    use fsps_interpolation, only: interpolate_linear, find_interval
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    implicit none
    private

    ! Public Interface
    ! public :: apply_dust_transmission_and_emission
    public :: compute_attenuation_curve
    ! public :: apply_agb_dust_screen
    ! public :: apply_agn_dust_emission

    ! ------------------------------------------------------------------------
    ! CONSTANTS: General Dust Parameters
    ! ------------------------------------------------------------------------
    real(sp), parameter :: V_BAND_ANGSTROMS  = 5500.0_sp
    real(sp), parameter :: UV_BUMP_CENTER    = 2175.0_sp
    real(sp), parameter :: CALZETTI_BREAK    = 6300.0_sp

    ! ------------------------------------------------------------------------
    ! CONSTANTS: Cardelli, Clayton, & Mathis (1989) Extinction Curve Parameters
    ! ------------------------------------------------------------------------

    ! Region Boundaries (in inverse microns, x = 1/lambda)
    real(sp), parameter :: CCM_X_IR_MIN  = 0.3_sp
    real(sp), parameter :: CCM_X_OPT_MIN = 1.1_sp
    real(sp), parameter :: CCM_X_NUV_MIN = 3.3_sp
    real(sp), parameter :: CCM_X_MUV_MIN = 5.9_sp
    real(sp), parameter :: CCM_X_FUV_MIN = 8.0_sp

    ! FUV Cutoff (Clamp values for x > 12.0, i.e., lambda < 833 A)
    real(sp), parameter :: CCM_X_CUTOFF = 12.0_sp
    
    ! Infrared Parameters (0.3 <= x < 1.1)
    real(sp), parameter :: CCM_IR_A_SCALE = 0.574_sp
    real(sp), parameter :: CCM_IR_B_SCALE = -0.527_sp
    real(sp), parameter :: CCM_IR_EXP     = 1.61_sp

    ! Optical Parameters (1.1 <= x < 3.3)
    ! Polynomial coefficients for y = x - 1.82 (Powers 0 through 7)
    real(sp), parameter :: CCM_OPT_Y_SHIFT = 1.82_sp
    real(sp), parameter :: CCM_OPT_A_COEFFS(0:7) = [ &
        1.0_sp,      0.17699_sp, -0.50447_sp, -0.02427_sp, &
        0.72085_sp,  0.01979_sp, -0.77530_sp,  0.32999_sp ]
    real(sp), parameter :: CCM_OPT_B_COEFFS(0:7) = [ &
        0.0_sp,      1.41338_sp,  2.28305_sp,  1.07233_sp, &
       -5.38434_sp, -0.62251_sp,  5.30260_sp, -2.09002_sp ]

    ! Near-UV Parameters (3.3 <= x < 5.9)
    ! Base Linear Terms: C1 + C2*x
    real(sp), parameter :: CCM_NUV_A_BASE(2) = [ 1.752_sp, -0.316_sp] 
    real(sp), parameter :: CCM_NUV_B_BASE(2) = [-3.09_sp,   1.825_sp]
    ! Drude Bump Terms: Scale / ((x-Pos)**2 + Width)
    real(sp), parameter :: CCM_NUV_A_BUMP(3) = [-0.104_sp, 4.67_sp, 0.341_sp]
    real(sp), parameter :: CCM_NUV_B_BUMP(3) = [ 1.206_sp, 4.62_sp, 0.263_sp]

    ! Mid-UV Parameters (5.9 <= x < 8.0) - Additions to NUV
    ! F(x) = C2*(x-5.9)**2 + C3*(x-5.9)**3
    real(sp), parameter :: CCM_MUV_A_POLY(2) = [-0.04473_sp, -0.009779_sp]
    real(sp), parameter :: CCM_MUV_B_POLY(2) = [ 0.2130_sp,   0.1207_sp]

    ! Far-UV Parameters (x >= 8.0)
    ! Polynomials in (x - 8.0)
    real(sp), parameter :: CCM_FUV_A_COEFFS(0:3) = [-1.073_sp, -0.628_sp, 0.137_sp, -0.070_sp]
    real(sp), parameter :: CCM_FUV_B_COEFFS(0:3) = [ 13.67_sp,  4.257_sp, -0.42_sp,  0.374_sp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Calzetti et al. (2000) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(sp), parameter :: CALZ_LAM_UV_MIN = 1200.0_sp  ! 0.12 microns
    real(sp), parameter :: CALZ_LAM_BREAK  = 6300.0_sp  ! 0.63 microns
    real(sp), parameter :: CALZ_LAM_IR_MAX = 22000.0_sp ! 2.20 microns

    ! General Parameters
    real(sp), parameter :: CALZ_R_V   = 4.05_sp
    real(sp), parameter :: CALZ_SCALE = 2.659_sp      ! Scaling factor k'

    ! UV/Optical Polynomial Coefficients (0.12 <= lambda < 0.63 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(sp), parameter :: CALZ_UV_COEFFS(0:3) = [-2.156_sp, 1.509_sp, -0.198_sp, 0.011_sp]

    ! Optical/NIR Linear Coefficients (0.63 <= lambda <= 2.2 um)
    ! Linear in (1/lambda_microns): c0 + c1*x
    real(sp), parameter :: CALZ_OPT_COEFFS(0:1) = [-1.857_sp, 1.040_sp]

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Kriek & Conroy (2013) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------
    
    ! Drude Profile Parameters
    real(sp), parameter :: KC13_BUMP_WIDTH = 350.0_sp  ! Delta lambda (Angstroms)
    
    ! Bump Amplitude Relationship: E_b = 0.85 - 1.9 * delta
    real(sp), parameter :: KC13_AMPL_INTERCEPT = 0.85_sp
    real(sp), parameter :: KC13_AMPL_SLOPE     = 1.9_sp

    ! The base model is tied to Calzetti's specific R_V
    real(sp), parameter :: KC13_R_V_BASE = 4.05_sp

    ! --------------------------------------------------------------------------
    ! CONSTANTS: Reddy et al. (2015) Attenuation Curve Parameters
    ! --------------------------------------------------------------------------

    ! Region Boundaries (Angstroms)
    real(sp), parameter :: REDDY_LAM_UV_MIN = 1500.0_sp
    real(sp), parameter :: REDDY_LAM_BREAK  = 6000.0_sp
    real(sp), parameter :: REDDY_LAM_IR_MAX = 28500.0_sp

    ! Region Boundary (Inverse Microns)
    real(sp), parameter :: REDDY_X_UV_MAX   = 1.0e4_sp / REDDY_LAM_UV_MIN

    ! General Parameters
    real(sp), parameter :: REDDY_R_V        = 2.505_sp
    real(sp), parameter :: REDDY_OFFSET     = 2.505_sp ! Base offset added to both curves

    ! UV Polynomial Coefficients (0.15 <= lambda < 0.60 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(sp), parameter :: REDDY_UV_COEFFS(0:3) = [-5.726_sp, 4.004_sp, -0.525_sp, 0.029_sp]
    
    ! Blueward Extrapolation Value (< 0.15 um)
    real(sp), parameter :: REDDY_UV_EXTRAP_VAL = 10.36_sp

    ! Optical/NIR Polynomial Coefficients (0.60 <= lambda < 2.85 um)
    ! Polynomial in (1/lambda_microns): c0 + c1*x + c2*x^2 + c3*x^3
    real(sp), parameter :: REDDY_OPT_COEFFS(0:3) = [-2.672_sp, -0.010_sp, 1.532_sp, -0.412_sp]
    
    ! Continuity Correction for Optical Range
    real(sp), parameter :: REDDY_OPT_CORRECTION = -0.036221981_sp

contains

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
        real(sp), dimension(:), intent(in) :: wavelengths
        integer, intent(in)                :: dust_type_id
        type(params), intent(in)           :: settings
        type(fsps_context_t), intent(in)   :: ctx
        real(sp), dimension(size(wavelengths)) :: attenuation_curve

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

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER ROUTINES
    ! ------------------------------------------------------------------------

    !> Implementation of Cardelli, Clayton, & Mathis (1989) extinction curve.
    !> Includes the "hack" for smooth transitions used in the original FSPS.
    pure function get_ccm89_curve(wavelengths, r_v, uv_bump_strength) result(curve)
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), intent(in) :: r_v, uv_bump_strength
        real(sp), dimension(size(wavelengths)) :: curve

        ! Array variables for the main calculation
        real(sp), dimension(size(wavelengths)) :: wavenumbers, wavenumber_term, poly_a, poly_b
        real(sp), dimension(size(wavelengths)) :: temp_curve, wavenumbers_clamped

        ! Scalar variables for calculating the Smoothing Hack (at x=3.3)
        real(sp) :: y_anchor, a_scalar, b_scalar
        real(sp) :: opt_val_at_break, nuv_val_at_break, continuity_correction

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
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), dimension(size(wavelengths)) :: curve
        real(sp), dimension(size(wavelengths)) :: wavenumbers, extinction_k

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
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), intent(in) :: tilt_index
        real(sp), dimension(size(wavelengths)) :: curve
        
        real(sp), dimension(size(wavelengths)) :: base_calzetti, drude_profile
        real(sp) :: bump_amplitude ! E_b in paper

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
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), dimension(size(wavelengths)) :: curve
        real(sp), dimension(size(wavelengths)) :: wavenumbers, extinction_k, wavenumbers_clamped
        
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
        real(sp), dimension(:), intent(in) :: wavelengths
        real(sp), dimension(size(wavelengths)) :: wavenumbers
        
        ! x = 1 / lambda_microns = 10000 / lambda_angstroms
        wavenumbers = 1.0e4_sp / wavelengths
    end function get_wavenumber

end module fsps_dust