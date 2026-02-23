module fsps_spectral_library
    !> @file fsps_spectral_library.f90
    !> @brief Core module for generating individual stellar spectra from theoretical and empirical libraries.
    !>
    !> @details
    !> This module serves as the "Spectral Engine" of FSPS. It replaces the legacy `GETSPEC` subroutine
    !> with a modernized, modular, and vectorized implementation. Its primary responsibility is to
    !> calculate the Spectral Energy Distribution (SED) for a single star (or isochrone point)
    !> defined by its physical parameters (Mass, Temperature, Luminosity, Gravity) and evolutionary phase.
    !>
    !> **Key Functionality:**
    !> 1. **Physical Consistency:** Enforces consistency between Mass, Luminosity, Temperature, and Gravity
    !>    by recalculating surface gravity and radius derived from the Stefan-Boltzmann law and Newton's laws.
    !>    This ensures that modifications to stellar properties (e.g., during the AGB phase) are physically 
    !>    reflected in the spectrum.
    !> 2. **Library Dispatch:** Automatically selects the appropriate spectral library based on:
    !>    - **Phase:** Main Sequence, Red Giant, TP-AGB, Post-AGB, Wolf-Rayet.
    !>    - **Temperature:** Hot star extensions (WMBasic) vs Cool star libraries.
    !>    - **Composition:** Carbon-rich vs Oxygen-rich AGB stars.
    !> 3. **High-Performance Interpolation:** Utilizes vectorized, whole-array arithmetic for
    !>    bilinear (2D) and linear (1D) interpolations, ensuring optimal performance on modern CPU architectures
    !>    by maintaining memory locality.
    !> 4. **Robustness:** Handles grid edge cases (partial overlaps) and numerical underflows using
    !>    safe floors and nearest-neighbor snapping where data is missing, preventing integration failures.
    !>
    !> **Supported Libraries:**
    !> - **Main:** BaSeL (theoretical) or MILES (empirical).
    !> - **Hot Stars:** WMBasic (O-stars), Rauch (Post-AGB).
    !> - **Wolf-Rayet:** CMFGEN (Smith et al. 2002).
    !> - **TP-AGB:** Lancon & Mouhcine (O-rich), Aringer or Lancon & Wood (C-rich).
    !>
    !> **Output Units:**
    !> Returns Luminosity Density in units of \f$ L_{\odot} / \mathrm{Hz} \f$. The output is automatically scaled:
    !> - For Surface Flux grids: Scaled by \f$ 4\pi R^2 \f$.
    !> - For Normalized SEDs: Scaled by \f$ L_{\mathrm{bol}} \f$.

    use fsps_precision, only: WP
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params
    use fsps_constants, only: SAFE_FLOOR, PI, YEAR_TO_SECOND, C_LIGHT, G_NEWTON, GRAVITY_L_M_T_COEFF, &
                              M_SOL, L_SOL, NDIM_WR, NDIM_PAGB, N_AGB_C, N_AGB_O, N_AGB_CAR, CSTAR_ARINGER
    use fsps_interpolation, only: find_interval
    use fsps_dust, only: apply_agb_dust_screen

    implicit none
    private

    public :: get_stellar_spectrum

    ! ------------------------------------------------------------------------
    ! CONSTANTS
    ! ------------------------------------------------------------------------
    integer, parameter :: SRC_MAIN_LIB = 0
    integer, parameter :: SRC_WMBASIC  = 1
    integer, parameter :: SRC_AGB_O    = 2
    integer, parameter :: SRC_AGB_C    = 3
    integer, parameter :: SRC_WR       = 4
    integer, parameter :: SRC_PAGB     = 5

    ! Library cuts
    real(WP), parameter :: LOGT_CUT_PAGB = 4.699_wp
    real(WP), parameter :: LOGT_CUT_AGB_COOL = 3.6_wp

    ! WR-specific constants
    real(WP), parameter :: V_INF_WR = 3.0e8_wp
    real(WP), parameter :: KAPPA_ES = 0.346_wp
    real(WP), parameter :: MIX_FAC  = 0.4_wp

contains

    !> @brief Main dispatcher: Calculates the spectrum for a single star.
    !>
    !> @details
    !> This routine replaces the legacy `GETSPEC`. It orchestrates the process:
    !> 1. Calculates consistent physical parameters (Radius, Gravity).
    !> 2. Determines the appropriate spectral library based on Phase/Temp.
    !> 3. Dispatches to the specific interpolation routine.
    !> 4. Applies physical scaling (Surface Flux -> Luminosity Density).
    !> 5. Applies circumstellar AGB dust (if applicable).
    !>
    !> @param[in]    ctx    The FSPS context (contains large lookup tables).
    !> @param[in]    pset   User parameters (contains dust/IMF settings).
    !> @param[in]    mact   Actual mass (M_sol).
    !> @param[in]    logt   Log(Teff).
    !> @param[in]    lbol   Luminosity (L_sol) (Linear units).
    !> @param[in]    logg   Log(g) (Input from isochrone).
    !> @param[in]    phase  Evolutionary phase flag.
    !> @param[in]    ffco   Composition flag (AGB C/O ratio, etc).
    !> @param[in]    lmdot  Log(Mass Loss Rate).
    !> @param[in]    wght   IMF Weight (used for context, not calculation).
    !> @param[out]   spec   The resulting spectrum (L_sun/Hz).
    subroutine get_stellar_spectrum(ctx, pset, mact, logt, lbol, logg, phase, ffco, lmdot, spec)
        !$acc routine seq
        type(fsps_context_t), intent(inout) :: ctx
        type(params),         intent(in)    :: pset
        real(WP),             intent(in)    :: mact, logt, lbol, logg
        real(WP),             intent(in)    :: phase, ffco, lmdot
        real(WP),             intent(out)   :: spec(:)

        ! Local variables
        real(WP) :: r2_cm, logg_calc, scale_factor
        integer  :: library_source
        integer  :: i_phase
        
        ! 1. Determine Spectral Library
        !    (Abstracts the complex phase/temperature decision tree)
        library_source = determine_library_source(ctx, phase, logt, ffco)
        i_phase = int(phase)

        ! 2. Compute physical quantities only when needed by the selected source.
        !    (Main library, WMBasic, and WR require consistent gravity/radius)
        select case (library_source)
        case (SRC_MAIN_LIB, SRC_WMBASIC, SRC_WR)
            call calculate_physical_parameters(mact, lbol, logt, logg, r2_cm, logg_calc)
        case default
            r2_cm = 0.0_wp
            logg_calc = logg
        end select

        ! 3. Calculate scale factor pre-dispatch
        select case (library_source)
        case (SRC_MAIN_LIB)
            ! Surface Flux -> Luminosity: Scale by Surface Area (4*pi*R^2) * 4pi (Eddington)
            ! Note: dividing by L_SOL to get solar units
            scale_factor = (16.0_wp * PI * PI * r2_cm) / L_SOL
        
        case default
            ! Normalized Libraries -> Luminosity: Scale by Lbol
            scale_factor = lbol
        end select

        ! 4. Dispatch to Interpolation Routines
        select case (library_source)
        
        case (SRC_PAGB)
            ! Post-AGB (Rauch)
            call get_pagb_spectrum(ctx, pset, logt, scale_factor, spec)
        
        case (SRC_WR)
            ! Wolf-Rayet (CMFGEN)
            call get_wr_spectrum(ctx, pset, logt, lmdot, r2_cm, ffco, scale_factor, spec)
            
        case (SRC_AGB_O)
            ! Oxygen-rich TP-AGB (Lancon & Mouhcine)
            call get_agb_o_spectrum(ctx, pset, logt, scale_factor, spec)
            
        case (SRC_AGB_C)
            ! Carbon-rich TP-AGB (Aringer or Lancon & Wood)
            call get_agb_c_spectrum(ctx, logt, scale_factor, spec)
            
        case (SRC_WMBASIC)
            ! Hot Main Sequence (WMBasic)
            ! Note: Uses the CALCULATED logg_calc for physical consistency
            call get_wmbasic_spectrum(ctx, pset, logt, logg_calc, scale_factor, spec)
            
        case default ! SRC_MAIN_LIB
            ! Standard Library (BaSeL/MILES/etc)
            ! Note: Uses the CALCULATED logg_calc
            call get_main_lib_spectrum(ctx, pset, logt, logg_calc, scale_factor, spec)
            
        end select

        ! 5. Apply AGB Circumstellar Dust (If applicable)
        !    Only for AGB phases (4=RGB/E-AGB, 5=TP-AGB) if model is enabled.
        if ( (i_phase == 4 .or. i_phase == 5) .and. &
             ctx%add_agb_dust_model_val == 1 .and. &
             pset%agb_dust > SAFE_FLOOR ) then
             
             ! Note: Passing LOG10(lbol) as expected by the legacy interface.
             ! Note: Passing input 'logg' to match legacy behavior exactly.
             call apply_agb_dust_screen(ctx, pset%agb_dust, spec, mact, &
                                        logt, log10(lbol), logg, ffco, lmdot)
        end if

        ! 6. Final safety clamp after all optional post-processing.
        spec = max(spec, SAFE_FLOOR)

    end subroutine get_stellar_spectrum

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPER PROCEDURES
    ! ------------------------------------------------------------------------

    !> @brief Calculates derived physical parameters (Radius, Gravity).
    !>
    !> @details
    !> Re-computes the surface gravity (log g) and Radius squared based on the 
    !> input Mass, Luminosity, and Temperature. 
    !>
    !> **Physical Consistency:**
    !> This routine intentionally ignores the `logg_in` (from the isochrone) 
    !> and recalculates it. This is strictly necessary because FSPS modifies 
    !> Log(L) and Log(T) for various evolutionary phases (e.g., AGB dust, 
    !> changes in Teff). Using the consistent parameters ensures the 
    !> Stefan-Boltzmann law and Newton's Gravity are satisfied.
    !>
    !> **Formulae:**
    !> 1. \f$ g = \frac{4 \pi G \sigma M T^4}{L} \f$
    !> 2. \f$ R^2 = \frac{G M}{g} \f$
    !>
    !> @param[in]  mact      Actual Mass [M_sol].
    !> @param[in]  lbol      Bolometric Luminosity [L_sol] (Linear, not Log).
    !> @param[in]  logt      Log(Temperature) [K].
    !> @param[in]  logg_in   Input Log(g) from isochrone (Unused for calculation).
    !> @param[out] r2        Stellar Radius Squared [cm^2].
    !> @param[out] logg_out  Calculated consistent Log(g) [cgs].
    pure subroutine calculate_physical_parameters(mact, lbol, logt, logg_in, r2, logg_out)
        !$acc routine seq
        real(WP), intent(in)  :: mact, lbol, logt, logg_in
        real(WP), intent(out) :: r2, logg_out

        ! Local variables
        real(WP) :: gravity_cgs

        ! 1. Calculate Consistent Surface Gravity (logg_out)
        !    logg = log10( C * M / L ) + 4 * logT
        !    GRAVITY_L_M_T_COEFF encapsulates (4 * pi * G * sigma) and unit conversions.
        
        if (lbol > tiny(0.0_wp)) then
            logg_out = log10(GRAVITY_L_M_T_COEFF * mact / lbol) + (4.0_wp * logt)
        else
            ! Fallback for zero-luminosity objects (e.g., dark remnants) to avoid NaN.
            ! We default to the input isochrone gravity if available.
            logg_out = logg_in 
        end if

        ! 2. Calculate Radius Squared (cm^2)
        !    R^2 = (G * M) / g
        !    Note: We calculate g = 10^logg_out
        
        gravity_cgs = 10.0_wp**logg_out
        
        ! Protect against division by zero if gravity is extremely small (unlikely)
        if (gravity_cgs > tiny(0.0_wp)) then
            r2 = (mact * M_SOL * G_NEWTON) / gravity_cgs
        else
            r2 = 0.0_wp
        end if

    end subroutine calculate_physical_parameters

    !> @brief Determines which spectral library to use for a given star.
    !>
    !> @details
    !> Implements the FSPS decision tree based on evolutionary phase, 
    !> temperature, composition, and user settings.
    !>
    !> **Decision Logic:**
    !> 1. **Post-AGB (Phase 6):** If LogT >= 4.699.
    !> 2. **Wolf-Rayet (Phase 9):** If enabled by user.
    !> 3. **TP-AGB (Phase 5):** If LogT < 3.6 (Cool Giants).
    !>    - O-rich if ffco <= 1.0
    !>    - C-rich if ffco > 1.0
    !> 4. **Hot Main Sequence:** If LogT > Cutoff (WMBasic).
    !> 5. **Default:** Main Library (BaSeL/MILES/etc).
    !>
    !> @param[in] ctx    FSPS context (for settings and grid bounds).
    !> @param[in] phase  Evolutionary phase (0=MS, 5=AGB, 6=PAGB, 9=WR, etc).
    !> @param[in] logt   Log(Temperature).
    !> @param[in] ffco   Composition flag (C/O ratio).
    !> 
    !> @return Integer ID of the library source (SRC_* constants).
    pure function determine_library_source(ctx, phase, logt, ffco) result(src)
        !$acc routine seq
        type(fsps_context_t), intent(in) :: ctx
        real(WP),             intent(in) :: phase, logt, ffco
        integer                          :: src
        
        ! Local variables
        real(WP) :: logt_cut_wmb
        integer  :: i_phase

        ! Convert phase to integer for cleaner comparisons
        i_phase = int(phase)
        
        ! Default to Main Library
        src = SRC_MAIN_LIB

        ! 1. Check Post-AGB (Phase 6)
        !    Only used for hot Post-AGB stars (LogT >= 4.699)
        if (i_phase == 6 .and. logt >= LOGT_CUT_PAGB) then
            src = SRC_PAGB
            return
        end if

        ! 2. Check Wolf-Rayet (Phase 9)
        if (i_phase == 9 .and. ctx%use_wr_spectra_val == 1) then
            src = SRC_WR
            return
        end if

        ! 3. Check TP-AGB (Phase 5) - Cool Giants
        !    Specific cutoff at LogT < 3.6
        if (i_phase == 5 .and. logt < LOGT_CUT_AGB_COOL) then
            if (ffco <= 1.0_wp) then
                src = SRC_AGB_O
            else
                src = SRC_AGB_C
            end if
            return
        end if

        ! 4. Check Hot Main Sequence (WMBasic)
        !    Usually applies to Phase 0 (Main Sequence) or generic phases.
        !    Calculate the switch-over temperature logic.
        !    Cutoff is MAX of the library start and the user setting.
        logt_cut_wmb = max(ctx%state%wmb_logt(1), ctx%logt_wmb_hot_val)

        if (logt > logt_cut_wmb) then
            ! Note: Legacy code checks phase=0, but often allows others if hot enough?
            ! Legacy strict check: IF (phase.EQ.0.0.AND.logt.GT.logt_cut)
            if (i_phase == 0) then
                src = SRC_WMBASIC
                return
            end if
        end if

        ! 5. Fallback: Main Library
        !    Used for MS, RGB, Horizontal Branch, and anything not caught above.
        src = SRC_MAIN_LIB

    end function determine_library_source

    ! ------------------------------------------------------------------------
    ! SPECIFIC LIBRARY RETRIEVALS
    ! ------------------------------------------------------------------------
    
    !> @brief Interpolates the Post-AGB spectral library (Rauch 2003).
    !>
    !> @details
    !> Used for high-temperature stars in the Post-AGB phase (Phase 6).
    !> Performs 1D Linear Interpolation in Log(Teff).
    !>
    !> **Metallicity Handling:**
    !> The library only supports two metallicities:
    !> * **Solar (Index 2):** Used if Z > 0.5 Z_sol.
    !> * **0.1 Solar (Index 1):** Used if Z <= 0.5 Z_sol.
    !>
    !> @param[in]    ctx            FSPS context.
    !> @param[in]    pset           Parameters (provides metallicity index).
    !> @param[in]    logt           Log(Temperature).
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated spectrum (Normalized).
    pure subroutine get_pagb_spectrum(ctx, pset, logt, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        type(params),         intent(in)  :: pset
        real(WP),             intent(in)  :: logt
        real(WP),             intent(in)  :: scale_factor
        real(WP),             intent(out) :: spec(:)

        ! Local variables
        integer  :: jlo, k_z
        real(WP) :: t, w1, w2, metallicity_ratio

        ! 1. Determine Metallicity Grid (Binary Selection)
        !    Compare current Z against Solar Z.
        metallicity_ratio = ctx%state%zlegend(pset%zmet) / ctx%state%zsol
        
        if (metallicity_ratio > 0.5_wp) then
            k_z = 2 ! Use Solar models
        else
            k_z = 1 ! Use 0.1 Solar models
        end if

        ! 2. Find Indices (Temperature)
        jlo = find_interval(ctx%state%pagb_logt, logt)
        jlo = max(1, min(jlo, NDIM_PAGB - 1))

        ! 3. Calculate Weight
        t = (logt - ctx%state%pagb_logt(jlo)) / &
            (ctx%state%pagb_logt(jlo+1) - ctx%state%pagb_logt(jlo))
        
        t = max(0.0_wp, min(t, 1.0_wp)) ! Clamp

        ! 4. Compute Scaled weights
        w1 = (1.0_wp - t) * scale_factor
        w2 = (t)          * scale_factor
        
        ! 5. Interpolate
        !    pagb_spec is (Lambda, Temp, Metallicity)
        spec = w1 * ctx%state%pagb_spec(:, jlo,   k_z) + &
               w2 * ctx%state%pagb_spec(:, jlo+1, k_z)

    end subroutine get_pagb_spectrum

    !> @brief Interpolates Wolf-Rayet spectral libraries (WN/WC).
    !>
    !> @details
    !> Source: Smith et al. (2002) (CMFGEN).
    !> 
    !> **Physics (MIST Isochrones):**
    !> If using MIST, this routine calculates a modified "Wind Temperature" (`twr`)
    !> accounting for the optical depth of the stellar wind (Maeder 1990).
    !> 
    !> **Library Selection:**
    !> Determined by `ffco` (composition flag):
    !> * `ffco < 10`: WN (Nitrogen-rich).
    !> * `ffco >= 10`: WC (Carbon-rich).
    !>
    !> @param[in]    ctx            FSPS context.
    !> @param[in]    pset           Parameters (provides metallicity index).
    !> @param[in]    logt           Log(Temperature) (Photospheric).
    !> @param[in]    lmdot          Log(Mass Loss Rate) [M_sol/yr].
    !> @param[in]    r2             Stellar Radius Squared [cm^2].
    !> @param[in]    ffco           Composition flag.
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated spectrum (Normalized).
    pure subroutine get_wr_spectrum(ctx, pset, logt, lmdot, r2, ffco, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        type(params),         intent(in)  :: pset
        real(WP),             intent(in)  :: logt, lmdot, r2, ffco, scale_factor
        real(WP),             intent(out) :: spec(:)

        ! Local variables
        integer  :: jlo, nz
        real(WP) :: t, w1, w2, twr
        real(WP) :: r_phot, r_wr, mdot_cgs

        ! 1. Calculate Effective "Wind Temperature" (twr)
        if (ctx%state%isoc_type == 'mist') then
            ! Convert Mdot to g/s
            mdot_cgs = (10.0_wp**lmdot) * M_SOL / YEAR_TO_SECOND
            
            ! Calculate Radius of the pseudo-photosphere in the wind
            r_phot = sqrt(r2)
            r_wr   = r_phot + (3.0_wp * KAPPA_ES * abs(mdot_cgs)) / &
                              (8.0_wp * PI * V_INF_WR)

            ! Calculate Temperature at R_wr
            twr = (10.0_wp**logt) * sqrt(r_phot / r_wr)

            ! Empirical mixing from Smith et al. 2002
            twr = log10((1.0_wp - MIX_FAC) * (10.0_wp**logt) + MIX_FAC * twr)
        else
            ! Standard Isochrones: Use photospheric temperature
            twr = logt
        end if

        nz = pset%zmet

        ! 2. Select Grid and Interpolate (Directly, no pointers)
        if (ffco < 10.0_wp) then
            ! --- WN Spectra ---
            jlo = find_interval(ctx%state%wrn_logt, twr)
            jlo = max(1, min(jlo, NDIM_WR - 1))

            t = (twr - ctx%state%wrn_logt(jlo)) / &
                (ctx%state%wrn_logt(jlo+1) - ctx%state%wrn_logt(jlo))
            t = max(0.0_wp, min(t, 1.0_wp))

            w1 = (1.0_wp - t) * scale_factor
            w2 = (t)          * scale_factor

            spec = w1 * ctx%state%wrn_spec(:, jlo, nz) + &
                   w2 * ctx%state%wrn_spec(:, jlo+1, nz)
        else
            ! --- WC Spectra ---
            jlo = find_interval(ctx%state%wrc_logt, twr)
            jlo = max(1, min(jlo, NDIM_WR - 1))

            t = (twr - ctx%state%wrc_logt(jlo)) / &
                (ctx%state%wrc_logt(jlo+1) - ctx%state%wrc_logt(jlo))
            t = max(0.0_wp, min(t, 1.0_wp))

            w1 = (1.0_wp - t) * scale_factor
            w2 = (t)          * scale_factor

            spec = w1 * ctx%state%wrc_spec(:, jlo, nz) + &
                   w2 * ctx%state%wrc_spec(:, jlo+1, nz)
        end if

    end subroutine get_wr_spectrum

    !> @brief Interpolates the Oxygen-rich TP-AGB spectral library.
    !>
    !> @details
    !> Source: Lancon & Mouhcine (2002) Empirical Spectra.
    !> Performs 1D Linear Interpolation in Log(Teff).
    !>
    !> **Metallicity Handling:**
    !> The temperature grid `agb_logt_o` varies with metallicity (Z), so the 
    !> lookup depends on `pset%zmet`.
    !>
    !> **Unit Conversion:**
    !> Source units are F_lambda. Converted to F_nu internal units using 
    !> (lambda^2 / c). The result is normalized to L_bol = 1.
    !>
    !> @param[in]    ctx            FSPS context.
    !> @param[in]    pset           Parameters (provides metallicity index).
    !> @param[in]    logt           Log(Temperature).
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated spectrum (Normalized L_nu).
    pure subroutine get_agb_o_spectrum(ctx, pset, logt, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        type(params),         intent(in)  :: pset
        real(WP),             intent(in)  :: logt
        real(WP),             intent(in)  :: scale_factor
        real(WP),             intent(out) :: spec(:)

        integer  :: jlo, nz
        real(WP) :: t, w1, w2
        
        nz = pset%zmet

        ! Search directly on the strided row to avoid per-call copies.
        jlo = find_interval_row_2d(ctx%state%agb_logt_o, nz, N_AGB_O, logt)
        jlo = max(1, min(jlo, N_AGB_O - 1))

        t = (logt - ctx%state%agb_logt_o(nz, jlo)) / &
            (ctx%state%agb_logt_o(nz, jlo+1) - ctx%state%agb_logt_o(nz, jlo))
        
        t = max(0.0_wp, min(t, 1.0_wp)) 

        w1 = (1.0_wp - t)
        w2 = t

        spec = (scale_factor / C_LIGHT) * (ctx%state%spec_lambda**2) * &
               ( w1 * ctx%state%agb_spec_o(:, jlo) + &
                 w2 * ctx%state%agb_spec_o(:, jlo+1) )

    end subroutine get_agb_o_spectrum

    !> @brief Find interpolation interval on a strided row in a 2D array.
    !>
    !> @details
    !> Specialized helper for lookups like `array2d(row, :)` where row slices
    !> are non-contiguous in memory (column-major layout). This avoids creating
    !> temporary contiguous copies on every call.
    pure function find_interval_row_2d(array2d, row, ncol, value) result(idx)
        !$acc routine seq
        real(WP), intent(in) :: array2d(:, :)
        integer,  intent(in) :: row, ncol
        real(WP), intent(in) :: value
        integer :: idx

        integer :: lower, upper, mid
        logical :: is_ascending

        if (ncol < 2) then
            idx = 1
            return
        end if

        is_ascending = (array2d(row, ncol) >= array2d(row, 1))

        lower = 0
        upper = ncol + 1

        do while (upper - lower > 1)
            mid = (upper + lower) / 2
            if (is_ascending .eqv. (value >= array2d(row, mid))) then
                lower = mid
            else
                upper = mid
            end if
        end do

        if (value == array2d(row, 1)) then
            idx = 1
        else if (value == array2d(row, ncol)) then
            idx = ncol - 1
        else
            idx = lower
        end if

    end function find_interval_row_2d

    !> @brief Interpolates the Carbon-rich TP-AGB spectral library.
    !>
    !> @details
    !> Handles two possible source libraries based on the configuration:
    !> 1. **Aringer et al. (2009):** Synthetic spectra (Default if CSTAR_ARINGER=1).
    !>    Native units are F_nu / L_bol.
    !> 2. **Lancon & Wood (2002):** Empirical spectra.
    !>    Native units are F_lambda. This routine converts them to F_nu
    !>    using (lambda^2 / c) so the output units are consistent.
    !>
    !> Performs 1D Linear Interpolation in Log(Teff).
    !>
    !> @param[in]    ctx            FSPS context containing AGB grids.
    !> @param[in]    logt           Log(Temperature).
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated spectrum (Normalized L_nu).
    pure subroutine get_agb_c_spectrum(ctx, logt, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        real(WP),             intent(in)  :: logt
        real(WP),             intent(in)  :: scale_factor
        real(WP),             intent(out) :: spec(:)

        ! Local variables
        integer  :: jlo, n_grid
        real(WP) :: t, w1, w2
        
        ! 1. Select the Library Source
        if (CSTAR_ARINGER == 1) then
            ! --- Aringer et al. (2009) ---
            ! Point to Aringer Data
            n_grid = N_AGB_CAR
            
            ! 2. Find Indices (1D Interpolation)
            jlo = find_interval(ctx%state%agb_logt_car, logt)
            jlo = max(1, min(jlo, n_grid - 1))

            ! 3. Calculate Weight
            t = (logt - ctx%state%agb_logt_car(jlo)) / &
                (ctx%state%agb_logt_car(jlo+1) - ctx%state%agb_logt_car(jlo))
            t = max(0.0_wp, min(t, 1.0_wp)) ! Clamp

            w1 = (1.0_wp - t) * scale_factor
            w2 = (t)          * scale_factor

            ! 4. Interpolate (Units are already F_nu)
            spec = w1 * ctx%state%agb_spec_car(:, jlo) + &
                   w2 * ctx%state%agb_spec_car(:, jlo+1)

        else
            ! --- Lancon & Wood (2002) ---
            ! Point to LW02 Data
            n_grid = N_AGB_C

            ! 2. Find Indices
            jlo = find_interval(ctx%state%agb_logt_c, logt)
            jlo = max(1, min(jlo, n_grid - 1))

            ! 3. Calculate Weight
            t = (logt - ctx%state%agb_logt_c(jlo)) / &
                (ctx%state%agb_logt_c(jlo+1) - ctx%state%agb_logt_c(jlo))
            t = max(0.0_wp, min(t, 1.0_wp)) ! Clamp

            w1 = (1.0_wp - t)
            w2 = t

            ! 4. Interpolate and Convert Units (F_lambda -> F_nu)
            !    F_nu ~ lambda^2 * F_lambda
            spec = (scale_factor / C_LIGHT) * (ctx%state%spec_lambda**2) * &
               ( w1 * ctx%state%agb_spec_c(:, jlo) + &
                 w2 * ctx%state%agb_spec_c(:, jlo+1) )

        end if

    end subroutine get_agb_c_spectrum

    !> @brief Interpolates the WMBasic spectral library (Hot Stars / O-stars).
    !>
    !> @details
    !> Used for high-temperature Main Sequence stars (typically > 25,000K or > 50,000K).
    !> Performs a 2D Bilinear Interpolation in (logT, logg) for a fixed Metallicity.
    !> The resulting spectrum is normalized to unity (L_bol = 1), to be scaled later.
    !>
    !> @param[in]    ctx            FSPS context containing the `wmb_spec` grid.
    !> @param[in]    pset           Parameters (provides metallicity index).
    !> @param[in]    logt           Log(Temperature) of the star.
    !> @param[in]    logg           Log(Gravity) of the star.
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated spectrum (Normalized L_nu).
    pure subroutine get_wmbasic_spectrum(ctx, pset, logt, logg, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        type(params),         intent(in)  :: pset
        real(WP),             intent(in)  :: logt, logg
        real(WP),             intent(in)  :: scale_factor
        real(WP),             intent(out) :: spec(:)

        ! Local variables
        integer  :: jlo, klo   ! Grid indices
        real(WP) :: t, u       ! Interpolation weights
        real(WP) :: w1, w2, w3, w4
        integer  :: nz, nt, ng

        ! 1. Setup Dimensions
        nz = pset%zmet
        nt = size(ctx%state%wmb_logt)
        ng = size(ctx%state%wmb_logg)

        ! 2. Find Indices & Clamp
        !    Ensure we stay within bounds [1, N-1] for valid interpolation
        jlo = find_interval(ctx%state%wmb_logt, logt)
        jlo = max(1, min(jlo, nt - 1))

        klo = find_interval(ctx%state%wmb_logg, logg)
        klo = max(1, min(klo, ng - 1))

        ! 3. Calculate Weights (Bilinear)
        !    t = fraction along logT, u = fraction along logg
        t = (logt - ctx%state%wmb_logt(jlo)) / &
            (ctx%state%wmb_logt(jlo+1) - ctx%state%wmb_logt(jlo))
        
        u = (logg - ctx%state%wmb_logg(klo)) / &
            (ctx%state%wmb_logg(klo+1) - ctx%state%wmb_logg(klo))

        ! 4. Strict Clamping (No Extrapolation)
        !    If logT > max_grid_T, we clamp to the hottest model.
        t = max(0.0_wp, min(t, 1.0_wp))
        u = max(0.0_wp, min(u, 1.0_wp))

        ! 5. Pre-calculate Scalar Coefficients
        w1 = scale_factor * (1.0_wp - t) * (1.0_wp - u)  ! Corner (j, k)
        w2 = scale_factor * t * (1.0_wp - u)             ! Corner (j+1, k)
        w3 = scale_factor * t * u                        ! Corner (j+1, k+1)
        w4 = scale_factor * (1.0_wp - t) * u             ! Corner (j, k+1)

        ! 6. Vectorized Interpolation
        !    wmb_spec is usually rank 4: (lambda, Z, logT, logg)
        spec = w1 * ctx%state%wmb_spec(:, nz, jlo,   klo)   + &
               w2 * ctx%state%wmb_spec(:, nz, jlo+1, klo)   + &
               w3 * ctx%state%wmb_spec(:, nz, jlo+1, klo+1) + &
               w4 * ctx%state%wmb_spec(:, nz, jlo,   klo+1)

    end subroutine get_wmbasic_spectrum

    !> @brief Interpolates the main spectral library (BaSeL, MILES, etc.).
    !>
    !> @details
    !> Performs a 2D Bilinear Interpolation in (logT, logg) for a fixed Metallicity (Z).
    !> 
    !> **Grid Boundary Handling:**
    !> It checks the flux at a reference wavelength (typically 5000A) to ensure 
    !> the grid points contain valid data. 
    !> 1. If all 4 corners are valid: Standard Bilinear Interpolation.
    !> 2. If all 4 corners are invalid: Returns SAFE_FLOOR.
    !> 3. If some corners are invalid (partial overlap): Snaps to the *closest* !>    valid neighbor (highest weight) to avoid interpolating zeros.
    !>
    !> @param[in]    ctx            FSPS context containing the main `speclib` (4D grid).
    !> @param[in]    pset           Parameters (provides metallicity index).
    !> @param[in]    logt           Log(Temperature) of the star.
    !> @param[in]    logg           Log(Gravity) of the star.
    !> @param[in]    scale_factor   Scaling factor to convert normalized spectrum to Luminosity Density.
    !> @param[out]   spec           Interpolated Surface Flux (normalized units).
    pure subroutine get_main_lib_spectrum(ctx, pset, logt, logg, scale_factor, spec)
        !$acc routine seq
        type(fsps_context_t), intent(in)  :: ctx
        type(params),         intent(in)  :: pset
        real(WP),             intent(in)  :: logt, logg
        real(WP),             intent(in)  :: scale_factor
        real(WP),             intent(out) :: spec(:)

        ! Local variables
        integer  :: jlo, klo   ! Grid indices for logT, logg
        real(WP) :: t, u       ! Interpolation weights (0..1)
        real(WP) :: w1, w2, w3, w4 ! Final bilinear weights
        integer  :: nz, nt, ng ! Grid dimensions
        integer  :: idx_check  ! Wavelength index for validity check
        
        ! Validity flags
        logical  :: valid(4)
        real(WP) :: check_flux

        ! 1. Setup & Dimensions
        nz = pset%zmet
        nt = size(ctx%state%speclib_logt)
        ng = size(ctx%state%speclib_logg)
        idx_check = ctx%state%whlam5000

        ! 2. Find Indices & Clamp
        !    We clamp to (1, N-1) to ensure we always have a valid upper neighbor (j+1).
        jlo = find_interval(ctx%state%speclib_logt, logt)
        jlo = max(1, min(jlo, nt - 1))

        klo = find_interval(ctx%state%speclib_logg, logg)
        klo = max(1, min(klo, ng - 1))

        ! 3. Calculate Linear Weights
        !    t = fraction along logT, u = fraction along logg
        t = (logt - ctx%state%speclib_logt(jlo)) / &
            (ctx%state%speclib_logt(jlo+1) - ctx%state%speclib_logt(jlo))
        t = max(0.0_wp, min(t, 1.0_wp)) ! Strict clamping to avoid extrapolation

        u = (logg - ctx%state%speclib_logg(klo)) / &
            (ctx%state%speclib_logg(klo+1) - ctx%state%speclib_logg(klo))
        u = max(0.0_wp, min(u, 1.0_wp))

        ! 4. Check Grid Validity (Off-Grid Detection)
        !    Check the flux at the reference wavelength for all 4 corners.
        !    Corner mapping: 
        !    1: (j, k),  2: (j+1, k),  3: (j+1, k+1),  4: (j, k+1)
        
        ! Corner 1: (jlo, klo)
        check_flux = ctx%state%speclib(idx_check, nz, jlo, klo)
        valid(1)   = (check_flux > SAFE_FLOOR)

        ! Corner 2: (jlo+1, klo)
        check_flux = ctx%state%speclib(idx_check, nz, jlo+1, klo)
        valid(2)   = (check_flux > SAFE_FLOOR)

        ! Corner 3: (jlo+1, klo+1)
        check_flux = ctx%state%speclib(idx_check, nz, jlo+1, klo+1)
        valid(3)   = (check_flux > SAFE_FLOOR)

        ! Corner 4: (jlo, klo+1)
        check_flux = ctx%state%speclib(idx_check, nz, jlo, klo+1)
        valid(4)   = (check_flux > SAFE_FLOOR)

        ! 5. Interpolation Dispatch
        if (all(valid)) then
            ! --- Standard Case: All corners valid ---
            
            ! Pre-calculate scalar weights
            w1 = scale_factor * (1.0_wp - t) * (1.0_wp - u)
            w2 = scale_factor * t * (1.0_wp - u)
            w3 = scale_factor * t * u
            w4 = scale_factor * (1.0_wp - t) * u

            ! Vectorized accumulation
            ! Note: Accessing speclib is contiguous in the first dimension (wavelength).
            spec = w1 * ctx%state%speclib(:, nz, jlo,   klo)   + &
                   w2 * ctx%state%speclib(:, nz, jlo+1, klo)   + &
                   w3 * ctx%state%speclib(:, nz, jlo+1, klo+1) + &
                   w4 * ctx%state%speclib(:, nz, jlo,   klo+1)

        else if (.not. any(valid)) then
            ! --- Failure Case: All corners invalid ---
            spec = SAFE_FLOOR

        else
            ! --- Edge Case: Partial Overlap ---
            ! The legacy code blindly overwrote spectra based on loop order.
            ! We improve this by snapping to the valid corner with the HIGHEST weight.
            ! This minimizes the jump when crossing from a valid region to a partial hole.

            w1 = (1.0_wp - t) * (1.0_wp - u)
            w2 = t * (1.0_wp - u)
            w3 = t * u
            w4 = (1.0_wp - t) * u
            
            ! Zero out weights for invalid corners
            if (.not. valid(1)) w1 = -1.0_wp
            if (.not. valid(2)) w2 = -1.0_wp
            if (.not. valid(3)) w3 = -1.0_wp
            if (.not. valid(4)) w4 = -1.0_wp

            ! Find max weight
            if (w1 >= w2 .and. w1 >= w3 .and. w1 >= w4) then
                spec = ctx%state%speclib(:, nz, jlo, klo) * scale_factor
            else if (w2 >= w1 .and. w2 >= w3 .and. w2 >= w4) then
                spec = ctx%state%speclib(:, nz, jlo+1, klo) * scale_factor
            else if (w3 >= w1 .and. w3 >= w2 .and. w3 >= w4) then
                spec = ctx%state%speclib(:, nz, jlo+1, klo+1) * scale_factor
            else
                spec = ctx%state%speclib(:, nz, jlo, klo+1) * scale_factor
            end if

        end if

    end subroutine get_main_lib_spectrum

end module fsps_spectral_library