module fsps_io
    !> @brief 
    !> Unified module for FSPS Input/Output operations.
    !>
    !> @details
    !> This module consolidates all disk I/O functionality, separating these concerns from the 
    !> core physics engines. It serves two primary roles:
    !>
    !> **1. Input Handling (Tabular Data):**
    !> Responsible for reading and parsing external data files required for flexible star 
    !> formation histories. Specifically, it loads user-supplied tabular SFHs (`sfh_tab`)
    !> from disk, handling error checking and unit conversion (Gyr -> Years).
    !>
    !> **2. Output Formatting (Results):**
    !> Handles the writing of scientific results to standard FSPS output formats.
    !> - `.mags`: Photometric magnitudes and physical properties (Mass, SFR, Lbol).
    !> - `.spec`: Full spectral energy distributions (SEDs).
    !> - `.indx`: Spectral indices (e.g., Lick indices, D4000).
    !>
    !> **Design Philosophy:**
    !> This separation allows the physics modules (`fsps_csp`, `fsps_ssp`) to remain 
    !> "pure" calculation engines that return structured data types (`compspout`), 
    !> making them easier to wrap for Python/C++ interfaces without side effects (file handles).
    !>
    !> @note
    !> Legacy compatibility: The output file headers and column formats strictly match 
    !> the legacy FSPS behavior to ensure compatibility with existing analysis pipelines.

    use fsps_precision, only: WP
    use fsps_constants, only: SAFE_FLOOR, NTABMAX, NEMLINE, VERBOSE
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params, compspout
    use fsps_cosmology, only: vacuum_to_air

    implicit none
    private

    public :: load_tabular_sfh
    public :: write_csp_output_files

    ! Output Mode Constants (legacy write_compsp flags)
    integer, parameter :: OUTPUT_MAGS = 1
    integer, parameter :: OUTPUT_SPEC = 2
    integer, parameter :: OUTPUT_BOTH = 3
    integer, parameter :: OUTPUT_INDX = 4

contains

    ! ========================================================================
    ! INPUT ROUTINES
    ! ========================================================================

    !> @brief Loads a tabular SFH from file or validates pre-loaded arrays.
    !>
    !> @details
    !> Handles the setup for arbitrary star formation histories:
    !> * **SFH=2:** Reads the file specified by `sfh_filename` (default: 'sfh.dat') 
    !>   from the data directory. Expects columns: Time(Gyr), SFR, [Metallicity].
    !>   Converts time to Years and populates `ctx%state%sfh_tab`.
    !> * **SFH=3:** Validates that `sfh_tab` has been populated via the API.
    !>
    !> Also clips SFR values to `SAFE_FLOOR` to prevent log-domain errors.
    !>
    !> @param[inout] ctx  The FSPS context (updates sfh_tab, ntabsfh).
    !> @param[in]    pset User parameters.
    !> @param[in]    nzin Current metallicity index (legacy logic for Z-column reading).
    subroutine load_tabular_sfh(ctx, pset, nzin)
        type(fsps_context_t), intent(inout) :: ctx
        type(params), intent(in) :: pset
        integer, intent(in) :: nzin

        integer :: io_stat, file_unit, i_row
        character(len=256) :: file_path
        logical :: file_exists

        ! --------------------------------------------------------------------
        ! CASE 1: Read from File (SFH = 2)
        ! --------------------------------------------------------------------
        if (pset%sfh == 2) then

            ! Constraint: Tabular SFHs in FSPS start at T=0.
            if (pset%sf_start > SAFE_FLOOR) then
                print *, 'Error: [load_tabular_sfh] sf_start must be 0 for tabular SFH (type 2).'
                stop
            end if

            ! Construct File Path
            if (len_trim(pset%sfh_filename) == 0) then
                file_path = trim(ctx%sps_home) // '/data/sfh.dat'
            else
                file_path = trim(ctx%sps_home) // '/data/' // trim(pset%sfh_filename)
            end if

            ! Check Existence
            inquire(file=trim(file_path), exist=file_exists)
            if (.not. file_exists) then
                print *, 'Error: [load_tabular_sfh] File not found: ', trim(file_path)
                stop
            end if

            ! Open File
            open(newunit=file_unit, file=trim(file_path), status='old', &
                 action='read', iostat=io_stat)
            
            if (io_stat /= 0) then
                print *, 'Error: [load_tabular_sfh] Could not open file: ', trim(file_path)
                stop
            end if

            ! Read Loop
            do i_row = 1, NTABMAX
                ! Logic: Only read the 3rd column (Metallicity) if we are on the 
                ! specific Z index matching the total count (legacy behavior).
                ! Otherwise, read 2 columns and set Z=0.
                if (nzin == ctx%state%nz) then
                    read(file_unit, *, iostat=io_stat) ctx%state%sfh_tab(1, i_row), &
                                                       ctx%state%sfh_tab(2, i_row), &
                                                       ctx%state%sfh_tab(3, i_row)
                else
                    read(file_unit, *, iostat=io_stat) ctx%state%sfh_tab(1, i_row), &
                                                       ctx%state%sfh_tab(2, i_row)
                    ctx%state%sfh_tab(3, i_row) = 0.0_wp
                end if

                ! Check for End of File or Error
                if (io_stat < 0) exit ! EOF
                if (io_stat > 0) then
                    print *, 'Error: [load_tabular_sfh] Read error at line ', i_row
                    close(file_unit)
                    stop
                end if
            end do

            close(file_unit)

            ! Check if we hit the array limit without finding EOF
            if (io_stat == 0 .and. i_row > NTABMAX) then
                print *, 'Error: [load_tabular_sfh] SFH file exceeds NTABMAX rows.'
                print *, '       Increase NTABMAX in fsps_constants.'
                stop
            end if

            ! Store number of points (i_row is now the index of the first BAD read, so subtract 1)
            ctx%state%ntabsfh = i_row - 1

            ! Convert Time from Gyr to Years
            ctx%state%sfh_tab(1, 1:ctx%state%ntabsfh) = &
                ctx%state%sfh_tab(1, 1:ctx%state%ntabsfh) * 1.0e9_wp

        ! --------------------------------------------------------------------
        ! CASE 2: Pre-filled Array (SFH = 3)
        ! --------------------------------------------------------------------
        else if (pset%sfh == 3) then
            
            if (ctx%state%ntabsfh == 0) then
                print *, 'Error: [load_tabular_sfh] SFH=3 selected but sfh_tab is empty.'
                stop
            end if

        end if

        ! --------------------------------------------------------------------
        ! CLEANUP: Clip SFRs
        ! --------------------------------------------------------------------
        ! Ensure no negative or zero values cause log errors later.
        ! We use max(val, SAFE_FLOOR).
        if (pset%sfh == 2 .or. pset%sfh == 3) then
            where (ctx%state%sfh_tab(2, 1:ctx%state%ntabsfh) < SAFE_FLOOR)
                ctx%state%sfh_tab(2, 1:ctx%state%ntabsfh) = SAFE_FLOOR
            end where
        end if

    end subroutine load_tabular_sfh

    ! ========================================================================
    ! OUTPUT ROUTINES
    ! ========================================================================

    !> @brief Writes the computed CSP results to disk files.
    !>
    !> @details
    !> Matches legacy functionality to write:
    !> * `.mags`: Magnitudes and physical properties.
    !> * `.spec`: Full spectra (optionally converted to air wavelengths).
    !> * `.indx`: Spectral indices (if mode=4).
    !>
    !> @param[in] ctx        Context.
    !> @param[in] pset       User parameters.
    !> @param[in] results    Array of computed CSP structures.
    !> @param[in] outfile    Base filename (without extension).
    !> @param[in] write_mode Mode flag (1=Mags, 2=Spec, 3=Both, 4=Indx).
    subroutine write_csp_output_files(ctx, pset, results, outfile, write_mode)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in) :: pset
        type(compspout), intent(in) :: results(:)
        character(len=*), intent(in) :: outfile
        integer, intent(in) :: write_mode

        integer :: u_mag, u_spec, u_indx
        integer :: i, n_out
        character(len=1024) :: filepath
        character(len=32) :: fmt_mags
        integer :: nbands

        if (write_mode == 0) return

        n_out = size(results)

        ! 1. OPEN FILES AND WRITE HEADERS
        ! -------------------------------

        ! Mags File (.mags)
        if (write_mode == OUTPUT_MAGS .or. write_mode == OUTPUT_BOTH) then
            filepath = trim(ctx%output_home) // '/OUTPUTS/' // trim(outfile) // '.mags'
            open(newunit=u_mag, file=trim(filepath), status='replace', action='write')
            
            call write_header(ctx, pset, u_mag)
            
            if (pset%sfh == 0) then
                write(u_mag, '("#   Processing SSP", /,"#")')
                write(u_mag, '("#   log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)")')
            else
                write(u_mag, '("#")')
                write(u_mag, '("#   log(age) log(mass) Log(lbol) log(SFR) mags (see FILTER_LIST)")')
            end if
        end if

        ! Spec File (.spec)
        if (write_mode == OUTPUT_SPEC .or. write_mode == OUTPUT_BOTH) then
            filepath = trim(ctx%output_home) // '/OUTPUTS/' // trim(outfile) // '.spec'
            open(newunit=u_spec, file=trim(filepath), status='replace', action='write')
            
            call write_header(ctx, pset, u_spec)
            
            if (pset%sfh == 0) then
                write(u_spec, '("#   Processing SSP", /,"#")')
                write(u_spec, '("#   log(age) log(mass) Log(lbol) log(SFR) spectra")')
                if (n_out > 1) then
                    write(u_spec, '(I3,1x,I6)') ctx%state%ntfull, ctx%state%nspec
                else
                    write(u_spec, '(I3,1x,I6)') 1, ctx%state%nspec
                end if
            else
                write(u_spec, '("#")')
                write(u_spec, '("#   log(age) log(mass) Log(lbol) log(SFR) spectra")')
                if (n_out > 1) then
                    write(u_spec, '(I3,1x,I6)') ctx%state%ntfull, ctx%state%nspec
                else
                    write(u_spec, '(I3,1x,I6)') 1, ctx%state%nspec
                end if
            end if
            
            ! Write Wavelength Grid (Vacuum or Air)
            if (ctx%vactoair_flag_val == 0) then
                write(u_spec, '(50000(F15.4))') ctx%state%spec_lambda
            else
                write(u_spec, '(50000(F15.4))') vacuum_to_air(ctx%state%spec_lambda)
            end if
        end if

        ! Indices File (.indx)
        if (write_mode == OUTPUT_INDX) then
            filepath = trim(ctx%output_home) // '/OUTPUTS/' // trim(outfile) // '.indx'
            open(newunit=u_indx, file=trim(filepath), status='replace', action='write')
            call write_header(ctx, pset, u_indx)
            write(u_indx, '("#")')
            write(u_indx, '("#   log(age) indices (see allindices.dat)")')
        end if

        ! 2. WRITE DATA LOOPS
        ! -------------------

        ! Pre-format for mags: (Time, Mass, Lbol, SFR, [Mags...])
        nbands = size(results(1)%mags)
        write(fmt_mags, '("(F7.4,1x,3(F8.4,1x),", I3, "(F7.3,1x))")') nbands

        do i = 1, n_out
            
            ! Write Mags
            if (write_mode == OUTPUT_MAGS .or. write_mode == OUTPUT_BOTH) then
                write(u_mag, fmt_mags) &
                     results(i)%age, &
                     log10(results(i)%mass_csp + SAFE_FLOOR), &
                     results(i)%lbol_csp, &
                     log10(results(i)%sfr + SAFE_FLOOR), &
                     results(i)%mags
            end if

            ! Write Spec
            if (write_mode == OUTPUT_SPEC .or. write_mode == OUTPUT_BOTH) then
                write(u_spec, '(4(F8.4,1x))') &
                     results(i)%age, &
                     log10(results(i)%mass_csp + SAFE_FLOOR), &
                     results(i)%lbol_csp, &
                     log10(results(i)%sfr + SAFE_FLOOR)
                
                write(u_spec, '(50000(E14.6))') max(results(i)%spec, SAFE_FLOOR)
            end if

            ! Write Indx
            if (write_mode == OUTPUT_INDX) then
                write(u_indx, '(F8.4,99(F7.3,1x))') results(i)%age, results(i)%indx
            end if

        end do

        ! 3. CLOSE FILES
        ! --------------
        if (write_mode == OUTPUT_MAGS .or. write_mode == OUTPUT_BOTH) close(u_mag)
        if (write_mode == OUTPUT_SPEC .or. write_mode == OUTPUT_BOTH) close(u_spec)
        if (write_mode == OUTPUT_INDX) close(u_indx)

    end subroutine write_csp_output_files

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPERS
    ! ------------------------------------------------------------------------

    !> @brief Writes standard FSPS header information.
    subroutine write_header(ctx, pset, unit)
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in) :: pset
        integer, intent(in) :: unit
        
        real(WP) :: writeage

        ! Metallicity Info
        if (pset%sfh /= 2) then
             write(unit,'("#   Log(Z/Zsol): ",F6.3)') log10(ctx%state%zlegend(pset%zmet)/ctx%state%zsol)
        else
             write(unit,'("#   Log(Z/Zsol): tabulated")')
        end if

        ! Stellar Physics Info
        write(unit,'("#   Fraction of blue HB stars: ",F6.3,"; Ratio of BS to HB stars: ",F6.3)') &
             pset%fbhb, pset%sbss
        write(unit,'("#   Shift to TP-AGB [log(Teff),log(Lbol)]: ",F5.2,1x,F5.2)') &
             pset%delt, pset%dell

        ! IMF Info
        if (ctx%imf_type_val == 2) then
             write(unit,'("#   IMF: ",I1,", slopes= ",3F4.1)') ctx%imf_type_val, pset%imf1, pset%imf2, pset%imf3
        else if (ctx%imf_type_val == 3) then
             write(unit,'("#   IMF: ",I1,", cut-off= ",F4.2)') ctx%imf_type_val, pset%vdmc
        else
             write(unit,'("#   IMF: ",I1)') ctx%imf_type_val
        end if

        ! Mag System
        if (ctx%compute_vega_mags_val == 1) then
             write(unit,'("#   Mag Zero Point: Vega (not relevant for spec/indx files)")')
        else
             write(unit,'("#   Mag Zero Point: AB (not relevant for spec/indx files)")')
        end if

        ! SFH Parameter Info
        if (pset%sfh == 2 .or. pset%sfh == 3) then
            write(unit, '("#   SFH: tabulated input, dust=(",F6.2,",",F6.2,")")') &
                 pset%dust1, pset%dust2
        elseif (pset%sfh /= 0) then
            if (pset%tage > SAFE_FLOOR) then
                 writeage = pset%tage
            else
                 writeage = 10.0_wp**ctx%state%time_full(ctx%state%ntfull)/1.0e9_wp
            end if
            
            write(unit, '("#   SFH: Tage=",F6.2," Gyr, log(tau/Gyr)= ",F6.3, &
                         &", const= ",F6.3,", fb= ",F6.3,", tb= ",F6.2, &
                         &" Gyr, sf_start= ",F6.3,", dust=(",F6.2,",",F6.2,")")') &
                 writeage, log10(pset%tau), pset%const, pset%fburst, &
                 pset%tburst, pset%sf_start, pset%dust1, pset%dust2
        end if

    end subroutine write_header

end module fsps_io