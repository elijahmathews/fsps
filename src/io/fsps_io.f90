module fsps_io
    !> @brief Unified Input/Output abstraction layer for the Flexible Stellar Population Synthesis (FSPS) code.
    !>
    !> @details
    !> The `fsps_io` module serves as the exclusive interface between FSPS internal data structures 
    !> (stored in `fsps_context_t`) and the disk file system. It encapsulates all logic required to 
    !> parse, load, and write the heterogeneous file formats used by FSPS.
    !>
    !> ### Functional Scope
    !>
    !> **1. System Initialization (Setup Phase):**
    !> Responsible for populating the core physics arrays during startup.
    !> - **Stellar Models:** Parsing ASCII isochrone tables (Padova, MIST, BaSTI, etc.) and 
    !>   handling BPASS binary models.
    !> - **Spectral Libraries:** Loading pre-compiled binary spectral cubes (MILES, BaSeL) and 
    !>   compiling raw ASCII libraries (C3K, CKC14).
    !> - **Auxiliary Physics:** !>     - Wolf-Rayet (CMFGEN) and O-Star (WMBasic) libraries.
    !>     - AGB (Lancon & Wood, Aringer) and Post-AGB (Rauch) spectra.
    !>     - Nebular emission grids (Cloudy) and Dust emission templates (Draine & Li, THEMIS).
    !>     - Dust attenuation curves (Calzetti, SMC, Witt & Gordon).
    !>     - Line Spread Functions (LSF) for instrument smoothing.
    !> - **Calibration:** Loading filter transmission curves and standard SEDs (Vega, Sun).
    !>
    !> **2. Runtime Operations:**
    !> - **Tabular SFHs:** Efficient parsing of user-supplied star formation history files (`sfh_tab`).
    !>
    !> **3. Data Export (Results):**
    !> - **CSP Outputs:** Writing computed spectral energy distributions (`.spec`), photometry (`.mags`), 
    !>   and spectral indices (`.indx`).
    !> - **Isochrone Data:** Dumping color-magnitude diagrams (`.cmd`) for specific populations.
    !>
    !> ### Design Philosophy
    !> - **Separation of Concerns:** This module handles *how* data is read (file formats, parsing, 
    !>   precision promotion, byte-swapping), while the physics modules (`fsps_initialization`, `fsps_ssp`) 
    !>   determine *what* the data represents physically (e.g. applying gravity corrections or 
    !>   interpolating weights).
    !> - **Legacy Compatibility:** Strictly adheres to the binary and ASCII file formats established 
    !>   in FSPS v3.2 to ensure compatibility with existing data directories.
    !> - **Memory Safety:** Utilizes modern Fortran `allocatable` arrays for large temporary buffers 
    !>   to prevent stack overflow errors common in the legacy codebase.

    use fsps_precision, only: WP
    use fsps_constants, only: VERBOSE, SAFE_FLOOR, C_LIGHT, &
                              NM, NLINES, NDIM_LOGT, NDIM_LOGG, NTABMAX, &
                              NLAM_NEBCONT, NEBNZ, NEBNAGE, NEBNIP, &
                              NDIM_WR, BASEL_STR, N_AGB_C, N_AGB_O, N_AGB_CAR, &
                              NDIM_PAGB, NDIM_WMB_LOGT, NDIM_WMB_LOGG, &
                              NAGNDUST, NAGNDUST_SPEC
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params, compspout
    use fsps_cosmology, only: vacuum_to_air
    use fsps_integration, only: integrate_trapezoid_array
    use fsps_interpolation, only: find_interval, interpolate_linear
    use, intrinsic :: iso_fortran_env, only: real32, iostat_end, error_unit, file_storage_size

    implicit none
    private

    ! --- Public API: Input (Setup) ---
    public :: load_zlegend_file
    public :: load_wavelength_grid
    public :: load_spectral_resolution
    public :: read_isochrone_database
    public :: read_spectral_binary
    public :: read_bpass_data
    public :: load_filter_definitions
    public :: load_nebular_grid
    public :: load_dust_emission_table
    public :: load_attenuation_curves
    public :: apply_legacy_filter_norm
    public :: load_standard_sed
    public :: load_index_definitions
    public :: load_wr_spectra
    public :: load_agb_spectra
    public :: load_post_agb_spectra
    public :: load_wmbasic_spectra
    public :: load_lsf_data
    
    ! --- Public API: Input (Runtime) ---
    public :: load_tabular_sfh

    ! --- Public API: Output ---
    public :: write_csp_output_files
    public :: write_isochrone_cmd
    public :: write_binary_spectral_lib

    ! Output Mode Constants (legacy write_compsp flags)
    integer, parameter :: OUTPUT_MAGS = 1
    integer, parameter :: OUTPUT_SPEC = 2
    integer, parameter :: OUTPUT_BOTH = 3
    integer, parameter :: OUTPUT_INDX = 4

    ! Dust emission table constants
    real(WP), parameter :: QPAH_ARR_DL07(7) = [ &
        0.47_wp, 1.12_wp, 1.77_wp, 2.50_wp, 3.19_wp, 3.90_wp, 4.58_wp &
    ]
    real(WP), parameter :: UMIN_ARR_DL07(22) = [ &
        0.1_wp,  0.15_wp, 0.2_wp,  0.3_wp,  0.4_wp, 0.5_wp, &
        0.7_wp,  0.8_wp,  1.0_wp,  1.2_wp,  1.5_wp, 2.0_wp, &
        2.5_wp,  3.0_wp,  4.0_wp,  5.0_wp,  7.0_wp, 8.0_wp, &
        12.0_wp, 15.0_wp, 20.0_wp, 25.0_wp &
    ]
    real(WP), parameter :: QPAH_ARR_THEMIS(11) = (100.0_wp / 2.2_wp) * [ &
        0.02_wp, 0.06_wp, 0.10_wp, 0.14_wp, 0.17_wp, 0.20_wp, &
        0.24_wp, 0.28_wp, 0.32_wp, 0.36_wp, 0.40_wp &
    ]
    real(WP), parameter :: UMIN_ARR_THEMIS(37) = [ &
        0.1_wp,  0.12_wp, 0.15_wp, 0.17_wp, 0.2_wp,  0.25_wp, &
        0.3_wp,  0.35_wp, 0.4_wp,  0.5_wp,  0.6_wp,  0.7_wp, &
        0.8_wp,  1.0_wp,  1.2_wp,  1.5_wp,  1.7_wp,  2.0_wp, &
        2.5_wp,  3.0_wp,  3.5_wp,  4.0_wp,  5.0_wp,  6.0_wp, &
        7.0_wp,  8.0_wp,  10.0_wp, 12.0_wp, 15.0_wp, 17.0_wp, &
        20.0_wp, 25.0_wp, 30.0_wp, 35.0_wp, 40.0_wp, 50.0_wp, &
        80.0_wp &
    ]

    ! Filter hardcoded constants (legacy indices and wavelengths for normalization)
    integer, parameter :: FILTER_IND_IR(14)   = [53, 54, 55, 56, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104]
    integer, parameter :: FILTER_IND_MIPS(3)  = [90, 91, 92]
    real(WP), parameter :: FILTER_LAM_IR(14)  = 1.0e4_wp * [ &
        3.550_wp, 4.493_wp, 5.731_wp, 7.872_wp, 70.0_wp, 100.0_wp, &
        160.0_wp, 250.0_wp, 350.0_wp, 500.0_wp, 12.0_wp, 25.0_wp, &
        60.0_wp, 100.0_wp &
    ]
    real(WP), parameter :: FILTER_LAM_MIPS(3) = 1.0e4_wp * [23.68_wp, 71.42_wp, 155.9_wp]

    ! Other constants
    integer, parameter :: MAX_POINTS = 50000
    integer, parameter :: N_SED_POINTS = 1221 ! Matches ntlam in legacy
    integer, parameter :: NLAMWR = 1963
    integer, parameter :: NSPEC_PAGB = 9281
    integer, parameter :: NSPEC_AGB = 6146
    integer, parameter :: NSPEC_ARINGER = 9032
    integer, parameter :: NZWMB = 12
    integer, parameter :: NSPEC_WMB = 5508

contains

    ! ========================================================================
    ! ISOCHRONE I/O
    ! ========================================================================

    !> @brief Reads the isochrone database for a specific metallicity index.
    !>
    !> @details
    !> Parses the ASCII isochrone files (Padova, MIST, BaSTI, etc.) for a single metallicity 
    !> slice. The routine identifies tracks (separated by '#') and steps within tracks 
    !> to populate the global isochrone arrays.
    !>
    !> **Optimization Note:**
    !> Uses internal string reading ("internal files") to avoid expensive `BACKSPACE` operations 
    !> found in the legacy code. It reads an entire line into a buffer and then parses 
    !> based on content, robustly handling headers and comments.
    !>
    !> **MIST vs Standard Format:**
    !> - **Standard (Padova/BaSTI):** Expects 8 columns. `lmdot` is not provided (set to -99.0).
    !> - **MIST:** Expects 9 columns. The 9th column is explicitly read as `lmdot`.
    !>
    !> @param[inout] ctx       The FSPS context. Populates `*_isoc` arrays at index `z_idx`.
    !> @param[in]    isoc_type String identifier (e.g., 'mist', 'pdva', 'bsti').
    !> @param[in]    z_idx     The metallicity index to read (1..nz).
    subroutine read_isochrone_database(ctx, isoc_type, z_idx)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: isoc_type
        integer, intent(in) :: z_idx

        integer :: u_file, io_stat, i_line
        integer :: current_track, current_step
        character(len=1024) :: file_path
        character(len=4096) :: line_buf
        logical :: is_mist
        
        ! Temporary variables for reading
        real(WP) :: r_logage, r_mini, r_mact, r_logl, r_logt, r_logg
        real(WP) :: r_ffco, r_phase, r_lmdot

        ! 1. Setup
        call get_isochrone_filename(ctx, isoc_type, z_idx, file_path)
        is_mist = (trim(isoc_type) == 'mist')
        
        if (VERBOSE > 1) print *, 'Reading Isochrone: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)

        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error: Cannot open isochrone file: ", A)') trim(file_path)
            error stop
        end if

        ! 2. Parse Loop
        current_track = 0
        current_step  = 0
        
        line_loop: do i_line = 1, NLINES
            
            ! Read entire line into buffer (Fast)
            read(u_file, '(A)', iostat=io_stat) line_buf
            if (io_stat == iostat_end) exit line_loop
            if (io_stat /= 0) then
                write(error_unit, '("[FSPS_IO] Error reading line ", I0, " of ", A)') i_line, trim(file_path)
                error stop
            end if

            ! Handle indented comments/headers
            line_buf = adjustl(line_buf)

            ! Skip empty lines
            if (len_trim(line_buf) == 0) cycle

            ! Check for Track Delimiter ('#')
            if (line_buf(1:1) == '#') then
                ! Just reset the step counter. Do NOT increment track yet.
                ! We wait for the first data line to confirm the track started.
                current_step = 0
                cycle 
            end if

            ! Parse Data Line
            ! ---------------
            ! We attempt to read FIRST. Only if successful do we update counters.
            if (is_mist) then
                ! MIST Format (9 columns)
                read(line_buf, *, iostat=io_stat) &
                    r_logage, r_mini, r_mact, r_logl, r_logt, r_logg, &
                    r_ffco, r_phase, r_lmdot
            else
                ! Standard Format (8 columns)
                read(line_buf, *, iostat=io_stat) &
                    r_logage, r_mini, r_mact, r_logl, r_logt, r_logg, &
                    r_ffco, r_phase
            end if

            ! If the line doesn't parse as numeric data, treat it as a delimiter/header.
            if (io_stat /= 0) then
                current_step = 0
                cycle line_loop
            end if

            ! Counters updated here, after verifying we have valid data.
            ! If current_step is 0, this is the first data line after a delimiter.
            if (current_step == 0) then
                current_track = current_track + 1
                
                if (current_track > ctx%state%nt) then
                      write(error_unit, '("[FSPS_IO] Error: Tracks > NT in ", A)') trim(file_path)
                      error stop
                end if
            end if
            
            current_step = current_step + 1
            if (current_step > NM) then
                write(error_unit, '("[FSPS_IO] Error: Steps > NM in ", A)') trim(file_path)
                error stop
            end if

            if (is_mist) then
                ! Store directly
                ctx%state%lmdot_isoc(z_idx, current_track, current_step) = r_lmdot
            else
                ! Default unavailable value
                ctx%state%lmdot_isoc(z_idx, current_track, current_step) = -99.0_wp
            end if

            ! Store Common Properties
            ctx%state%mini_isoc(z_idx, current_track, current_step) = r_mini
            ctx%state%mact_isoc(z_idx, current_track, current_step) = r_mact
            ctx%state%logl_isoc(z_idx, current_track, current_step) = r_logl
            ctx%state%logt_isoc(z_idx, current_track, current_step) = r_logt
            ctx%state%logg_isoc(z_idx, current_track, current_step) = r_logg
            ctx%state%ffco_isoc(z_idx, current_track, current_step) = r_ffco
            ctx%state%phase_isoc(z_idx, current_track, current_step) = r_phase

            ! Store Age (only needed on the first step of the track effectively)
            if (current_step == 1) then
                ctx%state%timestep_isoc(z_idx, current_track) = r_logage
            end if
            
            ! Update count
            ctx%state%nmass_isoc(z_idx, current_track) = current_step

        end do line_loop

        close(u_file)

        ! Validation
        if (current_track /= ctx%state%nt) then
            write(error_unit, '("[FSPS_IO] Warning: Expected ", I0, " tracks, found ", I0, " in ", A)') &
                ctx%state%nt, current_track, trim(file_path)
        end if

    end subroutine read_isochrone_database

    ! ========================================================================
    ! SPECTRAL LIBRARY I/O
    ! ========================================================================

    !> @brief Reads a pre-compiled binary spectral library for a specific metallicity.
    !>
    !> @details
    !> Loads a single metallicity slab of the stellar spectral library (e.g., BaSeL, MILES) 
    !> from disk. These files are pre-compiled binaries created by `write_binary_spectral_lib` 
    !> (formerly the `spec_bin` program).
    !>
    !> **Data Format:**
    !> - The files are **Unformatted Direct Access**.
    !> - Data is strictly stored as **32-bit floats** (`real32`) to minimize disk usage and 
    !>   maintain backward compatibility with legacy FSPS datasets.
    !> - Structure: A single record containing the full `(Nspec, NlogT, NlogG)` cube.
    !>
    !> **Precision Handling:**
    !> This routine reads the raw data into a temporary 32-bit buffer and explicitly promotes 
    !> the values to the working precision (`WP`) of the context (typically 64-bit) during 
    !> the transfer to `array_out`.
    !>
    !> @param[in]  ctx       The FSPS context (provides paths and dimensions).
    !> @param[in]  spec_type String identifier for the library (e.g., 'miles', 'basel').
    !> @param[in]  z_idx     The metallicity index to read (1..nzinit).
    !> @param[out] array_out The output spectral cube with dimensions (nspec, NDIM_LOGT, NDIM_LOGG).
    subroutine read_spectral_binary(ctx, spec_type, z_idx, array_out)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        integer, intent(in) :: z_idx
        real(WP), dimension(:,:,:), intent(out), contiguous :: array_out

        integer :: u_file, io_stat, rec_len, file_unit_size
        character(len=1024) :: file_path
        
        ! Temporary buffer matching the LEGACY DISK FORMAT (32-bit reals)
        real(real32), allocatable, dimension(:,:,:) :: buffer_32

        ! 1. Setup
        ! --------
        call get_spectral_filename(ctx, spec_type, z_idx, file_path)

        if (VERBOSE > 1) print *, 'Reading Spectrum: ', trim(file_path)

        ! Get the unit size (in bits) for file storage
        ! usually 8 (bytes) or 32 (words)
        file_unit_size = file_storage_size / 8 
        
        ! Adjust record length calculation
        rec_len = (ctx%state%nspec * NDIM_LOGT * NDIM_LOGG * 4) / file_unit_size

        allocate(buffer_32(ctx%state%nspec, NDIM_LOGT, NDIM_LOGG))

        ! 2. Read Binary
        ! --------------
        open(newunit=u_file, file=trim(file_path), status='old', &
             access='direct', recl=rec_len, form='unformatted', &
             action='read', iostat=io_stat)

        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error: Cannot open binary spec file: ", A)') trim(file_path)
            error stop
        end if

        ! Read Record 1 (The entire cube for this metallicity)
        read(u_file, rec=1, iostat=io_stat) buffer_32

        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error reading binary record in ", A)') trim(file_path)
            error stop
        end if

        close(u_file)

        ! 3. Promotion to WP
        ! ------------------
        ! This explicit cast ensures we handle the 32-bit -> 64-bit conversion safely
        array_out = real(buffer_32, kind=WP)
        
        deallocate(buffer_32)

    end subroutine read_spectral_binary

    !> @brief Reads all data required for BPASS models (wavelengths, mass, spectra).
    !>
    !> @details
    !> BPASS (Binary Population and Spectral Synthesis) models are handled distinctly 
    !> from the standard isochrone/spectral library architecture. This routine loads 
    !> the specific BPASS datasets:
    !>
    !> 1. **Wavelengths:** Reads `bpass.lambda` into `spec_lambda`.
    !> 2. **Mass/Age:** Reads `bpass.mass` into `time_full` and `bpass_mass_ssp`.
    !> 3. **Spectra:** Reads the compiled binary SSP cube `bpass_v2.2_salpeter100.ssp.bin` 
    !>    into `bpass_spec_ssp`.
    !>
    !> **Format Note:**
    !> The BPASS binary file stores data as 64-bit doubles (matching `real(WP)`), 
    !> unlike the standard spectral libraries which use 32-bit floats.
    !> The data cube dimensions are (NSPEC, NT, NZ).
    !>
    !> @param[inout] ctx The FSPS context. Populates `spec_lambda`, `time_full`, 
    !>                   `bpass_mass_ssp`, and `bpass_spec_ssp`.
    subroutine read_bpass_data(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i, rec_len
        character(len=1024) :: file_path

        if (VERBOSE > 0) print *, 'Loading BPASS Data...'

        file_path = trim(ctx%sps_home) // '/data/isochrones/BPASS/bpass.lambda'
        
        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        ! [cite_start]Legacy code loop [cite: 61-62]
        do i = 1, ctx%state%nspec
            read(u_file, *, iostat=io_stat) ctx%state%spec_lambda(i)
            if (io_stat /= 0) call handle_read_error(file_path, i)
        end do
        close(u_file)

        ! Update frequency array immediately
        where (ctx%state%spec_lambda > SAFE_FLOOR)
            ctx%state%spec_nu = C_LIGHT / ctx%state%spec_lambda
        elsewhere
            ctx%state%spec_nu = 0.0_wp
        end where

        file_path = trim(ctx%sps_home) // '/data/isochrones/BPASS/bpass.mass'
        
        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        do i = 1, ctx%state%nt
            read(u_file, *, iostat=io_stat) ctx%state%time_full(i), &
                                            ctx%state%bpass_mass_ssp(i, :)
            if (io_stat /= 0) call handle_read_error(file_path, i)
        end do
        close(u_file)

        file_path = trim(ctx%sps_home) // &
                    '/data/isochrones/BPASS/bpass_v2.2_salpeter100.ssp.bin'

        rec_len = ctx%state%nspec * ctx%state%nt * ctx%state%nz * 8

        open(newunit=u_file, file=trim(file_path), status='old', &
             access='direct', recl=rec_len, form='unformatted', &
             action='read', iostat=io_stat)
        
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening BPASS binary: ", A)') trim(file_path)
            error stop
        end if

        read(u_file, rec=1, iostat=io_stat) ctx%state%bpass_spec_ssp
        
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error reading BPASS binary record.")')
            error stop
        end if
        close(u_file)

    end subroutine read_bpass_data

    ! ========================================================================
    ! METADATA LOADERS (LEGENDS & GRIDS)
    ! ========================================================================

    !> @brief Loads the metallicity legend for Isochrones or Spectra.
    !>
    !> @details
    !> Reads the `zlegend.dat` file which defines the available metallicity grid points 
    !> for the selected isochrone set or spectral library. This file maps an integer 
    !> index (1..N) to a specific metallicity value (Z).
    !>
    !> **MIST Special Handling:**
    !> MIST isochrones use a string-based format in `zlegend.dat` (e.g., "m1.50", "p0.25") 
    !> representing [Fe/H]. This routine parses these strings:
    !> - 'm' denotes minus (negative [Fe/H]).
    !> - 'p' denotes plus (positive [Fe/H]).
    !> - Values are converted from log(Z/Zsol) to absolute Z using `ctx%state%zsol`.
    !>
    !> **Standard Handling:**
    !> For all other libraries (Padova, BaSTI, MILES, etc.), the file is assumed to contain 
    !> a simple column of floating-point Z values.
    !>
    !> @param[inout] ctx      The FSPS context. Populates `zlegend` or `zlegendinit`.
    !> @param[in]    lib_type String identifier for the library (e.g., 'mist', 'miles').
    !> @param[in]    is_spec  Logical flag:
    !>                        - `.true.`: Loading spectral library legend (populates `zlegendinit`).
    !>                        - `.false.`: Loading isochrone legend (populates `zlegend`).
    subroutine load_zlegend_file(ctx, lib_type, is_spec)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: lib_type
        logical, intent(in) :: is_spec

        integer :: u_file, io_stat, i, n_max
        character(len=1024) :: file_path
        character(len=5) :: mist_str
        real(WP) :: val_read
        
        ! Determine target array and size
        if (is_spec) then
            n_max = ctx%state%nzinit
        else
            n_max = ctx%state%nz
        end if

        ! 1. Construct Path
        call get_zlegend_filename(ctx, lib_type, is_spec, file_path)
        
        if (VERBOSE > 1) print *, 'Reading Z-Legend: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
             
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening zlegend: ", A)') trim(file_path)
            error stop
        end if

        ! 2. Parse Loop
        do i = 1, n_max
            
            if (.not. is_spec .and. trim(lib_type) == 'mist') then
                ! -- MIST SPECIAL CASE --
                ! MIST zlegend.dat contains strings: "m1.50", "p0.25"
                ! m = minus, p = plus. Values are log(Z/Zsol).
                read(u_file, '(A5)', iostat=io_stat) mist_str
                
                if (io_stat /= 0) call handle_read_error(file_path, i)

                read(mist_str(2:5), '(F4.2)', iostat=io_stat) val_read
                if (io_stat /= 0) call handle_read_error(file_path, i)

                if (mist_str(1:1) == 'm') then
                    val_read = 10.0_wp**(-val_read) * ctx%state%zsol
                else
                    val_read = 10.0_wp**(val_read) * ctx%state%zsol
                end if

            else
                ! -- STANDARD CASE --
                ! Just a float value
                read(u_file, *, iostat=io_stat) val_read
                if (io_stat /= 0) call handle_read_error(file_path, i)
            end if

            ! 3. Store Result
            if (is_spec) then
                ctx%state%zlegendinit(i) = val_read
            else
                ctx%state%zlegend(i) = val_read
            end if

        end do

        close(u_file)

        if (env_flag_true('FSPS_DEBUG_IO')) then
            if (is_spec) then
                call debug_stats_1d('[DEBUG_IO] zlegendinit', ctx%state%zlegendinit)
            else
                call debug_stats_1d('[DEBUG_IO] zlegend', ctx%state%zlegend)
            end if
        end if

    end subroutine load_zlegend_file

    !> @brief Loads the master wavelength grid for the active spectral library.
    !>
    !> @details
    !> Reads the `*.lambda` file associated with the chosen spectral library (e.g., 
    !> `miles.lambda` or `basel.lambda`). This grid defines the common wavelength 
    !> points (in Angstroms) for all subsequent spectral operations.
    !>
    !> **Optimization:**
    !> Immediately calculates and caches the corresponding frequency array (`spec_nu`) 
    !> to avoid repeated divisions by lambda during flux unit conversions later in the pipeline.
    !>
    !> @param[inout] ctx       The FSPS context. Populates `ctx%state%spec_lambda` and `spec_nu`.
    !> @param[in]    spec_type String identifier for the library (e.g., 'miles').
    subroutine load_wavelength_grid(ctx, spec_type)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: spec_type
        
        integer :: u_file, io_stat, i
        character(len=1024) :: file_path
        
        call get_lambda_filename(ctx, spec_type, file_path)
        
        if (VERBOSE > 1) print *, 'Reading Wavelengths: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
             
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening lambda file: ", A)') trim(file_path)
            error stop
        end if

        do i = 1, ctx%state%nspec
            read(u_file, *, iostat=io_stat) ctx%state%spec_lambda(i)
            if (io_stat /= 0) call handle_read_error(file_path, i)
        end do

        close(u_file)

        if (env_flag_true('FSPS_DEBUG_IO')) then
            call debug_stats_1d('[DEBUG_IO] spec_lambda', ctx%state%spec_lambda)
            call debug_stats_1d('[DEBUG_IO] spec_nu', ctx%state%spec_nu)
        end if
        
        ! Vectorized operation
        where (ctx%state%spec_lambda > SAFE_FLOOR)
            ctx%state%spec_nu = C_LIGHT / ctx%state%spec_lambda
        elsewhere
            ctx%state%spec_nu = 0.0_wp
        end where

    end subroutine load_wavelength_grid

    !> @brief Loads the spectral resolution (FWHM) for the current library.
    !>
    !> @details
    !> Reads the `*.res` file associated with the chosen spectral library (e.g., `miles.res`, 
    !> `basel.res`). This file defines the instrumental resolution (FWHM in Angstroms) 
    !> corresponding to each point on the master wavelength grid.
    !>
    !> **Purpose:**
    !> This array is used by smoothing routines (velocity dispersion or LSF application) 
    !> to determine the intrinsic broadening already present in the templates.
    !>
    !> @param[inout] ctx       The FSPS context. Populates `ctx%state%spec_res`.
    !> @param[in]    spec_type String identifier for the library (e.g., 'miles').
    subroutine load_spectral_resolution(ctx, spec_type)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: spec_type
        
        integer :: u_file, io_stat, i
        character(len=1024) :: file_path
        
        call get_res_filename(ctx, spec_type, file_path)
        
        if (VERBOSE > 1) print *, 'Reading Resolution: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
             
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening res file: ", A)') trim(file_path)
            error stop
        end if

        do i = 1, ctx%state%nspec
            read(u_file, *, iostat=io_stat) ctx%state%spec_res(i)
            if (io_stat /= 0) call handle_read_error(file_path, i)
        end do

        close(u_file)

        if (env_flag_true('FSPS_DEBUG_IO')) then
            call debug_stats_1d('[DEBUG_IO] spec_res', ctx%state%spec_res)
        end if

    end subroutine load_spectral_resolution

    ! ========================================================================
    ! NEBULAR & DUST I/O (Load + Interpolate)
    ! ========================================================================

    !> @brief Loads pre-computed Nebular Continuum and Line grids.
    !>
    !> @details
    !> Reads the Cloudy-based nebular emission models stored in `data/nebular/`. 
    !> The grids are defined over 3 dimensions: Metallicity (Z), Age, and Ionization Parameter (U), 
    !> with sizes defined by `NEBNZ`, `NEBNAGE`, and `NEBNIP`.
    !>
    !> **File Selection:**
    !> - If `use_cloudy=.true.`: Loads `ZAU_WD_...` (With Dust).
    !> - If `use_cloudy=.false.`: Loads `ZAU_ND_...` (No Dust).
    !>
    !> **Processing Steps:**
    !> 1. **Continuum:** Reads the high-res nebular continuum. 
    !>    - Interpolates from the native Cloudy wavelength grid (`NLAM_NEBCONT`) onto the 
    !>      master FSPS grid (`spec_lambda`). 
    !>    - Interpolation is performed on log10(flux) to preserve dynamic range.
    !> 2. **Lines:** Reads the specific line luminosities for major emission lines.
    !> 3. **Metadata:** Stores the grid axes (LogZ, LogU, Age).
    !>
    !> **Post-Processing:**
    !> - Converts the `nebem_age` grid to Log10(Age).
    !> - Converts line luminosities to Log10 units.
    !>
    !> @param[inout] ctx        The FSPS context. Populates `nebem_cont` and `nebem_line`.
    !> @param[in]    isoc_type  Isochrone identifier (e.g. 'mist') to match the stellar grid.
    !> @param[in]    use_cloudy Logical flag to select dust-attenuated or dust-free models.
    subroutine load_nebular_grid(ctx, isoc_type, use_cloudy)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: isoc_type
        logical, intent(in) :: use_cloudy

        integer :: u_file, io_stat, i, j, k, i_spec
        character(len=1024) :: file_path_cont, file_path_lines
        character(len=32) :: suffix
        
        ! Temporary arrays for the raw grid reading
        real(WP), allocatable :: raw_lam(:), raw_spec(:), raw_spec_log(:)
        real(WP) :: val_logz, val_age, val_logu
        real(WP), parameter :: NEBULAR_FLOOR = 10.0_wp**(-95.0_wp)

        ! 1. Construct File Paths
        ! -----------------------
        ! "WD" = With Dust, "ND" = No Dust
        if (use_cloudy) then
            suffix = '_WD_'
        else
            suffix = '_ND_'
        end if
        
        file_path_cont = trim(ctx%sps_home) // '/data/nebular/ZAU' // trim(suffix) // &
                         trim(isoc_type) // '.cont'
        
        file_path_lines = trim(ctx%sps_home) // '/data/nebular/ZAU' // trim(suffix) // &
                          trim(isoc_type) // '.lines'

        ! 2. Read Continuum (.cont)
        ! -------------------------
        if (VERBOSE > 1) print *, 'Reading Nebular Cont: ', trim(file_path_cont)

        open(newunit=u_file, file=trim(file_path_cont), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path_cont, 0)

        ! Header skip
        read(u_file, *) 
        
        ! Allocate using the Constant from fsps_constants
        allocate(raw_lam(NLAM_NEBCONT))
        allocate(raw_spec(NLAM_NEBCONT))
        allocate(raw_spec_log(NLAM_NEBCONT))

        read(u_file, *) raw_lam

        ! Nested Loops: Z -> Age -> U
        ! Uses constants NEBNZ, NEBNAGE, NEBNIP directly
        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    
                    ! Read Metadata (LogZ, Age, LogU) and Spectrum
                    read(u_file, *, iostat=io_stat) val_logz, val_age, val_logu
                    if (io_stat /= 0) call handle_read_error(file_path_cont, i)
                    
                    read(u_file, *, iostat=io_stat) raw_spec
                    if (io_stat /= 0) call handle_read_error(file_path_cont, i)

                    ! Store Metadata
                    ctx%state%nebem_logz(i) = val_logz
                    ctx%state%nebem_age(j)  = val_age  ! Will be logged later if needed
                    ctx%state%nebem_logu(k) = val_logu

                    ! INTERPOLATION
                    ! Legacy logic: raw_spec can be 0.0, set floor 1e-95, log10, interpolate.
                    raw_spec_log = log10(raw_spec + NEBULAR_FLOOR)
                    do i_spec = 1, ctx%state%nspec
                        ctx%state%nebem_cont(i_spec, i, j, k) = interpolate_linear( &
                            raw_lam, &
                            raw_spec_log, &
                            ctx%state%spec_lambda(i_spec))
                    end do

                end do
            end do
        end do
        close(u_file)
        
        if (env_flag_true('FSPS_DEBUG_IO')) then
            call debug_stats_1d('[DEBUG_IO] nebem_logz', ctx%state%nebem_logz)
            call debug_stats_1d('[DEBUG_IO] nebem_logu', ctx%state%nebem_logu)
            call debug_stats_1d('[DEBUG_IO] nebem_age(log)', ctx%state%nebem_age)
            write(*,'(A,1x,ES13.5)') '[DEBUG_IO] nebem_cont(1,1,1,1)=', ctx%state%nebem_cont(1,1,1,1)
            write(*,'(A,1x,ES13.5)') '[DEBUG_IO] nebem_line(1,1,1,1)=', ctx%state%nebem_line(1,1,1,1)
        end if

        deallocate(raw_lam, raw_spec, raw_spec_log)


        ! 3. Read Lines (.lines)
        ! ----------------------
        ! Only available for certain isochrones in legacy code check, 
        ! but we assume caller handles that logic.
        if (VERBOSE > 1) print *, 'Reading Nebular Lines: ', trim(file_path_lines)

        open(newunit=u_file, file=trim(file_path_lines), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path_lines, 0)

        read(u_file, *) ! Header
        
        ! Read Line Wavelengths
        read(u_file, *) ctx%state%nebem_line_pos

        ! Nested Loops
        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    ! Metadata (redundant overwrite, but consistent with file structure)
                    read(u_file, *) val_logz, val_age, val_logu

                    ! Store the age, as legacy code did
                    ctx%state%nebem_age(j) = val_age
                    
                    ! Read Line Luminosities directly into state
                    read(u_file, *) ctx%state%nebem_line(:, i, j, k)
                end do
            end do
        end do
        close(u_file)

        ! Post-processing (matches legacy setup)
        ! Convert Age array to Log10
        ctx%state%nebem_age = log10(ctx%state%nebem_age)
        ! Convert Line Luminosities to Log10
        ctx%state%nebem_line = log10(ctx%state%nebem_line + NEBULAR_FLOOR)

    end subroutine load_nebular_grid

    !> @brief Loads Draine & Li (2007) or THEMIS dust emission templates.
    !>
    !> @details
    !> Initializes the dust emission physics by loading the appropriate template library 
    !> based on `dust_type`. This routine performs several key setup tasks:
    !>
    !> 1. **Model Configuration:** Sets the dimensions (`nqpah`, `numin`) and populates 
    !>    the available parameter grids (`qpaharr`, `uminarr`) for the chosen model.
    !> 2. **File Loading:** Iterates through files corresponding to different PAH mass 
    !>    fractions (e.g., `_MW3.1_X0.dat`).
    !> 3. **Unit Conversion:** Converts input wavelengths from microns to Angstroms.
    !> 4. **Interpolation:** Interpolates the raw emission SEDs onto the master 
    !>    `spec_lambda` grid for use during synthesis.
    !>
    !> **Supported Models:**
    !> - **'DL07':** Draine & Li (2007). 7 PAH fractions, 22 U_min values.
    !> - **'THEMIS':** Jones et al. (2013, 2017). 11 PAH fractions, 35+ U_min values.
    !>
    !> @param[inout] ctx       The FSPS context. Populates `dustem2_dustem`, `qpaharr`, and `uminarr`.
    !> @param[in]    dust_type String identifier for the dust model ('DL07' or 'THEMIS').
    subroutine load_dust_emission_table(ctx, dust_type)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: dust_type

        integer :: u_file, io_stat, k, i_spec, i_lam
        character(len=1024) :: file_path
        character(len=1) :: qpah_char
        real(WP), allocatable :: raw_lam(:), raw_spec(:,:)
        integer :: ndim_dust, numin_dust, nqpah
        integer :: start_idx

        ! Setup dimensions based on type (logic moved to fsps_initialization)
        if (trim(dust_type) == 'THEMIS') then
            ndim_dust  = 576
            numin_dust = 37
            nqpah      = 11
            ctx%state%qpaharr = QPAH_ARR_THEMIS
            ctx%state%uminarr = UMIN_ARR_THEMIS
        else ! DL07
            ndim_dust  = 1001
            numin_dust = 22
            nqpah      = 7
            ctx%state%qpaharr = QPAH_ARR_DL07
            ctx%state%uminarr = UMIN_ARR_DL07
        end if

        allocate(raw_lam(ndim_dust))
        allocate(raw_spec(ndim_dust, numin_dust * 2)) ! *2 because file has multiple cols per Umin

        ! Loop over PAH fractions (qpah)
        ! The files are numbered 0..N-1 in the filename char
        do k = 1, nqpah
            
            ! Construct Filename "MW3.1_X0.dat"
            if (k - 1 == 10) then
                file_path = trim(ctx%sps_home) // '/data/dust/dustem/' // &
                            trim(dust_type) // '_MW3.1_100.dat'
            else
                write(qpah_char, '(I1)') k - 1
                file_path = trim(ctx%sps_home) // '/data/dust/dustem/' // &
                            trim(dust_type) // '_MW3.1_' // qpah_char // '0.dat'
            end if

            if (VERBOSE > 1) print *, 'Reading Dust Em: ', trim(file_path)

            open(newunit=u_file, file=trim(file_path), status='old', &
                 action='read', iostat=io_stat)
            if (io_stat /= 0) call handle_read_error(file_path, 0)

            ! Burn Header (2 lines)
            read(u_file, *)
            read(u_file, *)

            ! Read Data Table
            ! Columns: Lambda, [Spectrum for Umin 1..N]...
            do i_spec = 1, ndim_dust
                read(u_file, *) raw_lam(i_spec), raw_spec(i_spec, :)
            end do
            close(u_file)

            ! Conversion: microns -> Angstroms
            raw_lam = raw_lam * 1.0e4_wp
            ! Interpolate every U_min column for this Q_PAH in linear space.
            ! Match legacy behavior: only populate wavelengths >= 1 micron
            ! and allow linear extrapolation beyond the native dust table range.
            start_idx = max(find_interval(ctx%state%spec_lambda, 1.0e4_wp), 1)
            do i_spec = 1, numin_dust * 2
                ctx%state%dustem2_dustem(:, k, i_spec) = 0.0_wp
                do i_lam = start_idx, ctx%state%nspec
                    ctx%state%dustem2_dustem(i_lam, k, i_spec) = interpolate_linear( &
                        raw_lam, raw_spec(:, i_spec), ctx%state%spec_lambda(i_lam))
                end do
            end do

        end do

        if (env_flag_true('FSPS_DEBUG_IO')) then
            write(*,'(A,1x,I0,1x,I0)') '[DEBUG_IO] dustem dims (qpah,numin*2)=', nqpah, numin_dust*2
            call debug_stats_1d('[DEBUG_IO] dustem raw_lam (um->A)', raw_lam)
            write(*,'(A,1x,ES13.5)') '[DEBUG_IO] dustem2(1,1,1)=', ctx%state%dustem2_dustem(1,1,1)
        end if

        deallocate(raw_lam, raw_spec)

    end subroutine load_dust_emission_table

    !> @brief Loads Nenkova et al. (2008) AGN dust torus models and interpolates onto the master grid.
    !>
    !> @details
    !> This routine reads the theoretical AGN torus emission templates from the `data/dust/` 
    !> directory. It populates the AGN dust state within the context, which is used to add 
    !> the obscured AGN component to the integrated SED when `fagn > 0`.
    !>
    !> **Processing Steps:**
    !> 1. **File Parsing:** Opens the hardcoded Nenkova08 torus model file and skips the 
    !>    3-line metadata header.
    !> 2. **Optical Depth Grid:** Reads the available optical depth values (`agndust_tau`) 
    !>    into the context state.
    !> 3. **Raw Spectral Load:** Reads the raw AGN spectra (Wavelength vs Flux density) 
    !>    for each model in the library.
    !> 4. **Wavelength Clipping:** Identifies the overlap interval between the AGN model 
    !>    coverage and the master FSPS wavelength grid (`spec_lambda`) to avoid out-of-bounds 
    !>    interpolation.
    !> 5. **Log-Space Interpolation:** Interpolates the raw spectra onto the master grid. 
    !>    Interpolation is performed in `log10(Flux + SAFE_FLOOR)` to maintain numerical 
    !>    stability across high-contrast emission features and prevent negative flux artifacts.
    !>
    !> **Note on Units:**
    !> The resulting `agndust_spec` is stored in the context's working precision (`WP`) 
    !> and is ready for scaling by the `fagn` parameter during CSP synthesis.
    !>
    !> @param[inout] ctx The FSPS context containing the system state and path configuration. 
    !>                   Populates `ctx%state%agndust_tau` and `ctx%state%agndust_spec`.
    subroutine load_agn_dust_models(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        
        integer :: u_file, io_stat, i, i1, i2
        character(len=1024) :: file_path
        real(WP), dimension(NAGNDUST_SPEC) :: agn_lam
        real(WP), dimension(NAGNDUST_SPEC, NAGNDUST) :: agn_raw
        
        file_path = trim(ctx%sps_home) // '/data/dust/Nenkova08_y010_torusg_n10_q2.0.dat'
        
        if (VERBOSE > 1) print *, 'Reading AGN Dust: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        ! Burn header (3 lines)
        do i = 1, 3
            read(u_file, *)
        end do

        ! Read Optical Depths (agndust_tau is in context state)
        read(u_file, *) ctx%state%agndust_tau

        ! Read Spectra
        do i = 1, NAGNDUST_SPEC
            read(u_file, *) agn_lam(i), agn_raw(i, :)
        end do
        close(u_file)

        ! Interpolate onto master grid
        ! Models are usually defined over a subset of the FSPS range
        i1 = max(find_interval(ctx%state%spec_lambda, agn_lam(1)), 1)
        i2 = max(find_interval(ctx%state%spec_lambda, agn_lam(NAGNDUST_SPEC)), 1)

        ctx%state%agndust_spec = 0.0_wp

        do i = 1, NAGNDUST
            ctx%state%agndust_spec(i1:i2, i) = 10.0_wp**interpolate_linear( &
                log10(agn_lam), &
                log10(agn_raw(:, i) + SAFE_FLOOR), &
                log10(ctx%state%spec_lambda(i1:i2)) ) - SAFE_FLOOR
        end do

    end subroutine load_agn_dust_models

    ! ========================================================================
    ! PHOTOMETRY & CALIBRATION
    ! ========================================================================

    !> @brief Loads filter transmission curves from the database.
    !>
    !> @details
    !> Reads photometric filter transmission curves from `allfilters.dat` (or a custom file 
    !> specified by `alt_filename`). This file contains ~160+ filters stacked sequentially. 
    !> The routine iterates through the file, identifying each filter block, reading its 
    !> (lambda, transmission) pairs, and interpolating them onto the master `spec_lambda` grid.
    !>
    !> **Normalization:**
    !> Filters are normalized such that the integral of (Transmission / lambda) over all 
    !> wavelengths equals 1. This ensures standard magnitude integration logic (photon counting) 
    !> works correctly.
    !>
    !> **Effective Wavelength:**
    !> Calculates the "pivot wavelength" (effective lambda) for each filter, defined as:
    !> \f[ \lambda_{eff} = \sqrt{ \frac{\int \lambda \cdot T(\lambda) d\lambda}{\int (T(\lambda)/\lambda) d\lambda} } \f]
    !>
    !> @param[inout] ctx          The FSPS context. Populates `bands` and `filter_leff`.
    !> @param[in]    alt_filename Optional custom filter filename (relative to `data/`).
    subroutine load_filter_definitions(ctx, alt_filename)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in), optional :: alt_filename

        integer :: u_file, io_stat, i_filt, i_point, i_grid, n_points, start_idx, end_idx
        character(len=256) :: line_buf
        character(len=1024) :: file_path
        real(WP), allocatable :: raw_lam(:), raw_trans(:)
        real(WP) :: val_lam, val_trans, norm_fac
        
        ! 1. Determine File Path
        if (present(alt_filename) .and. len_trim(alt_filename) > 0) then
            file_path = trim(ctx%sps_home) // '/data/' // trim(alt_filename)
        else
            file_path = trim(ctx%sps_home) // '/data/allfilters.dat'
        end if

        if (VERBOSE > 1) print *, 'Reading Filters: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        ! Header Skip
        read(u_file, *) 

        ! Allocate temporary buffers for raw filter shapes
        allocate(raw_lam(MAX_POINTS))
        allocate(raw_trans(MAX_POINTS))

        ! 2. Filter Loop
        ! --------------
        ! Iterate over the number of bands defined in the context (nbands)
        do i_filt = 1, ctx%state%nbands
            
            n_points = 0

            block_read: do i_point = 1, MAX_POINTS
                ! Read line as text first
                read(u_file, '(A)', iostat=io_stat) line_buf
                
                ! Check EOF
                if (io_stat < 0) then
                    if (n_points > 0 .and. i_filt == ctx%state%nbands) exit block_read
                    call handle_read_error(file_path, 0)
                end if

                ! Handle indented comments/headers
                line_buf = adjustl(line_buf)

                ! Skip blank lines (matches list-directed read behavior)
                if (len_trim(line_buf) == 0) cycle block_read
                
                ! Check for Comment/Header (Legacy filters start with #)
                if (line_buf(1:1) == '#') then
                    if (n_points > 0) then
                        ! We finished a block and hit the header of the next one
                        ! Backspace so the next i_filt loop can read this header/comment if needed
                        ! (Or essentially, we are just done with this filter)
                        backspace(u_file)
                        exit block_read
                    else
                        ! We are likely sitting on the header of the current filter
                        cycle block_read
                    end if
                end if

                ! Parse data from string
                read(line_buf, *, iostat=io_stat) val_lam, val_trans
                if (io_stat /= 0) then
                    ! If we fail to parse numbers here, it's a true error
                    call handle_read_error(file_path, i_point)
                end if

                if (n_points > 0) then
                    if (abs(val_lam - raw_lam(n_points)) < SAFE_FLOOR) then
                        ! Duplicate point, skip
                        cycle block_read
                    end if
                end if

                n_points = n_points + 1
                raw_lam(n_points)   = val_lam
                raw_trans(n_points) = max(val_trans, 0.0_wp)
            end do block_read

            ! 3. Interpolate onto Master Grid
            ! -------------------------------
            if (n_points > 1) then
                ! Initialize to 0.0 to ensure no "phantom" throughput outside defined range
                ctx%state%bands(:, i_filt) = 0.0_wp

                ! Only interpolate within the valid range of the raw filter data
                ! This preserves Legacy behavior and prevents linear extrapolation artifacts
                call find_bounds(ctx%state%spec_lambda, raw_lam(1), raw_lam(n_points), start_idx, end_idx)
                
                ! Restore strict legacy check (start_idx < end_idx)
                if (end_idx > start_idx) then
                    do i_grid = start_idx, end_idx
                        ctx%state%bands(i_grid, i_filt) = interpolate_linear( &
                            raw_lam(1:n_points), &
                            raw_trans(1:n_points), &
                            ctx%state%spec_lambda(i_grid))
                    end do
                end if
            end if

            ! 4a. Normalize (Integral T/lambda dlambda = 1)
            ! --------------------------------------------
            ! We integrate (T / lambda)
            ! Avoid division by zero with SAFE_FLOOR
            norm_fac = integrate_trapezoid_array( &
                ctx%state%spec_lambda, &
                ctx%state%bands(:, i_filt) / (ctx%state%spec_lambda + SAFE_FLOOR) )
            
            ! Normalize
            if (norm_fac > SAFE_FLOOR) then
                ctx%state%bands(:, i_filt) = ctx%state%bands(:, i_filt) / norm_fac
            else
                ctx%state%bands(:, i_filt) = 0.0_wp ! Dead filter
            end if

            ! Clamp to zero after normalization
            where (ctx%state%bands(:, i_filt) < 0.0_wp)
                ctx%state%bands(:, i_filt) = 0.0_wp
            end where

        end do

        close(u_file)
        deallocate(raw_lam, raw_trans)

        ! ! 4b. IR Filter Renormalization (Legacy Compatibility)
        ! ! ----------------------------------------------------
        ! ! Only apply if using the standard filter list (not custom)
        ! if (.not. present(alt_filename) .or. (present(alt_filename) .and. len_trim(alt_filename) == 0)) then
        !     call apply_legacy_filter_norm(ctx)
        ! end if

        ! 5. Calculate Effective Wavelengths (Pivot Lambda)
        ! -------------------------------------------------
        
        do i_filt = 1, ctx%state%nbands
             val_lam = integrate_trapezoid_array( &
                 ctx%state%spec_lambda, &
                 ctx%state%spec_lambda * ctx%state%bands(:, i_filt) )
             
             ! Calculate denominator (integral of T/lambda)
             ! Even though bands are normalized, re-calculating ensures consistency
             ! and matches legacy behavior for re-normalized bands.
             norm_fac = integrate_trapezoid_array( &
                 ctx%state%spec_lambda, &
                 ctx%state%bands(:, i_filt) / (ctx%state%spec_lambda + SAFE_FLOOR) )
             
             if (norm_fac > SAFE_FLOOR) then
                 ctx%state%filter_leff(i_filt) = sqrt(val_lam / norm_fac)
             else
                 ctx%state%filter_leff(i_filt) = 0.0_wp
             end if
        end do

    end subroutine load_filter_definitions

    !> @brief Applies frequency-dependent normalization to Infrared (IR) filter transmission curves.
    !>
    !> @details
    !> This routine overrides the default photon-counting normalization for specific IR instruments 
    !> to maintain consistency with legacy FSPS flux calibration and specific telescope standards. 
    !> It applies two distinct types of re-normalization:
    !>
    !> **1. Constant $\nu F_\nu$ Normalization:**
    !> Applied to IRAC, PACS, SPIRE, and IRAS filters. The transmission curve $T(\lambda)$ is 
    !> normalized such that the integral of $(T(\lambda) / \lambda) \cdot (\lambda / \lambda_0)^{-1}$ 
    !> equals unity. This matches the calibration assumption that the source spectrum is 
    !> flat in $\nu F_\nu$ (i.e., $F_\nu \propto \nu^{-1}$).
    !>
    !> **2. Blackbody $\beta=2$ Normalization:**
    !> Applied specifically to MIPS filters. The transmission is normalized against a 
    !> modified blackbody slope where the integral includes a $(\lambda / \lambda_0)^{-2}$ 
    !> weighting factor.
    !>
    !> **Numerical Implementation:**
    !> - Uses `integrate_trapezoid_array` to perform the integration over the master wavelength grid.
    !> - Implements defensive checks against `SAFE_FLOOR` to prevent division by zero for 
    !>   unloaded or "dead" filter indices.
    !> - Utilizes the `FILTER_LAM_IR` and `FILTER_LAM_MIPS` parameter arrays from `fsps_constants` 
    !>   as the reference pivot wavelengths ($\lambda_0$).
    !>
    !> @param[inout] ctx The FSPS context containing the populated `bands` array. Only filters 
    !>                   matching the indices in `FILTER_IND_IR` and `FILTER_IND_MIPS` are modified.
    subroutine apply_legacy_filter_norm(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        
        integer :: j, i_filt
        real(WP) :: d

        ! 1. Normalize nu*fnu = const (IRAC, PACS, SPIRE, IRAS)
        do j = 1, 14
            i_filt = FILTER_IND_IR(j)
            if (i_filt > ctx%state%nbands) cycle

            ! Integrate (spec_lambda / lam_0)^(-1) * T / lambda
            d = integrate_trapezoid_array(ctx%state%spec_lambda, &
                (ctx%state%spec_lambda / FILTER_LAM_IR(j))**(-1.0_wp) * &
                ctx%state%bands(:, i_filt) / (ctx%state%spec_lambda + SAFE_FLOOR))
            
            if (d > SAFE_FLOOR) then
                ctx%state%bands(:, i_filt) = ctx%state%bands(:, i_filt) / d
            end if
        end do

        ! 2. Normalize Blackbody beta=2 (MIPS)
        do j = 1, 3
            i_filt = FILTER_IND_MIPS(j)
            if (i_filt > ctx%state%nbands) cycle

            ! Integrate (spec_lambda / lam_0)^(-2) * T / lambda
            d = integrate_trapezoid_array(ctx%state%spec_lambda, &
                (ctx%state%spec_lambda / FILTER_LAM_MIPS(j))**(-2.0_wp) * &
                ctx%state%bands(:, i_filt) / (ctx%state%spec_lambda + SAFE_FLOOR))

            if (d > SAFE_FLOOR) then
                ctx%state%bands(:, i_filt) = ctx%state%bands(:, i_filt) / d
            end if
        end do

    end subroutine apply_legacy_filter_norm

    !> @brief Loads standard SEDs (Vega, Sun) for photometric calibration.
    !>
    !> @details
    !> Reads the reference spectra for Vega (`A0V_KURUCZ_92.SED`) and the Sun (`SUN_STScI.SED`) 
    !> from the data directory. These are used to establish the photometric zero points.
    !>
    !> **Processing Steps:**
    !> 1. **Vega:** Reads raw data, interpolates to `spec_lambda`. Converts units from 
    !>    F_lambda to F_nu (multiplying by lambda^2).
    !> 2. **Sun:** Reads raw data, interpolates to `spec_lambda`.
    !> 3. **Zero Points:** Calculates the AB magnitude of both Vega and the Sun in every 
    !>    loaded filter band. This pre-calculation allows for fast magnitude computation 
    !>    during population synthesis.
    !>
    !> **Prerequisite:**
    !> `load_filter_definitions` must be called BEFORE this routine, as it relies on 
    !> `ctx%state%nbands` and `ctx%state%bands` being populated.
    !>
    !> @param[inout] ctx The FSPS context. Populates `vega_spec`, `sun_spec`, `magvega`, `magsun`.
    subroutine load_standard_sed(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i
        character(len=1024) :: file_path
        real(WP), allocatable :: raw_lam(:), raw_flux(:)

        allocate(raw_lam(N_SED_POINTS))
        allocate(raw_flux(N_SED_POINTS))

        ! ! Add defensive check
        ! if (all(ctx%state%filter_leff < tiny(0.0_wp))) then
        !     write(error_unit, '("[FSPS_IO] Error: Filters must be loaded before Standard SEDs.")')
        !     error stop
        ! end if

        ! --------------------------------------------------------------------
        ! 1. VEGA (Kurucz Model)
        ! --------------------------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/spectra/A0V_KURUCZ_92.SED'
        if (VERBOSE > 1) print *, 'Reading Vega SED: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        read(u_file, *) ! Header

        do i = 1, N_SED_POINTS
            read(u_file, *) raw_lam(i), raw_flux(i)
        end do
        close(u_file)
        
        ctx%state%vega_spec = 10.0_wp**interpolate_linear( &
             log10(raw_lam), log10(raw_flux + SAFE_FLOOR), log10(ctx%state%spec_lambda) )
        
        ! Apply F_lambda -> F_nu scaling factor (lambda^2)
        ctx%state%vega_spec = ctx%state%vega_spec * (ctx%state%spec_lambda**2)

        ! Only clamp the upper tail (> maxval).
        ! Legacy code allows interpolation to extrapolate blueward (lower indices).
        where (ctx%state%spec_lambda > maxval(raw_lam))
            ctx%state%vega_spec = SAFE_FLOOR
        end where

        ! --------------------------------------------------------------------
        ! 2. SUN (STScI)
        ! --------------------------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/spectra/SUN_STScI.SED'
        if (VERBOSE > 1) print *, 'Reading Sun SED: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        do i = 1, N_SED_POINTS
            read(u_file, *) raw_lam(i), raw_flux(i)
        end do
        close(u_file)
        
        ctx%state%sun_spec = 10.0_wp**interpolate_linear( &
             log10(raw_lam), log10(raw_flux + SAFE_FLOOR), log10(ctx%state%spec_lambda) )

        ! ! Apply F_lambda -> F_nu scaling factor (lambda^2) to Sun spectrum
        ! ctx%state%sun_spec = ctx%state%sun_spec * (ctx%state%spec_lambda**2)

        ! Only clamp the upper tail.
        where (ctx%state%spec_lambda > maxval(raw_lam))
            ctx%state%sun_spec = SAFE_FLOOR
        end where


        ! --------------------------------------------------------------------
        ! 3. PRE-CALCULATE ZERO POINTS
        ! --------------------------------------------------------------------
        ! We can calculate the magnitudes of Vega and Sun in every filter now.
        ! This saves time during runtime.
        
        do i = 1, ctx%state%nbands
            ! Mag = -2.5 * log10( Integral(F_nu * T / lam) ) - 48.60
            
            ! Vega Zero Point
            ctx%state%magvega(i) = integrate_trapezoid_array( &
                ctx%state%spec_lambda, &
                ctx%state%vega_spec * ctx%state%bands(:, i) / ctx%state%spec_lambda )
            
            if (ctx%state%magvega(i) > SAFE_FLOOR) then
                ctx%state%magvega(i) = -2.5_wp * log10(ctx%state%magvega(i)) - 48.60_wp
            else
                ctx%state%magvega(i) = 99.0_wp
            end if

            ! Sun Magnitude
            ctx%state%magsun(i) = integrate_trapezoid_array( &
                ctx%state%spec_lambda, &
                ctx%state%sun_spec * ctx%state%bands(:, i) / ctx%state%spec_lambda )

            if (ctx%state%magsun(i) > SAFE_FLOOR) then
                ctx%state%magsun(i) = -2.5_wp * log10(ctx%state%magsun(i)) - 48.60_wp
            else
                ctx%state%magsun(i) = 99.0_wp
            end if
        end do

        if (env_flag_true('FSPS_DEBUG_IO')) then
            call debug_stats_1d('[DEBUG_IO] vega_spec', ctx%state%vega_spec)
            call debug_stats_1d('[DEBUG_IO] sun_spec', ctx%state%sun_spec)
            call debug_stats_1d('[DEBUG_IO] magvega', ctx%state%magvega)
            call debug_stats_1d('[DEBUG_IO] magsun', ctx%state%magsun)
        end if
        
        deallocate(raw_lam, raw_flux)

    end subroutine load_standard_sed

    !> @brief Loads spectral index definitions from disk.
    !>
    !> @details
    !> Reads the `allindices.dat` file which defines the wavelength windows for 
    !> various spectral indices (e.g., Lick indices, D4000, Balmer lines).
    !>
    !> **Air-to-Vacuum Conversion:**
    !> The first 25 indices in `allindices.dat` correspond to the standard Lick system, 
    !> which is defined in air wavelengths. This routine automatically converts these 
    !> definitions to vacuum wavelengths to match the internal FSPS spectral grid. 
    !> Indices > 25 are assumed to be defined in vacuum already.
    !>
    !> @param[inout] ctx The FSPS context. Populates `ctx%state%indexdefined`.
    subroutine load_index_definitions(ctx)
        use fsps_cosmology, only: air_to_vacuum
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i
        character(len=1024) :: file_path

        file_path = trim(ctx%sps_home) // '/data/allindices.dat'
        
        if (VERBOSE > 1) print *, 'Reading Indices: ', trim(file_path)

        open(newunit=u_file, file=trim(file_path), status='old', &
             action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        ! Burn header (4 lines) 
        do i = 1, 4
            read(u_file, *)
        end do

        do i = 1, ctx%state%nindx
            read(u_file, *, iostat=io_stat) ctx%state%indexdefined(:, i)
            if (io_stat /= 0) call handle_read_error(file_path, i)

            ! Convert Lick indices from air to vacuum [cite: 162]
            if (i <= 25) then
                ctx%state%indexdefined(1:6, i) = air_to_vacuum(ctx%state%indexdefined(1:6, i))
            end if
        end do

        close(u_file)
    end subroutine load_index_definitions

    ! ========================================================================
    ! AUXILIARY SPECTRA I/O
    ! ========================================================================

    !> @brief Loads Wolf-Rayet (WR) spectra and aligns them to the context's Z grid.
    !> 
    !> @details
    !> Reads the CMFGEN model grids for Nitrogen-rich (WN) and Carbon-rich (WC) stars.
    !> The raw data comes on a specific, coarse metallicity grid (typically 5 points).
    !> This routine reads that data and immediately interpolates it onto:
    !> 1. The master wavelength grid (`spec_lambda`).
    !> 2. The target isochrone metallicity grid (`zlegend`).
    !>
    !> **Processing Steps:**
    !> 1. Reads effective temperature grids from `CMFGEN_WN.teff` and `CMFGEN_WC.teff`.
    !> 2. Reads raw spectra from `CMFGEN_WN_Zall.spec` and `CMFGEN_WC_Zall.spec`.
    !> 3. Interpolates the raw spectra (log-space flux) from the input Z-grid to the 
    !>    simulation's Z-grid defined in `zlegend`.
    !> 4. Interpolates the result onto the master wavelength grid.
    !>
    !> @param[inout] ctx The FSPS context. Populates `wrn_spec`, `wrc_spec`, and their logT arrays.
    subroutine load_wr_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i, j, i1
        character(len=1024) :: file_path
        real(WP) :: d1, dz
        
        ! Raw data buffers (Hardcoded sizes from legacy)
        real(WP), dimension(5) :: twrzmet
        real(WP), dimension(NLAMWR) :: tlamwr
        real(WP), dimension(NLAMWR) :: tspecwr
        ! Heap allocation for large array
        real(WP), allocatable :: twr_raw(:,:,:)
        
        allocate(twr_raw(NLAMWR, NDIM_WR, 5))

        ! 1. Load Teff Arrays
        ! -------------------
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/CMFGEN_WN.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        read(u_file, *) ctx%state%wrn_logt
        close(u_file)

        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/CMFGEN_WC.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        read(u_file, *) ctx%state%wrc_logt
        close(u_file)

        ! 2. Load and Process WN Spectra
        ! ------------------------------
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/CMFGEN_WN_Zall.spec'
        if (VERBOSE > 1) print *, 'Reading/Aligning WR-N: ', trim(file_path)
        
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        read(u_file, *) tlamwr
        do j = 1, 5
            do i = 1, NDIM_WR
                read(u_file, *) d1, twrzmet(j)
                read(u_file, *) twr_raw(:, i, j)
            end do
        end do
        close(u_file)

        ! Convert raw Z to log10(Z/Zsol) for interpolation
        twrzmet = log10(twrzmet / ctx%state%zsol_spec)

        ! Interpolate to target Z grid (ctx%state%zlegend) and Wavelengths
        do j = 1, ctx%state%nz
            ! Find Z interval
            i1 = min(max(find_interval(twrzmet, log10(ctx%state%zlegend(j)/ctx%state%zsol_spec)), 1), 4)
            dz = (log10(ctx%state%zlegend(j)/ctx%state%zsol_spec) - twrzmet(i1)) / (twrzmet(i1+1) - twrzmet(i1))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            do i = 1, NDIM_WR
                ! Interpolate in Z (log-space flux)
                tspecwr = (1.0_wp - dz) * log10(twr_raw(:, i, i1) + SAFE_FLOOR) + &
                          dz * log10(twr_raw(:, i, i1+1) + SAFE_FLOOR)
                
                ! Interpolate in Wavelength
                ctx%state%wrn_spec(:, i, j) = 10.0_wp**interpolate_linear( &
                    log10(tlamwr), tspecwr, log10(ctx%state%spec_lambda)) - SAFE_FLOOR
            end do
        end do

        ! 3. Load and Process WC Spectra (Identical Logic)
        ! ------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/CMFGEN_WC_Zall.spec'
        if (VERBOSE > 1) print *, 'Reading/Aligning WR-C: ', trim(file_path)
        
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        read(u_file, *) tlamwr
        do j = 1, 5
            do i = 1, NDIM_WR
                read(u_file, *) d1, twrzmet(j)
                read(u_file, *) twr_raw(:, i, j)
            end do
        end do
        close(u_file)

        twrzmet = log10(twrzmet / ctx%state%zsol_spec)

        do j = 1, ctx%state%nz
            i1 = min(max(find_interval(twrzmet, log10(ctx%state%zlegend(j)/ctx%state%zsol_spec)), 1), 4)
            dz = (log10(ctx%state%zlegend(j)/ctx%state%zsol_spec) - twrzmet(i1)) / (twrzmet(i1+1) - twrzmet(i1))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            do i = 1, NDIM_WR
                tspecwr = (1.0_wp - dz) * log10(twr_raw(:, i, i1) + SAFE_FLOOR) + &
                          dz * log10(twr_raw(:, i, i1+1) + SAFE_FLOOR)
                
                ctx%state%wrc_spec(:, i, j) = 10.0_wp**interpolate_linear( &
                    log10(tlamwr), tspecwr, log10(ctx%state%spec_lambda)) - SAFE_FLOOR
            end do
        end do

        deallocate(twr_raw)

    end subroutine load_wr_spectra

    !> @brief Loads AGB (Asymptotic Giant Branch) spectral libraries.
    !>
    !> @details
    !> Loads distinct datasets required for accurate modeling of the TP-AGB phase.
    !> This routine handles three specific libraries:
    !> 1. **O-rich AGB stars:** (Empirical, Lancon & Wood). Reads `Orich.teff` and `Orich.spec`.
    !>    - Interpolates the Teff grid to the target isochrone metallicities.
    !>    - Interpolates spectra onto the master wavelength grid.
    !> 2. **C-rich AGB stars:** (Empirical, Lancon & Wood). Reads `Crich.teff` and `Crich.spec`.
    !>    - Similar processing to O-rich.
    !> 3. **Aringer C-rich models:** (Theoretical). Reads `Crich_Aringer.teff` and `Crich_Aringer.spec`.
    !>
    !> **Interpolation:**
    !> All raw spectra are immediately interpolated onto the global `spec_lambda` grid 
    !> upon loading to ensure compatibility with the rest of the synthesis engine.
    !>
    !> @param[inout] ctx The FSPS context. Populates `agb_spec_o`, `agb_spec_c`, 
    !>                   `agb_spec_car`, and their associated `logt` arrays.
    subroutine load_agb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i, i1, i_spec
        character(len=1024) :: file_path
        character(len=1) :: char_dummy
        real(WP) :: dumr1, dz
        
        ! Constants
        integer, parameter :: NSPEC_AGB = 6146
        integer, parameter :: NSPEC_ARINGER = 9032
        
        ! Allocatable buffers for large arrays
        real(WP), allocatable :: agb_lam(:)
        real(WP), allocatable :: temp_row_o(:), temp_row_c(:)
        real(WP), allocatable :: raw_spec_o(:,:), raw_spec_c(:,:)
        
        real(WP), allocatable :: aringer_lam(:)
        real(WP), allocatable :: temp_row_car(:)
        real(WP), allocatable :: raw_spec_car(:,:)
        
        real(WP), dimension(22, N_AGB_O) :: tagb_logt_o_raw
        real(WP), dimension(22) :: tagb_logz_o

        ! Allocate heap memory
        allocate(agb_lam(NSPEC_AGB), aringer_lam(NSPEC_ARINGER))
        allocate(temp_row_o(N_AGB_O), temp_row_c(N_AGB_C), temp_row_car(N_AGB_CAR))
        allocate(raw_spec_o(NSPEC_AGB, N_AGB_O))
        allocate(raw_spec_c(NSPEC_AGB, N_AGB_C))
        allocate(raw_spec_car(NSPEC_ARINGER, N_AGB_CAR))

        ! 1. O-rich Teff Grid
        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Orich.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        read(u_file, *) char_dummy
        read(u_file, *) dumr1, tagb_logz_o
        do i = 1, N_AGB_O
            read(u_file, *) dumr1, tagb_logt_o_raw(:, i)
        end do
        close(u_file)

        do i = 1, ctx%state%nz
            i1 = min(max(find_interval(tagb_logz_o, log10(ctx%state%zlegend(i)/ctx%state%zsol_spec)), 1), 21)
            dz = (log10(ctx%state%zlegend(i)/ctx%state%zsol_spec) - tagb_logz_o(i1)) / (tagb_logz_o(i1+1) - tagb_logz_o(i1))
            ctx%state%agb_logt_o(i, :) = (1.0_wp - dz) * tagb_logt_o_raw(i1, :) + dz * tagb_logt_o_raw(i1+1, :)
        end do
        ctx%state%agb_logt_o = log10(ctx%state%agb_logt_o)

        ! 2. C-rich Teff Grid
        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Crich.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        read(u_file, *) char_dummy
        do i = 1, N_AGB_C
            read(u_file, *) dumr1, ctx%state%agb_logt_c(i)
        end do
        close(u_file)
        ctx%state%agb_logt_c = log10(ctx%state%agb_logt_c)

        ! 3. O-rich Spectra (Read & Interpolate)
        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Orich.spec'
        if (VERBOSE > 1) print *, 'Reading AGB Orich: ', trim(file_path)
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        do i = 1, NSPEC_AGB
            read(u_file, *) agb_lam(i), temp_row_o
            raw_spec_o(i, :) = temp_row_o
        end do
        close(u_file)
        
        do i = 1, N_AGB_O
            do i_spec = 1, ctx%state%nspec
                ctx%state%agb_spec_o(i_spec, i) = max(interpolate_linear(agb_lam, raw_spec_o(:, i), &
                                                   ctx%state%spec_lambda(i_spec)), SAFE_FLOOR)
            end do
        end do

        ! 4. C-rich Spectra (Read & Interpolate)
        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Crich.spec'
        if (VERBOSE > 1) print *, 'Reading AGB Crich: ', trim(file_path)
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        do i = 1, NSPEC_AGB
            read(u_file, *) agb_lam(i), temp_row_c
            raw_spec_c(i, :) = temp_row_c
        end do
        close(u_file)
        
        do i = 1, N_AGB_C
            do i_spec = 1, ctx%state%nspec
                ctx%state%agb_spec_c(i_spec, i) = max(interpolate_linear(agb_lam, raw_spec_c(:, i), &
                                                   ctx%state%spec_lambda(i_spec)), SAFE_FLOOR)
            end do
        end do

        ! 5. Aringer C-rich (Read & Interpolate)
        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Crich_Aringer.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        read(u_file, *) ctx%state%agb_logt_car
        close(u_file)
        ctx%state%agb_logt_car = log10(ctx%state%agb_logt_car)

        file_path = trim(ctx%sps_home) // '/data/spectra/AGB_spectra/Crich_Aringer.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        do i = 1, NSPEC_ARINGER
            read(u_file, *) aringer_lam(i), temp_row_car
            raw_spec_car(i, :) = temp_row_car
        end do
        close(u_file)
        
        do i = 1, N_AGB_CAR
            do i_spec = 1, ctx%state%nspec
                ctx%state%agb_spec_car(i_spec, i) = max(interpolate_linear(aringer_lam, raw_spec_car(:, i), &
                                                     ctx%state%spec_lambda(i_spec)), SAFE_FLOOR)
            end do
        end do

        ! Deallocate temporary buffers
        deallocate(agb_lam, aringer_lam)
        deallocate(temp_row_o, temp_row_c, temp_row_car)
        deallocate(raw_spec_o, raw_spec_c, raw_spec_car)

    end subroutine load_agb_spectra

    !> @brief Loads Post-AGB spectral libraries (Rauch 2003).
    !>
    !> @details
    !> Reads the H-deficient Post-AGB spectral models provided by Rauch (2003).
    !> These are used to model the hot, evolved phases of low-to-intermediate mass stars.
    !>
    !> **Processing Steps:**
    !> 1. Reads the effective temperature grid from `ipagb.teff`.
    !> 2. Reads the "Solar" metallicity spectra from `ipagb_solar.spec`.
    !> 3. Reads the "Halo" metallicity spectra from `ipagb_halo.spec`.
    !> 4. Linearly interpolates both grids onto the master wavelength grid (`spec_lambda`).
    !>
    !> **Storage:**
    !> The internal `pagb_spec` array is stored as (Lambda, Teff, Metallicity), where:
    !> - Metallicity index 1 = Halo
    !> - Metallicity index 2 = Solar
    !>
    !> @param[inout] ctx The FSPS context. Populates `ctx%state%pagb_spec` and `pagb_logt`.
    subroutine load_post_agb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        
        integer :: u_file, io_stat, i, j, i_spec
        character(len=1024) :: file_path
        integer, parameter :: NSPEC_PAGB = 9281
        
        ! Allocatable buffers
        real(WP), allocatable :: pagb_lam(:)
        real(WP), allocatable :: pagb_raw(:,:,:) ! 1=Halo, 2=Solar

        allocate(pagb_lam(NSPEC_PAGB))
        allocate(pagb_raw(NSPEC_PAGB, NDIM_PAGB, 2))
        
        ! 1. Teff
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/ipagb.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        do i = 1, NDIM_PAGB
            read(u_file, *) ctx%state%pagb_logt(i)
        end do
        close(u_file)
        ctx%state%pagb_logt = log10(ctx%state%pagb_logt)

        ! 2. Solar Spec
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/ipagb_solar.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        do i = 1, NSPEC_PAGB
            read(u_file, *) pagb_lam(i), pagb_raw(i, :, 2)
        end do
        close(u_file)

        ! 3. Halo Spec
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/ipagb_halo.spec'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        do i = 1, NSPEC_PAGB
            read(u_file, *) pagb_lam(i), pagb_raw(i, :, 1)
        end do
        close(u_file)

        ! 4. Interpolate
        do j = 1, 2
            do i = 1, NDIM_PAGB
                do i_spec = 1, ctx%state%nspec
                    ctx%state%pagb_spec(i_spec, i, j) = max(interpolate_linear(pagb_lam, pagb_raw(:, i, j), &
                                                        ctx%state%spec_lambda(i_spec)), SAFE_FLOOR)
                end do
            end do
        end do

        deallocate(pagb_lam, pagb_raw)

    end subroutine load_post_agb_spectra

    !> @brief Loads dust attenuation curves from Witt & Gordon (2000) and Gordon et al. (2003).
    !>
    !> @details
    !> This routine populates the dust attenuation arrays required for the complex dust models 
    !> (dust_type=1 and dust_type=4). It performs two main operations:
    !>
    !> **1. Witt & Gordon (2000) "Dirty" Dust:**
    !> Reads `alldirty_h.dat` (Shell geometry) and `alldirty_c.dat` (Dusty geometry). 
    !> The raw data is read, permuted, and interpolated onto the master wavelength grid.
    !> - **Input File Format:** (Lambda, Geom, Type, Tau)
    !> - **Internal Storage:** `ctx%state%wgdust` (Lambda, Geom, Tau, Type)
    !>
    !> **2. SMC Bar Extinction (Gordon et al. 2003):**
    !> Reads `Gordon03_table4.dat` which defines the SMC Bar extinction curve.
    !> Interpolates onto `spec_lambda` to populate `ctx%state%g03smcextn`.
    !>
    !> **Interpolation Behavior:**
    !> - Uses linear interpolation for all curves.
    !> - Wavelengths falling outside the range defined in the data files are padded with 0.0.
    !>
    !> @param[inout] ctx The FSPS context. Populates `wgdust` and `g03smcextn`.
    subroutine load_attenuation_curves(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i, j, k, n
        character(len=1024) :: file_path
        real(WP) :: d1
        
        ! Constants for WG00 Data structure
        integer, parameter :: NWG_LAM = 25
        integer, parameter :: NWG_GEOM = 18
        integer, parameter :: NWG_TAU = 6
        integer, parameter :: NSMC_LAM = 30

        ! Temporary buffers
        real(WP), dimension(NWG_LAM) :: wglam
        ! Raw: (Lambda, Geom, Type[H/C], Tau)
        real(WP), dimension(NWG_LAM, NWG_GEOM, 2, NWG_TAU) :: wgtmp
        
        real(WP), dimension(NSMC_LAM) :: g03lam, g03smc

        ! --------------------------------------------------------------------
        ! 1. Witt & Gordon (2000) - Shell (H)
        ! --------------------------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/dust/alldirty_h.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        ! Burn header (3 lines)
        read(u_file, *)
        read(u_file, *)
        read(u_file, *)
        
        do i = 1, NWG_GEOM
            do j = 1, NWG_LAM
                ! File cols: Lambda, Albedo(?), [6 Tau values]
                read(u_file, *) wglam(j), d1, wgtmp(j, i, 1, :)
            end do
        end do
        close(u_file)

        ! --------------------------------------------------------------------
        ! 2. Witt & Gordon (2000) - Dust (C)
        ! --------------------------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/dust/alldirty_c.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        
        ! Burn header (3 lines)
        read(u_file, *)
        read(u_file, *)
        read(u_file, *)
        
        do i = 1, NWG_GEOM
            do j = 1, NWG_LAM
                read(u_file, *) wglam(j), d1, wgtmp(j, i, 2, :)
            end do
        end do
        close(u_file)

        ! --------------------------------------------------------------------
        ! 3. Interpolate WG00 to Master Grid
        ! --------------------------------------------------------------------
        ! wgdust dimensions: (nspec, 18, 6, 2) => (Lam, Geom, Tau, Type)
        
        do k = 1, 2 ! Type (1=H, 2=C)
            do i = 1, NWG_GEOM
                do j = 1, NWG_TAU
                    do n = 1, ctx%state%nspec
                        
                        if (ctx%state%spec_lambda(n) > wglam(NWG_LAM)) then
                            ctx%state%wgdust(n, i, j, k) = 0.0_wp
                        else if (ctx%state%spec_lambda(n) < wglam(1)) then
                            ctx%state%wgdust(n, i, j, k) = wgtmp(1, i, k, j)
                        else
                            ctx%state%wgdust(n, i, j, k) = interpolate_linear( &
                                wglam, wgtmp(:, i, k, j), ctx%state%spec_lambda(n))
                        end if
                        
                    end do
                end do
            end do
        end do

        ! --------------------------------------------------------------------
        ! 4. Gordon et al. (2003) SMC Bar
        ! --------------------------------------------------------------------
        file_path = trim(ctx%sps_home) // '/data/dust/Gordon03_table4.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        do i = 1, NSMC_LAM
            ! Legacy reads in reverse order (30 down to 1)
            read(u_file, *) g03lam(NSMC_LAM - i + 1), d1, g03smc(NSMC_LAM - i + 1)
        end do
        close(u_file)

        ! Convert microns to Angstroms
        g03lam = g03lam * 1.0e4_wp

        ! Interpolate SMC
        do n = 1, ctx%state%nspec
            if (ctx%state%spec_lambda(n) > g03lam(NSMC_LAM)) then
                ctx%state%g03smcextn(n) = 0.0_wp
            else if (ctx%state%spec_lambda(n) < g03lam(1)) then
                ctx%state%g03smcextn(n) = g03smc(1)
            else
                ctx%state%g03smcextn(n) = interpolate_linear( &
                    g03lam, g03smc, ctx%state%spec_lambda(n))
            end if
        end do

    end subroutine load_attenuation_curves

    !> @brief Loads the WMBasic spectral library for hot O-stars.
    !>
    !> @details
    !> Reads the theoretical O-star spectral grid from Eldridge et al. (based on Pauldrach 
    !> et al. 2001 models). This library is essential for accurate UV modeling of young 
    !> stellar populations.
    !>
    !> **Processing Steps:**
    !> 1. Reads the Teff grid from `WMBASIC.teff`.
    !> 2. Sets a hardcoded surface gravity grid (log g = 3.5, 4.0, 4.5).
    !> 3. Reads spectra for all metallicities defined in `WMBASIC_zlegend.dat`.
    !> 4. Interpolates the raw spectra onto:
    !>    - The master wavelength grid (`spec_lambda`).
    !>    - The target isochrone metallicity grid (`zlegend`).
    !>
    !> @note
    !> The internal storage `wmb_spec` ends up being 4D: (Lambda, Z, LogT, LogG).
    !>
    !> @param[inout] ctx The FSPS context. Populates `ctx%state%wmb_spec`.
    subroutine load_wmbasic_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: u_file, io_stat, i, j, z, i1, i_spec
        character(len=1024) :: file_path
        character(len=6) :: zstype
        real(WP) :: dz, log_spec_val
        
        real(WP), dimension(NZWMB) :: zwmb
        
        ! Allocatable buffers for raw data handling
        real(WP), allocatable :: wmb_lam(:)
        real(WP), allocatable :: wmb_spec_raw(:,:,:) ! (nspec, logt, logg)
        real(WP), allocatable :: wmbsi(:,:,:,:)      ! (nspec_master, nz, logt, logg)

        allocate(wmb_lam(NSPEC_WMB))
        allocate(wmb_spec_raw(NSPEC_WMB, NDIM_WMB_LOGT, NDIM_WMB_LOGG))
        allocate(wmbsi(ctx%state%nspec, NZWMB, NDIM_WMB_LOGT, NDIM_WMB_LOGG))

        ! 1. Read Teff
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/WMBASIC.teff'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)
        do i = 1, NDIM_WMB_LOGT
            read(u_file, *) ctx%state%wmb_logt(i)
        end do
        close(u_file)

        ! Hardcoded Logg grid from legacy
        ctx%state%wmb_logg = [3.5_wp, 4.0_wp, 4.5_wp]

        ! 2. Read Z-Legend and Spectra
        file_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/WMBASIC_zlegend.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        do z = 1, NZWMB
            read(u_file, *) zwmb(z)
            write(zstype, '(F6.4)') zwmb(z)
            
            ! Read raw file for this Z
            call load_wmb_raw_z(ctx, zstype, wmb_lam, wmb_spec_raw)
            
            ! Interpolate onto master wavelength grid immediately
            do j = 1, NDIM_WMB_LOGG
                do i = 1, NDIM_WMB_LOGT
                    do i_spec = 1, ctx%state%nspec
                        wmbsi(i_spec, z, i, j) = max(interpolate_linear(wmb_lam, wmb_spec_raw(:, i, j), &
                                                     ctx%state%spec_lambda(i_spec)), SAFE_FLOOR)
                    end do
                end do
            end do
        end do
        close(u_file)

        ! 3. Interpolate to Isochrone Z grid
        do z = 1, ctx%state%nz
            i1 = min(max(find_interval(log10(zwmb/ctx%state%zsol_spec), &
                                       log10(ctx%state%zlegend(z)/ctx%state%zsol)), 1), NZWMB-1)
            
            dz = (log10(ctx%state%zlegend(z)/ctx%state%zsol) - log10(zwmb(i1)/ctx%state%zsol_spec)) / &
                 (log10(zwmb(i1+1)/ctx%state%zsol_spec) - log10(zwmb(i1)/ctx%state%zsol_spec))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            ! Interpolate in metallicity with scalar loops (ifx-safe path).
            do j = 1, NDIM_WMB_LOGG
                do i = 1, NDIM_WMB_LOGT
                    do i_spec = 1, ctx%state%nspec
                        log_spec_val = (1.0_wp - dz) * log10(wmbsi(i_spec, i1, i, j) + SAFE_FLOOR) + &
                                       dz * log10(wmbsi(i_spec, i1 + 1, i, j) + SAFE_FLOOR)
                        ctx%state%wmb_spec(i_spec, z, i, j) = 10.0_wp**log_spec_val
                    end do
                end do
            end do
        end do

        deallocate(wmb_lam, wmb_spec_raw, wmbsi)
    end subroutine load_wmbasic_spectra

    !> @brief Loads the Line Spread Function (LSF) for spectral smoothing.
    !>
    !> @details
    !> Reads the user-supplied LSF curve from `data/lsf.dat`. This file defines the 
    !> instrumental or physical broadening profile (wavelength vs sigma) to be applied 
    !> to the SSPs if smoothing is enabled.
    !>
    !> **Processing Steps:**
    !> 1. Reads raw (Wavelength, Sigma) pairs from the ASCII file.
    !> 2. Identifies the valid wavelength range (`minlam`, `maxlam`).
    !> 3. Linearly interpolates the LSF sigma onto the master wavelength grid (`spec_lambda`).
    !>
    !> **Boundary Conditions:**
    !> Master wavelengths falling outside the range defined in `lsf.dat` are assigned 
    !> an LSF sigma of 0.0 (no smoothing).
    !>
    !> @param[inout] ctx The FSPS context. Populates `ctx%state%lsfinfo`.
    subroutine load_lsf_data(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        
        integer :: u_file, io_stat, i
        character(len=1024) :: file_path
        ! Use allocatable arrays to avoid stack overflow (NTABMAX can be large)
        real(WP), allocatable :: lsflam(:), lsfsig(:)
        integer :: n_points

        allocate(lsflam(NTABMAX), lsfsig(NTABMAX))

        ! Initialize to zero to prevent garbage data
        lsflam = 0.0_wp
        lsfsig = 0.0_wp

        file_path = trim(ctx%sps_home) // '/data/lsf.dat'
        open(newunit=u_file, file=trim(file_path), status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_read_error(file_path, 0)

        n_points = 0
        do i = 1, NTABMAX
            read(u_file, *, iostat=io_stat) lsflam(i), lsfsig(i)
            if (io_stat /= 0) exit
            n_points = i
        end do
        close(u_file)

        if (n_points == 0) then
             write(error_unit, '("[FSPS_IO] Error: lsf.dat is empty or invalid")')
             error stop
        end if

        ctx%state%lsfinfo%minlam = lsflam(1)
        ctx%state%lsfinfo%maxlam = lsflam(n_points)

        ! Interpolate onto master grid
        do i = 1, ctx%state%nspec
            if (ctx%state%spec_lambda(i) >= lsflam(1) .and. &
                ctx%state%spec_lambda(i) <= lsflam(n_points)) then
                
                ctx%state%lsfinfo%lsf(i) = interpolate_linear( &
                    lsflam(1:n_points), lsfsig(1:n_points), ctx%state%spec_lambda(i))
            else
                ctx%state%lsfinfo%lsf(i) = 0.0_wp
            end if
        end do

        deallocate(lsflam, lsfsig)

    end subroutine load_lsf_data

    ! ========================================================================
    ! TABULAR SFH
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
        character(len=128) :: fmt_mags
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

    !> @brief Writes the Isochrone Color-Magnitude Diagram (CMD) to a file.
    !>
    !> @details
    !> Dumps the detailed isochrone properties for a specific metallicity to a disk file 
    !> (e.g., `OUTPUTS/mysps.cmd`). This matches the legacy output format of `write_isochrone.f90`.
    !>
    !> **Design Philosophy:**
    !> This routine is strictly an I/O handler. It expects that all physical modifications 
    !> (e.g., Blue Stragglers, Horizontal Branch shifts, TP-AGB corrections) and magnitude 
    !> calculations have already been performed by the caller. It simply formats and writes 
    !> the provided arrays to disk.
    !>
    !> **Output Format:**
    !> - **Header:** Single line describing columns.
    !> - **Data:** One row per mass point per isochrone time step.
    !> - **Columns:** !>   1. Log(Age) (Years)
    !>   2. Log(Z/Zsol)
    !>   3. Initial Mass
    !>   4. Actual Mass (current)
    !>   5. Log(Luminosity/Lsol)
    !>   6. Log(Teff)
    !>   7. Log(g)
    !>   8. Evolutionary Phase Index
    !>   9. Composition Flag (C/O ratio for AGB)
    !>   10. Log(IMF Weight)
    !>   11. Log(Mass Loss Rate)
    !>   12+. Magnitudes in all active bands (nbands).
    !>
    !> @param[in] ctx        The FSPS context (provides output paths and metadata).
    !> @param[in] pset       User parameters (specifies the Z index to write).
    !> @param[in] outfile    Base filename for output (e.g., "mysps").
    !> @param[in] time_grid  Array of Log(Age) for each time step [nt].
    !> @param[in] nmass      Array containing the number of mass points per isochrone [nt].
    !> @param[in] mini       Initial Mass grid [nt, NM].
    !> @param[in] mact       Actual Mass grid [nt, NM].
    !> @param[in] logl       Log Luminosity grid [nt, NM].
    !> @param[in] logt       Log Temperature grid [nt, NM].
    !> @param[in] logg       Log Surface Gravity grid [nt, NM].
    !> @param[in] phase      Evolutionary Phase grid [nt, NM].
    !> @param[in] ffco       Composition flag grid [nt, NM].
    !> @param[in] lmdot      Log Mass Loss Rate grid [nt, NM].
    !> @param[in] weights    Integrated IMF weights [nt, NM].
    !> @param[in] mags       Computed magnitudes [nt, NM, nbands].
    subroutine write_isochrone_cmd(ctx, pset, outfile, time_grid, nmass, &
                                   mini, mact, logl, logt, logg, &
                                   phase, ffco, lmdot, weights, mags)
        
        type(fsps_context_t), intent(in) :: ctx
        type(params), intent(in) :: pset
        character(len=*), intent(in) :: outfile
        
        ! Data Arrays
        ! Note: We use assumed-shape arrays. The caller must ensure 
        ! dimensions match (nt, NM).
        real(WP), dimension(:), intent(in) :: time_grid
        integer, dimension(:), intent(in) :: nmass
        real(WP), dimension(:,:), intent(in) :: mini, mact, logl, logt, logg
        real(WP), dimension(:,:), intent(in) :: phase, ffco, lmdot, weights
        real(WP), dimension(:,:,:), intent(in) :: mags ! (nt, NM, nbands)
        
        ! Locals
        integer :: u_file, i_t, i_m, nt, z_idx
        character(len=1024) :: filepath
        character(len=128) :: fmt_str
        real(WP) :: z_log
        
        nt = size(time_grid)
        z_idx = pset%zmet
        
        ! Calculate Log(Z/Zsol) for the header column
        ! Matches legacy behavior of writing log(Z) not Z_index
        if (ctx%state%zlegend(z_idx) > 0.0_wp) then
            z_log = log10(ctx%state%zlegend(z_idx))
        else
            z_log = -99.0_wp
        end if

        ! Construct File Path
        filepath = trim(ctx%output_home) // '/OUTPUTS/' // trim(outfile) // '.cmd'
        
        if (VERBOSE > 0) print *, 'Writing Isochrone CMD: ', trim(filepath)

        open(newunit=u_file, file=trim(filepath), status='replace', action='write')
        
        ! Header
        write(u_file, '(A)') '# age log(Z) mini mact logl logt logg phase composition log(weight) log(mdot) mags'
        
        ! Construct Format String (Dynamic based on nbands)
        ! Legacy format string construction: 
        ! '(F7.4,1x,F8.4,1x,F14.9,1x,F14.9,1x,7(F8.4,1x),000(F7.3,1x))'
        write(fmt_str, '("(F7.4,1x,F8.4,1x,F14.9,1x,F14.9,1x,7(F8.4,1x),", I0, "(F7.3,1x))")') ctx%state%nbands

        ! Loop over time steps
        do i_t = 1, nt
            ! Loop over mass points in this isochrone
            do i_m = 1, nmass(i_t)
                
                write(u_file, fmt_str) &
                    time_grid(i_t), &
                    z_log, &
                    mini(i_t, i_m), &
                    mact(i_t, i_m), &
                    logl(i_t, i_m), &
                    logt(i_t, i_m), &
                    logg(i_t, i_m), &
                    phase(i_t, i_m), &
                    ffco(i_t, i_m), &
                    log10(max(weights(i_t, i_m), SAFE_FLOOR)), &
                    lmdot(i_t, i_m), &
                    mags(i_t, i_m, :)
                    
            end do
        end do
        
        close(u_file)

    end subroutine write_isochrone_cmd

    !> @brief Compiles a raw ASCII spectral library into the FSPS binary format.
    !>
    !> @details
    !> Reads a raw ASCII spectral file (e.g., from BaSeL, MILES, or C3K) and compiles it 
    !> into an unformatted direct-access binary file for efficient runtime loading. 
    !> This routine replicates the functionality of the legacy `spec_bin` program.
    !>
    !> **Input Format (ASCII):**
    !> - Expected filename pattern matches `get_ascii_spectral_filename`.
    !> - Structure: Nested loops over Gravity (outer) and Temperature (inner).
    !> - Each block has a header line followed by the spectral flux array.
    !>
    !> **Output Format (Binary):**
    !> - Filename pattern matches `get_spectral_filename`.
    !> - Format: Unformatted Direct Access.
    !> - Record Length: `NSPEC * NDIM_LOGT * NDIM_LOGG * 4` bytes.
    !> - Precision: **32-bit floats** (`real32`). Data is strictly stored as single precision 
    !>   to minimize disk footprint and maintain backward compatibility.
    !>
    !> @param[in] ctx       The FSPS context (provides dimensions and paths).
    !> @param[in] spec_type String identifier for the library (e.g., 'miles', 'basel').
    !> @param[in] z_idx     Metallicity index to process (1..nzinit).
    subroutine write_binary_spectral_lib(ctx, spec_type, z_idx)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        integer, intent(in) :: z_idx

        integer :: u_in, u_out, io_stat, i_g, i_t, rec_len, file_unit_size
        character(len=1024) :: file_ascii, file_bin
        
        ! Buffers
        ! Only allocate the 32-bit output cube (Legacy format requirements)
        real(real32), allocatable :: spec_cube_32(:,:,:) 
        ! Tiny 1D buffer for reading a single spectrum line
        real(WP), dimension(:), allocatable :: temp_spec_row

        ! 1. Setup Paths
        ! --------------
        call get_ascii_spectral_filename(ctx, spec_type, z_idx, file_ascii)
        ! Reuse the existing binary path getter
        call get_spectral_filename(ctx, spec_type, z_idx, file_bin) 

        if (VERBOSE > 0) then
            print *, 'Compiling Spectrum: ', trim(spec_type)
            print *, '   Input:  ', trim(file_ascii)
            print *, '   Output: ', trim(file_bin)
        end if

        ! Allocate huge output buffer
        allocate(spec_cube_32(ctx%state%nspec, NDIM_LOGT, NDIM_LOGG))
        ! Allocate tiny read buffer
        allocate(temp_spec_row(ctx%state%nspec))

        ! 2. Read ASCII Input
        ! -------------------
        open(newunit=u_in, file=trim(file_ascii), status='old', &
             action='read', iostat=io_stat)
        
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening ASCII spec: ", A)') trim(file_ascii)
            error stop
        end if

        ! The ASCII format is strictly nested loops matching sps_program_spec_bin.f90:
        ! Outer Loop: logg (1..NDIM_LOGG)
        ! Inner Loop: logt (1..NDIM_LOGT)
        ! Header line: index, T, g, etc (read as dummy)
        ! Body: Spectrum array
        
        do i_g = 1, NDIM_LOGG
            do i_t = 1, NDIM_LOGT
                
                ! Skip metadata header
                read(u_in, *, iostat=io_stat) 
                if (io_stat /= 0) call handle_read_error(file_ascii, 0)
                
                ! Read into small temporary buffer
                read(u_in, *, iostat=io_stat) temp_spec_row
                if (io_stat /= 0) call handle_read_error(file_ascii, 0)
                
                ! Convert and store immediately
                spec_cube_32(:, i_t, i_g) = real(temp_spec_row, kind=real32)
                
            end do
        end do

        deallocate(temp_spec_row) ! Free small buffer
        close(u_in)

        ! 3. Write Binary Output
        ! ----------------------
        ! Get the unit size (in bits) for file storage
        ! usually 8 (bytes) or 32 (words)
        file_unit_size = file_storage_size / 8
        ! Adjust record length calculation
        rec_len = (ctx%state%nspec * NDIM_LOGT * NDIM_LOGG * 4) / file_unit_size
        
        open(newunit=u_out, file=trim(file_bin), status='replace', &
             access='direct', recl=rec_len, form='unformatted', &
             action='write', iostat=io_stat)
             
        if (io_stat /= 0) then
            write(error_unit, '("[FSPS_IO] Error opening Binary output: ", A)') trim(file_bin)
            error stop
        end if

        ! Write the entire cube as a single record
        write(u_out, rec=1) spec_cube_32
        
        close(u_out)
        
        deallocate(spec_cube_32)

    end subroutine write_binary_spectral_lib

    ! ------------------------------------------------------------------------
    ! PRIVATE HELPERS
    ! ------------------------------------------------------------------------

    subroutine handle_read_error(file_path, line_num)
        character(len=*), intent(in) :: file_path
        integer, intent(in) :: line_num
        write(error_unit, '("[FSPS_IO] Parsing error at line ", I0, " of file: ", A)') &
            line_num, trim(file_path)
        error stop
    end subroutine handle_read_error

    !> @brief Constructs the path for isochrone files based on type and metallicity.
    subroutine get_isochrone_filename(ctx, isoc_type, z_idx, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: isoc_type
        integer, intent(in) :: z_idx
        character(len=*), intent(out) :: file_path
        
        character(len=6) :: z_str
        real(WP) :: val_log, z_val
        
        ! Get the numerical Z value from the state
        z_val = ctx%state%zlegend(z_idx)

        ! Handle MIST's specific "mX.XX" / "pX.XX" naming convention
        if (trim(isoc_type) == 'mist') then
             ! Guard against log(0)
             if (z_val < tiny(0.0_wp)) then
                 val_log = -99.0_wp 
             else
                 val_log = log10(z_val / ctx%state%zsol)
             end if

             ! Construct the string (e.g., "m1.50" or "p0.25")
             ! Note: We add a small epsilon for rounding safety before formatting
             if (val_log < -0.001_wp) then
                 write(z_str, '("m", F4.2)') abs(val_log)
             else
                 write(z_str, '("p", F4.2)') abs(val_log)
             end if
             
             file_path = trim(ctx%sps_home) // '/data/isochrones/MIST/isoc_z' // &
                         trim(z_str) // '.dat'

        else
             ! Standard format (F6.4) e.g., "0.0190"
             write(z_str, '(F6.4)') z_val
             
             select case (trim(isoc_type))
             case ('pdva')
                 file_path = trim(ctx%sps_home) // '/data/isochrones/Padova/Padova2007/isoc_z' // z_str // '.dat'
             case ('prsc')
                 file_path = trim(ctx%sps_home) // '/data/isochrones/PARSEC/isoc_z' // z_str // '.dat'
             case ('bsti')
                 file_path = trim(ctx%sps_home) // '/data/isochrones/BaSTI/isoc_z' // z_str // '.dat'
             case ('gnva')
                 file_path = trim(ctx%sps_home) // '/data/isochrones/Geneva/isoc_z' // z_str // '.dat'
             case default
                 write(error_unit, '("[FSPS_IO] Unknown isochrone type: ", A)') trim(isoc_type)
                 error stop
             end select
        end if
    end subroutine get_isochrone_filename

    !> @brief Constructs the path for spectral library binaries.
    subroutine get_spectral_filename(ctx, spec_type, z_idx, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        integer, intent(in) :: z_idx
        character(len=*), intent(out) :: file_path
        
        character(len=6) :: z_str
        
        ! Note: Spectral libraries use zlegendinit for their Z values
        write(z_str, '(F6.4)') ctx%state%zlegendinit(z_idx)

        select case (trim(spec_type))
        case ('basel')
            file_path = trim(ctx%sps_home) // '/data/spectra/BaSeL3.1/basel_' // &
                        trim(BASEL_STR) // '_z' // z_str // '.spectra.bin'
        
        case ('miles')
            file_path = trim(ctx%sps_home) // '/data/spectra/MILES/imiles_z' // &
                        z_str // '.spectra.bin'
        
        case default
            if (spec_type(1:3) == 'c3k') then
                file_path = trim(ctx%sps_home) // '/data/spectra/C3K/' // &
                            trim(spec_type) // '_z' // z_str // '.spectra.bin'
            else if (trim(spec_type) == 'ckc14') then
                ! CKC14 support
                file_path = trim(ctx%sps_home) // '/data/spectra/CKC14/' // &
                            trim(spec_type) // '_z' // z_str // '.spectra.bin'
            else
                write(error_unit, '("[FSPS_IO] Unknown spec_type: ", A)') trim(spec_type)
                error stop
            end if
        end select
    end subroutine get_spectral_filename

    subroutine get_zlegend_filename(ctx, lib_type, is_spec, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: lib_type
        logical, intent(in) :: is_spec
        character(len=*), intent(out) :: file_path

        if (is_spec) then
            select case (trim(lib_type))
            case ('basel')
                file_path = trim(ctx%sps_home) // '/data/spectra/BaSeL3.1/zlegend.dat'
            case ('miles')
                file_path = trim(ctx%sps_home) // '/data/spectra/MILES/zlegend.dat'
            case default
                if (lib_type(1:3) == 'c3k') then
                    file_path = trim(ctx%sps_home) // '/data/spectra/C3K/zlegend.dat'
                else
                    write(error_unit, '("[FSPS_IO] Unknown spec type for zlegend: ", A)') trim(lib_type)
                    error stop
                end if
            end select
        else
            select case (trim(lib_type))
            case ('mist')
                file_path = trim(ctx%sps_home) // '/data/isochrones/MIST/zlegend.dat'
            case ('pdva')
                file_path = trim(ctx%sps_home) // '/data/isochrones/Padova/Padova2007/zlegend.dat'
            case ('prsc')
                file_path = trim(ctx%sps_home) // '/data/isochrones/PARSEC/zlegend.dat'
            case ('bsti')
                file_path = trim(ctx%sps_home) // '/data/isochrones/BaSTI/zlegend.dat'
            case ('gnva')
                file_path = trim(ctx%sps_home) // '/data/isochrones/Geneva/zlegend.dat'
            case ('bpss')
                file_path = trim(ctx%sps_home) // '/data/isochrones/BPASS/zlegend.dat'
            case default
                write(error_unit, '("[FSPS_IO] Unknown isochrone type for zlegend: ", A)') trim(lib_type)
                error stop
            end select
        end if
    end subroutine get_zlegend_filename

    subroutine get_lambda_filename(ctx, spec_type, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        character(len=*), intent(out) :: file_path
        
        select case (trim(spec_type))
        case ('basel')
            file_path = trim(ctx%sps_home) // '/data/spectra/BaSeL3.1/basel.lambda'
        case ('miles')
            file_path = trim(ctx%sps_home) // '/data/spectra/MILES/miles.lambda'
        case default
            if (spec_type(1:3) == 'c3k') then
                file_path = trim(ctx%sps_home) // '/data/spectra/C3K/' // trim(spec_type) // '.lambda'
            else
                write(error_unit, '("[FSPS_IO] Unknown spec type for lambda: ", A)') trim(spec_type)
                error stop
            end if
        end select
    end subroutine get_lambda_filename

    subroutine get_res_filename(ctx, spec_type, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        character(len=*), intent(out) :: file_path
        
        select case (trim(spec_type))
        case ('basel')
            file_path = trim(ctx%sps_home) // '/data/spectra/BaSeL3.1/basel.res'
        case ('miles')
            file_path = trim(ctx%sps_home) // '/data/spectra/MILES/miles.res'
        case default
            if (spec_type(1:3) == 'c3k') then
                file_path = trim(ctx%sps_home) // '/data/spectra/C3K/' // trim(spec_type) // '.res'
            else
                write(error_unit, '("[FSPS_IO] Unknown spec type for res: ", A)') trim(spec_type)
                error stop
            end if
        end select
    end subroutine get_res_filename
    
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

    !> @brief Helper to construct the filename for the RAW ASCII spectral files.
    subroutine get_ascii_spectral_filename(ctx, spec_type, z_idx, file_path)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: spec_type
        integer, intent(in) :: z_idx
        character(len=*), intent(out) :: file_path
        
        character(len=6) :: z_str
        
        write(z_str, '(F6.4)') ctx%state%zlegendinit(z_idx)

        select case (trim(spec_type))
        case ('basel')
            file_path = trim(ctx%sps_home) // '/data/spectra/BaSeL3.1/basel_' // &
                        trim(BASEL_STR) // '_z' // z_str // '.spectra'
        
        case ('miles')
            file_path = trim(ctx%sps_home) // '/data/spectra/MILES/imiles_z' // &
                        z_str // '.spectra'
        
        case default
            if (spec_type(1:3) == 'c3k') then
                file_path = trim(ctx%sps_home) // '/data/spectra/C3K/' // &
                            trim(spec_type) // '_z' // z_str // '.spectra'
            else if (trim(spec_type) == 'ckc14') then
                 ! [cite: 175] CKC14 support
                file_path = trim(ctx%sps_home) // '/data/spectra/CKC14/' // &
                            trim(spec_type) // '_z' // z_str // '.spectra'
            else
                write(error_unit, '("[FSPS_IO] Unknown spec_type: ", A)') trim(spec_type)
                error stop
            end if
        end select
    end subroutine get_ascii_spectral_filename

    !> @brief Helper to read the specific WMBasic format
    subroutine load_wmb_raw_z(ctx, z_str, w_lam, w_out)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(in) :: z_str
        real(WP), dimension(:), intent(out) :: w_lam
        real(WP), dimension(:,:,:), intent(out) :: w_out ! (nspec, logt, logg)
        
        integer :: u_spec, io, k
        character(len=1024) :: f_path

        f_path = trim(ctx%sps_home) // '/data/spectra/Hot_spectra/WMBASIC_z' // trim(z_str) // '.spec'
        open(newunit=u_spec, file=trim(f_path), status='old', action='read', iostat=io)
        if (io /= 0) call handle_read_error(f_path, 0)

        do k = 1, NSPEC_WMB
            read(u_spec, *) w_lam(k), w_out(k, :, 1), w_out(k, :, 2), w_out(k, :, 3)
        end do
        close(u_spec)
    end subroutine load_wmb_raw_z

    !> @brief Finds the start and end indices in the sorted array `grid` 
    !> that bracket the range [val_min, val_max].
    subroutine find_bounds(grid, val_min, val_max, start_idx, end_idx)
        real(WP), dimension(:), intent(in) :: grid
        real(WP), intent(in) :: val_min, val_max
        integer, intent(out) :: start_idx, end_idx
        
        integer :: n
        logical :: is_ascending

        n = size(grid)
        is_ascending = (grid(n) >= grid(1))
        
        ! Find first point. find_interval returns i where grid(i) <= val.
        ! We start loosely at the bracket lower bound.
        if (is_ascending) then
            if (val_min <= grid(1)) then
                start_idx = 1
            else
                start_idx = max(find_interval(grid, val_min), 1)
            end if

            if (val_max >= grid(n)) then
                end_idx = n
            else
                end_idx = min(find_interval(grid, val_max), n)
            end if
        else
            if (val_min >= grid(1)) then
                start_idx = 1
            else
                start_idx = max(find_interval(grid, val_min), 1)
            end if

            if (val_max <= grid(n)) then
                end_idx = n
            else
                end_idx = min(find_interval(grid, val_max), n)
            end if
        end if
        
        ! Ensure we don't cross over if the grid is coarse
        if (start_idx > end_idx) then 
           start_idx = 1
           end_idx = 0
        end if
        
    end subroutine find_bounds

    !> @brief Returns .true. if the environment variable is set to a truthy value.
    logical function env_flag_true(var_name)
        character(len=*), intent(in) :: var_name
        character(len=32) :: buf
        integer :: stat

        env_flag_true = .false.
        call get_environment_variable(var_name, value=buf, status=stat)
        if (stat /= 0) return
        if (len_trim(buf) == 0) return
        select case (buf(1:1))
        case ('1','t','T','y','Y')
            env_flag_true = .true.
        case default
            env_flag_true = .false.
        end select
    end function env_flag_true

    !> @brief Prints basic statistics for a 1D array.
    subroutine debug_stats_1d(label, arr)
        character(len=*), intent(in) :: label
        real(WP), dimension(:), intent(in) :: arr
        real(WP) :: amin, amax
        integer :: n

        n = size(arr)
        if (n <= 0) return
        amin = minval(arr)
        amax = maxval(arr)
        write(*,'(A,1x,I0,1x,ES13.5,1x,ES13.5,1x,ES13.5,1x,ES13.5)') &
            trim(label), n, amin, amax, arr(1), arr(n)
    end subroutine debug_stats_1d

end module fsps_io