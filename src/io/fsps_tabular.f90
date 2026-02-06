module fsps_tabular

    use fsps_precision, only: WP
    use fsps_constants, only: SAFE_FLOOR, NTABMAX
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params

    implicit none
    private

    public :: load_tabular_sfh

contains

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

end module fsps_tabular