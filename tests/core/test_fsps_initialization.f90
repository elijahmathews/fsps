module test_fsps_initialization_mod
    use fsps_precision, only: WP
    use fsps_context_types, only: fsps_context_t
    use fsps_types, only: params
    use fsps_api, only: fsps_create, fsps_setup, fsps_destroy
    use test_utils_mod, only: print_group, print_summary_line, print_minor_header, assert_true, assert_int_equals

    implicit none

    integer :: total_failures = 0
    integer :: total_tests = 0

    public :: run_fsps_initialization_tests, total_failures, total_tests

contains

    subroutine run_fsps_initialization_tests()
        call print_minor_header("fsps_initialization")

        call test_fsps_setup_mist_miles_hdf5()
        call test_fsps_setup_pdva_basel_hdf5()
        call test_fsps_setup_name_normalization_hdf5()
        call test_fsps_setup_default_aliases_hdf5()

        call print_summary_line("Module Summary", total_tests - total_failures, total_tests)
    end subroutine run_fsps_initialization_tests

    subroutine test_fsps_setup_mist_miles_hdf5()
        type(fsps_context_t) :: ctx
        logical :: have_hdf5
        character(len=1024) :: h5_path
        integer :: st

        call print_group("fsps_setup smoke: mist+miles+DL07")

        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=h5_path, status=st)
        have_hdf5 = (st == 0 .and. len_trim(h5_path) > 0)
        if (have_hdf5) then
            inquire(file=trim(h5_path), exist=have_hdf5)
        end if

        if (.not. have_hdf5) then
            call assert_true(.true., 'FSPS_HDF5_DATA_PATH not set/found; skipping setup smoke test', total_tests, total_failures)
            return
        end if

        call fsps_create(ctx)
        call fsps_setup(ctx, 10, isoc_type_in='mist', spec_type_in='miles', dust_type_in='DL07')

        call assert_true(ctx%initialized, 'context marked initialized', total_tests, total_failures)
        call assert_int_equals(1, ctx%state%check_sps_setup, 'check_sps_setup set', total_tests, total_failures)
        call assert_true(trim(ctx%state%isoc_type) == 'mist', 'isoc_type resolved to mist', total_tests, total_failures)
        call assert_true(trim(ctx%state%spec_type) == 'miles', 'spec_type resolved to miles', total_tests, total_failures)
        call assert_true(trim(ctx%state%str_dustem) == 'DL07', 'dust model resolved to DL07', total_tests, total_failures)
        call assert_true(ctx%state%nspec > 0, 'nspec populated', total_tests, total_failures)
        call assert_true(ctx%state%nt > 0, 'nt populated', total_tests, total_failures)
        call assert_true(ctx%state%nz > 0, 'nz populated', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_fsps_setup_mist_miles_hdf5

    subroutine test_fsps_setup_pdva_basel_hdf5()
        type(fsps_context_t) :: ctx
        logical :: have_hdf5
        character(len=1024) :: h5_path
        integer :: st

        call print_group("fsps_setup smoke: pdva+basel+THEMIS")

        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=h5_path, status=st)
        have_hdf5 = (st == 0 .and. len_trim(h5_path) > 0)
        if (have_hdf5) then
            inquire(file=trim(h5_path), exist=have_hdf5)
        end if

        if (.not. have_hdf5) then
            call assert_true(.true., 'FSPS_HDF5_DATA_PATH not set/found; skipping setup smoke test', total_tests, total_failures)
            return
        end if

        call fsps_create(ctx)
        call fsps_setup(ctx, 6, isoc_type_in='pdva', spec_type_in='basel', dust_type_in='THEMIS')

        call assert_true(ctx%initialized, 'context marked initialized', total_tests, total_failures)
        call assert_int_equals(1, ctx%state%check_sps_setup, 'check_sps_setup set', total_tests, total_failures)
        call assert_true(trim(ctx%state%isoc_type) == 'pdva', 'isoc_type resolved to pdva', total_tests, total_failures)
        call assert_true(trim(ctx%state%spec_type) == 'basel', 'spec_type resolved to basel', total_tests, total_failures)
        call assert_true(trim(ctx%state%str_dustem) == 'THEMIS', 'dust model resolved to THEMIS', total_tests, total_failures)
        call assert_true(ctx%state%nspec > 0, 'nspec populated', total_tests, total_failures)
        call assert_true(ctx%state%nt > 0, 'nt populated', total_tests, total_failures)
        call assert_true(ctx%state%nz > 0, 'nz populated', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_fsps_setup_pdva_basel_hdf5

    subroutine test_fsps_setup_name_normalization_hdf5()
        type(fsps_context_t) :: ctx
        logical :: have_hdf5
        character(len=1024) :: h5_path
        integer :: st

        call print_group("fsps_setup canonicalization: aliases + casing")

        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=h5_path, status=st)
        have_hdf5 = (st == 0 .and. len_trim(h5_path) > 0)
        if (have_hdf5) then
            inquire(file=trim(h5_path), exist=have_hdf5)
        end if

        if (.not. have_hdf5) then
            call assert_true(.true., 'FSPS_HDF5_DATA_PATH not set/found; skipping normalization test', total_tests, total_failures)
            return
        end if

        call fsps_create(ctx)
        call fsps_setup(ctx, 6, isoc_type_in='padova', spec_type_in='BaSeL3.1', dust_type_in='themis')

        call assert_true(ctx%initialized, 'context initialized after canonicalized setup', total_tests, total_failures)
        call assert_true(trim(ctx%state%isoc_type) == 'pdva', 'padova alias normalized to pdva', total_tests, total_failures)
        call assert_true(trim(ctx%state%spec_type) == 'basel', 'BaSeL3.1 normalized to basel', total_tests, total_failures)
        call assert_true(trim(ctx%state%str_dustem) == 'THEMIS', &
                 'lowercase themis normalized to THEMIS', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_fsps_setup_name_normalization_hdf5

    subroutine test_fsps_setup_default_aliases_hdf5()
        type(fsps_context_t) :: ctx
        logical :: have_hdf5
        character(len=1024) :: h5_path
        integer :: st

        call print_group("fsps_setup canonicalization: default/padova2007 + basel2.2 + dust fallback")

        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=h5_path, status=st)
        have_hdf5 = (st == 0 .and. len_trim(h5_path) > 0)
        if (have_hdf5) then
            inquire(file=trim(h5_path), exist=have_hdf5)
        end if

        if (.not. have_hdf5) then
            call assert_true(.true., 'FSPS_HDF5_DATA_PATH not set/found; skipping alias fallback test', total_tests, total_failures)
            return
        end if

        call fsps_create(ctx)
        call fsps_setup(ctx, 6, isoc_type_in='padova2007', spec_type_in='basel2.2', dust_type_in='invalid_dust')

        call assert_true(ctx%initialized, 'context initialized for alias/fallback setup', total_tests, total_failures)
        call assert_true(trim(ctx%state%isoc_type) == 'pdva', 'padova2007 normalized to pdva', total_tests, total_failures)
        call assert_true(trim(ctx%state%spec_type) == 'basel', 'basel2.2 normalized to basel', total_tests, total_failures)
        call assert_true(trim(ctx%state%str_dustem) == 'DL07', 'invalid dust value falls back to DL07', total_tests, total_failures)
        call assert_int_equals(1001, ctx%state%ndim_dustem, 'DL07 ndim_dustem selected', total_tests, total_failures)
        call assert_int_equals(22, ctx%state%numin_dustem, 'DL07 numin_dustem selected', total_tests, total_failures)
        call assert_int_equals(7, ctx%state%nqpah_dustem, 'DL07 nqpah_dustem selected', total_tests, total_failures)

        call fsps_destroy(ctx)
    end subroutine test_fsps_setup_default_aliases_hdf5

end module test_fsps_initialization_mod
