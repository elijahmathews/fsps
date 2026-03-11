module test_fsps_hdf5_fixture_support_mod
    use fsps_precision, only: WP
    use hdf5
    implicit none

    public :: create_test_hdf5_fixture
    public :: remove_test_hdf5_fixture
    public :: make_backend_uri

contains

    subroutine create_test_hdf5_fixture(file_path, ok)
        character(len=:), allocatable, intent(out) :: file_path
        logical, intent(out) :: ok

        integer(HID_T) :: file_id
        integer :: hdferr
        integer :: clk
        real(WP) :: r
        integer :: suffix
        character(len=256) :: tmp_path

        real(WP), dimension(3) :: axis_lambda
        real(WP), dimension(2) :: axis_z
        real(WP), dimension(2) :: axis_logt
        real(WP), dimension(2) :: axis_logg
        real(WP), dimension(2) :: axis_afe
        real(WP), dimension(2) :: axis_age
        real(WP), dimension(3) :: axis_nebular_lambda

        real(WP), dimension(3) :: r1
        real(WP), dimension(2,3) :: r2
        real(WP), dimension(2,2,2) :: r3
        real(WP), dimension(2,2,2,2) :: r4

        integer, dimension(3) :: i1
        integer, dimension(2,2) :: i2
        integer, dimension(2,2,2) :: i3

        real(WP), dimension(3,2,2,2,2) :: spec5
        real(WP), dimension(3,2,2,2) :: spec4
        integer, dimension(2,2,2,2) :: mask4
        integer, dimension(2,2,2) :: mask3
        integer :: il, iz, ia, it, ig

        ok = .true.

        axis_lambda = [100.0_wp, 200.0_wp, 300.0_wp]
        axis_z = [0.01_wp, 0.02_wp]
        axis_logt = [3.5_wp, 3.6_wp]
        axis_logg = [4.0_wp, 4.5_wp]
        axis_afe = [0.0_wp, 0.2_wp]
        axis_age = [1.0_wp, 5.0_wp]
        axis_nebular_lambda = [95.0_wp, 205.0_wp, 305.0_wp]

        r1 = [1.0_wp, 2.0_wp, 3.0_wp]
        r2 = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp], [2,3])
        r3 = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp], [2,2,2])
        r4 = reshape([1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, &
                      9.0_wp, 10.0_wp, 11.0_wp, 12.0_wp, 13.0_wp, 14.0_wp, 15.0_wp, 16.0_wp], [2,2,2,2])

        i1 = [4, 5, 6]
        i2 = reshape([1, 2, 3, 4], [2,2])
        i3 = reshape([1, 2, 3, 4, 5, 6, 7, 8], [2,2,2])

        do ig = 1, 2
            do it = 1, 2
                do ia = 1, 2
                    do iz = 1, 2
                        do il = 1, 3
                            spec5(il, iz, ia, it, ig) = real(10000*il + 1000*iz + 100*ia + 10*it + ig, WP)
                        end do
                    end do
                end do
            end do
        end do

        do ig = 1, 2
            do it = 1, 2
                do iz = 1, 2
                    do il = 1, 3
                        spec4(il, iz, it, ig) = real(10000*il + 1000*iz + 10*it + ig, WP)
                    end do
                end do
            end do
        end do

        mask4 = 1
        mask4(2, 1, 2, 2) = 0
        mask3 = 1
        mask3(2, 2, 1) = 0

        call system_clock(count=clk)
        call random_number(r)
        suffix = int(r * 1000000.0_wp) + abs(clk)

        write(tmp_path, '(A,I0,A)') '/tmp/fsps_hdf5_backend_fixture_', suffix, '.h5'
        file_path = trim(tmp_path)

        call h5open_f(hdferr)
        if (hdferr < 0) then
            ok = .false.
            return
        end if

        call h5fcreate_f(trim(file_path), H5F_ACC_TRUNC_F, file_id, hdferr)
        if (hdferr < 0) then
            ok = .false.
            call h5close_f(hdferr)
            return
        end if

        call write_real_1d(file_id, '/libraries/spectra/miles/axes/lambda', axis_lambda, ok)
        call write_real_1d(file_id, '/libraries/spectra/miles/axes/z', axis_z, ok)
        call write_real_1d(file_id, '/libraries/spectra/miles/axes/logt', axis_logt, ok)
        call write_real_1d(file_id, '/libraries/spectra/miles/axes/logg', axis_logg, ok)
        call write_real_1d(file_id, '/libraries/spectra/miles/axes/afe', axis_afe, ok)
        call write_real_1d(file_id, '/libraries/spectra/miles/axes/alpha_fe', axis_afe, ok)

        call write_real_1d(file_id, '/axes/lambda', axis_lambda, ok)
        call write_real_1d(file_id, '/axes/z', axis_z, ok)
        call write_real_1d(file_id, '/axes/logt', axis_logt, ok)
        call write_real_1d(file_id, '/axes/logg', axis_logg, ok)
        call write_real_1d(file_id, '/axes/age', axis_age, ok)
        call write_real_1d(file_id, '/axes/nebular_lambda', axis_nebular_lambda, ok)

        call write_real_1d(file_id, '/libraries/nebular/mist/WD/line_pos', [4861.0_wp, 6563.0_wp], ok)
        call write_real_1d(file_id, '/libraries/nebular/WD/line_pos', [4861.0_wp, 6563.0_wp], ok)

        call write_real_1d(file_id, '/tests/r1', r1, ok)
        call write_real_2d(file_id, '/tests/r2', r2, ok)
        call write_real_3d(file_id, '/tests/r3', r3, ok)
        call write_real_4d(file_id, '/tests/r4', r4, ok)

        call write_int_1d(file_id, '/tests/i1', i1, ok)
        call write_int_2d(file_id, '/tests/i2', i2, ok)
        call write_int_3d(file_id, '/tests/i3', i3, ok)

        call write_real_5d(file_id, '/tests/spec5', spec5, ok)
        call write_real_4d(file_id, '/tests/spec4', spec4, ok)
        call write_int_4d(file_id, '/tests/spec5_valid', mask4, ok)
        call write_int_3d(file_id, '/tests/spec4_valid', mask3, ok)

        call h5fclose_f(file_id, hdferr)
        if (hdferr /= 0) ok = .false.

        call h5close_f(hdferr)
        if (hdferr /= 0) ok = .false.
    end subroutine create_test_hdf5_fixture

    subroutine remove_test_hdf5_fixture(file_path)
        character(len=*), intent(in) :: file_path
        integer :: u, ios

        open(newunit=u, file=trim(file_path), status='old', action='readwrite', iostat=ios)
        if (ios == 0) then
            close(u, status='delete', iostat=ios)
        end if
    end subroutine remove_test_hdf5_fixture

    function make_backend_uri(file_path, isoc_type, spec_type, dust_type) result(uri)
        character(len=*), intent(in) :: file_path
        character(len=*), intent(in) :: isoc_type
        character(len=*), intent(in) :: spec_type
        character(len=*), intent(in) :: dust_type
        character(len=:), allocatable :: uri

        uri = trim(file_path)//'|'//trim(isoc_type)//'|'//trim(spec_type)//'|'//trim(dust_type)
    end function make_backend_uri

    subroutine write_real_1d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        real(WP), intent(in) :: values(:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values), HSIZE_T)]
        call h5screate_simple_f(1, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), get_h5_real_type(), space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, get_h5_real_type(), values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_1d

    subroutine write_real_2d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        real(WP), intent(in) :: values(:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T)]
        call h5screate_simple_f(2, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), get_h5_real_type(), space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, get_h5_real_type(), values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_2d

    subroutine write_real_3d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        real(WP), intent(in) :: values(:,:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(3) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T), int(size(values,3), HSIZE_T)]
        call h5screate_simple_f(3, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), get_h5_real_type(), space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, get_h5_real_type(), values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_3d

    subroutine write_real_4d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        real(WP), intent(in) :: values(:,:,:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(4) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T), &
            int(size(values,3), HSIZE_T), int(size(values,4), HSIZE_T)]
        call h5screate_simple_f(4, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), get_h5_real_type(), space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, get_h5_real_type(), values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_4d

    subroutine write_real_5d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        real(WP), intent(in) :: values(:,:,:,:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(5) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T), int(size(values,3), HSIZE_T), &
                int(size(values,4), HSIZE_T), int(size(values,5), HSIZE_T)]
        call h5screate_simple_f(5, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), get_h5_real_type(), space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, get_h5_real_type(), values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_5d

    subroutine write_int_1d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        integer, intent(in) :: values(:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values), HSIZE_T)]
        call h5screate_simple_f(1, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_1d

    subroutine write_int_2d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        integer, intent(in) :: values(:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T)]
        call h5screate_simple_f(2, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_2d

    subroutine write_int_3d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        integer, intent(in) :: values(:,:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(3) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T), int(size(values,3), HSIZE_T)]
        call h5screate_simple_f(3, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_3d

    subroutine write_int_4d(file_id, path, values, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        integer, intent(in) :: values(:,:,:,:)
        logical, intent(inout) :: ok

        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(4) :: dims
        integer :: hdferr

        call ensure_parent_groups(file_id, path, ok)
        if (.not. ok) return

        dims = [int(size(values,1), HSIZE_T), int(size(values,2), HSIZE_T), &
            int(size(values,3), HSIZE_T), int(size(values,4), HSIZE_T)]
        call h5screate_simple_f(4, dims, space_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            return
        end if
        call h5dcreate_f(file_id, trim(path), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        if (hdferr /= 0) then
            ok = .false.
            call h5sclose_f(space_id, hdferr)
            return
        end if
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, values, dims, hdferr)
        if (hdferr /= 0) ok = .false.
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_4d

    subroutine ensure_parent_groups(file_id, path, ok)
        integer(HID_T), intent(in) :: file_id
        character(len=*), intent(in) :: path
        logical, intent(inout) :: ok

        integer :: n, p_last, p_start, p_rel
        logical :: exists
        character(len=:), allocatable :: parent, token, current
        integer(HID_T) :: group_id
        integer :: hdferr

        if (.not. ok) return

        n = len_trim(path)
        if (n <= 1) return

        p_last = scan(trim(path), '/', back=.true.)
        if (p_last <= 1) return

        parent = path(1:p_last-1)
        p_start = 2
        current = ''

        do while (p_start <= len_trim(parent))
            p_rel = index(parent(p_start:), '/')
            if (p_rel == 0) then
                token = trim(parent(p_start:len_trim(parent)))
                p_start = len_trim(parent) + 1
            else
                token = trim(parent(p_start:p_start+p_rel-2))
                p_start = p_start + p_rel
            end if

            if (len_trim(token) == 0) cycle
            current = trim(current)//'/'//trim(token)

            call h5lexists_f(file_id, trim(current), exists, hdferr)
            if (hdferr /= 0) then
                ok = .false.
                return
            end if

            if (.not. exists) then
                call h5gcreate_f(file_id, trim(current), group_id, hdferr)
                if (hdferr /= 0) then
                    ok = .false.
                    return
                end if
                call h5gclose_f(group_id, hdferr)
                if (hdferr /= 0) then
                    ok = .false.
                    return
                end if
            end if
        end do
    end subroutine ensure_parent_groups

    pure function get_h5_real_type() result(h5_type)
        integer(HID_T) :: h5_type

        if (kind(1.0_wp) == kind(1.0d0)) then
            h5_type = H5T_NATIVE_DOUBLE
        else
            h5_type = H5T_NATIVE_REAL
        end if
    end function get_h5_real_type

end module test_fsps_hdf5_fixture_support_mod