#include "fsps_build_config.h"

module fsps_initialization
    !> @brief
    !> Handles the initialization, memory allocation, and data loading for FSPS contexts.
    !>
    !> @details
    !> This module replaces the legacy `sps_setup.f90`. It orchestrates the setup process:
    !> 1. Determines array dimensions based on selected physics models (MIST vs Padova, MILES vs BaSeL).
    !> 2. Manages the shared data cache (allocating new memory or linking to existing caches).
    !> 3. Calls `fsps_io` to load data from disk if the cache is cold.
    !> 4. Pre-calculates derived physical quantities (redshift grids, cosmology).

    use iso_fortran_env, only: error_unit
    use fsps_precision, only: WP
    use fsps_constants, only: VERBOSE, TIME_RES_INCR, SAFE_FLOOR, C_LIGHT, L_SOL, PI, &
                              NM, NDIM_LOGT, NDIM_LOGG, NDIM_WMB_LOGT, NDIM_WMB_LOGG, &
                              NDIM_PAGB, NDIM_WR, N_AGB_O, N_AGB_C, N_AGB_CAR, &
                              NTAU_DAGB, NTEFF_DAGB, NEMLINE, NEBNZ, NEBNAGE, NEBNIP, &
                              NAGNDUST, NTABMAX, NT_XRB, NZ_XRB, NSPEC_XRB, &
                              NT_MIST, NZ_MIST, NT_PADOVA, NZ_PADOVA, NT_PARSEC, NZ_PARSEC, &
                              NT_BASTI, NZ_BASTI, NT_GENEVA, NZ_GENEVA
    use fsps_strings, only: to_lower
    use fsps_context_types, only: fsps_context_t
    use fsps_cache, only: fsps_setup_cache_t, fsps_cache_get_setup
    use fsps_environment, only: fsps_resolve_paths, fsps_cleanup
    use fsps_io, only: load_zlegend_file, load_wavelength_grid, load_spectral_resolution, &
                       read_bpass_data, read_isochrone_database_legacy => read_isochrone_database, &
                       read_spectral_binary_legacy => read_spectral_binary, &
                       load_nebular_grid_legacy => load_nebular_grid, &
                       load_dust_emission_table_legacy => load_dust_emission_table, &
                       load_agn_dust_models_legacy => load_agn_dust_models, &
                       load_filter_definitions, &
                       load_standard_sed, load_index_definitions, &
                       load_lsf_data, &
                       load_attenuation_curves_legacy => load_attenuation_curves, &
                       load_wmbasic_spectra_legacy => load_wmbasic_spectra, &
                       load_agb_spectra_legacy => load_agb_spectra, &
                       load_post_agb_spectra_legacy => load_post_agb_spectra, &
                       load_wr_spectra_legacy => load_wr_spectra, &
                       apply_legacy_filter_norm
    use fsps_data_backend, only: backend_status_t, data_backend_t, backend_status_ok
    use fsps_data_registry, only: create_data_backend
    use fsps_data_mapper, only: fsps_data_mapper_t
    use fsps_data_schema, only: spectral_grid_t, isochrone_grid_t, nebular_grid_t, aux_wmbasic_t, aux_pagb_t, aux_wr_t, aux_agb_t, &
                                dust_emission_t, agn_dust_t, dust_attenuation_t, xrb_spectra_t, library_manifest_t, dataset_desc_t
    use fsps_cosmology, only: get_universe_age, get_luminosity_distance
    use fsps_interpolation, only: find_interval, interpolate_linear
    use fsps_integration, only: integrate_trapezoid_array
    

    implicit none
    private

    ! Only the main driver is public
    public :: fsps_initialize_data

    ! ---------------------------------------------------------------------
    ! Module constants
    ! ---------------------------------------------------------------------
    integer, parameter :: NZWMB = 12
    class(data_backend_t), allocatable, save :: fsps_backend
    type(fsps_data_mapper_t), save :: fsps_mapper
    type(library_manifest_t), save :: fsps_manifest
    logical, save :: fsps_backend_open = .false.

contains

    !> @brief Main entry point for initializing an FSPS context.
    !>
    !> @details
    !> This routine coordinates the entire setup lifecycle. It is the only routine
    !> that `fsps_api` needs to call.
    !>
    !> @param[inout] ctx             The context to initialize.
    !> @param[in]    zin             Specific metallicity index to load (-1 for all).
    !> @param[in]    isoc_type_in    (Optional) Isochrone library string.
    !> @param[in]    spec_type_in    (Optional) Spectral library string.
    !> @param[in]    dust_type_in    (Optional) Dust emission model string.
    subroutine fsps_initialize_data(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: zin
        character(len=*), intent(in), optional :: isoc_type_in
        character(len=*), intent(in), optional :: spec_type_in
        character(len=*), intent(in), optional :: dust_type_in
        logical :: cache_is_new
        character(len=512) :: cache_key
        type(fsps_setup_cache_t), pointer :: cache_ptr
        character(len=64) :: isoc_name, spec_name, dust_name

        call fsps_cleanup(ctx)
        ctx%zin = zin

        call normalize_library_names(isoc_type_in, spec_type_in, dust_type_in, isoc_name, spec_name, dust_name)
        ctx%isoc_type_name = trim(isoc_name)
        ctx%spec_type_name = trim(spec_name)
        ctx%dust_type_name = trim(dust_name)
        ctx%state%isoc_type = trim(isoc_name)
        ctx%state%spec_type = trim(spec_name)
        ctx%state%str_dustem = trim(dust_name)

        call set_simulation_dimensions(ctx)
        call configure_dust_model(ctx)

        if (zin > ctx%state%nz) then
            write (error_unit, '(A,1x,I0,1x,I0)') '[FSPS_INIT] Error: zin > nz', zin, ctx%state%nz
            error stop 1
        end if

        call fsps_resolve_paths(ctx)

        call build_cache_key(ctx, zin, cache_key)
        call fsps_cache_get_setup(cache_key, ctx%setup_cache, cache_is_new)
        cache_ptr => ctx%setup_cache

        if (cache_is_new) then
            call populate_cache_metadata(ctx, cache_ptr)
            call allocate_shared_memory(cache_ptr, ctx)
        end if

        call link_context_to_cache(ctx, cache_ptr)
        call allocate_local_workspace(ctx)

        if (cache_is_new) then
            call initialize_cache_arrays(ctx)
        end if

        call load_stellar_data(ctx, zin)
        call load_interstellar_physics(ctx)
        call load_photometry_data(ctx)
        call compute_filter_leff(ctx)
        call load_agn_dust_models(ctx)
        call load_index_definitions(ctx)
        call build_time_grid(ctx)
        if (ctx%smooth_lsf_val == 1) call load_lsf_data(ctx)

        call compute_local_properties(ctx)

        ctx%state%check_sps_setup = 1
        ctx%initialized = .true.

        if (VERBOSE == 1) then
            write (*, '(A)') '      ...done'
        end if

    end subroutine fsps_initialize_data

    ! ========================================================================
    ! PRIVATE HELPER SUBROUTINES
    ! ========================================================================

    !> @brief Determines grid dimensions (NT, NZ, NSPEC) based on library selection.
    !>
    !> @details
    !> Sets the fundamental array shapes for the simulation based on the selected
    !> Isochrone (`isoc_type`) and Spectral Library (`spec_type`).
    !>
    !> Logic corresponds to legacy `sps_setup` lines ~749-756.
    !>
    !> **Dimensions Calculated:**
    !> - `nt`: Number of age steps in the isochrone.
    !> - `nz`: Number of metallicity steps in the isochrone.
    !> - `nspec`: Number of wavelength points in the spectral library.
    !> - `ntfull`: The expanded time grid size (including interpolation).
    !>
    !> @param[inout] ctx The context to configure.
    subroutine set_simulation_dimensions(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=64) :: iso_str, spec_str

        iso_str = trim(ctx%state%isoc_type)
        spec_str = trim(ctx%state%spec_type)

        select case (iso_str)
        case ('mist')
            ctx%state%zsol = 0.0142_wp
            ctx%state%nt = NT_MIST
            ctx%state%nz = NZ_MIST
        case ('pdva')
            ctx%state%zsol = 0.019_wp
            ctx%state%nt = NT_PADOVA
            ctx%state%nz = NZ_PADOVA
        case ('prsc')
            ctx%state%zsol = 0.01524_wp
            ctx%state%nt = NT_PARSEC
            ctx%state%nz = NZ_PARSEC
        case ('bsti')
            ctx%state%zsol = 0.020_wp
            ctx%state%nt = NT_BASTI
            ctx%state%nz = NZ_BASTI
        case ('gnva')
            ctx%state%zsol = 0.020_wp
            ctx%state%nt = NT_GENEVA
            ctx%state%nz = NZ_GENEVA
        case ('bpss')
            ctx%state%zsol = 0.020_wp
            ctx%state%nt = 43
            ctx%state%nz = 12
        case default
            write (error_unit, '(A,1x,A)') '[FSPS_INIT] Error: Unknown isoc_type', trim(iso_str)
            error stop 1
        end select

        select case (spec_str)
        case ('miles')
            ctx%state%zsol_spec = 0.019_wp
            ctx%state%nzinit = 5
            ctx%state%nspec = 5994
        case ('basel')
            ctx%state%zsol_spec = 0.020_wp
            ctx%state%nzinit = 6
            ctx%state%nspec = 1963
        case ('bpass')
            ctx%state%zsol_spec = 0.020_wp
            ctx%state%nzinit = 1
            ctx%state%nspec = 15000
        case default
            if (len_trim(spec_str) >= 3 .and. spec_str(1:3) == 'c3k') then
                ctx%state%zsol_spec = 0.0134_wp
                ctx%state%nzinit = 11
                ctx%state%nspec = 11149
            else
                write (error_unit, '(A,1x,A)') '[FSPS_INIT] Error: Unknown spec_type', trim(spec_str)
                error stop 1
            end if
        end select

        ctx%state%nbands = 159
        ctx%state%nindx = 30
        ctx%state%ntfull = TIME_RES_INCR*ctx%state%nt

        ctx%state%nspec_xrb = NSPEC_XRB
        ctx%state%nt_xrb = NT_XRB
        ctx%state%nz_xrb = NZ_XRB

        if (ctx%state%isoc_type == 'bpss' .and. TIME_RES_INCR /= 1) then
            write (error_unit, '(A)') '[FSPS_INIT] Error: TIME_RES_INCR>1 not supported for BPASS.'
            error stop 1
        end if

    end subroutine set_simulation_dimensions

    !> @brief Allocates the massive shared arrays within the cache object.
    !>
    !> @details
    !> Uses dimensions stored in `ctx%state` (NT, NZ, NSPEC) and global constants
    !> to allocate the shared physics grids.
    !>
    !> **Allocation Strategy:**
    !> 1. **Isochrones:** Allocated using `NT` and `NZ`.
    !> 2. **Auxiliary Stars:** Allocated using `NSPEC` and `fsps_constants` (WMB, AGB, WR).
    !> 3. **Interstellar Physics:** Allocated using `fsps_constants` (Nebular, Dust).
    !> 4. **Spectral Library:** The main `speclib` is *not* allocated here; it is
    !>    deferred to `read_spectral_binary` which reads the grid dimensions (T, g)
    !>    directly from the file header.
    !>
    !> @param[inout] cache The cache object to allocate.
    !> @param[in]    ctx   The context containing dimension info (nt, nz, nspec).
    subroutine allocate_shared_memory(cache, ctx)
        type(fsps_setup_cache_t), intent(inout) :: cache
        type(fsps_context_t), intent(in) :: ctx

        integer :: nt, nz, nspec, ntfull
        integer :: nspec_xrb, nt_xrb, nz_xrb
        integer :: ndim_dustem, numin_dustem, nqpah_dustem

        nt = ctx%state%nt
        nz = ctx%state%nz
        nspec = ctx%state%nspec
        ntfull = ctx%state%ntfull
        nspec_xrb = ctx%state%nspec_xrb
        nt_xrb = ctx%state%nt_xrb
        nz_xrb = ctx%state%nz_xrb
        ndim_dustem = ctx%state%ndim_dustem
        numin_dustem = ctx%state%numin_dustem
        nqpah_dustem = ctx%state%nqpah_dustem

        allocate (cache%indexdefined(7, ctx%state%nindx))
        allocate (cache%wgdust(nspec, 18, 6, 2))
        allocate (cache%g03smcextn(nspec))
        allocate (cache%bands(nspec, ctx%state%nbands))
        allocate (cache%magsun(ctx%state%nbands), cache%magvega(ctx%state%nbands), &
                  cache%filter_leff(ctx%state%nbands))
        allocate (cache%vega_spec(nspec), cache%sun_spec(nspec))
        allocate (cache%spec_lambda(nspec), cache%spec_nu(nspec))
        allocate (cache%spec_res(nspec))
        allocate (cache%speclib(nspec, nz, NDIM_LOGT, NDIM_LOGG))
        allocate (cache%wmb_spec(nspec, nz, NDIM_WMB_LOGT, NDIM_WMB_LOGG))
        allocate (cache%agb_spec_o(nspec, N_AGB_O))
        allocate (cache%agb_logt_o(nz, N_AGB_O))
        allocate (cache%agb_spec_c(nspec, N_AGB_C))
        allocate (cache%agb_logt_c(N_AGB_C))
        allocate (cache%agb_spec_car(nspec, N_AGB_CAR))
        allocate (cache%pagb_spec(nspec, NDIM_PAGB, 2))
        allocate (cache%wrn_spec(nspec, NDIM_WR, nz), cache%wrc_spec(nspec, NDIM_WR, nz))

        allocate (cache%qpaharr(nqpah_dustem))
        allocate (cache%uminarr(numin_dustem))
        allocate (cache%lambda_dustem(ndim_dustem))
        allocate (cache%dustem_dustem(ndim_dustem, numin_dustem*2))
        allocate (cache%dustem2_dustem(nspec, nqpah_dustem, numin_dustem*2))
        allocate (cache%flux_dagb(nspec, 2, NTEFF_DAGB, NTAU_DAGB))
        allocate (cache%nebem_cont(nspec, NEBNZ, NEBNAGE, NEBNIP), &
                  cache%xnebem_cont(nspec, NEBNZ, NEBNAGE, NEBNIP))
        allocate (cache%neb_res_min(nspec))
        allocate (cache%gaussnebarr(nspec, NEMLINE))
        allocate (cache%agndust_spec(nspec, NAGNDUST))
        allocate (cache%mact_isoc(nz, nt, NM), cache%logl_isoc(nz, nt, NM), &
                  cache%logt_isoc(nz, nt, NM), cache%logg_isoc(nz, nt, NM))
        allocate (cache%ffco_isoc(nz, nt, NM), cache%phase_isoc(nz, nt, NM), &
                  cache%mini_isoc(nz, nt, NM), cache%lmdot_isoc(nz, nt, NM))
        allocate (cache%nmass_isoc(nz, nt))
        allocate (cache%timestep_isoc(nz, nt))
        allocate (cache%zlegend(nz))
        allocate (cache%zlegendinit(ctx%state%nzinit))
        allocate (cache%time_full(ntfull))
        allocate (cache%bpass_spec_ssp(nspec, nt, nz))
        allocate (cache%bpass_mass_ssp(nt, nz))
        allocate (cache%lam_xrb(nspec_xrb))
        allocate (cache%spec_xrb(nspec, nt_xrb, nz_xrb))
        allocate (cache%ages_xrb(nt_xrb))
        allocate (cache%zmet_xrb(nz_xrb))

    end subroutine allocate_shared_memory

    !> @brief Orchestrates I/O for Isochrones, Spectral Libraries, and Auxiliary Stars.
    !>
    !> @details
    !> This is the core loading routine for stellar physics. It handles:
    !> 1. **Isochrones:** Evolution tracks (MIST, Padova, etc.) via backend dataset loading.
    !> 2. **Spectral Libraries:** Base stellar spectra (MILES, BaSeL, C3K) via backend dataset loading.
    !> 3. **Auxiliary Spectra:** Special handling for O-stars (WMBasic),
    !>    Wolf-Rayet stars (CMFGEN), and AGB stars (Lancon & Wood).
    !>
    !> @param[inout] ctx   The context containing path configuration and model selections.
    !> @param[inout] cache The cache object to populate.
    subroutine load_stellar_data(ctx, zin)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: zin
        integer :: z, zmin, zmax, i, i1, j, k, nzinit, nm_data
        real(WP) :: dz, log_spec_val
        real(WP), allocatable :: speclibinit(:, :, :, :), speclib_slice(:, :, :)
        type(backend_status_t) :: io_status
        type(spectral_grid_t) :: base_grid
        type(isochrone_grid_t) :: iso_grid
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri

        call load_zlegend_file(ctx, ctx%state%isoc_type, .false.)

        if (zin <= 0) then
            zmin = 1
            zmax = ctx%state%nz
        else
            zmin = zin
            zmax = zin
        end if

        if (ctx%state%isoc_type == 'bpss') then
            call read_bpass_data(ctx)
            return
        end if

        call load_wavelength_grid(ctx, ctx%state%spec_type)
        call load_spectral_resolution(ctx, ctx%state%spec_type)
        call load_zlegend_file(ctx, ctx%state%spec_type, .true.)

        open (91, file=trim(ctx%sps_home)//'/data/spectra/BaSeL3.1/basel_logt.dat', status='old', action='read')
        do i = 1, NDIM_LOGT
            read (91, *) ctx%state%speclib_logt(i)
        end do
        close (91)

        open (91, file=trim(ctx%sps_home)//'/data/spectra/BaSeL3.1/basel_logg.dat', status='old', action='read')
        do i = 1, NDIM_LOGG
            read (91, *) ctx%state%speclib_logg(i)
        end do
        close (91)

        nzinit = ctx%state%nzinit
        allocate (speclibinit(ctx%state%nspec, nzinit, NDIM_LOGT, NDIM_LOGG))
        speclibinit = 0.0_wp

        call resolve_data_backend_uri(ctx, backend_mode, uri)

        if (trim(backend_mode) == 'legacy') then
            allocate(speclib_slice(ctx%state%nspec, NDIM_LOGT, NDIM_LOGG))
            do z = 1, nzinit
                call read_spectral_binary_legacy(ctx, ctx%state%spec_type, z, speclib_slice)
                speclibinit(:, z, :, :) = speclib_slice
            end do
            deallocate(speclib_slice)

            do z = 1, ctx%state%nz
                i1 = min(max(find_interval(log10(ctx%state%zlegendinit/ctx%state%zsol_spec), &
                                           log10(ctx%state%zlegend(z)/ctx%state%zsol)), 1), nzinit - 1)
                dz = (log10(ctx%state%zlegend(z)/ctx%state%zsol) - &
                      log10(ctx%state%zlegendinit(i1)/ctx%state%zsol_spec))/ &
                     (log10(ctx%state%zlegendinit(i1 + 1)/ctx%state%zsol_spec) - &
                      log10(ctx%state%zlegendinit(i1)/ctx%state%zsol_spec))
                dz = min(max(dz, 0.0_wp), 1.0_wp)

                do k = 1, NDIM_LOGG
                    do j = 1, NDIM_LOGT
                        do i = 1, ctx%state%nspec
                            log_spec_val = (1.0_wp - dz) * log10(speclibinit(i, i1, j, k) + SAFE_FLOOR) + &
                                           dz * log10(speclibinit(i, i1 + 1, j, k) + SAFE_FLOOR)
                            ctx%state%speclib(i, z, j, k) = 10.0_wp**log_spec_val
                        end do
                    end do
                end do
            end do

            deallocate(speclibinit)

            call load_wmbasic_spectra_legacy(ctx)
            call load_agb_spectra_legacy(ctx)
            call load_post_agb_spectra_legacy(ctx)
            call load_wr_spectra_legacy(ctx)

            do z = zmin, zmax
                call read_isochrone_database_legacy(ctx, ctx%state%isoc_type, z)
            end do

            if (ctx%state%isoc_type == 'gnva') then
                ctx%state%imf_lower_bound = minval(ctx%state%mini_isoc(zmin, 1, 1:ctx%state%nmass_isoc(zmin, 1)))*0.99_wp
            else
                ctx%state%imf_lower_bound = ctx%state%imf_lower_limit
            end if

            return
        end if

        call fsps_data_open(uri, backend_mode, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_open failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call fsps_data_load_spectral_library('spectral_base', base_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_spectral_library failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(base_grid%flux)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: base spectral grid is unallocated.'
            error stop 1
        end if

        if (size(base_grid%flux, 1) /= ctx%state%nspec .or. size(base_grid%flux, 2) /= nzinit .or. &
            size(base_grid%flux, 3) < 1 .or. size(base_grid%flux, 4) /= NDIM_LOGT .or. &
            size(base_grid%flux, 5) /= NDIM_LOGG) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: base spectral grid shape mismatch.'
            error stop 1
        end if

        speclibinit(:, :, :, :) = base_grid%flux(:, :, 1, :, :)

        do z = 1, ctx%state%nz
            i1 = min(max(find_interval(log10(ctx%state%zlegendinit/ctx%state%zsol_spec), &
                                       log10(ctx%state%zlegend(z)/ctx%state%zsol)), 1), nzinit - 1)
            dz = (log10(ctx%state%zlegend(z)/ctx%state%zsol) - &
                  log10(ctx%state%zlegendinit(i1)/ctx%state%zsol_spec))/ &
                 (log10(ctx%state%zlegendinit(i1 + 1)/ctx%state%zsol_spec) - &
                  log10(ctx%state%zlegendinit(i1)/ctx%state%zsol_spec))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            do k = 1, NDIM_LOGG
                do j = 1, NDIM_LOGT
                    do i = 1, ctx%state%nspec
                        log_spec_val = (1.0_wp - dz) * log10(speclibinit(i, i1, j, k) + SAFE_FLOOR) + &
                                       dz * log10(speclibinit(i, i1 + 1, j, k) + SAFE_FLOOR)
                        ctx%state%speclib(i, z, j, k) = 10.0_wp**log_spec_val
                    end do
                end do
            end do
        end do

        deallocate (speclibinit)

        call load_wmbasic_spectra_legacy(ctx)
        call load_agb_spectra_legacy(ctx)
        call load_post_agb_spectra_legacy(ctx)
        call load_wr_spectra_legacy(ctx)

        call fsps_data_load_isochrones(iso_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_isochrones failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(iso_grid%nmass) .or. .not. allocated(iso_grid%timestep_logyr) .or. &
            .not. allocated(iso_grid%mini) .or. .not. allocated(iso_grid%mact) .or. &
            .not. allocated(iso_grid%logl) .or. .not. allocated(iso_grid%logt) .or. &
            .not. allocated(iso_grid%logg) .or. .not. allocated(iso_grid%phase) .or. &
            .not. allocated(iso_grid%ffco) .or. .not. allocated(iso_grid%lmdot)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: isochrone grid is not fully allocated.'
            error stop 1
        end if

        if (size(iso_grid%nmass, 1) /= ctx%state%nt .or. &
            size(iso_grid%nmass, 2) /= ctx%state%nz .or. &
            size(iso_grid%timestep_logyr, 1) /= ctx%state%nt .or. &
            size(iso_grid%timestep_logyr, 2) /= ctx%state%nz .or. &
            size(iso_grid%mini, 1) > NM .or. &
            size(iso_grid%mini, 2) /= ctx%state%nt .or. &
            size(iso_grid%mini, 3) /= ctx%state%nz) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: isochrone grid shape mismatch.'
            error stop 1
        end if

        nm_data = size(iso_grid%mini, 1)

        do z = zmin, zmax
            ctx%state%nmass_isoc(z, :) = iso_grid%nmass(:, z)
            ctx%state%timestep_isoc(z, :) = real(iso_grid%timestep_logyr(:, z), kind(ctx%state%timestep_isoc))
            ctx%state%mini_isoc(z, :, :) = 0.0
            ctx%state%mact_isoc(z, :, :) = 0.0
            ctx%state%logl_isoc(z, :, :) = 0.0
            ctx%state%logt_isoc(z, :, :) = 0.0
            ctx%state%logg_isoc(z, :, :) = 0.0
            ctx%state%phase_isoc(z, :, :) = 0.0
            ctx%state%ffco_isoc(z, :, :) = 0.0
            ctx%state%lmdot_isoc(z, :, :) = -99.0
            ctx%state%mini_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%mini(:, :, z)), kind(ctx%state%mini_isoc))
            ctx%state%mact_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%mact(:, :, z)), kind(ctx%state%mact_isoc))
            ctx%state%logl_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%logl(:, :, z)), kind(ctx%state%logl_isoc))
            ctx%state%logt_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%logt(:, :, z)), kind(ctx%state%logt_isoc))
            ctx%state%logg_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%logg(:, :, z)), kind(ctx%state%logg_isoc))
            ctx%state%phase_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%phase(:, :, z)), kind(ctx%state%phase_isoc))
            ctx%state%ffco_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%ffco(:, :, z)), kind(ctx%state%ffco_isoc))
            ctx%state%lmdot_isoc(z, :, 1:nm_data) = real(transpose(iso_grid%lmdot(:, :, z)), kind(ctx%state%lmdot_isoc))
        end do

        call fsps_data_close(io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_close failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call base_grid%clear()
        call iso_grid%clear()

        if (ctx%state%isoc_type == 'gnva') then
            ctx%state%imf_lower_bound = minval(ctx%state%mini_isoc(zmin, 1, 1:ctx%state%nmass_isoc(zmin, 1)))*0.99_wp
        else
            ctx%state%imf_lower_bound = ctx%state%imf_lower_limit
        end if

    end subroutine load_stellar_data

    !> @brief Orchestrates I/O for Dust Emission, Attenuation, and Nebular grids.
    !>
    !> @details
    !> Calls `fsps_io` routines to populate the shared cache with:
    !> 1. **Dust Emission:** Draine & Li (2007) or THEMIS IR models.
    !> 2. **Nebular Emission:** CLOUDY model grids for HII regions (lines and continuum).
    !> 3. **AGN Dust:** Nenkova et al. (2008) CLUMPY torus models.
    !> 4. **Attenuation:** Tabulated extinction curves (Witt & Gordon, SMC).
    !> 5. **X-Ray Binaries:** Theoretical XRB spectral grids.
    !>
    !> @param[inout] ctx   The context containing path configuration and model selections.
    !> @param[inout] cache The cache object to populate.
    subroutine load_interstellar_physics(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        type(backend_status_t) :: io_status
        type(nebular_grid_t) :: neb_grid
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri

            call load_dust_emission_table_legacy(ctx, ctx%state%str_dustem)
            call load_attenuation_curves_legacy(ctx)
            call load_dusty_agb_spectra(ctx)

        if (ctx%state%isoc_type == 'mist' .or. ctx%state%isoc_type == 'pdva' .or. &
            ctx%state%isoc_type == 'prsc' .or. ctx%state%isoc_type == 'bpss') then
            call resolve_data_backend_uri(ctx, backend_mode, uri)

            if (trim(backend_mode) == 'legacy') then
                call load_nebular_grid_legacy(ctx, ctx%state%isoc_type, ctx%cloudy_dust_val == 1)
                call compute_nebular_kernels(ctx)
            else

                call fsps_data_open(uri, backend_mode, io_status)
                if (io_status%code /= 0) then
                    write(error_unit, '(A,1x,I0,1x,A)') &
                        '[FSPS_INIT] Error: fsps_data_open failed for nebular load', io_status%code, trim(io_status%message)
                    error stop 1
                end if

                if (ctx%cloudy_dust_val == 1) then
                    call fsps_data_load_nebular('WD', neb_grid, io_status)
                else
                    call fsps_data_load_nebular('ND', neb_grid, io_status)
                end if
                if (io_status%code /= 0) then
                    write(error_unit, '(A,1x,I0,1x,A)') &
                        '[FSPS_INIT] Error: fsps_data_load_nebular failed', io_status%code, trim(io_status%message)
                    call fsps_data_close(io_status)
                    error stop 1
                end if

                if (allocated(neb_grid%cont)) ctx%state%nebem_cont = neb_grid%cont
                if (allocated(neb_grid%line)) ctx%state%nebem_line = neb_grid%line
                if (allocated(neb_grid%line_pos)) ctx%state%nebem_line_pos = neb_grid%line_pos
                if (allocated(neb_grid%logz)) ctx%state%nebem_logz = neb_grid%logz
                if (allocated(neb_grid%age)) ctx%state%nebem_age = neb_grid%age
                if (allocated(neb_grid%logu)) ctx%state%nebem_logu = neb_grid%logu

                call fsps_data_close(io_status)
                if (io_status%code /= 0) then
                    write(error_unit, '(A,1x,I0,1x,A)') &
                        '[FSPS_INIT] Error: fsps_data_close failed for nebular load', io_status%code, trim(io_status%message)
                    error stop 1
                end if

                call neb_grid%clear()
                call compute_nebular_kernels(ctx)
            end if
        end if

        if (ctx%state%isoc_type == 'bpss') then
            call load_xray_nebular_grid(ctx)
        end if

        call load_xrb_spectra(ctx)

    end subroutine load_interstellar_physics

    subroutine load_dust_emission_table(ctx, dust_type)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: dust_type

        type(dust_emission_t) :: em_grid
        type(backend_status_t) :: io_status
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri
        integer :: i_spec, k, start_idx, nqpah, numin_cols

        call resolve_data_backend_uri(ctx, backend_mode, uri, dust_type)

        if (trim(backend_mode) == 'legacy') then
            call load_dust_emission_table_legacy(ctx, dust_type)
            return
        end if

        call fsps_data_open(uri, backend_mode, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_open failed for dust emission load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call fsps_data_load_dust_emission(em_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_dust_emission failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(em_grid%qpah) .or. .not. allocated(em_grid%umin) .or. &
            .not. allocated(em_grid%lam) .or. .not. allocated(em_grid%spec)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: dust emission grid is not fully allocated.'
            error stop 1
        end if

        nqpah = size(em_grid%qpah)
        numin_cols = size(em_grid%spec, 3)

        if (size(em_grid%spec, 2) /= nqpah .or. numin_cols /= 2*size(em_grid%umin)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: dust emission grid shape mismatch.'
            error stop 1
        end if

        if (nqpah > size(ctx%state%dustem2_dustem, 2) .or. numin_cols > size(ctx%state%dustem2_dustem, 3)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: dust emission output array is too small for loaded grid.'
            error stop 1
        end if

        ctx%state%qpaharr(1:nqpah) = em_grid%qpah
        ctx%state%uminarr(1:size(em_grid%umin)) = em_grid%umin

        ctx%state%dustem2_dustem = 0.0_wp
        start_idx = max(find_interval(ctx%state%spec_lambda, 1.0e4_wp), 1)
        do k = 1, nqpah
            do i_spec = 1, numin_cols
                ctx%state%dustem2_dustem(start_idx:ctx%state%nspec, k, i_spec) = interpolate_linear( &
                    em_grid%lam, em_grid%spec(:, k, i_spec), ctx%state%spec_lambda(start_idx:ctx%state%nspec))
            end do
        end do

        call fsps_data_close(io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_close failed for dust emission load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call em_grid%clear()
    end subroutine load_dust_emission_table

    subroutine load_dust_attenuation_curves(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(dust_attenuation_t) :: att_grid
        type(backend_status_t) :: io_status
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri
        integer :: n, i, j, k

        call resolve_data_backend_uri(ctx, backend_mode, uri)

        if (trim(backend_mode) == 'legacy') then
            call load_attenuation_curves_legacy(ctx)
            return
        end if

        call fsps_data_open(uri, backend_mode, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_open failed for attenuation load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call fsps_data_load_dust_attenuation(att_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_dust_attenuation failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(att_grid%wg_lam) .or. .not. allocated(att_grid%wg_spec) .or. &
            .not. allocated(att_grid%smc_lam) .or. .not. allocated(att_grid%smc_ext)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: dust attenuation grid is not fully allocated.'
            error stop 1
        end if

        if (size(att_grid%wg_spec, 2) /= 18 .or. size(att_grid%wg_spec, 3) /= 6 .or. size(att_grid%wg_spec, 4) /= 2) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: dust attenuation WG grid shape mismatch.'
            error stop 1
        end if

        do k = 1, 2
            do i = 1, 18
                do j = 1, 6
                    do n = 1, ctx%state%nspec
                        if (ctx%state%spec_lambda(n) > att_grid%wg_lam(size(att_grid%wg_lam))) then
                            ctx%state%wgdust(n, i, j, k) = 0.0_wp
                        else if (ctx%state%spec_lambda(n) < att_grid%wg_lam(1)) then
                            ctx%state%wgdust(n, i, j, k) = att_grid%wg_spec(1, i, j, k)
                        else
                            ctx%state%wgdust(n, i, j, k) = interpolate_linear( &
                                att_grid%wg_lam, att_grid%wg_spec(:, i, j, k), ctx%state%spec_lambda(n))
                        end if
                    end do
                end do
            end do
        end do

        do n = 1, ctx%state%nspec
            if (ctx%state%spec_lambda(n) > att_grid%smc_lam(size(att_grid%smc_lam))) then
                ctx%state%g03smcextn(n) = 0.0_wp
            else if (ctx%state%spec_lambda(n) < att_grid%smc_lam(1)) then
                ctx%state%g03smcextn(n) = att_grid%smc_ext(1)
            else
                ctx%state%g03smcextn(n) = interpolate_linear(att_grid%smc_lam, att_grid%smc_ext, ctx%state%spec_lambda(n))
            end if
        end do

        call fsps_data_close(io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_close failed for attenuation load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call att_grid%clear()
    end subroutine load_dust_attenuation_curves

    subroutine load_agn_dust_models(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(agn_dust_t) :: agn_grid
        type(backend_status_t) :: io_status
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri
        integer :: i, i1, i2

        call resolve_data_backend_uri(ctx, backend_mode, uri)

        if (trim(backend_mode) == 'legacy') then
            call load_agn_dust_models_legacy(ctx)
            return
        end if

        call fsps_data_open(uri, backend_mode, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_open failed for AGN dust load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call fsps_data_load_agn_dust(agn_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_agn_dust failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(agn_grid%tau) .or. .not. allocated(agn_grid%lam) .or. .not. allocated(agn_grid%spec)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGN dust grid is not fully allocated.'
            error stop 1
        end if

        if (size(agn_grid%spec, 2) /= size(agn_grid%tau)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGN dust grid shape mismatch.'
            error stop 1
        end if

        if (size(agn_grid%tau) > size(ctx%state%agndust_tau) .or. size(agn_grid%spec, 2) > size(ctx%state%agndust_spec, 2)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGN dust output arrays are too small for loaded grid.'
            error stop 1
        end if

        ctx%state%agndust_tau(1:size(agn_grid%tau)) = agn_grid%tau
        ctx%state%agndust_spec = 0.0_wp

        i1 = max(find_interval(ctx%state%spec_lambda, agn_grid%lam(1)), 1)
        i2 = max(find_interval(ctx%state%spec_lambda, agn_grid%lam(size(agn_grid%lam))), 1)
        if (i2 < i1) then
            i = i1
            i1 = i2
            i2 = i
        end if

        do i = 1, size(agn_grid%tau)
            ctx%state%agndust_spec(i1:i2, i) = 10.0_wp**interpolate_linear( &
                log10(agn_grid%lam), log10(agn_grid%spec(:, i) + SAFE_FLOOR), &
                log10(ctx%state%spec_lambda(i1:i2))) - SAFE_FLOOR
        end do

        call fsps_data_close(io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_close failed for AGN dust load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call agn_grid%clear()
    end subroutine load_agn_dust_models

    subroutine load_wmbasic_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(aux_wmbasic_t) :: wmb_grid
        type(backend_status_t) :: io_status
        integer :: z, i, j, i1, nzwmb
        real(WP) :: dz
        real(WP), allocatable :: wmbsi(:, :, :, :)

        call fsps_data_load_wmbasic(wmb_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_wmbasic failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        if (.not. allocated(wmb_grid%lam) .or. .not. allocated(wmb_grid%logt) .or. &
            .not. allocated(wmb_grid%z) .or. .not. allocated(wmb_grid%spec)) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: WMBasic grid is not fully allocated.'
            error stop 1
        end if

        nzwmb = size(wmb_grid%z)
        if (size(wmb_grid%logt) /= NDIM_WMB_LOGT .or. size(wmb_grid%spec, 2) /= NDIM_WMB_LOGT .or. &
            size(wmb_grid%spec, 3) /= NDIM_WMB_LOGG .or. size(wmb_grid%spec, 4) /= nzwmb) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: WMBasic grid shape mismatch.'
            error stop 1
        end if

        ctx%state%wmb_logt = wmb_grid%logt
        ctx%state%wmb_logg = [3.5_wp, 4.0_wp, 4.5_wp]

        allocate(wmbsi(ctx%state%nspec, nzwmb, NDIM_WMB_LOGT, NDIM_WMB_LOGG))

        do z = 1, nzwmb
            do j = 1, NDIM_WMB_LOGG
                do i = 1, NDIM_WMB_LOGT
                    wmbsi(:, z, i, j) = max(interpolate_linear(wmb_grid%lam, wmb_grid%spec(:, i, j, z), &
                                                              ctx%state%spec_lambda), SAFE_FLOOR)
                end do
            end do
        end do

        do z = 1, ctx%state%nz
            i1 = min(max(find_interval(log10(wmb_grid%z/ctx%state%zsol_spec), &
                                       log10(ctx%state%zlegend(z)/ctx%state%zsol)), 1), nzwmb - 1)

            dz = (log10(ctx%state%zlegend(z)/ctx%state%zsol) - log10(wmb_grid%z(i1)/ctx%state%zsol_spec)) / &
                 (log10(wmb_grid%z(i1 + 1)/ctx%state%zsol_spec) - log10(wmb_grid%z(i1)/ctx%state%zsol_spec))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            ctx%state%wmb_spec(:, z, :, :) = real(10.0_wp**((1.0_wp - dz)*log10(wmbsi(:, i1, :, :) + SAFE_FLOOR) + &
                                                             dz*log10(wmbsi(:, i1 + 1, :, :) + SAFE_FLOOR)), &
                                                  kind=kind(ctx%state%wmb_spec))
        end do

        deallocate(wmbsi)
        call wmb_grid%clear()
    end subroutine load_wmbasic_spectra

    subroutine load_post_agb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(aux_pagb_t) :: pagb_grid
        type(backend_status_t) :: io_status
        integer :: i, j

        call fsps_data_load_pagb(pagb_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_pagb failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        if (.not. allocated(pagb_grid%lam) .or. .not. allocated(pagb_grid%logt) .or. .not. allocated(pagb_grid%spec)) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: Post-AGB grid is not fully allocated.'
            error stop 1
        end if

        if (size(pagb_grid%logt) /= NDIM_PAGB .or. size(pagb_grid%spec, 2) /= NDIM_PAGB .or. size(pagb_grid%spec, 3) /= 2) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: Post-AGB grid shape mismatch.'
            error stop 1
        end if

        ctx%state%pagb_logt = pagb_grid%logt

        do j = 1, 2
            do i = 1, NDIM_PAGB
                ctx%state%pagb_spec(:, i, j) = max(interpolate_linear(pagb_grid%lam, pagb_grid%spec(:, i, j), &
                                                                      ctx%state%spec_lambda), SAFE_FLOOR)
            end do
        end do

        call pagb_grid%clear()
    end subroutine load_post_agb_spectra

    subroutine load_wr_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(aux_wr_t) :: wr_grid
        type(backend_status_t) :: io_status
        integer :: i, j, i1, nz_wr, nlam_wr
        real(WP) :: dz
        real(WP) :: target_logz
        real(WP), allocatable :: wrn_interp(:, :, :), wrc_interp(:, :, :)
        real(WP), allocatable :: wr_z_log(:), wr_lam_log(:), target_lam_log(:)

        call fsps_data_load_wr(wr_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_wr failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        if (.not. allocated(wr_grid%logt_wn) .or. .not. allocated(wr_grid%logt_wc) .or. .not. allocated(wr_grid%z) .or. &
            .not. allocated(wr_grid%lam) .or. .not. allocated(wr_grid%spec_wn) .or. .not. allocated(wr_grid%spec_wc)) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: WR grid is not fully allocated.'
            error stop 1
        end if

        nz_wr = size(wr_grid%z)
        nlam_wr = size(wr_grid%lam)
        if (size(wr_grid%logt_wn) /= NDIM_WR .or. size(wr_grid%logt_wc) /= NDIM_WR .or. &
            size(wr_grid%spec_wn, 1) /= nlam_wr .or. size(wr_grid%spec_wn, 2) /= NDIM_WR .or. &
            size(wr_grid%spec_wn, 3) /= nz_wr .or. size(wr_grid%spec_wc, 1) /= nlam_wr .or. &
            size(wr_grid%spec_wc, 2) /= NDIM_WR .or. size(wr_grid%spec_wc, 3) /= nz_wr) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: WR grid shape mismatch.'
            error stop 1
        end if

        ctx%state%wrn_logt = wr_grid%logt_wn
        ctx%state%wrc_logt = wr_grid%logt_wc

        allocate(wrn_interp(ctx%state%nspec, NDIM_WR, nz_wr))
        allocate(wrc_interp(ctx%state%nspec, NDIM_WR, nz_wr))
        allocate(wr_z_log(nz_wr), wr_lam_log(nlam_wr), target_lam_log(ctx%state%nspec))

        wr_z_log = log10(wr_grid%z / ctx%state%zsol_spec)
        wr_lam_log = log10(wr_grid%lam)
        target_lam_log = log10(ctx%state%spec_lambda)

        do j = 1, nz_wr
            do i = 1, NDIM_WR
                wrn_interp(:, i, j) = 10.0_wp**interpolate_linear( &
                    wr_lam_log, log10(wr_grid%spec_wn(:, i, j) + SAFE_FLOOR), target_lam_log) - SAFE_FLOOR
                wrc_interp(:, i, j) = 10.0_wp**interpolate_linear( &
                    wr_lam_log, log10(wr_grid%spec_wc(:, i, j) + SAFE_FLOOR), target_lam_log) - SAFE_FLOOR
            end do
        end do

        do j = 1, ctx%state%nz
            target_logz = log10(ctx%state%zlegend(j) / ctx%state%zsol_spec)
            i1 = min(max(find_interval(wr_z_log, target_logz), 1), nz_wr - 1)
            dz = (target_logz - wr_z_log(i1)) / (wr_z_log(i1 + 1) - wr_z_log(i1))
            dz = min(max(dz, 0.0_wp), 1.0_wp)

            ctx%state%wrn_spec(:, :, j) = real(10.0_wp**((1.0_wp - dz)*log10(wrn_interp(:, :, i1) + SAFE_FLOOR) + &
                                                          dz*log10(wrn_interp(:, :, i1 + 1) + SAFE_FLOOR)), &
                                               kind=kind(ctx%state%wrn_spec))
            ctx%state%wrc_spec(:, :, j) = real(10.0_wp**((1.0_wp - dz)*log10(wrc_interp(:, :, i1) + SAFE_FLOOR) + &
                                                          dz*log10(wrc_interp(:, :, i1 + 1) + SAFE_FLOOR)), &
                                               kind=kind(ctx%state%wrc_spec))
        end do

        deallocate(wrn_interp, wrc_interp, wr_z_log, wr_lam_log, target_lam_log)
        call wr_grid%clear()
    end subroutine load_wr_spectra

    subroutine load_agb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(aux_agb_t) :: agb_grid
        type(backend_status_t) :: io_status
        integer :: i, iz, i1, nz_o
        real(WP) :: dz
        real(WP), allocatable :: z_o_log(:)

        call fsps_data_load_agb(agb_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_agb failed', io_status%code, trim(io_status%message)
            error stop 1
        end if

        if (.not. allocated(agb_grid%z_o) .or. .not. allocated(agb_grid%logt_o) .or. .not. allocated(agb_grid%logt_c) .or. &
            .not. allocated(agb_grid%logt_car) .or. .not. allocated(agb_grid%lam_o) .or. .not. allocated(agb_grid%lam_c) .or. &
            .not. allocated(agb_grid%lam_car) .or. .not. allocated(agb_grid%spec_o) .or. .not. allocated(agb_grid%spec_c) .or. &
            .not. allocated(agb_grid%spec_car)) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGB grid is not fully allocated.'
            error stop 1
        end if

        if (size(agb_grid%logt_o, 2) /= N_AGB_O .or. size(agb_grid%logt_c) /= N_AGB_C .or. &
            size(agb_grid%logt_car) /= N_AGB_CAR .or. size(agb_grid%spec_o, 2) /= N_AGB_O .or. &
            size(agb_grid%spec_c, 2) /= N_AGB_C .or. size(agb_grid%spec_car, 2) /= N_AGB_CAR) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGB grid shape mismatch.'
            error stop 1
        end if

        nz_o = size(agb_grid%z_o)
        if (size(agb_grid%logt_o, 1) /= nz_o) then
            write(error_unit, '(A)') '[FSPS_INIT] Error: AGB O-rich logT grid mismatch.'
            error stop 1
        end if

        allocate(z_o_log(nz_o))
        z_o_log = agb_grid%z_o

        do iz = 1, ctx%state%nz
              i1 = min(max(find_interval(z_o_log, log10(ctx%state%zlegend(iz)/ctx%state%zsol_spec)), 1), nz_o - 1)
              dz = (log10(ctx%state%zlegend(iz)/ctx%state%zsol_spec) - z_o_log(i1)) / &
                  (z_o_log(i1+1) - z_o_log(i1))
            dz = min(max(dz, 0.0_wp), 1.0_wp)
              ctx%state%agb_logt_o(iz, :) = (1.0_wp - dz) * agb_grid%logt_o(i1, :) + &
                                      dz * agb_grid%logt_o(i1+1, :)
        end do

        ctx%state%agb_logt_c = agb_grid%logt_c
        ctx%state%agb_logt_car = agb_grid%logt_car

        do i = 1, N_AGB_O
            ctx%state%agb_spec_o(:, i) = max(interpolate_linear(agb_grid%lam_o, agb_grid%spec_o(:, i), &
                                                                 ctx%state%spec_lambda), SAFE_FLOOR)
        end do

        do i = 1, N_AGB_C
            ctx%state%agb_spec_c(:, i) = max(interpolate_linear(agb_grid%lam_c, agb_grid%spec_c(:, i), &
                                                                 ctx%state%spec_lambda), SAFE_FLOOR)
        end do

        do i = 1, N_AGB_CAR
            ctx%state%agb_spec_car(:, i) = max(interpolate_linear(agb_grid%lam_car, agb_grid%spec_car(:, i), &
                                                                   ctx%state%spec_lambda), SAFE_FLOOR)
        end do

        deallocate(z_o_log)
        call agb_grid%clear()
    end subroutine load_agb_spectra

    !> @brief Loads Filter Curves, Solar/Vega SEDs, and Spectral Indices.
    !>
    !> @details
    !> Orchestrates the loading of photometric and calibration data via `fsps_io`.
    !>
    !> **Data Dependency Note:**
    !> The Solar and Vega spectra are loaded *before* the filters. This allows
    !> `load_filter_definitions` (in `fsps_io`) to integrate the transmission curves
    !> against these standards immediately upon loading, populating `magsun` and
    !> `magvega` without requiring permanent storage of the full transmission curves
    !> in the cache.
    !>
    !> @param[inout] ctx   The context containing path configuration.
    !> @param[inout] cache The cache object to populate.
    subroutine load_photometry_data(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        if (len_trim(ctx%state%alt_filter_file) == 0 .or. trim(ctx%state%alt_filter_file) == 'allfilters.dat') then
            call load_filter_definitions(ctx)
        else
            call load_filter_definitions(ctx, trim(ctx%state%alt_filter_file))
        end if

        call load_standard_sed(ctx)

        if (ctx%compute_vega_mags_val == 1) then
            call convert_to_vega_system(ctx)
        end if

        if (len_trim(ctx%state%alt_filter_file) == 0 .or. trim(ctx%state%alt_filter_file) == 'allfilters.dat') then
            call apply_legacy_filter_norm(ctx)
        end if

    end subroutine load_photometry_data

    !> @brief Computes data that depends on loaded physics but isn't read from files.
    !>
    !> @details
    !> Populates shared cache arrays that are derived from the primary loaded data.
    !>
    !> 1. **Frequency Grid (`spec_nu`):** Calculated from `spec_lambda`.
    !> 2. **Expanded Time Grid (`time_full`):** Interpolates the isochrone age grid
    !>    if `TIME_RES_INCR > 1`.
    !>
    !> @note
    !> This routine populates *Cached* data. Context-specific derived data (like
    !> `cosmospl` or `mwdindex`) which are not pointers in the context state must
    !> be computed in a separate local initialization step.
    !>
    !> @param[inout] ctx   The context (used for accessing constants/state).
    !> @param[inout] cache The cache object to populate.
    subroutine compute_derived_quantities(ctx, cache)
        type(fsps_context_t), intent(inout) :: ctx
        type(fsps_setup_cache_t), intent(inout) :: cache

        integer :: i, j, k
        integer :: nt_iso, nt_full
        real(WP) :: dt, t_start

        ! --------------------------------------------------------------------
        ! 1. COMPUTE FREQUENCY GRID (spec_nu)
        ! --------------------------------------------------------------------
        ! c / lambda (converted to appropriate units if necessary, but FSPS usually
        ! keeps standard units). Lambda is in Angstroms.
        ! C_LIGHT is typically in Ang/s (2.99792458e18).

        where (cache%spec_lambda > SAFE_FLOOR)
            cache%spec_nu = C_LIGHT/cache%spec_lambda
        elsewhere
            cache%spec_nu = 0.0_wp
        end where

        ! --------------------------------------------------------------------
        ! 2. COMPUTE EXTENDED TIME GRID (time_full)
        ! --------------------------------------------------------------------
        ! If TIME_RES_INCR > 1, we interpolate between the native isochrone
        ! time steps to create a finer grid for the SSP generation.

        nt_iso = ctx%state%nt
        nt_full = ctx%state%ntfull

        ! Use the first available metallicity (z=1) as the reference time grid.
        ! Isochrone libraries (MIST, Padova) enforce consistent age grids across Z.

        if (TIME_RES_INCR == 1) then
            ! Simple copy if no oversampling
            cache%time_full(1:nt_iso) = cache%timestep_isoc(1, 1:nt_iso)
        else
            ! Interpolate
            do i = 1, nt_iso - 1

                ! The start time of this isochrone interval
                t_start = cache%timestep_isoc(1, i)

                ! The time step size (delta t) between oversampled points
                dt = (cache%timestep_isoc(1, i + 1) - t_start)/real(TIME_RES_INCR, kind=WP)

                ! Fill the sub-steps
                do j = 0, TIME_RES_INCR - 1
                    k = (i - 1)*TIME_RES_INCR + j + 1
                    if (k > nt_full) exit
                    cache%time_full(k) = t_start + real(j, kind=WP)*dt
                end do

            end do

            ! Handle the final point
            cache%time_full(nt_full) = cache%timestep_isoc(1, nt_iso)
        end if

    end subroutine compute_derived_quantities

    !> @brief Computes local physical properties specific to this context instance.
    !>
    !> @details
    !> These quantities are stored in `fsps_context_state_t` but are *not* pointers,
    !> meaning they must be re-calculated for every new context (they cannot be
    !> linked from the cache).
    !>
    !> Calculations:
    !> 1. **MW Dust Indices (`mwdindex`):** Defines wavelength breakpoints for the
    !>    Cardelli, Clayton, & Mathis (1989) extinction curve.
    !> 2. **Cosmology Grid (`cosmospl`):** A lookup table for (Redshift, Age, Distance)
    !>    used for fast interpolation during runtime.
    !> 3. **Universe Age (`tuniv`):** Age at z=0.
    !> 4. **Key Indices:** `whlam5000` (5000A) and `whlylim` (912A).
    !>
    !> @param[inout] ctx The context to populate.
    subroutine compute_local_properties(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: i, j
        real(WP) :: x_inv_micron, a_scale
        real(WP), dimension(:), pointer :: lam

        ! Alias for readability
        lam => ctx%state%spec_lambda

        ! --------------------------------------------------------------------
        ! 1. MW Dust Extinction Indices (CCM89 Breakpoints)
        ! --------------------------------------------------------------------
        ! x = 10000 / lambda(A)
        ! Breaks at x = [0.1, 1.1, 3.3, 5.9, 8.0, 12.0]

        ctx%state%mwdindex = 0

        do j = 1, ctx%state%nspec
            if (lam(j) > SAFE_FLOOR) then
                x_inv_micron = 1.0e4_wp/lam(j)

                ! Legacy logic finds the *last* index where the condition holds
                ! assuming lambda is increasing (and x is decreasing).
                if (x_inv_micron >= 0.1_wp) ctx%state%mwdindex(1) = j
                if (x_inv_micron >= 1.1_wp) ctx%state%mwdindex(2) = j
                if (x_inv_micron >= 3.3_wp) ctx%state%mwdindex(3) = j
                if (x_inv_micron >= 5.9_wp) ctx%state%mwdindex(4) = j
                if (x_inv_micron >= 8.0_wp) ctx%state%mwdindex(5) = j
                if (x_inv_micron > 12.0_wp) ctx%state%mwdindex(6) = j
            end if
        end do

        ! --------------------------------------------------------------------
        ! 2. Cosmology Lookup Table (cosmospl)
        ! --------------------------------------------------------------------
        ! Generates a grid of 500 points scaling factor 'a' from ~0.001 to 1.0.
        ! Stores: [Redshift, Age(Gyr), Luminosity Distance(pc)]

        do i = 1, 500
            ! Scale factor 'a' distribution (legacy formula)
            a_scale = real(i - 1, kind=WP)/499.0_wp*(1.0_wp - 1.0_wp/1001.0_wp) + &
                      1.0_wp/1001.0_wp

            ! Col 1: Redshift (z = 1/a - 1)
            ctx%state%cosmospl(i, 1) = (1.0_wp/a_scale) - 1.0_wp

            ! Col 2: Age of Universe at z
            ctx%state%cosmospl(i, 2) = get_universe_age(ctx, ctx%state%cosmospl(i, 1))

            ! Col 3: Luminosity Distance at z
            ctx%state%cosmospl(i, 3) = get_luminosity_distance(ctx, ctx%state%cosmospl(i, 1))
        end do

        ! --------------------------------------------------------------------
        ! 3. Age of the Universe (z=0)
        ! --------------------------------------------------------------------
        ctx%state%tuniv = get_universe_age(ctx, 0.0_wp)

        ! --------------------------------------------------------------------
        ! 4. Key Wavelength Indices
        ! --------------------------------------------------------------------
        ! Find index for 5000 Angstroms
        ctx%state%whlam5000 = find_interval(lam, 5000.0_wp)

        ! Find index for Lyman Limit (912 Angstroms)
        ctx%state%whlylim = find_interval(lam, 912.0_wp)

        ! Recompute frequency grid
        ctx%state%spec_nu = 0.0_wp
        where (lam > SAFE_FLOOR)
            ctx%state%spec_nu = C_LIGHT/lam
        end where

        ctx%state%sfh_tab = 0.0_wp
        ctx%state%ntabsfh = 0

    end subroutine compute_local_properties

    !> @brief Maps the pointers in `ctx%state` to the arrays in `cache`.
    !>
    !> @details
    !> Performs two types of linking:
    !> 1. **Pointer Association (`=>`):** For large dynamic arrays (spectra, isochrones).
    !>    This allows multiple contexts to share read-only memory.
    !> 2. **Value Copy (`=`):** For small fixed-size arrays (grid axes like `logt`, `logg`).
    !>    This ensures that contexts initialized from the cache still have valid
    !>    axis definitions, fixing a potential issue in the legacy architecture.
    !>
    !> @param[inout] ctx   The context to populate.
    !> @param[in]    cache The shared cache object containing the loaded data.
    subroutine link_context_to_cache(ctx, cache)
        type(fsps_context_t), intent(inout) :: ctx
        type(fsps_setup_cache_t), target, intent(in) :: cache

        ctx%state%indexdefined => cache%indexdefined
        ctx%state%wgdust => cache%wgdust
        ctx%state%g03smcextn => cache%g03smcextn
        ctx%state%bands => cache%bands
        ctx%state%magsun => cache%magsun
        ctx%state%magvega => cache%magvega
        ctx%state%filter_leff => cache%filter_leff
        ctx%state%vega_spec => cache%vega_spec
        ctx%state%sun_spec => cache%sun_spec
        ctx%state%spec_lambda => cache%spec_lambda
        ctx%state%spec_nu => cache%spec_nu
        ctx%state%spec_res => cache%spec_res
        ctx%state%speclib => cache%speclib
        ctx%state%wmb_spec => cache%wmb_spec
        ctx%state%agb_spec_o => cache%agb_spec_o
        ctx%state%agb_logt_o => cache%agb_logt_o
        ctx%state%agb_spec_c => cache%agb_spec_c
        ctx%state%agb_logt_c => cache%agb_logt_c
        ctx%state%agb_spec_car => cache%agb_spec_car
        ctx%state%pagb_spec => cache%pagb_spec
        ctx%state%wrn_spec => cache%wrn_spec
        ctx%state%wrc_spec => cache%wrc_spec
        ctx%state%qpaharr => cache%qpaharr
        ctx%state%uminarr => cache%uminarr
        ctx%state%lambda_dustem => cache%lambda_dustem
        ctx%state%dustem_dustem => cache%dustem_dustem
        ctx%state%dustem2_dustem => cache%dustem2_dustem
        ctx%state%flux_dagb => cache%flux_dagb
        ctx%state%nebem_cont => cache%nebem_cont
        ctx%state%xnebem_cont => cache%xnebem_cont
        ctx%state%neb_res_min => cache%neb_res_min
        ctx%state%gaussnebarr => cache%gaussnebarr
        ctx%state%agndust_spec => cache%agndust_spec
        ctx%state%mact_isoc => cache%mact_isoc
        ctx%state%logl_isoc => cache%logl_isoc
        ctx%state%logt_isoc => cache%logt_isoc
        ctx%state%logg_isoc => cache%logg_isoc
        ctx%state%ffco_isoc => cache%ffco_isoc
        ctx%state%phase_isoc => cache%phase_isoc
        ctx%state%mini_isoc => cache%mini_isoc
        ctx%state%lmdot_isoc => cache%lmdot_isoc
        ctx%state%nmass_isoc => cache%nmass_isoc
        ctx%state%timestep_isoc => cache%timestep_isoc
        ctx%state%zlegend => cache%zlegend
        ctx%state%zlegendinit => cache%zlegendinit
        ctx%state%time_full => cache%time_full
        ctx%state%bpass_spec_ssp => cache%bpass_spec_ssp
        ctx%state%bpass_mass_ssp => cache%bpass_mass_ssp
        ctx%state%lam_xrb => cache%lam_xrb
        ctx%state%spec_xrb => cache%spec_xrb
        ctx%state%ages_xrb => cache%ages_xrb
        ctx%state%zmet_xrb => cache%zmet_xrb

    end subroutine link_context_to_cache

    !> @brief Allocates arrays that are specific to this context instance.
    !>
    !> @details
    !> Allocates working buffers for SSP generation, spectral smoothing, and intermediate results.
    !> Unlike the cached physics arrays, these arrays are mutable and specific to the
    !> current simulation state.
    !>
    !> Allocates:
    !> - `spec_ssp_zz`: The computed SSP grid [nspec, ntfull, nz].
    !> - `mass_ssp_zz`, `lbol_ssp_zz`: Stellar mass and luminosity grids.
    !> - `weight_ssp`: Weights for interpolation.
    !> - `spec_young`, `spec_old`: Buffers for separating young/old stellar populations.
    !> - `lsfinfo%lsf`: Buffer for the Line Spread Function.
    !>
    !> @param[inout] ctx The FSPS context to allocate.
    subroutine allocate_local_workspace(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        ! Clean up existing allocations if this context is being re-initialized
        if (allocated(ctx%state%spec_ssp_zz)) deallocate (ctx%state%spec_ssp_zz)
        if (allocated(ctx%state%mass_ssp_zz)) deallocate (ctx%state%mass_ssp_zz)
        if (allocated(ctx%state%lbol_ssp_zz)) deallocate (ctx%state%lbol_ssp_zz)
        if (allocated(ctx%state%weight_ssp)) deallocate (ctx%state%weight_ssp)
        if (allocated(ctx%state%spec_young)) deallocate (ctx%state%spec_young)
        if (allocated(ctx%state%spec_old)) deallocate (ctx%state%spec_old)
        if (allocated(ctx%state%lsfinfo%lsf)) deallocate (ctx%state%lsfinfo%lsf)

        ! Allocate SSP Grids
        ! Dimensions: (Wavelengths, Time Steps, Metallicity)
        allocate (ctx%state%spec_ssp_zz(ctx%state%nspec, ctx%state%ntfull, ctx%state%nz))
        allocate (ctx%state%mass_ssp_zz(ctx%state%ntfull, ctx%state%nz))
        allocate (ctx%state%lbol_ssp_zz(ctx%state%ntfull, ctx%state%nz))

        ! Allocate Weighting Buffer
        allocate (ctx%state%weight_ssp(ctx%state%ntfull, ctx%state%nz))

        ! Allocate Component Spectral Buffers
        allocate (ctx%state%spec_young(ctx%state%nspec, ctx%state%ntfull))
        allocate (ctx%state%spec_old(ctx%state%nspec, ctx%state%ntfull))

        ! Allocate Line Spread Function Buffer
        allocate (ctx%state%lsfinfo%lsf(ctx%state%nspec))

        ! Initialize to zero to prevent garbage data
        ctx%state%spec_ssp_zz = 0.0_wp
        ctx%state%mass_ssp_zz = 0.0_wp
        ctx%state%lbol_ssp_zz = 0.0_wp
        ctx%state%weight_ssp = 0.0_wp
        ctx%state%spec_young = 0.0_wp
        ctx%state%spec_old = 0.0_wp
        ctx%state%lsfinfo%lsf = 0.0_wp

    end subroutine allocate_local_workspace

    ! =====================================================================
    ! SUPPORT HELPERS
    ! =====================================================================

    subroutine normalize_library_names(isoc_in, spec_in, dust_in, isoc_out, spec_out, dust_out)
        character(len=*), intent(in), optional :: isoc_in
        character(len=*), intent(in), optional :: spec_in
        character(len=*), intent(in), optional :: dust_in
        character(len=64), intent(out) :: isoc_out
        character(len=64), intent(out) :: spec_out
        character(len=64), intent(out) :: dust_out
        character(len=64) :: iso, spec, dust

        iso = 'mist'
        spec = 'miles'
        dust = 'DL07'

        if (present(isoc_in)) then
            if (len_trim(isoc_in) > 0) iso = to_lower(trim(isoc_in))
        end if
        if (present(spec_in)) then
            if (len_trim(spec_in) > 0) spec = to_lower(trim(spec_in))
        end if
        if (present(dust_in)) then
            if (len_trim(dust_in) > 0) dust = trim(dust_in)
        end if

        select case (iso)
        case ('padova', 'padova2007', 'default')
            iso = 'pdva'
        case ('parsec')
            iso = 'prsc'
        case ('basti')
            iso = 'bsti'
        case ('geneva')
            iso = 'gnva'
        case ('bpass', 'bpss')
            iso = 'bpss'
        case default
            ! keep as-is
        end select

        if (iso == 'bpss') then
            spec = 'bpass'
        end if

        select case (spec)
        case ('basel', 'basel2.2', 'basel3.1')
            spec = 'basel'
        case ('miles')
            spec = 'miles'
        case ('bpass')
            spec = 'bpass'
        case default
            if (len_trim(spec) >= 3) then
                if (spec(1:3) == 'c3k') spec = trim(spec)
            end if
        end select

        if (trim(dust) == 'themis') dust = 'THEMIS'
        if (trim(dust) /= 'THEMIS') dust = 'DL07'

        isoc_out = iso
        spec_out = spec
        dust_out = dust
    end subroutine normalize_library_names

    subroutine configure_dust_model(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        if (trim(ctx%state%str_dustem) == 'THEMIS') then
            ctx%state%ndim_dustem = 576
            ctx%state%numin_dustem = 37
            ctx%state%nqpah_dustem = 11
        else
            ctx%state%str_dustem = 'DL07'
            ctx%state%ndim_dustem = 1001
            ctx%state%numin_dustem = 22
            ctx%state%nqpah_dustem = 7
        end if
    end subroutine configure_dust_model

    subroutine build_cache_key(ctx, zin, cache_key)
        type(fsps_context_t), intent(in) :: ctx
        integer, intent(in) :: zin
        character(len=*), intent(out) :: cache_key

        write (cache_key, '(A,"|",A,"|",A,"|",A,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",I0,"|",A)') &
            trim(ctx%sps_home), trim(ctx%state%isoc_type), trim(ctx%state%spec_type), trim(ctx%state%str_dustem), &
            zin, ctx%smooth_velocity_val, ctx%setup_nebular_gaussians_val, ctx%add_neb_emission_val, &
            ctx%add_neb_continuum_val, ctx%add_dust_emission_val, ctx%add_agn_dust_val, &
            ctx%add_xrb_emission_val, ctx%add_agb_dust_model_val, trim(ctx%state%alt_filter_file)
    end subroutine build_cache_key

    subroutine populate_cache_metadata(ctx, cache)
        type(fsps_context_t), intent(in) :: ctx
        type(fsps_setup_cache_t), intent(inout) :: cache

        cache%isoc_type = trim(ctx%state%isoc_type)
        cache%spec_type = trim(ctx%state%spec_type)
        cache%dust_type = trim(ctx%state%str_dustem)
        cache%alt_filter_file = ctx%state%alt_filter_file
        cache%nz = ctx%state%nz
        cache%nt = ctx%state%nt
        cache%nspec = ctx%state%nspec
        cache%nzinit = ctx%state%nzinit
        cache%nbands = ctx%state%nbands
        cache%nindx = ctx%state%nindx
        cache%ntfull = ctx%state%ntfull
        cache%nspec_xrb = ctx%state%nspec_xrb
        cache%nt_xrb = ctx%state%nt_xrb
        cache%nz_xrb = ctx%state%nz_xrb
        cache%smooth_velocity = ctx%smooth_velocity_val
        cache%setup_nebular_gaussians = ctx%setup_nebular_gaussians_val
        cache%add_neb_emission = ctx%add_neb_emission_val
        cache%add_neb_continuum = ctx%add_neb_continuum_val
        cache%add_dust_emission = ctx%add_dust_emission_val
        cache%add_agn_dust = ctx%add_agn_dust_val
        cache%add_xrb_emission = ctx%add_xrb_emission_val
        cache%add_agb_dust_model = ctx%add_agb_dust_model_val
        cache%use_wr_spectra = ctx%use_wr_spectra_val
    end subroutine populate_cache_metadata

    subroutine initialize_cache_arrays(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        ctx%state%bands = 0.0_wp
        ctx%state%magsun = 0.0_wp
        ctx%state%magvega = 0.0_wp
        ctx%state%filter_leff = 0.0_wp
        ctx%state%vega_spec = 0.0_wp
        ctx%state%sun_spec = 0.0_wp
        ctx%state%spec_lambda = 0.0_wp
        ctx%state%spec_nu = 0.0_wp
        ctx%state%spec_res = 0.0_wp
        ctx%state%speclib = 0.0
        ctx%state%wmb_spec = 0.0
        ctx%state%agb_spec_o = 0.0_wp
        ctx%state%agb_logt_o = 0.0_wp
        ctx%state%agb_spec_c = 0.0_wp
        ctx%state%agb_logt_c = 0.0_wp
        ctx%state%agb_spec_car = 0.0_wp
        ctx%state%pagb_spec = 0.0_wp
        ctx%state%wrn_spec = 0.0_wp
        ctx%state%wrc_spec = 0.0_wp

        if (trim(ctx%state%str_dustem) == 'THEMIS') then
         ctx%state%qpaharr = (/0.02_wp, 0.06_wp, 0.10_wp, 0.14_wp, 0.17_wp, 0.20_wp, 0.24_wp, 0.28_wp, 0.32_wp, 0.36_wp, 0.40_wp/) &
                                /2.2_wp*100.0_wp
        ctx%state%uminarr = (/0.1_wp, 0.12_wp, 0.15_wp, 0.17_wp, 0.2_wp, 0.25_wp, 0.3_wp, 0.35_wp, 0.4_wp, 0.5_wp, 0.6_wp, 0.7_wp, &
                           0.8_wp, 1.0_wp, 1.2_wp, 1.5_wp, 1.7_wp, 2.0_wp, 2.5_wp, 3.0_wp, 3.5_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, &
                         8.0_wp, 10.0_wp, 12.0_wp, 15.0_wp, 17.0_wp, 20.0_wp, 25.0_wp, 30.0_wp, 35.0_wp, 40.0_wp, 50.0_wp, 80.0_wp/)
        else
            ctx%state%qpaharr = (/0.47_wp, 1.12_wp, 1.77_wp, 2.50_wp, 3.19_wp, 3.90_wp, 4.58_wp/)
            ctx%state%uminarr = (/0.1_wp, 0.15_wp, 0.2_wp, 0.3_wp, 0.4_wp, 0.5_wp, 0.7_wp, 0.8_wp, 1.0_wp, 1.2_wp, 1.5_wp, 2.0_wp, &
                                  2.5_wp, 3.0_wp, 4.0_wp, 5.0_wp, 7.0_wp, 8.0_wp, 12.0_wp, 15.0_wp, 20.0_wp, 25.0_wp/)
        end if

        ctx%state%lambda_dustem = 0.0_wp
        ctx%state%dustem_dustem = 0.0_wp
        ctx%state%dustem2_dustem = 0.0_wp
        ctx%state%flux_dagb = 0.0_wp
        ctx%state%nebem_cont = 0.0_wp
        ctx%state%xnebem_cont = 0.0_wp
        ctx%state%neb_res_min = 0.0_wp
        ctx%state%gaussnebarr = 0.0_wp
        ctx%state%agndust_spec = 0.0_wp
        ctx%state%mact_isoc = 0.0_wp
        ctx%state%logl_isoc = 0.0_wp
        ctx%state%logt_isoc = 0.0_wp
        ctx%state%logg_isoc = 0.0_wp
        ctx%state%ffco_isoc = 0.0_wp
        ctx%state%phase_isoc = 0.0_wp
        ctx%state%mini_isoc = 0.0_wp
        ctx%state%lmdot_isoc = 0.0_wp
        ctx%state%nmass_isoc = 0
        ctx%state%timestep_isoc = 0.0_wp
        ctx%state%zlegend = -99.0_wp
        ctx%state%zlegendinit = -99.0_wp
        ctx%state%time_full = 0.0_wp
        ctx%state%bpass_spec_ssp = 0.0_wp
        ctx%state%bpass_mass_ssp = 0.0_wp
        ctx%state%lam_xrb = 0.0_wp
        ctx%state%spec_xrb = 0.0_wp
        ctx%state%ages_xrb = 0.0_wp
        ctx%state%zmet_xrb = 0.0_wp
        ctx%state%indexdefined = 0.0_wp
        ctx%state%wgdust = 0.0_wp
        ctx%state%g03smcextn = 0.0_wp
        ctx%state%sfh_tab = 0.0_wp
        ctx%state%ntabsfh = 0
    end subroutine initialize_cache_arrays

    subroutine compute_filter_leff(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: i
        real(WP) :: d

        do i = 1, ctx%state%nbands
            ctx%state%filter_leff(i) = integrate_trapezoid_array(ctx%state%spec_lambda, &
                                                                 ctx%state%spec_lambda*ctx%state%bands(:, i))
            d = integrate_trapezoid_array(ctx%state%spec_lambda, ctx%state%bands(:, i)/ctx%state%spec_lambda)
            if (ctx%state%filter_leff(i) /= ctx%state%filter_leff(i) .or. d /= d .or. d <= SAFE_FLOOR) then
                ctx%state%filter_leff(i) = 0.0_wp
            else
                ctx%state%filter_leff(i) = sqrt(ctx%state%filter_leff(i)/d)
            end if
        end do
    end subroutine compute_filter_leff

    subroutine build_time_grid(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: i
        real(WP) :: d1
        integer :: zref

        if (ctx%state%isoc_type == 'bpss') return

        zref = 1
        if (ctx%zin >= 1 .and. ctx%zin <= ctx%state%nz) then
            zref = ctx%zin
        end if
        do i = 1, ctx%state%ntfull
            if (mod(i - 1, TIME_RES_INCR) == 0) then
                ctx%state%time_full(i) = ctx%state%timestep_isoc(zref, (i - 1)/TIME_RES_INCR + 1)
            else
                if ((i - 1)/TIME_RES_INCR + 2 < ctx%state%nt) then
                    d1 = (ctx%state%timestep_isoc(zref, (i - 1)/TIME_RES_INCR + 2) - &
                          ctx%state%timestep_isoc(zref, (i - 1)/TIME_RES_INCR + 1))/TIME_RES_INCR
                end if
                ctx%state%time_full(i) = ctx%state%time_full(i - 1) + d1
            end if
        end do
    end subroutine build_time_grid

    subroutine compute_nebular_kernels(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: i, j
        real(WP) :: dlam

        do i = 1, NEMLINE
            j = min(max(find_interval(ctx%state%spec_lambda, ctx%state%nebem_line_pos(i)), 1), ctx%state%nspec - 1)
            ctx%state%neb_res_min(i) = ctx%state%spec_lambda(j + 1) - ctx%state%spec_lambda(j)
        end do

        if (ctx%setup_nebular_gaussians_val == 1) then
            do i = 1, NEMLINE
                if (ctx%smooth_velocity_val == 1) then
                    dlam = ctx%state%nebem_line_pos(i)*ctx%nebular_smooth_init_val/C_LIGHT*1.0e13_wp
                else
                    dlam = ctx%nebular_smooth_init_val
                end if
                dlam = max(dlam, ctx%state%neb_res_min(i)*2.0_wp)
                ctx%state%gaussnebarr(:, i) = 1.0_wp/sqrt(2.0_wp*PI)/dlam* &
                                              exp(-(ctx%state%spec_lambda - ctx%state%nebem_line_pos(i))**2/(2.0_wp*dlam**2))/ &
                                              C_LIGHT*ctx%state%nebem_line_pos(i)**2
            end do
        end if
    end subroutine compute_nebular_kernels

    subroutine load_dusty_agb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: i, j, nlam, stat
        integer :: jj, i_spec
        real(WP), allocatable :: lambda_dagb(:), fluxin_dagb(:)

        open (99, file=trim(ctx%sps_home)//'/data/dust/dusty/Orich_dusty.spec', status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (error_unit, '(A)') '[FSPS_INIT] Error: dusty O-rich models missing'
            error stop 1
        end if

        read (99, *) nlam
        allocate (lambda_dagb(nlam), fluxin_dagb(nlam))
        read (99, *) lambda_dagb

        do i = 1, NTEFF_DAGB
            do j = 1, NTAU_DAGB
                read (99, *) ctx%state%teff_dagb(1, i), ctx%state%tau1_dagb(1, j)
                read (99, *) fluxin_dagb
                jj = max(find_interval(ctx%state%spec_lambda, lambda_dagb(1)), 1)
                do i_spec = jj, ctx%state%nspec
                    ctx%state%flux_dagb(i_spec, 1, i, j) = interpolate_linear(lambda_dagb, fluxin_dagb, &
                                                                              ctx%state%spec_lambda(i_spec))
                end do
            end do
        end do
        close (99)

        deallocate (lambda_dagb, fluxin_dagb)

        open (99, file=trim(ctx%sps_home)//'/data/dust/dusty/Crich_dusty.spec', status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (error_unit, '(A)') '[FSPS_INIT] Error: dusty C-rich models missing'
            error stop 1
        end if

        read (99, *) nlam
        allocate (lambda_dagb(nlam), fluxin_dagb(nlam))
        read (99, *) lambda_dagb

        do i = 1, NTEFF_DAGB
            do j = 1, NTAU_DAGB
                read (99, *) ctx%state%teff_dagb(2, i), ctx%state%tau1_dagb(2, j)
                read (99, *) fluxin_dagb
                jj = max(find_interval(ctx%state%spec_lambda, lambda_dagb(1)), 1)
                do i_spec = jj, ctx%state%nspec
                    ctx%state%flux_dagb(i_spec, 2, i, j) = interpolate_linear(lambda_dagb, fluxin_dagb, &
                                                                              ctx%state%spec_lambda(i_spec))
                end do
            end do
        end do
        close (99)

        deallocate (lambda_dagb, fluxin_dagb)

    end subroutine load_dusty_agb_spectra

    subroutine load_xray_nebular_grid(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        integer :: i, j, k, stat, i_spec
        real(WP), dimension(NEBNIP) :: readcontneb
        real(WP), dimension(NEBNIP) :: raw_spec_log
        real(WP), dimension(NEBNAGE) :: tmp_age
        real(WP), dimension(NEBNZ) :: tmp_logz
        real(WP), dimension(NEBNIP) :: tmp_logu
        real(WP), dimension(:), allocatable :: readlambneb
        character(len=1024) :: file_path

        allocate (readlambneb(ctx%state%nspec))

        if (ctx%cloudy_dust_val == 1) then
            file_path = trim(ctx%sps_home)//'/data/nebular/ZAU_WX_WD_'//trim(ctx%state%isoc_type)//'.cont'
        else
            file_path = trim(ctx%sps_home)//'/data/nebular/ZAU_WX_ND_'//trim(ctx%state%isoc_type)//'.cont'
        end if

        open (99, file=trim(file_path), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (error_unit, '(A)') '[FSPS_INIT] Error: X-ray nebular cont file missing'
            error stop 1
        end if
        read (99, *)
        read (99, *) readlambneb
        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    read (99, *, iostat=stat) tmp_logz(i), tmp_age(j), tmp_logu(k)
                    read (99, *, iostat=stat) readcontneb
                    raw_spec_log = log10(readcontneb + 10.0_wp**(-95.0_wp))
                    do i_spec = 1, ctx%state%nspec
                        ctx%state%xnebem_cont(i_spec, i, j, k) = interpolate_linear(readlambneb, raw_spec_log, &
                                                                                      ctx%state%spec_lambda(i_spec))
                    end do
                end do
            end do
        end do
        close (99)

        if (ctx%cloudy_dust_val == 1) then
            file_path = trim(ctx%sps_home)//'/data/nebular/ZAU_WX_WD_'//trim(ctx%state%isoc_type)//'.lines'
        else
            file_path = trim(ctx%sps_home)//'/data/nebular/ZAU_WX_ND_'//trim(ctx%state%isoc_type)//'.lines'
        end if

        open (99, file=trim(file_path), status='old', action='read', iostat=stat)
        if (stat /= 0) then
            write (error_unit, '(A)') '[FSPS_INIT] Error: X-ray nebular line file missing'
            error stop 1
        end if
        read (99, *)
        read (99, *) ctx%state%nebem_line_pos
        do i = 1, NEBNZ
            do j = 1, NEBNAGE
                do k = 1, NEBNIP
                    read (99, *, iostat=stat) ctx%state%nebem_logz(i), ctx%state%nebem_age(j), ctx%state%nebem_logu(k)
                    read (99, *, iostat=stat) ctx%state%xnebem_line(:, i, j, k)
                end do
            end do
        end do
        close (99)

        ctx%state%nebem_age = log10(ctx%state%nebem_age)
        ctx%state%xnebem_line = log10(ctx%state%xnebem_line + 10.0_wp**(-95.0_wp))

        deallocate (readlambneb)
    end subroutine load_xray_nebular_grid

    subroutine load_xrb_spectra(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        type(xrb_spectra_t) :: xrb_grid
        type(backend_status_t) :: io_status
        character(len=32) :: backend_mode
        character(len=:), allocatable :: uri
        integer :: i, j

        call resolve_data_backend_uri(ctx, backend_mode, uri)

        call fsps_data_open(uri, backend_mode, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_open failed for XRB load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call fsps_data_load_xrb(xrb_grid, io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_load_xrb failed', io_status%code, trim(io_status%message)
            call fsps_data_close(io_status)
            error stop 1
        end if

        if (.not. allocated(xrb_grid%lam) .or. .not. allocated(xrb_grid%age) .or. &
            .not. allocated(xrb_grid%z) .or. .not. allocated(xrb_grid%spec)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: XRB grid is not fully allocated.'
            error stop 1
        end if

        if (size(xrb_grid%spec, 1) /= size(xrb_grid%lam) .or. size(xrb_grid%spec, 2) /= size(xrb_grid%age) .or. &
            size(xrb_grid%spec, 3) /= size(xrb_grid%z)) then
            call fsps_data_close(io_status)
            write(error_unit, '(A)') '[FSPS_INIT] Error: XRB grid shape mismatch.'
            error stop 1
        end if

        ctx%state%lam_xrb = xrb_grid%lam
        ctx%state%ages_xrb = xrb_grid%age
        ctx%state%zmet_xrb = xrb_grid%z

        do j = 1, size(xrb_grid%z)
            do i = 1, size(xrb_grid%age)
                ctx%state%spec_xrb(:, i, j) = max(interpolate_linear( &
                    xrb_grid%lam, xrb_grid%spec(:, i, j), ctx%state%spec_lambda), SAFE_FLOOR)
            end do
        end do

        ctx%state%spec_xrb = ctx%state%spec_xrb * L_SOL

        call fsps_data_close(io_status)
        if (io_status%code /= 0) then
            write(error_unit, '(A,1x,I0,1x,A)') &
                '[FSPS_INIT] Error: fsps_data_close failed for XRB load', io_status%code, trim(io_status%message)
            error stop 1
        end if

        call xrb_grid%clear()
    end subroutine load_xrb_spectra

    subroutine resolve_data_backend_uri(ctx, backend_mode, uri, dust_name)
        type(fsps_context_t), intent(in) :: ctx
        character(len=32), intent(out) :: backend_mode
        character(len=:), allocatable, intent(out) :: uri
        character(len=*), intent(in), optional :: dust_name

        character(len=64) :: backend_mode_env
        character(len=1024) :: hdf5_file_path
        character(len=1024) :: hdf5_file_path_env
        character(len=64) :: dust_part
        integer :: env_stat
        logical :: use_hdf5_uri

        backend_mode = 'legacy'
        backend_mode_env = ''
        call get_environment_variable('FSPS_DATA_BACKEND', value=backend_mode_env, status=env_stat)
        if (env_stat == 0 .and. len_trim(backend_mode_env) > 0) then
            backend_mode = trim(to_lower(trim(backend_mode_env)))
        end if

        select case (trim(backend_mode))
        case ('fsds_hdf5')
            backend_mode = 'hdf5'
        case ('fsds_legacy')
            backend_mode = 'legacy'
        case ('fsds_auto')
            backend_mode = 'auto'
        end select

        use_hdf5_uri = .false.
        select case (trim(backend_mode))
        case ('hdf5')
            use_hdf5_uri = .true.
        case ('auto')
#if FSPS_HAS_HDF5 == 1
            use_hdf5_uri = .true.
            backend_mode = 'hdf5'
#else
            backend_mode = 'legacy'
#endif
        end select

        hdf5_file_path = trim(ctx%sps_home)//'/data/fsps_data_v1.h5'
        hdf5_file_path_env = ''
        call get_environment_variable('FSPS_HDF5_DATA_PATH', value=hdf5_file_path_env, status=env_stat)
        if (env_stat == 0 .and. len_trim(hdf5_file_path_env) > 0) then
            hdf5_file_path = trim(hdf5_file_path_env)
        end if

        if (present(dust_name)) then
            dust_part = trim(dust_name)
        else
            dust_part = trim(ctx%state%str_dustem)
        end if

        if (use_hdf5_uri) then
            uri = trim(hdf5_file_path)//'|'//trim(ctx%state%isoc_type)//'|'// &
                  trim(ctx%state%spec_type)//'|'//trim(dust_part)
        else
            uri = trim(ctx%sps_home)//'|'//trim(ctx%state%isoc_type)//'|'// &
                  trim(ctx%state%spec_type)//'|'//trim(dust_part)
        end if
    end subroutine resolve_data_backend_uri

    subroutine fsps_data_open(uri, backend_mode, status)
        character(len=*), intent(in) :: uri
        character(len=*), intent(in) :: backend_mode
        type(backend_status_t), intent(out) :: status
        type(backend_status_t) :: close_status
        character(len=:), allocatable :: backend_uri

        call status%set_ok()

        if (fsps_backend_open) then
            call fsps_data_close(close_status)
            if (close_status%code /= 0) then
                call status%set_error(close_status%code, trim(close_status%message))
                return
            end if
        end if

        call create_data_backend(backend_mode, fsps_backend, status)
        if (.not. backend_status_ok(status)) return

        backend_uri = trim(uri)
        call normalize_backend_uri(backend_mode, backend_uri)

        call fsps_backend%open(backend_uri, status)
        if (.not. backend_status_ok(status)) then
            if (allocated(fsps_backend)) deallocate(fsps_backend)
            return
        end if

        call fsps_manifest%clear()
        call fsps_backend%read_manifest(fsps_manifest, status)
        if (.not. backend_status_ok(status)) then
            call fsps_backend%close(close_status)
            if (allocated(fsps_backend)) deallocate(fsps_backend)
            call fsps_manifest%clear()
            return
        end if

        fsps_backend_open = .true.
    end subroutine fsps_data_open

    subroutine fsps_data_close(status)
        type(backend_status_t), intent(out) :: status
        type(backend_status_t) :: close_status

        call status%set_ok()

        if (allocated(fsps_backend)) then
            call fsps_backend%close(close_status)
            if (.not. backend_status_ok(close_status)) then
                call status%set_error(close_status%code, trim(close_status%message))
            end if
            deallocate(fsps_backend)
        end if

        call fsps_manifest%clear()
        fsps_backend_open = .false.
    end subroutine fsps_data_close

    subroutine fsps_data_load_spectral_library(role, grid, status)
        character(len=*), intent(in) :: role
        type(spectral_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status
        type(dataset_desc_t) :: dataset

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call find_dataset_by_role(role, dataset, status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_spectral_grid(fsps_backend, dataset, grid, status)
    end subroutine fsps_data_load_spectral_library

    subroutine fsps_data_load_isochrones(grid, status)
        type(isochrone_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_isochrone_grid(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_isochrones

    subroutine fsps_data_load_nebular(component, grid, status)
        character(len=*), intent(in) :: component
        type(nebular_grid_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_nebular_grid(fsps_backend, fsps_manifest, component, grid, status)
    end subroutine fsps_data_load_nebular

    subroutine fsps_data_load_wmbasic(grid, status)
        type(aux_wmbasic_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_wmbasic(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_wmbasic

    subroutine fsps_data_load_pagb(grid, status)
        type(aux_pagb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_pagb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_pagb

    subroutine fsps_data_load_wr(grid, status)
        type(aux_wr_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_wr(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_wr

    subroutine fsps_data_load_agb(grid, status)
        type(aux_agb_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_agb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_agb

    subroutine fsps_data_load_dust_emission(grid, status)
        type(dust_emission_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_dust_emission(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_dust_emission

    subroutine fsps_data_load_agn_dust(grid, status)
        type(agn_dust_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_agn_dust(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_agn_dust

    subroutine fsps_data_load_dust_attenuation(grid, status)
        type(dust_attenuation_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_dust_attenuation(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_dust_attenuation

    subroutine fsps_data_load_xrb(grid, status)
        type(xrb_spectra_t), intent(inout) :: grid
        type(backend_status_t), intent(out) :: status

        call ensure_backend_ready(status)
        if (.not. backend_status_ok(status)) return

        call fsps_mapper%map_xrb(fsps_backend, fsps_manifest, grid, status)
    end subroutine fsps_data_load_xrb

    subroutine ensure_backend_ready(status)
        type(backend_status_t), intent(out) :: status

        call status%set_ok()

        if (.not. fsps_backend_open .or. .not. allocated(fsps_backend)) then
            call status%set_error(9101, 'FSPS data backend is not open.')
            return
        end if
    end subroutine ensure_backend_ready

    subroutine find_dataset_by_role(role, dataset, status)
        character(len=*), intent(in) :: role
        type(dataset_desc_t), intent(out) :: dataset
        type(backend_status_t), intent(out) :: status
        integer :: i

        call status%set_ok()
        call dataset%clear()

        if (.not. allocated(fsps_manifest%datasets)) then
            call status%set_error(9102, 'FSPS manifest has no datasets.')
            return
        end if

        do i = 1, size(fsps_manifest%datasets)
            if (.not. allocated(fsps_manifest%datasets(i)%role)) cycle
            if (trim(fsps_manifest%datasets(i)%role) /= trim(role)) cycle
            dataset = fsps_manifest%datasets(i)
            return
        end do

        call status%set_error(9103, 'Dataset role not found in manifest: '//trim(role))
    end subroutine find_dataset_by_role

    subroutine normalize_backend_uri(backend_mode, uri)
        character(len=*), intent(in) :: backend_mode
        character(len=:), allocatable, intent(inout) :: uri
        integer :: sep1, sep2_rel, sep3_rel, sep2, sep3
        character(len=:), allocatable :: home_part, isoc_part, spec_part, dust_part

        if (trim(to_lower(trim(backend_mode))) /= 'legacy') return

        sep1 = index(uri, '|')
        if (sep1 <= 1) return

        sep2_rel = index(uri(sep1 + 1:), '|')
        if (sep2_rel <= 1) return
        sep2 = sep1 + sep2_rel

        sep3_rel = index(uri(sep2 + 1:), '|')
        if (sep3_rel <= 1) return
        sep3 = sep2 + sep3_rel

        if (sep3 >= len_trim(uri)) return

        home_part = uri(1:sep1 - 1)
        isoc_part = uri(sep1 + 1:sep2 - 1)
        spec_part = uri(sep2 + 1:sep3 - 1)
        dust_part = uri(sep3 + 1:len_trim(uri))

        if (trim(to_lower(trim(isoc_part))) == 'bpss') then
            isoc_part = 'mist'
        end if

        if (trim(to_lower(trim(spec_part))) == 'bpass') then
            spec_part = 'miles'
        end if

        uri = trim(home_part)//'|'//trim(isoc_part)//'|'//trim(spec_part)//'|'//trim(dust_part)
    end subroutine normalize_backend_uri

    subroutine convert_to_vega_system(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: i

        do i = 1, ctx%state%nbands
            if (ctx%state%magsun(i) /= 99.0_wp) then
                ctx%state%magsun(i) = (ctx%state%magsun(i) - ctx%state%magsun(1)) - &
                                      (ctx%state%magvega(i) - ctx%state%magvega(1)) + ctx%state%magsun(1)
            end if
        end do
    end subroutine convert_to_vega_system

end module fsps_initialization
