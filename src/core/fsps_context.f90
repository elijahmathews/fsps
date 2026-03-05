module fsps_context
    !> @brief
    !> Context lifecycle and parameter management for FSPS.
    !>
    !> @details
    !> Provides creation, setup, teardown, and parameter mutation for
    !> `fsps_context_t`, along with high-level SSP/CSP execution helpers.

    use fsps_precision, only: WP
    use fsps_constants, only: NEMLINE
    use fsps_types, only: PARAMS, COMPSPOUT
    use fsps_context_types, only: fsps_context_t, fsps_context_state_destroy
    use fsps_environment, only: fsps_cleanup
    use fsps_initialization, only: fsps_initialize_data
    use fsps_ssp, only: generate_ssp_grid
    use fsps_csp, only: compute_csp_scenario
    use fsps_io, only: load_tabular_sfh, write_csp_output_files

    implicit none
    private

    public :: fsps_context_create
    public :: fsps_context_setup
    public :: fsps_context_ensure_setup
    public :: fsps_context_destroy
    public :: fsps_context_set_pset
    public :: fsps_context_get_pset
    public :: fsps_context_set_param_int
    public :: fsps_context_set_param_float
    public :: fsps_context_set_param_str
    public :: fsps_context_get_paths
    public :: fsps_context_prepare_pset
    public :: fsps_context_compute_ssp
    public :: fsps_context_compute_csp
    public :: fsps_context_update_ssp_basis
    public :: fsps_context_prepare_csp_workspace
    public :: fsps_context_move_to_device
    public :: fsps_context_remove_from_device

    ! ---------------------------------------------------------------------
    ! Module constants
    ! ---------------------------------------------------------------------

    integer, parameter :: FSPS_ERR_UNKNOWN_INT_PARAM = 101
    integer, parameter :: FSPS_ERR_UNKNOWN_FLOAT_PARAM = 102
    integer, parameter :: FSPS_ERR_UNKNOWN_STRING_PARAM = 103

contains

    !> @brief Create a new FSPS context with default values.
    !> @param[out] ctx Newly initialized context.
    subroutine fsps_context_create(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        ctx%initialized = .false.
        ctx%zin = 0
        ctx%isoc_type_name = ''
        ctx%spec_type_name = ''
        ctx%dust_type_name = ''
        ctx%sps_home = ''
        ctx%data_home = ''
        ctx%output_home = ''
        ctx%om0_val = 0.27
        ctx%ol0_val = 0.73
        ctx%H0_val = 72.0
        ctx%tpagb_norm_type_val = 2
        ctx%pzcon_val = 0
        ctx%interpolation_type_val = 0
        ctx%tiny_logt_val = 0.0
        ctx%compute_light_ages_val = 0
        ctx%add_dust_emission_val = 1
        ctx%add_agn_dust_val = 1
        ctx%add_agb_dust_model_val = 1
        ctx%use_wr_spectra_val = 1
        ctx%logt_wmb_hot_val = 0.0
        ctx%add_neb_emission_val = 0
        ctx%add_neb_continuum_val = 1
        ctx%cloudy_dust_val = 0
        ctx%add_igm_absorption_val = 0
        ctx%add_xrb_emission_val = 0
        ctx%add_stellar_remnants_val = 1
        ctx%smoothspec_fast_val = 1
        ctx%smooth_velocity_val = 1
        ctx%smooth_lsf_val = 0
        ctx%dust_type_val = 0
        ctx%imf_type_val = 2
        ctx%compute_vega_mags_val = 0
        ctx%vactoair_flag_val = 0
        ctx%redshift_colors_val = 0
        ctx%use_isoc_mdot_val = 0
        ctx%setup_nebular_gaussians_val = 0
        ctx%nebular_smooth_init_val = 100.0
        ctx%nebemlineinspec_val = 1
        ctx%imf_lower_limit_val = 0.08
        ctx%imf_upper_limit_val = 120.0
    end subroutine fsps_context_create

    !> @brief Initialize a context by loading data and binding cache resources.
    !> @param[inout] ctx Context to initialize.
    !> @param[in] zin Metallicity index (-1 for all).
    !> @param[in] isoc_type_in Optional isochrone library name.
    !> @param[in] spec_type_in Optional spectral library name.
    !> @param[in] dust_type_in Optional dust model name.
    subroutine fsps_context_setup(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: zin
        character(len=*), intent(in), optional :: isoc_type_in
        character(len=*), intent(in), optional :: spec_type_in
        character(len=*), intent(in), optional :: dust_type_in

        ctx%zin = zin
        if (present(isoc_type_in)) ctx%isoc_type_name = trim(isoc_type_in)
        if (present(spec_type_in)) ctx%spec_type_name = trim(spec_type_in)
        if (present(dust_type_in)) ctx%dust_type_name = trim(dust_type_in)

        !$omp critical(fsps_setup_initialize)
        call fsps_initialize_data(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)
        !$omp end critical(fsps_setup_initialize)

        ctx%initialized = .true.
        if (len_trim(ctx%isoc_type_name) == 0) then
            ctx%isoc_type_name = trim(ctx%state%isoc_type)
        end if
        if (len_trim(ctx%spec_type_name) == 0) then
            ctx%spec_type_name = trim(ctx%state%spec_type)
        end if
        if (len_trim(ctx%dust_type_name) == 0) then
            ctx%dust_type_name = trim(ctx%state%str_dustem)
        end if
        call fsps_context_prepare_pset(ctx)
    end subroutine fsps_context_setup

    !> @brief Ensure the context is initialized with the current library selections.
    !> @param[inout] ctx Context to check and initialize if needed.
    subroutine fsps_context_ensure_setup(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        logical :: need_setup
        character(len=64) :: isoc_name
        character(len=64) :: spec_name
        character(len=64) :: dust_name

        need_setup = .false.
        isoc_name = trim(ctx%isoc_type_name)
        spec_name = trim(ctx%spec_type_name)
        dust_name = trim(ctx%dust_type_name)

        if (len_trim(isoc_name) /= 0 .and. trim(ctx%state%isoc_type) /= isoc_name) then
            need_setup = .true.
        end if

        if (len_trim(spec_name) /= 0 .and. trim(ctx%state%spec_type) /= spec_name) then
            need_setup = .true.
        end if

        if (.not. need_setup) then
            if (len_trim(isoc_name) == 0) then
                ctx%isoc_type_name = trim(ctx%state%isoc_type)
            end if
            if (len_trim(spec_name) == 0) then
                ctx%spec_type_name = trim(ctx%state%spec_type)
            end if
            if (len_trim(dust_name) == 0) then
                ctx%dust_type_name = trim(ctx%state%str_dustem)
            end if
            return
        end if

        if (len_trim(isoc_name) == 0 .and. len_trim(spec_name) == 0 .and. len_trim(dust_name) == 0) then
            call fsps_initialize_data(ctx, ctx%zin)
        else if (len_trim(spec_name) == 0 .and. len_trim(dust_name) == 0) then
            call fsps_initialize_data(ctx, ctx%zin, isoc_type_in=isoc_name)
        else if (len_trim(dust_name) == 0) then
            call fsps_initialize_data(ctx, ctx%zin, isoc_type_in=isoc_name, spec_type_in=spec_name)
        else if (len_trim(spec_name) == 0) then
            call fsps_initialize_data(ctx, ctx%zin, isoc_type_in=isoc_name, dust_type_in=dust_name)
        else
            call fsps_initialize_data(ctx, ctx%zin, isoc_type_in=isoc_name, spec_type_in=spec_name, &
                                      dust_type_in=dust_name)
        end if

        ctx%isoc_type_name = trim(ctx%state%isoc_type)
        ctx%spec_type_name = trim(ctx%state%spec_type)
        ctx%dust_type_name = trim(ctx%state%str_dustem)
    end subroutine fsps_context_ensure_setup

    !> @brief Release all resources owned by a context.
    !> @param[inout] ctx Context to destroy.
    subroutine fsps_context_destroy(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        if (ctx%initialized) then
            call fsps_cleanup(ctx)
            ctx%initialized = .false.
        end if
        call fsps_context_state_destroy(ctx%state)
    end subroutine fsps_context_destroy

    !> @brief Replace the context parameter set.
    !> @param[inout] ctx Context to update.
    !> @param[in] pset_in New parameter set.
    subroutine fsps_context_set_pset(ctx, pset_in)
        type(fsps_context_t), intent(inout) :: ctx
        type(PARAMS), intent(in) :: pset_in

        ctx%pset = pset_in
    end subroutine fsps_context_set_pset

    !> @brief Retrieve the current context parameter set.
    !> @param[in] ctx Context to read.
    !> @param[out] pset_out Current parameter set.
    subroutine fsps_context_get_pset(ctx, pset_out)
        type(fsps_context_t), intent(in) :: ctx
        type(PARAMS), intent(out) :: pset_out

        pset_out = ctx%pset
    end subroutine fsps_context_get_pset

    !> @brief Set an integer parameter by key.
    !> @param[inout] ctx Context to update.
    !> @param[in] key Parameter name.
    !> @param[in] value Parameter value.
    !> @param[out] status Non-zero on unknown parameter.
    subroutine fsps_context_set_param_int(ctx, key, value, status)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: key
        integer, intent(in) :: value
        integer, intent(out) :: status

        status = 0
        select case (trim(key))
        case ('sfh')
            ctx%pset%sfh = value
        case ('zmet')
            ctx%pset%zmet = value
        case ('wgp1')
            ctx%pset%wgp1 = value
        case ('wgp2')
            ctx%pset%wgp2 = value
        case ('wgp3')
            ctx%pset%wgp3 = value
        case ('evtype')
            ctx%pset%evtype = value
        case ('imf_type')
            ctx%imf_type_val = value
        case ('tpagb_norm_type')
            ctx%tpagb_norm_type_val = value
        case ('pzcon')
            ctx%pzcon_val = value
        case ('interpolation_type')
            ctx%interpolation_type_val = value
        case ('add_agb_dust_model')
            ctx%add_agb_dust_model_val = value
        case ('dust_type')
            ctx%dust_type_val = value
        case ('add_dust_emission')
            ctx%add_dust_emission_val = value
        case ('compute_vega_mags')
            ctx%compute_vega_mags_val = value
        case ('vactoair_flag')
            ctx%vactoair_flag_val = value
        case ('add_agn_dust')
            ctx%add_agn_dust_val = value
        case ('use_wr_spectra')
            ctx%use_wr_spectra_val = value
        case ('add_neb_emission')
            ctx%add_neb_emission_val = value
        case ('add_neb_continuum')
            ctx%add_neb_continuum_val = value
        case ('cloudy_dust')
            ctx%cloudy_dust_val = value
        case ('add_igm_absorption')
            ctx%add_igm_absorption_val = value
        case ('nebemlineinspec')
            ctx%nebemlineinspec_val = value
        case ('add_xrb_emission')
            ctx%add_xrb_emission_val = value
        case ('add_stellar_remnants')
            ctx%add_stellar_remnants_val = value
        case ('smooth_velocity')
            ctx%smooth_velocity_val = value
        case ('smooth_lsf')
            ctx%smooth_lsf_val = value
        case ('smoothspec_fast')
            ctx%smoothspec_fast_val = value
        case ('redshift_colors')
            ctx%redshift_colors_val = value
        case ('compute_light_ages')
            ctx%compute_light_ages_val = value
        case ('use_isoc_mdot')
            ctx%use_isoc_mdot_val = value
        case ('setup_nebular_gaussians')
            ctx%setup_nebular_gaussians_val = value
        case ('compute_mags')
            ctx%pset%compute_mags = value
        case ('compute_indices')
            ctx%pset%compute_indices = value
        case default
            status = FSPS_ERR_UNKNOWN_INT_PARAM
        end select
    end subroutine fsps_context_set_param_int

    !> @brief Set a floating-point parameter by key.
    !> @param[inout] ctx Context to update.
    !> @param[in] key Parameter name.
    !> @param[in] value Parameter value.
    !> @param[out] status Non-zero on unknown parameter.
    subroutine fsps_context_set_param_float(ctx, key, value, status)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: key
        real(WP), intent(in) :: value
        integer, intent(out) :: status

        status = 0
        select case (trim(key))
        case ('om0')
            ctx%om0_val = value
        case ('ol0')
            ctx%ol0_val = value
        case ('H0')
            ctx%H0_val = value
        case ('tiny_logt')
            ctx%tiny_logt_val = value
        case ('imf_upper_limit')
            ctx%imf_upper_limit_val = value
        case ('imf_lower_limit')
            ctx%imf_lower_limit_val = value
        case ('logt_wmb_hot')
            ctx%logt_wmb_hot_val = value
        case ('nebular_smooth_init')
            ctx%nebular_smooth_init_val = value
        case ('imf1')
            ctx%pset%imf1 = value
        case ('imf2')
            ctx%pset%imf2 = value
        case ('imf3')
            ctx%pset%imf3 = value
        case ('vdmc')
            ctx%pset%vdmc = value
        case ('mdave')
            ctx%pset%mdave = value
        case ('dell')
            ctx%pset%dell = value
        case ('delt')
            ctx%pset%delt = value
        case ('sbss')
            ctx%pset%sbss = value
        case ('fbhb')
            ctx%pset%fbhb = value
        case ('pagb')
            ctx%pset%pagb = value
        case ('agb_dust')
            ctx%pset%agb_dust = value
        case ('redgb')
            ctx%pset%redgb = value
        case ('agb')
            ctx%pset%agb = value
        case ('masscut')
            ctx%pset%masscut = value
        case ('fcstar')
            ctx%pset%fcstar = value
        case ('frac_xrb')
            ctx%pset%frac_xrb = value
        case ('logzsol')
            ctx%pset%logzsol = value
        case ('tau')
            ctx%pset%tau = value
        case ('const')
            ctx%pset%const = value
        case ('tage')
            ctx%pset%tage = value
        case ('fburst')
            ctx%pset%fburst = value
        case ('tburst')
            ctx%pset%tburst = value
        case ('dust1')
            ctx%pset%dust1 = value
        case ('dust2')
            ctx%pset%dust2 = value
        case ('dust3')
            ctx%pset%dust3 = value
        case ('zred')
            ctx%pset%zred = value
        case ('pmetals')
            ctx%pset%pmetals = value
        case ('dust_clumps')
            ctx%pset%dust_clumps = value
        case ('frac_nodust')
            ctx%pset%frac_nodust = value
        case ('dust_index')
            ctx%pset%dust_index = value
        case ('dust_tesc')
            ctx%pset%dust_tesc = value
        case ('frac_obrun')
            ctx%pset%frac_obrun = value
        case ('uvb')
            ctx%pset%uvb = value
        case ('mwr')
            ctx%pset%mwr = value
        case ('dust1_index')
            ctx%pset%dust1_index = value
        case ('sf_start')
            ctx%pset%sf_start = value
        case ('sf_trunc')
            ctx%pset%sf_trunc = value
        case ('sf_slope')
            ctx%pset%sf_slope = value
        case ('duste_gamma')
            ctx%pset%duste_gamma = value
        case ('duste_umin')
            ctx%pset%duste_umin = value
        case ('duste_qpah')
            ctx%pset%duste_qpah = value
        case ('sigma_smooth')
            ctx%pset%sigma_smooth = value
        case ('min_wave_smooth')
            ctx%pset%min_wave_smooth = value
        case ('max_wave_smooth')
            ctx%pset%max_wave_smooth = value
        case ('gas_logu')
            ctx%pset%gas_logu = value
        case ('gas_logz')
            ctx%pset%gas_logz = value
        case ('igm_factor')
            ctx%pset%igm_factor = value
        case ('fagn')
            ctx%pset%fagn = value
        case ('agn_tau')
            ctx%pset%agn_tau = value
        case default
            status = FSPS_ERR_UNKNOWN_FLOAT_PARAM
        end select
    end subroutine fsps_context_set_param_float

    !> @brief Set a string parameter by key.
    !> @param[inout] ctx Context to update.
    !> @param[in] key Parameter name.
    !> @param[in] value Parameter value.
    !> @param[out] status Non-zero on unknown parameter.
    subroutine fsps_context_set_param_str(ctx, key, value, status)
        type(fsps_context_t), intent(inout) :: ctx
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value
        integer, intent(out) :: status

        status = 0
        select case (trim(key))
        case ('imf_filename')
            ctx%pset%imf_filename = trim(value)
        case ('sfh_filename')
            ctx%pset%sfh_filename = trim(value)
        case default
            status = FSPS_ERR_UNKNOWN_STRING_PARAM
        end select
    end subroutine fsps_context_set_param_str

    !> @brief Retrieve resolved FSPS data and output paths.
    !> @param[in] ctx Context to read.
    !> @param[out] sps_home_out FSPS root path.
    !> @param[out] data_home_out Data directory path.
    !> @param[out] output_home_out Output directory path.
    subroutine fsps_context_get_paths(ctx, sps_home_out, data_home_out, output_home_out)
        type(fsps_context_t), intent(in) :: ctx
        character(len=*), intent(out) :: sps_home_out
        character(len=*), intent(out) :: data_home_out
        character(len=*), intent(out) :: output_home_out

        sps_home_out = trim(ctx%sps_home)
        data_home_out = trim(ctx%data_home)
        output_home_out = trim(ctx%output_home)
    end subroutine fsps_context_get_paths

    !> @brief Compute an SSP grid for the current context and parameter set.
    !> @param[inout] ctx Context to use.
    !> @param[out] mass_ssp SSP mass history.
    !> @param[out] lbol_ssp SSP bolometric luminosity history.
    !> @param[out] spec_ssp SSP spectra.
    subroutine fsps_context_compute_ssp(ctx, mass_ssp, lbol_ssp, spec_ssp)
        type(fsps_context_t), intent(inout) :: ctx
        real(WP), dimension(:), intent(out) :: mass_ssp, lbol_ssp
        real(WP), dimension(:, :), intent(out) :: spec_ssp

        call fsps_context_ensure_setup(ctx)
        call fsps_context_prepare_pset(ctx)
        call generate_ssp_grid(ctx, ctx%pset, mass_ssp, lbol_ssp, spec_ssp)
    end subroutine fsps_context_compute_ssp

    !> @brief Compute CSP outputs for the current context.
    !> @param[inout] ctx Context to use.
    !> @param[in] write_compsp Output file flag.
    !> @param[in] nzin Number of metallicities.
    !> @param[in] outfile Output filename prefix.
    !> @param[in] mass_ssp SSP mass grid.
    !> @param[in] lbol_ssp SSP bolometric luminosity grid.
    !> @param[in] spec_ssp SSP spectra grid.
    !> @param[inout] ocompsp CSP output buffer.
    subroutine fsps_context_compute_csp(ctx, write_compsp, nzin, outfile, mass_ssp, lbol_ssp, spec_ssp, ocompsp)
        type(fsps_context_t), intent(inout) :: ctx
        integer, intent(in) :: write_compsp, nzin
        character(len=*), intent(in) :: outfile
        real(WP), dimension(:, :), intent(in) :: mass_ssp, lbol_ssp
        real(WP), dimension(:, :, :), intent(in) :: spec_ssp
        type(COMPSPOUT), dimension(:), intent(inout) :: ocompsp
        type(COMPSPOUT), allocatable :: results(:)
        integer :: i, n_out, status

        call fsps_context_ensure_setup(ctx)
        call fsps_context_prepare_pset(ctx)

        if (ctx%pset%sfh == 2 .or. ctx%pset%sfh == 3) then
            call load_tabular_sfh(ctx, ctx%pset, nzin)
        end if

        call fsps_context_update_ssp_basis(ctx, spec_ssp, mass_ssp, lbol_ssp, nzin)

        call fsps_context_prepare_csp_workspace(ctx)
        call compute_csp_scenario(ctx, ctx%pset, nzin, results, status)
        if (status /= 0) then
            return
        end if

        if (write_compsp > 0) then
            call write_csp_output_files(ctx, ctx%pset, results, outfile, write_compsp)
        end if

        n_out = min(size(ocompsp), size(results))
        do i = 1, n_out
            ocompsp(i) = results(i)
        end do
        deallocate (results)
    end subroutine fsps_context_compute_csp

    !> @brief Update persistent SSP basis arrays used by CSP integration.
    !> @param[inout] ctx Context to update.
    !> @param[in] spec_ssp SSP spectra [nspec, nt, nz_in].
    !> @param[in] mass_ssp SSP mass [nt, nz_in].
    !> @param[in] lbol_ssp SSP bolometric luminosity [nt, nz_in].
    !> @param[in] nzin Number of metallicity bins actively provided.
    subroutine fsps_context_update_ssp_basis(ctx, spec_ssp, mass_ssp, lbol_ssp, nzin)
        type(fsps_context_t), intent(inout) :: ctx
        real(WP), intent(in) :: spec_ssp(:, :, :)
        real(WP), intent(in) :: mass_ssp(:, :)
        real(WP), intent(in) :: lbol_ssp(:, :)
        integer, intent(in) :: nzin
        logical :: mapped_new, needs_realloc

        mapped_new = .false.
        needs_realloc = .false.

        if (nzin <= 0) return
        if (size(spec_ssp, 3) < nzin .or. size(mass_ssp, 2) < nzin .or. size(lbol_ssp, 2) < nzin) return

        if (.not. allocated(ctx%state%ssp_basis_spec) .or. &
            .not. allocated(ctx%state%ssp_basis_mass) .or. &
            .not. allocated(ctx%state%ssp_basis_lbol)) then
            needs_realloc = .true.
        else
            ! Safe to check shape because we know they are all allocated
            if (any(shape(ctx%state%ssp_basis_spec) /= shape(spec_ssp)) .or. &
                any(shape(ctx%state%ssp_basis_mass) /= shape(mass_ssp)) .or. &
                any(shape(ctx%state%ssp_basis_lbol) /= shape(lbol_ssp))) then
                needs_realloc = .true.
            end if
        end if

        if (needs_realloc) then
                if (allocated(ctx%state%ssp_basis_spec)) then
                    !$acc exit data delete(ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
                end if
                if (allocated(ctx%state%ssp_basis_spec)) deallocate(ctx%state%ssp_basis_spec)
                if (allocated(ctx%state%ssp_basis_mass)) deallocate(ctx%state%ssp_basis_mass)
                if (allocated(ctx%state%ssp_basis_lbol)) deallocate(ctx%state%ssp_basis_lbol)
        end if

        if (.not. allocated(ctx%state%ssp_basis_spec)) then
            allocate(ctx%state%ssp_basis_spec(size(spec_ssp, 1), size(spec_ssp, 2), size(spec_ssp, 3)))
            allocate(ctx%state%ssp_basis_mass(size(mass_ssp, 1), size(mass_ssp, 2)))
            allocate(ctx%state%ssp_basis_lbol(size(lbol_ssp, 1), size(lbol_ssp, 2)))
            mapped_new = .true.
        end if

        if (mapped_new) then
            !$acc enter data copyin(ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
            !$acc enter data attach(ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
        end if

        ctx%state%ssp_basis_spec(:, :, 1:nzin) = spec_ssp(:, :, 1:nzin)
        ctx%state%ssp_basis_mass(:, 1:nzin) = mass_ssp(:, 1:nzin)
        ctx%state%ssp_basis_lbol(:, 1:nzin) = lbol_ssp(:, 1:nzin)

        !$acc update device(ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
    end subroutine fsps_context_update_ssp_basis

    !> @brief Set up workspace for the CSP hot path.
    !> @param[inout] ctx Context to prepare CSP workspace for.
    subroutine fsps_context_prepare_csp_workspace(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: nspec, nt, nz_max
        
        nspec = ctx%state%nspec
        nt = ctx%state%ntfull
        nz_max = ctx%state%nz

        if (nspec == 0 .or. nt == 0 .or. nz_max == 0) return

        if (.not. allocated(ctx%state%csp_ssp_grid)) then
            allocate(ctx%state%csp_ssp_grid(nspec, nt, nz_max))
        end if

        if (.not. allocated(ctx%state%csp_emlin_grid)) then
            allocate(ctx%state%csp_emlin_grid(NEMLINE, nt, nz_max))
        end if

        if (.not. allocated(ctx%state%csp_ssp_lum_linear)) then
            allocate(ctx%state%csp_ssp_lum_linear(nt, nz_max))
        end if

        if (.not. allocated(ctx%state%csp_igm_transmission)) then
            allocate(ctx%state%csp_igm_transmission(nspec))
        end if

        if (.not. allocated(ctx%state%csp_spec_final)) then
            allocate(ctx%state%csp_spec_final(nspec))
        end if

        if (.not. allocated(ctx%state%csp_emlin_final)) then
            allocate(ctx%state%csp_emlin_final(NEMLINE))
        end if

        if (.not. allocated(ctx%state%csp_weights)) then
            allocate(ctx%state%csp_weights(nt, nz_max))
        end if

        if (.not. allocated(ctx%state%csp_emlin_young)) then
            allocate(ctx%state%csp_emlin_young(NEMLINE))
        end if

        if (.not. allocated(ctx%state%csp_emlin_old)) then
            allocate(ctx%state%csp_emlin_old(NEMLINE))
        end if

        if (.not. allocated(ctx%state%spec_young)) then
            allocate(ctx%state%spec_young(nspec))
        end if

        if (.not. allocated(ctx%state%spec_old)) then
            allocate(ctx%state%spec_old(nspec))
        end if
    end subroutine fsps_context_prepare_csp_workspace

    !> @brief Ensure allocatable parameter-set arrays are allocated.
    !> @param[inout] ctx Context to update.
    subroutine fsps_context_prepare_pset(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        if (.not. allocated(ctx%pset%mag_compute)) then
            allocate (ctx%pset%mag_compute(ctx%state%nbands))
            ctx%pset%mag_compute = 1
        end if
        if (.not. allocated(ctx%pset%ssp_gen_age)) then
            allocate (ctx%pset%ssp_gen_age(ctx%state%nt))
            ctx%pset%ssp_gen_age = 1
        end if
    end subroutine fsps_context_prepare_pset

    !> @brief Moves the context and all its associated data to the device.
    !> @details Performs a deep copy of the context state structure.
    !> @param[inout] ctx Context to move.
    subroutine fsps_context_move_to_device(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        ! Pre-allocate CSP workspace so it gets mapped correctly
        call fsps_context_prepare_csp_workspace(ctx)

        ! Copy the main structure
        !$acc enter data copyin(ctx)
        !$acc enter data copyin(ctx%state)
        
        ! Copy allocatable parameter arrays in pset
        if (allocated(ctx%pset%mag_compute)) then
            !$acc enter data copyin(ctx%pset%mag_compute)
            !$acc enter data attach(ctx%pset%mag_compute)
        end if
        if (allocated(ctx%pset%ssp_gen_age)) then
            !$acc enter data copyin(ctx%pset%ssp_gen_age)
            !$acc enter data attach(ctx%pset%ssp_gen_age)
        end if

        ! Copy pointer/allocatable components of state
        if (associated(ctx%state%indexdefined)) then
            !$acc enter data copyin(ctx%state%indexdefined)
            !$acc enter data attach(ctx%state%indexdefined)
        end if
        if (associated(ctx%state%wgdust)) then
            !$acc enter data copyin(ctx%state%wgdust)
            !$acc enter data attach(ctx%state%wgdust)
        end if
        if (associated(ctx%state%g03smcextn)) then
            !$acc enter data copyin(ctx%state%g03smcextn)
            !$acc enter data attach(ctx%state%g03smcextn)
        end if
        if (associated(ctx%state%bands)) then
            !$acc enter data copyin(ctx%state%bands)
            !$acc enter data attach(ctx%state%bands)
        end if
        if (associated(ctx%state%magsun)) then
            !$acc enter data copyin(ctx%state%magsun)
            !$acc enter data attach(ctx%state%magsun)
        end if
        if (associated(ctx%state%magvega)) then
            !$acc enter data copyin(ctx%state%magvega)
            !$acc enter data attach(ctx%state%magvega)
        end if
        if (associated(ctx%state%filter_leff)) then
            !$acc enter data copyin(ctx%state%filter_leff)
            !$acc enter data attach(ctx%state%filter_leff)
        end if
        if (associated(ctx%state%vega_spec)) then
            !$acc enter data copyin(ctx%state%vega_spec)
            !$acc enter data attach(ctx%state%vega_spec)
        end if
        if (associated(ctx%state%sun_spec)) then
            !$acc enter data copyin(ctx%state%sun_spec)
            !$acc enter data attach(ctx%state%sun_spec)
        end if
        if (associated(ctx%state%spec_lambda)) then
            !$acc enter data copyin(ctx%state%spec_lambda)
            !$acc enter data attach(ctx%state%spec_lambda)
        end if
        if (associated(ctx%state%spec_nu)) then
            !$acc enter data copyin(ctx%state%spec_nu)
            !$acc enter data attach(ctx%state%spec_nu)
        end if
        if (associated(ctx%state%spec_res)) then
            !$acc enter data copyin(ctx%state%spec_res)
            !$acc enter data attach(ctx%state%spec_res)
        end if
        if (associated(ctx%state%speclib)) then
            !$acc enter data copyin(ctx%state%speclib)
            !$acc enter data attach(ctx%state%speclib)
        end if
        if (associated(ctx%state%wmb_spec)) then
            !$acc enter data copyin(ctx%state%wmb_spec)
            !$acc enter data attach(ctx%state%wmb_spec)
        end if
        if (associated(ctx%state%agb_spec_o)) then
            !$acc enter data copyin(ctx%state%agb_spec_o)
            !$acc enter data attach(ctx%state%agb_spec_o)
        end if
        if (associated(ctx%state%agb_logt_o)) then
            !$acc enter data copyin(ctx%state%agb_logt_o)
            !$acc enter data attach(ctx%state%agb_logt_o)
        end if
        if (associated(ctx%state%agb_spec_c)) then
            !$acc enter data copyin(ctx%state%agb_spec_c)
            !$acc enter data attach(ctx%state%agb_spec_c)
        end if
        if (associated(ctx%state%agb_logt_c)) then
            !$acc enter data copyin(ctx%state%agb_logt_c)
            !$acc enter data attach(ctx%state%agb_logt_c)
        end if
        if (associated(ctx%state%agb_spec_car)) then
            !$acc enter data copyin(ctx%state%agb_spec_car)
            !$acc enter data attach(ctx%state%agb_spec_car)
        end if
        if (associated(ctx%state%pagb_spec)) then
            !$acc enter data copyin(ctx%state%pagb_spec)
            !$acc enter data attach(ctx%state%pagb_spec)
        end if
        if (associated(ctx%state%wrn_spec)) then
            !$acc enter data copyin(ctx%state%wrn_spec)
            !$acc enter data attach(ctx%state%wrn_spec)
        end if
        if (associated(ctx%state%wrc_spec)) then
            !$acc enter data copyin(ctx%state%wrc_spec)
            !$acc enter data attach(ctx%state%wrc_spec)
        end if
        if (associated(ctx%state%qpaharr)) then
            !$acc enter data copyin(ctx%state%qpaharr)
            !$acc enter data attach(ctx%state%qpaharr)
        end if
        if (associated(ctx%state%uminarr)) then
            !$acc enter data copyin(ctx%state%uminarr)
            !$acc enter data attach(ctx%state%uminarr)
        end if
        if (associated(ctx%state%lambda_dustem)) then
            !$acc enter data copyin(ctx%state%lambda_dustem)
            !$acc enter data attach(ctx%state%lambda_dustem)
        end if
        if (associated(ctx%state%dustem_dustem)) then
            !$acc enter data copyin(ctx%state%dustem_dustem)
            !$acc enter data attach(ctx%state%dustem_dustem)
        end if
        if (associated(ctx%state%dustem2_dustem)) then
            !$acc enter data copyin(ctx%state%dustem2_dustem)
            !$acc enter data attach(ctx%state%dustem2_dustem)
        end if
        if (associated(ctx%state%flux_dagb)) then
            !$acc enter data copyin(ctx%state%flux_dagb)
            !$acc enter data attach(ctx%state%flux_dagb)
        end if
        if (associated(ctx%state%nebem_cont)) then
            !$acc enter data copyin(ctx%state%nebem_cont)
            !$acc enter data attach(ctx%state%nebem_cont)
        end if
        if (associated(ctx%state%xnebem_cont)) then
            !$acc enter data copyin(ctx%state%xnebem_cont)
            !$acc enter data attach(ctx%state%xnebem_cont)
        end if
        if (associated(ctx%state%neb_res_min)) then
            !$acc enter data copyin(ctx%state%neb_res_min)
            !$acc enter data attach(ctx%state%neb_res_min)
        end if
        if (associated(ctx%state%gaussnebarr)) then
            !$acc enter data copyin(ctx%state%gaussnebarr)
            !$acc enter data attach(ctx%state%gaussnebarr)
        end if
        if (associated(ctx%state%agndust_spec)) then
            !$acc enter data copyin(ctx%state%agndust_spec)
            !$acc enter data attach(ctx%state%agndust_spec)
        end if
        if (associated(ctx%state%mact_isoc)) then
            !$acc enter data copyin(ctx%state%mact_isoc)
            !$acc enter data attach(ctx%state%mact_isoc)
        end if
        if (associated(ctx%state%logl_isoc)) then
            !$acc enter data copyin(ctx%state%logl_isoc)
            !$acc enter data attach(ctx%state%logl_isoc)
        end if
        if (associated(ctx%state%logt_isoc)) then
            !$acc enter data copyin(ctx%state%logt_isoc)
            !$acc enter data attach(ctx%state%logt_isoc)
        end if
        if (associated(ctx%state%logg_isoc)) then
            !$acc enter data copyin(ctx%state%logg_isoc)
            !$acc enter data attach(ctx%state%logg_isoc)
        end if
        if (associated(ctx%state%ffco_isoc)) then
            !$acc enter data copyin(ctx%state%ffco_isoc)
            !$acc enter data attach(ctx%state%ffco_isoc)
        end if
        if (associated(ctx%state%phase_isoc)) then
            !$acc enter data copyin(ctx%state%phase_isoc)
            !$acc enter data attach(ctx%state%phase_isoc)
        end if
        if (associated(ctx%state%mini_isoc)) then
            !$acc enter data copyin(ctx%state%mini_isoc)
            !$acc enter data attach(ctx%state%mini_isoc)
        end if
        if (associated(ctx%state%lmdot_isoc)) then
            !$acc enter data copyin(ctx%state%lmdot_isoc)
            !$acc enter data attach(ctx%state%lmdot_isoc)
        end if
        if (associated(ctx%state%nmass_isoc)) then
            !$acc enter data copyin(ctx%state%nmass_isoc)
            !$acc enter data attach(ctx%state%nmass_isoc)
        end if
        if (associated(ctx%state%timestep_isoc)) then
            !$acc enter data copyin(ctx%state%timestep_isoc)
            !$acc enter data attach(ctx%state%timestep_isoc)
        end if
        if (associated(ctx%state%zlegend)) then
            !$acc enter data copyin(ctx%state%zlegend)
            !$acc enter data attach(ctx%state%zlegend)
        end if
        if (associated(ctx%state%zlegendinit)) then
            !$acc enter data copyin(ctx%state%zlegendinit)
            !$acc enter data attach(ctx%state%zlegendinit)
        end if
        if (allocated(ctx%state%spec_ssp_zz)) then
            !$acc enter data copyin(ctx%state%spec_ssp_zz)
            !$acc enter data attach(ctx%state%spec_ssp_zz)
        end if
        if (allocated(ctx%state%mass_ssp_zz)) then
            !$acc enter data copyin(ctx%state%mass_ssp_zz)
            !$acc enter data attach(ctx%state%mass_ssp_zz)
        end if
        if (allocated(ctx%state%lbol_ssp_zz)) then
            !$acc enter data copyin(ctx%state%lbol_ssp_zz)
            !$acc enter data attach(ctx%state%lbol_ssp_zz)
        end if
        if (associated(ctx%state%time_full)) then
            !$acc enter data copyin(ctx%state%time_full)
            !$acc enter data attach(ctx%state%time_full)
        end if
        if (allocated(ctx%state%weight_ssp)) then
            !$acc enter data copyin(ctx%state%weight_ssp)
            !$acc enter data attach(ctx%state%weight_ssp)
        end if
        if (allocated(ctx%state%spec_young)) then
            !$acc enter data copyin(ctx%state%spec_young)
            !$acc enter data attach(ctx%state%spec_young)
        end if
        if (allocated(ctx%state%spec_old)) then
            !$acc enter data copyin(ctx%state%spec_old)
            !$acc enter data attach(ctx%state%spec_old)
        end if
        if (allocated(ctx%state%ssp_temp_grid)) then
            !$acc enter data copyin(ctx%state%ssp_temp_grid)
            !$acc enter data attach(ctx%state%ssp_temp_grid)
        end if
        if (allocated(ctx%state%ssp_active_idx)) then
            !$acc enter data copyin(ctx%state%ssp_active_idx)
            !$acc enter data attach(ctx%state%ssp_active_idx)
        end if
        if (allocated(ctx%state%ssp_active_w)) then
            !$acc enter data copyin(ctx%state%ssp_active_w)
            !$acc enter data attach(ctx%state%ssp_active_w)
        end if
        if (associated(ctx%state%bpass_spec_ssp)) then
            !$acc enter data copyin(ctx%state%bpass_spec_ssp)
            !$acc enter data attach(ctx%state%bpass_spec_ssp)
        end if
        if (associated(ctx%state%bpass_mass_ssp)) then
            !$acc enter data copyin(ctx%state%bpass_mass_ssp)
            !$acc enter data attach(ctx%state%bpass_mass_ssp)
        end if
        if (associated(ctx%state%lam_xrb)) then
            !$acc enter data copyin(ctx%state%lam_xrb)
            !$acc enter data attach(ctx%state%lam_xrb)
        end if
        if (associated(ctx%state%spec_xrb)) then
            !$acc enter data copyin(ctx%state%spec_xrb)
            !$acc enter data attach(ctx%state%spec_xrb)
        end if
        if (associated(ctx%state%ages_xrb)) then
            !$acc enter data copyin(ctx%state%ages_xrb)
            !$acc enter data attach(ctx%state%ages_xrb)
        end if
        if (associated(ctx%state%zmet_xrb)) then
            !$acc enter data copyin(ctx%state%zmet_xrb)
            !$acc enter data attach(ctx%state%zmet_xrb)
        end if
        if (allocated(ctx%state%lsfinfo%lsf)) then
            !$acc enter data copyin(ctx%state%lsfinfo%lsf)
            !$acc enter data attach(ctx%state%lsfinfo%lsf)
        end if
        ! Powell and Sedfit data are usually observation data, typically not needed for simulation,
        ! but we include them to be safe if they are present.
        if (allocated(ctx%state%powell_data%mags)) then
            !$acc enter data copyin(ctx%state%powell_data%mags)
            !$acc enter data attach(ctx%state%powell_data%mags)
        end if
        if (allocated(ctx%state%powell_data%magerr)) then
            !$acc enter data copyin(ctx%state%powell_data%magerr)
            !$acc enter data attach(ctx%state%powell_data%magerr)
        end if
        if (allocated(ctx%state%powell_data%spec)) then
            !$acc enter data copyin(ctx%state%powell_data%spec)
            !$acc enter data attach(ctx%state%powell_data%spec)
        end if
        if (allocated(ctx%state%powell_data%specerr)) then
            !$acc enter data copyin(ctx%state%powell_data%specerr)
            !$acc enter data attach(ctx%state%powell_data%specerr)
        end if
        if (allocated(ctx%state%sedfit_data%mags)) then
            !$acc enter data copyin(ctx%state%sedfit_data%mags)
            !$acc enter data attach(ctx%state%sedfit_data%mags)
        end if
        if (allocated(ctx%state%sedfit_data%magerr)) then
            !$acc enter data copyin(ctx%state%sedfit_data%magerr)
            !$acc enter data attach(ctx%state%sedfit_data%magerr)
        end if
        if (allocated(ctx%state%sedfit_data%spec)) then
            !$acc enter data copyin(ctx%state%sedfit_data%spec)
            !$acc enter data attach(ctx%state%sedfit_data%spec)
        end if
        if (allocated(ctx%state%sedfit_data%specerr)) then
            !$acc enter data copyin(ctx%state%sedfit_data%specerr)
            !$acc enter data attach(ctx%state%sedfit_data%specerr)
        end if

        ! --- Permanent CSP Workspace ---
        if (allocated(ctx%state%csp_ssp_grid)) then
            !$acc enter data copyin(ctx%state%csp_ssp_grid)
            !$acc enter data attach(ctx%state%csp_ssp_grid)
        end if
        if (allocated(ctx%state%csp_emlin_grid)) then
            !$acc enter data copyin(ctx%state%csp_emlin_grid)
            !$acc enter data attach(ctx%state%csp_emlin_grid)
        end if
        if (allocated(ctx%state%csp_ssp_lum_linear)) then
            !$acc enter data copyin(ctx%state%csp_ssp_lum_linear)
            !$acc enter data attach(ctx%state%csp_ssp_lum_linear)
        end if
        if (allocated(ctx%state%csp_igm_transmission)) then
            !$acc enter data copyin(ctx%state%csp_igm_transmission)
            !$acc enter data attach(ctx%state%csp_igm_transmission)
        end if
        if (allocated(ctx%state%csp_spec_final)) then
            !$acc enter data copyin(ctx%state%csp_spec_final)
            !$acc enter data attach(ctx%state%csp_spec_final)
        end if
        if (allocated(ctx%state%csp_emlin_final)) then
            !$acc enter data copyin(ctx%state%csp_emlin_final)
            !$acc enter data attach(ctx%state%csp_emlin_final)
        end if
        if (allocated(ctx%state%csp_weights)) then
            !$acc enter data copyin(ctx%state%csp_weights)
            !$acc enter data attach(ctx%state%csp_weights)
        end if
        if (allocated(ctx%state%csp_emlin_young)) then
            !$acc enter data copyin(ctx%state%csp_emlin_young)
            !$acc enter data attach(ctx%state%csp_emlin_young)
        end if
        if (allocated(ctx%state%csp_emlin_old)) then
            !$acc enter data copyin(ctx%state%csp_emlin_old)
            !$acc enter data attach(ctx%state%csp_emlin_old)
        end if
    end subroutine fsps_context_move_to_device

    !> @brief Removes the context and its data from the device.
    !> @param[inout] ctx Context to remove.
    subroutine fsps_context_remove_from_device(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        associate(s => ctx%state)
            if (allocated(s%sedfit_data%specerr)) then
                !$acc exit data delete(s%sedfit_data%specerr)
            end if
            if (allocated(s%sedfit_data%spec)) then
                !$acc exit data delete(s%sedfit_data%spec)
            end if
            if (allocated(s%sedfit_data%magerr)) then
                !$acc exit data delete(s%sedfit_data%magerr)
            end if
            if (allocated(s%sedfit_data%mags)) then
                !$acc exit data delete(s%sedfit_data%mags)
            end if
            if (allocated(s%powell_data%specerr)) then
                !$acc exit data delete(s%powell_data%specerr)
            end if
            if (allocated(s%powell_data%spec)) then
                !$acc exit data delete(s%powell_data%spec)
            end if
            if (allocated(s%powell_data%magerr)) then
                !$acc exit data delete(s%powell_data%magerr)
            end if
            if (allocated(s%powell_data%mags)) then
                !$acc exit data delete(s%powell_data%mags)
            end if
            if (allocated(s%lsfinfo%lsf)) then
                !$acc exit data delete(s%lsfinfo%lsf)
            end if
            if (associated(s%zmet_xrb)) then
                !$acc exit data delete(s%zmet_xrb)
            end if
            if (associated(s%ages_xrb)) then
                !$acc exit data delete(s%ages_xrb)
            end if
            if (associated(s%spec_xrb)) then
                !$acc exit data delete(s%spec_xrb)
            end if
            if (associated(s%lam_xrb)) then
                !$acc exit data delete(s%lam_xrb)
            end if
            if (associated(s%bpass_mass_ssp)) then
                !$acc exit data delete(s%bpass_mass_ssp)
            end if
            if (associated(s%bpass_spec_ssp)) then
                !$acc exit data delete(s%bpass_spec_ssp)
            end if
            if (allocated(s%spec_old)) then
                !$acc exit data delete(s%spec_old)
            end if
            if (allocated(s%ssp_temp_grid)) then
                !$acc exit data delete(s%ssp_temp_grid)
            end if
            if (allocated(s%ssp_active_idx)) then
                !$acc exit data delete(s%ssp_active_idx)
            end if
            if (allocated(s%ssp_active_w)) then
                !$acc exit data delete(s%ssp_active_w)
            end if
            if (allocated(s%spec_young)) then
                !$acc exit data delete(s%spec_young)
            end if
            if (allocated(s%weight_ssp)) then
                !$acc exit data delete(s%weight_ssp)
            end if
            if (associated(s%time_full)) then
                !$acc exit data delete(s%time_full)
            end if
            if (allocated(s%lbol_ssp_zz)) then
                !$acc exit data delete(s%lbol_ssp_zz)
            end if
            if (allocated(s%mass_ssp_zz)) then
                !$acc exit data delete(s%mass_ssp_zz)
            end if
            if (allocated(s%spec_ssp_zz)) then
                !$acc exit data delete(s%spec_ssp_zz)
            end if
            if (associated(s%zlegendinit)) then
                !$acc exit data delete(s%zlegendinit)
            end if
            if (associated(s%zlegend)) then
                !$acc exit data delete(s%zlegend)
            end if
            if (associated(s%timestep_isoc)) then
                !$acc exit data delete(s%timestep_isoc)
            end if
            if (associated(s%nmass_isoc)) then
                !$acc exit data delete(s%nmass_isoc)
            end if
            if (associated(s%lmdot_isoc)) then
                !$acc exit data delete(s%lmdot_isoc)
            end if
            if (associated(s%mini_isoc)) then
                !$acc exit data delete(s%mini_isoc)
            end if
            if (associated(s%phase_isoc)) then
                !$acc exit data delete(s%phase_isoc)
            end if
            if (associated(s%ffco_isoc)) then
                !$acc exit data delete(s%ffco_isoc)
            end if
            if (associated(s%logg_isoc)) then
                !$acc exit data delete(s%logg_isoc)
            end if
            if (associated(s%logt_isoc)) then
                !$acc exit data delete(s%logt_isoc)
            end if
            if (associated(s%logl_isoc)) then
                !$acc exit data delete(s%logl_isoc)
            end if
            if (associated(s%mact_isoc)) then
                !$acc exit data delete(s%mact_isoc)
            end if
            if (associated(s%agndust_spec)) then
                !$acc exit data delete(s%agndust_spec)
            end if
            if (associated(s%gaussnebarr)) then
                !$acc exit data delete(s%gaussnebarr)
            end if
            if (associated(s%neb_res_min)) then
                !$acc exit data delete(s%neb_res_min)
            end if
            if (associated(s%xnebem_cont)) then
                !$acc exit data delete(s%xnebem_cont)
            end if
            if (associated(s%nebem_cont)) then
                !$acc exit data delete(s%nebem_cont)
            end if
            if (associated(s%flux_dagb)) then
                !$acc exit data delete(s%flux_dagb)
            end if
            if (associated(s%dustem2_dustem)) then
                !$acc exit data delete(s%dustem2_dustem)
            end if
            if (associated(s%dustem_dustem)) then
                !$acc exit data delete(s%dustem_dustem)
            end if
            if (associated(s%lambda_dustem)) then
                !$acc exit data delete(s%lambda_dustem)
            end if
            if (associated(s%uminarr)) then
                !$acc exit data delete(s%uminarr)
            end if
            if (associated(s%qpaharr)) then
                !$acc exit data delete(s%qpaharr)
            end if
            if (associated(s%wrc_spec)) then
                !$acc exit data delete(s%wrc_spec)
            end if
            if (associated(s%wrn_spec)) then
                !$acc exit data delete(s%wrn_spec)
            end if
            if (associated(s%pagb_spec)) then
                !$acc exit data delete(s%pagb_spec)
            end if
            if (associated(s%agb_spec_car)) then
                !$acc exit data delete(s%agb_spec_car)
            end if
            if (associated(s%agb_logt_c)) then
                !$acc exit data delete(s%agb_logt_c)
            end if
            if (associated(s%agb_spec_c)) then
                !$acc exit data delete(s%agb_spec_c)
            end if
            if (associated(s%agb_logt_o)) then
                !$acc exit data delete(s%agb_logt_o)
            end if
            if (associated(s%agb_spec_o)) then
                !$acc exit data delete(s%agb_spec_o)
            end if
            if (associated(s%wmb_spec)) then
                !$acc exit data delete(s%wmb_spec)
            end if
            if (associated(s%speclib)) then
                !$acc exit data delete(s%speclib)
            end if
            if (associated(s%spec_res)) then
                !$acc exit data delete(s%spec_res)
            end if
            if (associated(s%spec_nu)) then
                !$acc exit data delete(s%spec_nu)
            end if
            if (associated(s%spec_lambda)) then
                !$acc exit data delete(s%spec_lambda)
            end if
            if (associated(s%sun_spec)) then
                !$acc exit data delete(s%sun_spec)
            end if
            if (associated(s%vega_spec)) then
                !$acc exit data delete(s%vega_spec)
            end if
            if (associated(s%filter_leff)) then
                !$acc exit data delete(s%filter_leff)
            end if
            if (associated(s%magvega)) then
                !$acc exit data delete(s%magvega)
            end if
            if (associated(s%magsun)) then
                !$acc exit data delete(s%magsun)
            end if
            if (associated(s%bands)) then
                !$acc exit data delete(s%bands)
            end if
            if (associated(s%g03smcextn)) then
                !$acc exit data delete(s%g03smcextn)
            end if
            if (associated(s%wgdust)) then
                !$acc exit data delete(s%wgdust)
            end if
            if (associated(s%indexdefined)) then
                !$acc exit data delete(s%indexdefined)
            end if

            ! --- Permanent CSP Workspace ---
            if (allocated(s%ssp_basis_spec)) then
                !$acc exit data delete(s%ssp_basis_spec)
            end if
            if (allocated(s%ssp_basis_mass)) then
                !$acc exit data delete(s%ssp_basis_mass)
            end if
            if (allocated(s%ssp_basis_lbol)) then
                !$acc exit data delete(s%ssp_basis_lbol)
            end if
            if (allocated(s%csp_ssp_grid)) then
                !$acc exit data delete(s%csp_ssp_grid)
            end if
            if (allocated(s%csp_emlin_grid)) then
                !$acc exit data delete(s%csp_emlin_grid)
            end if
            if (allocated(s%csp_ssp_lum_linear)) then
                !$acc exit data delete(s%csp_ssp_lum_linear)
            end if
            if (allocated(s%csp_igm_transmission)) then
                !$acc exit data delete(s%csp_igm_transmission)
            end if
            if (allocated(s%csp_spec_final)) then
                !$acc exit data delete(s%csp_spec_final)
            end if
            if (allocated(s%csp_emlin_final)) then
                !$acc exit data delete(s%csp_emlin_final)
            end if
            if (allocated(s%csp_weights)) then
                !$acc exit data delete(s%csp_weights)
            end if
            if (allocated(s%csp_emlin_young)) then
                !$acc exit data delete(s%csp_emlin_young)
            end if
            if (allocated(s%csp_emlin_old)) then
                !$acc exit data delete(s%csp_emlin_old)
            end if
        end associate

        if (allocated(ctx%pset%ssp_gen_age)) then
            !$acc exit data delete(ctx%pset%ssp_gen_age)
        end if
        if (allocated(ctx%pset%mag_compute)) then
            !$acc exit data delete(ctx%pset%mag_compute)
        end if

        !$acc exit data delete(ctx%state)
        !$acc exit data delete(ctx)
    end subroutine fsps_context_remove_from_device

end module fsps_context
