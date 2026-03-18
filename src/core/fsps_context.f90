module fsps_context
    !> @brief
    !> Context lifecycle and parameter management for FSPS.
    !>
    !> @details
    !> Provides creation, setup, teardown, and parameter mutation for
    !> `fsps_context_t`, along with high-level SSP/CSP execution helpers.

    use fsps_precision, only: WP
    use fsps_constants, only: NEMLINE, NEBNAGE, NTABMAX
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
    public :: fsps_context_set_fast_mode
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
        ctx%fast_mode = .false.
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
        ctx%state%ssp_basis_is_dirty = .true.
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

        if (pset_changes_ssp_basis(ctx%pset, pset_in)) then
            ctx%state%ssp_basis_is_dirty = .true.
        end if
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
        character(len=64) :: key_norm

        status = 0
        key_norm = trim(key)
        select case (key_norm)
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
        case ('fast_mode')
            call fsps_context_set_fast_mode(ctx, value /= 0)
        case default
            status = FSPS_ERR_UNKNOWN_INT_PARAM
        end select

        if (status == 0 .and. is_ssp_sensitive_int_param(key_norm)) then
            ctx%state%ssp_basis_is_dirty = .true.
        end if
    end subroutine fsps_context_set_param_int

    !> @brief Enable or disable CSP fast mode.
    !> @param[inout] ctx Context to update.
    !> @param[in] fast_mode_in Fast mode flag.
    subroutine fsps_context_set_fast_mode(ctx, fast_mode_in)
        type(fsps_context_t), intent(inout) :: ctx
        logical, intent(in) :: fast_mode_in

        ctx%fast_mode = fast_mode_in
    end subroutine fsps_context_set_fast_mode

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
        character(len=64) :: key_norm

        status = 0
        key_norm = trim(key)
        select case (key_norm)
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

        if (status == 0 .and. is_ssp_sensitive_float_param(key_norm)) then
            ctx%state%ssp_basis_is_dirty = .true.
        end if
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
        character(len=64) :: key_norm

        status = 0
        key_norm = trim(key)
        select case (key_norm)
        case ('imf_filename')
            ctx%pset%imf_filename = trim(value)
        case ('sfh_filename')
            ctx%pset%sfh_filename = trim(value)
        case default
            status = FSPS_ERR_UNKNOWN_STRING_PARAM
        end select

        if (status == 0 .and. is_ssp_sensitive_str_param(key_norm)) then
            ctx%state%ssp_basis_is_dirty = .true.
        end if
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

        if (write_compsp > 0 .and. .not. ctx%fast_mode) then
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
                    !$omp target exit data map(delete: ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, &
                    !$omp                              ctx%state%ssp_basis_lbol)
                end if
                if (allocated(ctx%state%ssp_basis_spec)) deallocate(ctx%state%ssp_basis_spec)
                if (allocated(ctx%state%ssp_basis_mass)) deallocate(ctx%state%ssp_basis_mass)
                if (allocated(ctx%state%ssp_basis_lbol)) deallocate(ctx%state%ssp_basis_lbol)
            ctx%state%ssp_basis_is_dirty = .true.
        end if

        if (.not. allocated(ctx%state%ssp_basis_spec)) then
            allocate(ctx%state%ssp_basis_spec(size(spec_ssp, 1), size(spec_ssp, 2), size(spec_ssp, 3)))
            allocate(ctx%state%ssp_basis_mass(size(mass_ssp, 1), size(mass_ssp, 2)))
            allocate(ctx%state%ssp_basis_lbol(size(lbol_ssp, 1), size(lbol_ssp, 2)))
            mapped_new = .true.
        end if

        if (mapped_new) then
            !$omp target enter data map(to: ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
            ctx%state%ssp_basis_is_dirty = .true.
        end if

        if (ctx%state%ssp_basis_is_dirty) then
            ctx%state%ssp_basis_spec(:, :, 1:nzin) = spec_ssp(:, :, 1:nzin)
            ctx%state%ssp_basis_mass(:, 1:nzin) = mass_ssp(:, 1:nzin)
            ctx%state%ssp_basis_lbol(:, 1:nzin) = lbol_ssp(:, 1:nzin)

            !$omp target update to(ctx%state%ssp_basis_spec, ctx%state%ssp_basis_mass, ctx%state%ssp_basis_lbol)
            ctx%state%ssp_basis_is_dirty = .false.
        end if
    end subroutine fsps_context_update_ssp_basis

    logical function is_ssp_sensitive_int_param(key) result(is_sensitive)
        character(len=*), intent(in) :: key

        is_sensitive = .false.
        select case (trim(key))
        case ('zmet', 'wgp1', 'wgp2', 'wgp3', 'evtype', 'imf_type', 'tpagb_norm_type', &
              'interpolation_type', 'use_wr_spectra', 'add_neb_emission', 'add_neb_continuum', &
              'add_xrb_emission', 'add_stellar_remnants', 'smooth_velocity', 'smooth_lsf', &
              'smoothspec_fast', 'vactoair_flag', 'use_isoc_mdot')
            is_sensitive = .true.
        end select
    end function is_ssp_sensitive_int_param

    logical function is_ssp_sensitive_float_param(key) result(is_sensitive)
        character(len=*), intent(in) :: key

        is_sensitive = .false.
        select case (trim(key))
        case ('tiny_logt', 'imf_upper_limit', 'imf_lower_limit', 'logt_wmb_hot', 'nebular_smooth_init', &
              'imf1', 'imf2', 'imf3', 'vdmc', 'mdave', 'dell', 'delt', 'sbss', 'fbhb', 'pagb', &
              'redgb', 'agb', 'masscut', 'fcstar', 'frac_xrb', 'sigma_smooth', 'min_wave_smooth', &
              'max_wave_smooth', 'gas_logu', 'gas_logz')
            is_sensitive = .true.
        end select
    end function is_ssp_sensitive_float_param

    logical function is_ssp_sensitive_str_param(key) result(is_sensitive)
        character(len=*), intent(in) :: key

        is_sensitive = .false.
        select case (trim(key))
        case ('imf_filename')
            is_sensitive = .true.
        end select
    end function is_ssp_sensitive_str_param

    logical function pset_changes_ssp_basis(old_pset, new_pset) result(changed)
        type(PARAMS), intent(in) :: old_pset
        type(PARAMS), intent(in) :: new_pset

        changed = .false.

        if (old_pset%zmet /= new_pset%zmet) changed = .true.
        if (old_pset%wgp1 /= new_pset%wgp1) changed = .true.
        if (old_pset%wgp2 /= new_pset%wgp2) changed = .true.
        if (old_pset%wgp3 /= new_pset%wgp3) changed = .true.
        if (old_pset%evtype /= new_pset%evtype) changed = .true.

        if (old_pset%imf1 /= new_pset%imf1) changed = .true.
        if (old_pset%imf2 /= new_pset%imf2) changed = .true.
        if (old_pset%imf3 /= new_pset%imf3) changed = .true.
        if (old_pset%vdmc /= new_pset%vdmc) changed = .true.
        if (old_pset%mdave /= new_pset%mdave) changed = .true.
        if (old_pset%dell /= new_pset%dell) changed = .true.
        if (old_pset%delt /= new_pset%delt) changed = .true.
        if (old_pset%sbss /= new_pset%sbss) changed = .true.
        if (old_pset%fbhb /= new_pset%fbhb) changed = .true.
        if (old_pset%pagb /= new_pset%pagb) changed = .true.
        if (old_pset%redgb /= new_pset%redgb) changed = .true.
        if (old_pset%agb /= new_pset%agb) changed = .true.
        if (old_pset%masscut /= new_pset%masscut) changed = .true.
        if (old_pset%fcstar /= new_pset%fcstar) changed = .true.
        if (old_pset%frac_xrb /= new_pset%frac_xrb) changed = .true.
        if (old_pset%sigma_smooth /= new_pset%sigma_smooth) changed = .true.
        if (old_pset%min_wave_smooth /= new_pset%min_wave_smooth) changed = .true.
        if (old_pset%max_wave_smooth /= new_pset%max_wave_smooth) changed = .true.
        if (old_pset%gas_logu /= new_pset%gas_logu) changed = .true.
        if (old_pset%gas_logz /= new_pset%gas_logz) changed = .true.
        if (trim(old_pset%imf_filename) /= trim(new_pset%imf_filename)) changed = .true.

        if (allocated(old_pset%ssp_gen_age) .neqv. allocated(new_pset%ssp_gen_age)) then
            changed = .true.
        else if (allocated(old_pset%ssp_gen_age)) then
            if (size(old_pset%ssp_gen_age) /= size(new_pset%ssp_gen_age)) then
                changed = .true.
            else if (any(old_pset%ssp_gen_age /= new_pset%ssp_gen_age)) then
                changed = .true.
            end if
        end if
    end function pset_changes_ssp_basis

    !> @brief Set up workspace for the CSP hot path.
    !> @param[inout] ctx Context to prepare CSP workspace for.
    subroutine fsps_context_prepare_csp_workspace(ctx)
        type(fsps_context_t), intent(inout) :: ctx
        integer :: nspec, nt, nz_max

        nspec = ctx%state%nspec
        nt = ctx%state%ntfull
        nz_max = ctx%state%nz

        ! Always allocate size-independent or constant-sized buffers
        if (.not. allocated(ctx%state%csp_emlin_old)) allocate(ctx%state%csp_emlin_old(NEMLINE, nt))
        if (.not. allocated(ctx%state%csp_emlin_young)) allocate(ctx%state%csp_emlin_young(NEMLINE, nt))
        if (.not. allocated(ctx%state%csp_emlin_final)) allocate(ctx%state%csp_emlin_final(NEMLINE))
        if (.not. allocated(ctx%state%gas_current_step_lines)) allocate(ctx%state%gas_current_step_lines(NEMLINE))
        if (.not. allocated(ctx%state%scalar_reductions)) allocate(ctx%state%scalar_reductions(2))
        if (.not. allocated(ctx%state%sfh_t_calc)) allocate(ctx%state%sfh_t_calc(NTABMAX))
        if (.not. allocated(ctx%state%sfh_sfr_calc)) allocate(ctx%state%sfh_sfr_calc(NTABMAX))
        if (.not. allocated(ctx%state%sfh_age_integrand)) allocate(ctx%state%sfh_age_integrand(NTABMAX))

        ! Allocate spectrum-dependent arrays
        if (nspec > 0) then
            ! Dust Workspace
            if (.not. allocated(ctx%state%dust_transmission_diffuse)) allocate(ctx%state%dust_transmission_diffuse(nspec))
            if (.not. allocated(ctx%state%dust_frequencies)) allocate(ctx%state%dust_frequencies(nspec))
            if (.not. allocated(ctx%state%dust_spec_total_work)) allocate(ctx%state%dust_spec_total_work(nspec, nt))
            if (.not. allocated(ctx%state%dust_emission_shape)) allocate(ctx%state%dust_emission_shape(nspec))
            if (.not. allocated(ctx%state%dust_emission_final)) allocate(ctx%state%dust_emission_final(nspec, nt))

            ! Gas Workspace
            if (.not. allocated(ctx%state%gas_neb_cont_reduced)) allocate(ctx%state%gas_neb_cont_reduced(nspec, NEBNAGE))
            if (.not. allocated(ctx%state%gas_neb_line_reduced)) allocate(ctx%state%gas_neb_line_reduced(NEMLINE, NEBNAGE))
            if (.not. allocated(ctx%state%gas_current_step_cont)) allocate(ctx%state%gas_current_step_cont(nspec))

            ! CSP Spectrum Workspace
            if (.not. allocated(ctx%state%csp_igm_transmission)) allocate(ctx%state%csp_igm_transmission(nspec))
            if (.not. allocated(ctx%state%csp_spec_final)) allocate(ctx%state%csp_spec_final(nspec))
            if (.not. allocated(ctx%state%spec_young)) allocate(ctx%state%spec_young(nspec, nt))
            if (.not. allocated(ctx%state%spec_old)) allocate(ctx%state%spec_old(nspec, nt))

            ! Fused Loop Output Buffers
            if (.not. allocated(ctx%state%out_csp_spec)) allocate(ctx%state%out_csp_spec(nspec, nt))
            if (.not. allocated(ctx%state%out_csp_emlin)) allocate(ctx%state%out_csp_emlin(NEMLINE, nt))
            if (.not. allocated(ctx%state%out_mass_csp)) allocate(ctx%state%out_mass_csp(nt))
            if (.not. allocated(ctx%state%out_lbol_csp)) allocate(ctx%state%out_lbol_csp(nt))
        end if

        ! Allocate age/metallicity dependent arrays
        if (nt > 0 .and. nz_max > 0) then
            if (.not. allocated(ctx%state%csp_weights)) allocate(ctx%state%csp_weights(nt, nz_max, nt))
            if (.not. allocated(ctx%state%sfh_w_tmp1)) allocate(ctx%state%sfh_w_tmp1(max(nt, 3), nt))
            if (.not. allocated(ctx%state%sfh_w_tmp2)) allocate(ctx%state%sfh_w_tmp2(max(nt, 3), nt))
            if (.not. allocated(ctx%state%csp_ssp_lum_linear)) allocate(ctx%state%csp_ssp_lum_linear(nt, nz_max))
            if (.not. allocated(ctx%state%csp_emlin_grid)) allocate(ctx%state%csp_emlin_grid(NEMLINE, nt, nz_max))

            if (nspec > 0) then
                if (.not. allocated(ctx%state%csp_ssp_grid)) allocate(ctx%state%csp_ssp_grid(nspec, nt, nz_max))
            end if
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
        !$omp target enter data map(to: ctx)
        !$omp target enter data map(to: ctx%state)
        
        ! Copy allocatable parameter arrays in pset
        if (allocated(ctx%pset%mag_compute)) then
            !$omp target enter data map(to: ctx%pset%mag_compute)
        end if
        if (allocated(ctx%pset%ssp_gen_age)) then
            !$omp target enter data map(to: ctx%pset%ssp_gen_age)
        end if

        ! Copy pointer/allocatable components of state
        if (associated(ctx%state%indexdefined)) then
            !$omp target enter data map(to: ctx%state%indexdefined)
        end if
        if (associated(ctx%state%wgdust)) then
            !$omp target enter data map(to: ctx%state%wgdust)
        end if
        if (associated(ctx%state%g03smcextn)) then
            !$omp target enter data map(to: ctx%state%g03smcextn)
        end if
        if (associated(ctx%state%bands)) then
            !$omp target enter data map(to: ctx%state%bands)
        end if
        if (associated(ctx%state%magsun)) then
            !$omp target enter data map(to: ctx%state%magsun)
        end if
        if (associated(ctx%state%magvega)) then
            !$omp target enter data map(to: ctx%state%magvega)
        end if
        if (associated(ctx%state%filter_leff)) then
            !$omp target enter data map(to: ctx%state%filter_leff)
        end if
        if (associated(ctx%state%vega_spec)) then
            !$omp target enter data map(to: ctx%state%vega_spec)
        end if
        if (associated(ctx%state%sun_spec)) then
            !$omp target enter data map(to: ctx%state%sun_spec)
        end if
        if (associated(ctx%state%spec_lambda)) then
            !$omp target enter data map(to: ctx%state%spec_lambda)
        end if
        if (associated(ctx%state%spec_nu)) then
            !$omp target enter data map(to: ctx%state%spec_nu)
        end if
        if (associated(ctx%state%spec_res)) then
            !$omp target enter data map(to: ctx%state%spec_res)
        end if
        if (associated(ctx%state%speclib)) then
            !$omp target enter data map(to: ctx%state%speclib)
        end if
        if (associated(ctx%state%wmb_spec)) then
            !$omp target enter data map(to: ctx%state%wmb_spec)
        end if
        if (associated(ctx%state%agb_spec_o)) then
            !$omp target enter data map(to: ctx%state%agb_spec_o)
        end if
        if (associated(ctx%state%agb_logt_o)) then
            !$omp target enter data map(to: ctx%state%agb_logt_o)
        end if
        if (associated(ctx%state%agb_spec_c)) then
            !$omp target enter data map(to: ctx%state%agb_spec_c)
        end if
        if (associated(ctx%state%agb_logt_c)) then
            !$omp target enter data map(to: ctx%state%agb_logt_c)
        end if
        if (associated(ctx%state%agb_spec_car)) then
            !$omp target enter data map(to: ctx%state%agb_spec_car)
        end if
        if (associated(ctx%state%pagb_spec)) then
            !$omp target enter data map(to: ctx%state%pagb_spec)
        end if
        if (associated(ctx%state%wrn_spec)) then
            !$omp target enter data map(to: ctx%state%wrn_spec)
        end if
        if (associated(ctx%state%wrc_spec)) then
            !$omp target enter data map(to: ctx%state%wrc_spec)
        end if
        if (associated(ctx%state%qpaharr)) then
            !$omp target enter data map(to: ctx%state%qpaharr)
        end if
        if (associated(ctx%state%uminarr)) then
            !$omp target enter data map(to: ctx%state%uminarr)
        end if
        if (associated(ctx%state%lambda_dustem)) then
            !$omp target enter data map(to: ctx%state%lambda_dustem)
        end if
        if (associated(ctx%state%dustem_dustem)) then
            !$omp target enter data map(to: ctx%state%dustem_dustem)
        end if
        if (associated(ctx%state%dustem2_dustem)) then
            !$omp target enter data map(to: ctx%state%dustem2_dustem)
        end if
        if (associated(ctx%state%flux_dagb)) then
            !$omp target enter data map(to: ctx%state%flux_dagb)
        end if
        if (associated(ctx%state%nebem_cont)) then
            !$omp target enter data map(to: ctx%state%nebem_cont)
        end if
        if (associated(ctx%state%xnebem_cont)) then
            !$omp target enter data map(to: ctx%state%xnebem_cont)
        end if
        if (associated(ctx%state%neb_res_min)) then
            !$omp target enter data map(to: ctx%state%neb_res_min)
        end if
        if (associated(ctx%state%gaussnebarr)) then
            !$omp target enter data map(to: ctx%state%gaussnebarr)
        end if
        if (associated(ctx%state%agndust_spec)) then
            !$omp target enter data map(to: ctx%state%agndust_spec)
        end if
        if (associated(ctx%state%mact_isoc)) then
            !$omp target enter data map(to: ctx%state%mact_isoc)
        end if
        if (associated(ctx%state%logl_isoc)) then
            !$omp target enter data map(to: ctx%state%logl_isoc)
        end if
        if (associated(ctx%state%logt_isoc)) then
            !$omp target enter data map(to: ctx%state%logt_isoc)
        end if
        if (associated(ctx%state%logg_isoc)) then
            !$omp target enter data map(to: ctx%state%logg_isoc)
        end if
        if (associated(ctx%state%ffco_isoc)) then
            !$omp target enter data map(to: ctx%state%ffco_isoc)
        end if
        if (associated(ctx%state%phase_isoc)) then
            !$omp target enter data map(to: ctx%state%phase_isoc)
        end if
        if (associated(ctx%state%mini_isoc)) then
            !$omp target enter data map(to: ctx%state%mini_isoc)
        end if
        if (associated(ctx%state%lmdot_isoc)) then
            !$omp target enter data map(to: ctx%state%lmdot_isoc)
        end if
        if (associated(ctx%state%nmass_isoc)) then
            !$omp target enter data map(to: ctx%state%nmass_isoc)
        end if
        if (associated(ctx%state%timestep_isoc)) then
            !$omp target enter data map(to: ctx%state%timestep_isoc)
        end if
        if (associated(ctx%state%zlegend)) then
            !$omp target enter data map(to: ctx%state%zlegend)
        end if
        if (associated(ctx%state%zlegendinit)) then
            !$omp target enter data map(to: ctx%state%zlegendinit)
        end if
        if (allocated(ctx%state%spec_ssp_zz)) then
            !$omp target enter data map(to: ctx%state%spec_ssp_zz)
        end if
        if (allocated(ctx%state%mass_ssp_zz)) then
            !$omp target enter data map(to: ctx%state%mass_ssp_zz)
        end if
        if (allocated(ctx%state%lbol_ssp_zz)) then
            !$omp target enter data map(to: ctx%state%lbol_ssp_zz)
        end if
        if (associated(ctx%state%time_full)) then
            !$omp target enter data map(to: ctx%state%time_full)
        end if
        if (allocated(ctx%state%weight_ssp)) then
            !$omp target enter data map(to: ctx%state%weight_ssp)
        end if
        if (allocated(ctx%state%spec_young)) then
            !$omp target enter data map(to: ctx%state%spec_young)
        end if
        if (allocated(ctx%state%spec_old)) then
            !$omp target enter data map(to: ctx%state%spec_old)
        end if
        if (allocated(ctx%state%ssp_temp_grid)) then
            !$omp target enter data map(to: ctx%state%ssp_temp_grid)
        end if
        if (allocated(ctx%state%ssp_active_idx)) then
            !$omp target enter data map(to: ctx%state%ssp_active_idx)
        end if
        if (allocated(ctx%state%ssp_active_w)) then
            !$omp target enter data map(to: ctx%state%ssp_active_w)
        end if
        if (associated(ctx%state%bpass_spec_ssp)) then
            !$omp target enter data map(to: ctx%state%bpass_spec_ssp)
        end if
        if (associated(ctx%state%bpass_mass_ssp)) then
            !$omp target enter data map(to: ctx%state%bpass_mass_ssp)
        end if
        if (associated(ctx%state%lam_xrb)) then
            !$omp target enter data map(to: ctx%state%lam_xrb)
        end if
        if (associated(ctx%state%spec_xrb)) then
            !$omp target enter data map(to: ctx%state%spec_xrb)
        end if
        if (associated(ctx%state%ages_xrb)) then
            !$omp target enter data map(to: ctx%state%ages_xrb)
        end if
        if (associated(ctx%state%zmet_xrb)) then
            !$omp target enter data map(to: ctx%state%zmet_xrb)
        end if
        if (allocated(ctx%state%lsfinfo%lsf)) then
            !$omp target enter data map(to: ctx%state%lsfinfo%lsf)
        end if
        ! Powell and Sedfit data are usually observation data, typically not needed for simulation,
        ! but we include them to be safe if they are present.
        if (allocated(ctx%state%powell_data%mags)) then
            !$omp target enter data map(to: ctx%state%powell_data%mags)
        end if
        if (allocated(ctx%state%powell_data%magerr)) then
            !$omp target enter data map(to: ctx%state%powell_data%magerr)
        end if
        if (allocated(ctx%state%powell_data%spec)) then
            !$omp target enter data map(to: ctx%state%powell_data%spec)
        end if
        if (allocated(ctx%state%powell_data%specerr)) then
            !$omp target enter data map(to: ctx%state%powell_data%specerr)
        end if
        if (allocated(ctx%state%sedfit_data%mags)) then
            !$omp target enter data map(to: ctx%state%sedfit_data%mags)
        end if
        if (allocated(ctx%state%sedfit_data%magerr)) then
            !$omp target enter data map(to: ctx%state%sedfit_data%magerr)
        end if
        if (allocated(ctx%state%sedfit_data%spec)) then
            !$omp target enter data map(to: ctx%state%sedfit_data%spec)
        end if
        if (allocated(ctx%state%sedfit_data%specerr)) then
            !$omp target enter data map(to: ctx%state%sedfit_data%specerr)
        end if

        ! --- Permanent CSP Workspace ---
        if (allocated(ctx%state%csp_ssp_grid)) then
            !$omp target enter data map(to: ctx%state%csp_ssp_grid)
        end if
        if (allocated(ctx%state%csp_emlin_grid)) then
            !$omp target enter data map(to: ctx%state%csp_emlin_grid)
        end if
        if (allocated(ctx%state%csp_ssp_lum_linear)) then
            !$omp target enter data map(to: ctx%state%csp_ssp_lum_linear)
        end if
        if (allocated(ctx%state%csp_igm_transmission)) then
            !$omp target enter data map(to: ctx%state%csp_igm_transmission)
        end if
        if (allocated(ctx%state%csp_spec_final)) then
            !$omp target enter data map(to: ctx%state%csp_spec_final)
        end if
        if (allocated(ctx%state%csp_emlin_final)) then
            !$omp target enter data map(to: ctx%state%csp_emlin_final)
        end if
        if (allocated(ctx%state%csp_weights)) then
            !$omp target enter data map(to: ctx%state%csp_weights)
        end if
        if (allocated(ctx%state%csp_emlin_young)) then
            !$omp target enter data map(to: ctx%state%csp_emlin_young)
        end if
        if (allocated(ctx%state%csp_emlin_old)) then
            !$omp target enter data map(to: ctx%state%csp_emlin_old)
        end if

        ! --- Persistent Physics Workspaces ---
        if (allocated(ctx%state%dust_transmission_diffuse)) then
            !$omp target enter data map(to: ctx%state%dust_transmission_diffuse)
        end if
        if (allocated(ctx%state%dust_frequencies)) then
            !$omp target enter data map(to: ctx%state%dust_frequencies)
        end if
        if (allocated(ctx%state%dust_spec_total_work)) then
            !$omp target enter data map(to: ctx%state%dust_spec_total_work)
        end if
        if (allocated(ctx%state%dust_emission_shape)) then
            !$omp target enter data map(to: ctx%state%dust_emission_shape)
        end if
        if (allocated(ctx%state%dust_emission_final)) then
            !$omp target enter data map(to: ctx%state%dust_emission_final)
        end if
        if (allocated(ctx%state%gas_neb_cont_reduced)) then
            !$omp target enter data map(to: ctx%state%gas_neb_cont_reduced)
        end if
        if (allocated(ctx%state%gas_neb_line_reduced)) then
            !$omp target enter data map(to: ctx%state%gas_neb_line_reduced)
        end if
        if (allocated(ctx%state%gas_current_step_cont)) then
            !$omp target enter data map(to: ctx%state%gas_current_step_cont)
        end if
        if (allocated(ctx%state%gas_current_step_lines)) then
            !$omp target enter data map(to: ctx%state%gas_current_step_lines)
        end if
        if (allocated(ctx%state%scalar_reductions)) then
            !$omp target enter data map(to: ctx%state%scalar_reductions)
        end if
        if (allocated(ctx%state%sfh_t_calc)) then
            !$omp target enter data map(to: ctx%state%sfh_t_calc)
        end if
        if (allocated(ctx%state%sfh_sfr_calc)) then
            !$omp target enter data map(to: ctx%state%sfh_sfr_calc)
        end if
        if (allocated(ctx%state%sfh_age_integrand)) then
            !$omp target enter data map(to: ctx%state%sfh_age_integrand)
        end if
        if (allocated(ctx%state%sfh_w_tmp1)) then
            !$omp target enter data map(to: ctx%state%sfh_w_tmp1)
        end if
        if (allocated(ctx%state%sfh_w_tmp2)) then
            !$omp target enter data map(to: ctx%state%sfh_w_tmp2)
        end if

        ! --- Fused Loop Output Buffers ---
        if (allocated(ctx%state%out_csp_spec)) then
            !$omp target enter data map(to: ctx%state%out_csp_spec)
        end if
        if (allocated(ctx%state%out_csp_emlin)) then
            !$omp target enter data map(to: ctx%state%out_csp_emlin)
        end if
        if (allocated(ctx%state%out_mass_csp)) then
            !$omp target enter data map(to: ctx%state%out_mass_csp)
        end if
        if (allocated(ctx%state%out_lbol_csp)) then
            !$omp target enter data map(to: ctx%state%out_lbol_csp)
        end if
    end subroutine fsps_context_move_to_device

    !> @brief Removes the context and its data from the device.
    !> @param[inout] ctx Context to remove.
    subroutine fsps_context_remove_from_device(ctx)
        type(fsps_context_t), intent(inout) :: ctx

        associate(s => ctx%state)
            if (allocated(s%sedfit_data%specerr)) then
                !$omp target exit data map(delete: s%sedfit_data%specerr)
            end if
            if (allocated(s%sedfit_data%spec)) then
                !$omp target exit data map(delete: s%sedfit_data%spec)
            end if
            if (allocated(s%sedfit_data%magerr)) then
                !$omp target exit data map(delete: s%sedfit_data%magerr)
            end if
            if (allocated(s%sedfit_data%mags)) then
                !$omp target exit data map(delete: s%sedfit_data%mags)
            end if
            if (allocated(s%powell_data%specerr)) then
                !$omp target exit data map(delete: s%powell_data%specerr)
            end if
            if (allocated(s%powell_data%spec)) then
                !$omp target exit data map(delete: s%powell_data%spec)
            end if
            if (allocated(s%powell_data%magerr)) then
                !$omp target exit data map(delete: s%powell_data%magerr)
            end if
            if (allocated(s%powell_data%mags)) then
                !$omp target exit data map(delete: s%powell_data%mags)
            end if
            if (allocated(s%lsfinfo%lsf)) then
                !$omp target exit data map(delete: s%lsfinfo%lsf)
            end if
            if (associated(s%zmet_xrb)) then
                !$omp target exit data map(delete: s%zmet_xrb)
            end if
            if (associated(s%ages_xrb)) then
                !$omp target exit data map(delete: s%ages_xrb)
            end if
            if (associated(s%spec_xrb)) then
                !$omp target exit data map(delete: s%spec_xrb)
            end if
            if (associated(s%lam_xrb)) then
                !$omp target exit data map(delete: s%lam_xrb)
            end if
            if (associated(s%bpass_mass_ssp)) then
                !$omp target exit data map(delete: s%bpass_mass_ssp)
            end if
            if (associated(s%bpass_spec_ssp)) then
                !$omp target exit data map(delete: s%bpass_spec_ssp)
            end if
            if (allocated(s%spec_old)) then
                !$omp target exit data map(delete: s%spec_old)
            end if
            if (allocated(s%ssp_temp_grid)) then
                !$omp target exit data map(delete: s%ssp_temp_grid)
            end if
            if (allocated(s%ssp_active_idx)) then
                !$omp target exit data map(delete: s%ssp_active_idx)
            end if
            if (allocated(s%ssp_active_w)) then
                !$omp target exit data map(delete: s%ssp_active_w)
            end if
            if (allocated(s%spec_young)) then
                !$omp target exit data map(delete: s%spec_young)
            end if
            if (allocated(s%weight_ssp)) then
                !$omp target exit data map(delete: s%weight_ssp)
            end if
            if (associated(s%time_full)) then
                !$omp target exit data map(delete: s%time_full)
            end if
            if (allocated(s%lbol_ssp_zz)) then
                !$omp target exit data map(delete: s%lbol_ssp_zz)
            end if
            if (allocated(s%mass_ssp_zz)) then
                !$omp target exit data map(delete: s%mass_ssp_zz)
            end if
            if (allocated(s%spec_ssp_zz)) then
                !$omp target exit data map(delete: s%spec_ssp_zz)
            end if
            if (associated(s%zlegendinit)) then
                !$omp target exit data map(delete: s%zlegendinit)
            end if
            if (associated(s%zlegend)) then
                !$omp target exit data map(delete: s%zlegend)
            end if
            if (associated(s%timestep_isoc)) then
                !$omp target exit data map(delete: s%timestep_isoc)
            end if
            if (associated(s%nmass_isoc)) then
                !$omp target exit data map(delete: s%nmass_isoc)
            end if
            if (associated(s%lmdot_isoc)) then
                !$omp target exit data map(delete: s%lmdot_isoc)
            end if
            if (associated(s%mini_isoc)) then
                !$omp target exit data map(delete: s%mini_isoc)
            end if
            if (associated(s%phase_isoc)) then
                !$omp target exit data map(delete: s%phase_isoc)
            end if
            if (associated(s%ffco_isoc)) then
                !$omp target exit data map(delete: s%ffco_isoc)
            end if
            if (associated(s%logg_isoc)) then
                !$omp target exit data map(delete: s%logg_isoc)
            end if
            if (associated(s%logt_isoc)) then
                !$omp target exit data map(delete: s%logt_isoc)
            end if
            if (associated(s%logl_isoc)) then
                !$omp target exit data map(delete: s%logl_isoc)
            end if
            if (associated(s%mact_isoc)) then
                !$omp target exit data map(delete: s%mact_isoc)
            end if
            if (associated(s%agndust_spec)) then
                !$omp target exit data map(delete: s%agndust_spec)
            end if
            if (associated(s%gaussnebarr)) then
                !$omp target exit data map(delete: s%gaussnebarr)
            end if
            if (associated(s%neb_res_min)) then
                !$omp target exit data map(delete: s%neb_res_min)
            end if
            if (associated(s%xnebem_cont)) then
                !$omp target exit data map(delete: s%xnebem_cont)
            end if
            if (associated(s%nebem_cont)) then
                !$omp target exit data map(delete: s%nebem_cont)
            end if
            if (associated(s%flux_dagb)) then
                !$omp target exit data map(delete: s%flux_dagb)
            end if
            if (associated(s%dustem2_dustem)) then
                !$omp target exit data map(delete: s%dustem2_dustem)
            end if
            if (associated(s%dustem_dustem)) then
                !$omp target exit data map(delete: s%dustem_dustem)
            end if
            if (associated(s%lambda_dustem)) then
                !$omp target exit data map(delete: s%lambda_dustem)
            end if
            if (associated(s%uminarr)) then
                !$omp target exit data map(delete: s%uminarr)
            end if
            if (associated(s%qpaharr)) then
                !$omp target exit data map(delete: s%qpaharr)
            end if
            if (associated(s%wrc_spec)) then
                !$omp target exit data map(delete: s%wrc_spec)
            end if
            if (associated(s%wrn_spec)) then
                !$omp target exit data map(delete: s%wrn_spec)
            end if
            if (associated(s%pagb_spec)) then
                !$omp target exit data map(delete: s%pagb_spec)
            end if
            if (associated(s%agb_spec_car)) then
                !$omp target exit data map(delete: s%agb_spec_car)
            end if
            if (associated(s%agb_logt_c)) then
                !$omp target exit data map(delete: s%agb_logt_c)
            end if
            if (associated(s%agb_spec_c)) then
                !$omp target exit data map(delete: s%agb_spec_c)
            end if
            if (associated(s%agb_logt_o)) then
                !$omp target exit data map(delete: s%agb_logt_o)
            end if
            if (associated(s%agb_spec_o)) then
                !$omp target exit data map(delete: s%agb_spec_o)
            end if
            if (associated(s%wmb_spec)) then
                !$omp target exit data map(delete: s%wmb_spec)
            end if
            if (associated(s%speclib)) then
                !$omp target exit data map(delete: s%speclib)
            end if
            if (associated(s%spec_res)) then
                !$omp target exit data map(delete: s%spec_res)
            end if
            if (associated(s%spec_nu)) then
                !$omp target exit data map(delete: s%spec_nu)
            end if
            if (associated(s%spec_lambda)) then
                !$omp target exit data map(delete: s%spec_lambda)
            end if
            if (associated(s%sun_spec)) then
                !$omp target exit data map(delete: s%sun_spec)
            end if
            if (associated(s%vega_spec)) then
                !$omp target exit data map(delete: s%vega_spec)
            end if
            if (associated(s%filter_leff)) then
                !$omp target exit data map(delete: s%filter_leff)
            end if
            if (associated(s%magvega)) then
                !$omp target exit data map(delete: s%magvega)
            end if
            if (associated(s%magsun)) then
                !$omp target exit data map(delete: s%magsun)
            end if
            if (associated(s%bands)) then
                !$omp target exit data map(delete: s%bands)
            end if
            if (associated(s%g03smcextn)) then
                !$omp target exit data map(delete: s%g03smcextn)
            end if
            if (associated(s%wgdust)) then
                !$omp target exit data map(delete: s%wgdust)
            end if
            if (associated(s%indexdefined)) then
                !$omp target exit data map(delete: s%indexdefined)
            end if

            ! --- Permanent CSP Workspace ---
            if (allocated(s%ssp_basis_spec)) then
                !$omp target exit data map(delete: s%ssp_basis_spec)
            end if
            if (allocated(s%ssp_basis_mass)) then
                !$omp target exit data map(delete: s%ssp_basis_mass)
            end if
            if (allocated(s%ssp_basis_lbol)) then
                !$omp target exit data map(delete: s%ssp_basis_lbol)
            end if
            if (allocated(s%csp_ssp_grid)) then
                !$omp target exit data map(delete: s%csp_ssp_grid)
            end if
            if (allocated(s%csp_emlin_grid)) then
                !$omp target exit data map(delete: s%csp_emlin_grid)
            end if
            if (allocated(s%csp_ssp_lum_linear)) then
                !$omp target exit data map(delete: s%csp_ssp_lum_linear)
            end if
            if (allocated(s%csp_igm_transmission)) then
                !$omp target exit data map(delete: s%csp_igm_transmission)
            end if
            if (allocated(s%csp_spec_final)) then
                !$omp target exit data map(delete: s%csp_spec_final)
            end if
            if (allocated(s%csp_emlin_final)) then
                !$omp target exit data map(delete: s%csp_emlin_final)
            end if
            if (allocated(s%csp_weights)) then
                !$omp target exit data map(delete: s%csp_weights)
            end if
            if (allocated(s%csp_emlin_young)) then
                !$omp target exit data map(delete: s%csp_emlin_young)
            end if
            if (allocated(s%csp_emlin_old)) then
                !$omp target exit data map(delete: s%csp_emlin_old)
            end if

            ! --- Persistent Physics Workspaces ---
            if (allocated(s%dust_transmission_diffuse)) then
                !$omp target exit data map(delete: s%dust_transmission_diffuse)
            end if
            if (allocated(s%dust_frequencies)) then
                !$omp target exit data map(delete: s%dust_frequencies)
            end if
            if (allocated(s%dust_spec_total_work)) then
                !$omp target exit data map(delete: s%dust_spec_total_work)
            end if
            if (allocated(s%dust_emission_shape)) then
                !$omp target exit data map(delete: s%dust_emission_shape)
            end if
            if (allocated(s%dust_emission_final)) then
                !$omp target exit data map(delete: s%dust_emission_final)
            end if
            if (allocated(s%gas_neb_cont_reduced)) then
                !$omp target exit data map(delete: s%gas_neb_cont_reduced)
            end if
            if (allocated(s%gas_neb_line_reduced)) then
                !$omp target exit data map(delete: s%gas_neb_line_reduced)
            end if
            if (allocated(s%gas_current_step_cont)) then
                !$omp target exit data map(delete: s%gas_current_step_cont)
            end if
            if (allocated(s%gas_current_step_lines)) then
                !$omp target exit data map(delete: s%gas_current_step_lines)
            end if
            if (allocated(s%scalar_reductions)) then
                !$omp target exit data map(delete: s%scalar_reductions)
            end if
            if (allocated(s%sfh_t_calc)) then
                !$omp target exit data map(delete: s%sfh_t_calc)
            end if
            if (allocated(s%sfh_sfr_calc)) then
                !$omp target exit data map(delete: s%sfh_sfr_calc)
            end if
            if (allocated(s%sfh_age_integrand)) then
                !$omp target exit data map(delete: s%sfh_age_integrand)
            end if
            if (allocated(s%sfh_w_tmp1)) then
                !$omp target exit data map(delete: s%sfh_w_tmp1)
            end if
            if (allocated(s%sfh_w_tmp2)) then
                !$omp target exit data map(delete: s%sfh_w_tmp2)
            end if

            ! --- Fused Loop Output Buffers ---
            if (allocated(s%out_csp_spec)) then
                !$omp target exit data map(delete: s%out_csp_spec)
            end if
            if (allocated(s%out_csp_emlin)) then
                !$omp target exit data map(delete: s%out_csp_emlin)
            end if
            if (allocated(s%out_mass_csp)) then
                !$omp target exit data map(delete: s%out_mass_csp)
            end if
            if (allocated(s%out_lbol_csp)) then
                !$omp target exit data map(delete: s%out_lbol_csp)
            end if
        end associate

        if (allocated(ctx%pset%ssp_gen_age)) then
            !$omp target exit data map(delete: ctx%pset%ssp_gen_age)
        end if
        if (allocated(ctx%pset%mag_compute)) then
            !$omp target exit data map(delete: ctx%pset%mag_compute)
        end if

        !$omp target exit data map(delete: ctx%state)
        !$omp target exit data map(delete: ctx)
    end subroutine fsps_context_remove_from_device

end module fsps_context
