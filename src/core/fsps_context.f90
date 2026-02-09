module fsps_context
    !> @brief
    !> Context lifecycle and parameter management for FSPS.
    !>
    !> @details
    !> Provides creation, setup, teardown, and parameter mutation for
    !> `fsps_context_t`, along with high-level SSP/CSP execution helpers.

    use fsps_precision, only: WP
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
        type(fsps_context_t), intent(out) :: ctx

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

        call fsps_initialize_data(ctx, zin, isoc_type_in, spec_type_in, dust_type_in)

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

        call compute_csp_scenario(ctx, ctx%pset, nzin, spec_ssp, mass_ssp, lbol_ssp, results, status)
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

end module fsps_context
