# This is essentially an older version of make_LUT_Moon() that used lime_tbx (Toledano+24) to find lunar brightness in mags
# It WORKS (though not copypastable due to refactored variables), but ultimately discarded because
# 1) lime_tbx is slow; Allen's empirical formula is way faster for LUT generation with negligible difference
# 2) lime_tbx is unreliable at phase angle > 90deg, which introduces lunar brightness discontinuity at that regime
#    as PREFACE switches to the Allen empirical model
# 3) lime_tbx cannot be installed with pip (as of June 2026), which makes PREFACE much harder to deploy in pip
# 4) Even if it's pip-deployable, lime_tbx EocfiPath and KernelPath files takes a lot of storage space
# Though, lime_tbx is still credited for creating empirical full moon mags per filter's effective wavelengths
# Installation guide of lime_tbx as Python code can be found at README_lime_tbx_installation_guide.txt


# Creates a lookup table (LUT) containing hourly precomputed lunar brightness in mags and Mie scattering parameters
# as a function of effective band wavelength (that is pulled from another lookup table).
def make_LUT_Moon(csvcorepath, Inst, ObsStart, ObsEnd, Moon_Noise,
                  latvalue, lonvalue, altvalue, LUT_FILEPATH_ALTAZS):
    LUT_FILEPATH_MOON = rf'{csvcorepath}/Lookup_tables/lunar_TOA_Mags_at_{Inst}_{ObsStart.strftime("%b %d %Y")}_to_{ObsEnd.strftime("%b %d %Y")}.csv'
    if Moon_Noise == 'Y_Moon' and not os.path.exists(LUT_FILEPATH_MOON):
        print('[MP_Wrapper_Datacenter] Creating lookup table for moon background metric - This may take a while!')
        start = time.time()
        from lime_tbx.application.simulation.lime_simulation import LimeSimulation
        from lime_tbx.application.simulation.moon_data_factory import MoonDataFactory
        from lime_tbx.common.logger import get_logger
        from lime_tbx.common.datatypes import KernelsPath, EocfiPath, SurfacePoint
        from lime_tbx.presentation.gui.settings import SettingsManager
        from lime_tbx.persistence.local_storage.appdata import get_appdata_folder
        logger = get_logger()
        appdata_folder = get_appdata_folder(logger)

        # Get effective wavelengths and zeropoint irradiance of filters
        df_filters = pd.read_csv(rf'{csvcorepath}/filter_information.csv')
        filter_names = df_filters['filter']
        effective_wavelengths = df_filters['effective_wavelength']
        zeropoint_irradiances = df_filters['f_lamdba_zeropoint_e-11']
        allen_full_moon_mags = df_filters['allen_fullmag']

        # Sample all hours from ObsStart to ObsEnd
        n_hours = int((ObsEnd - ObsStart).total_seconds() // 3600) + 1
        hourly_times = np.array([ObsStart.replace(hour=0, minute=0, second=0) + dt.timedelta(hours=int(i))
                                 for i in np.arange(n_hours)])

        # lime_tbx setups before simulation        
        eocfi = EocfiPath(main_eocfi_path=rf'{appdata_folder}\eocfi_data',
                          custom_eocfi_path=rf'{appdata_folder}\eocfi_data')
        kernels = KernelsPath(main_kernels_path=rf'{appdata_folder}\kernels',
                              custom_kernel_path=rf'{appdata_folder}\kernels')
        settings = SettingsManager()     # To get default coefficients & SRF for fit
        default_srf = settings.get_default_srf()
        locationPoints = SurfacePoint(latitude=latvalue,
                                      longitude=lonvalue,
                                      altitude=(altvalue * u.m).to(u.km).value,
                                      dt=hourly_times)
        
        # The LIME model only reliably works when |phase angle| <= 90 deg -- remove the rest! (Also massively improves time)
        moon_data = MoonDataFactory().get_md_from_surface(locationPoints, kernels)
        phase_angles = np.array([mdata.mpa_degrees for mdata in moon_data])
        phase_mask = (np.abs(phase_angles) <= 90)

        # In order to speed up simulation time even more, remove daytime because we don't use these
        SunAlt_table = QTable.read(LUT_FILEPATH_ALTAZS, format='ascii.ecsv')
        obstimes = SunAlt_table['obstime'].datetime
        minutes = np.array([dt.minute for dt in obstimes])
        SunAlt_table = SunAlt_table[(obstimes >= min(hourly_times)) & (obstimes <= max(hourly_times)) & (minutes == 0)]
        night_mask = (SunAlt_table['sun_alt'] <= 0) # This is intentionally not -18, because hourly precision

        # Specify times to begin simulation
        hourly_times_good = hourly_times[phase_mask & night_mask]
        locationPoints_good = SurfacePoint(latitude=latvalue,
                                           longitude=lonvalue,
                                           altitude=(altvalue * u.m).to(u.km).value,
                                           dt=hourly_times_good)
        
        # Begin simulation
        sim = LimeSimulation(eocfi_path=eocfi,
                            kernels_path=kernels,
                            settings_manager=settings,
                            verbose=False)
        sim.set_simulation_changed()  # Set invalid simulation; force update irradiance simulation afterwards
        sim.update_irradiance(srf=default_srf, # THE heavy function -- times minimized
                              signals_srf=default_srf,
                              point=locationPoints_good,
                              cimel_coeff=settings.get_cimel_coef())

        # Acquire irradiance for all wavelengths in filter_information
        irradiance_spectra = sim.get_elis()  # Get irradiance spectra at all midnights_at_telescope

        wavelength_mask = np.isin(irradiance_spectra[0].wlens, np.around(effective_wavelengths))
        irradiances = np.array([spectrum.data[wavelength_mask] for spectrum in irradiance_spectra])
        irradiances = (irradiances * u.W / u.m**2 / u.nm).to(1e-11 * u.erg / u.cm**2 / u.s / u.angstrom).value
        zeropoint_irradiances_to_divide = zeropoint_irradiances.to_numpy().reshape(1, -1)
        moon_mags = -2.5 * np.log10(irradiances / zeropoint_irradiances_to_divide)

        Vmag_allen = calcMoonMag(-12.73, phase_angles)

        # Create dataframe and merge calculated moon_mags
        df_final = pd.DataFrame({'time_UTC': hourly_times,
                                 'phase_angle_deg': phase_angles,
                                 'Vmag_allen': Vmag_allen,
                                 'LIME_used': 1})
        colnames = ['time_UTC'] + [f'{f}mag' for f in filter_names]
        moon_mags_and_time = np.hstack((hourly_times_good.reshape(-1, 1), moon_mags))
        df_good_times = pd.DataFrame(moon_mags_and_time, columns=colnames)
        
        df_final = pd.merge(df_final, df_good_times, how='left', on='time_UTC')

        # Extrapolate the missing magnitudes using the Allen (1976) model whose full mag values
        # are precalculated from LIME from Oct 25 - May 26 at TNT ULTRASPEC
        first_mag_col = 4
        df_final.loc[df_final.iloc[:, first_mag_col].isna(), 'LIME_used'] = 0  # Marks use of simplified model
        for i, col in enumerate(df_final.columns[first_mag_col:]):
            allen_mags = calcMoonMag(allen_full_moon_mags.iloc[i], df_final['phase_angle_deg'])
            df_final[col] = df_final[col].astype(float).fillna(allen_mags)

        print(f'>> Lunar hourly magnitudes computed. ({time.time()-start:.2f} s)')

        # Lunar mag calculations are done! Now onto Mie scattering AERONET parameters retrieval.
        start = time.time()
        df_Mie = pd.read_csv(rf'{csvcorepath}/AERONET_AOD+INV_Level2_Daily_V3_monthly-median.csv')
        # Sort by haversine distance
        df_Mie['distance_km'] = haversine_distances(latvalue, lonvalue,
                                                    df_Mie['Site_Latitude(Degrees)'].values,
                                                    df_Mie['Site_Longitude(Degrees)'].values)
        df_Mie = df_Mie.sort_values(by=['distance_km', 'Month']).reset_index(drop=True)

        # Acquire inter/extrapolation parameters based on band effective wavelength as list of Series
        condlist = [effective_wavelengths < 440,
                    (effective_wavelengths >= 440) & (effective_wavelengths < 675),
                    effective_wavelengths >= 675]
        
        angstrom_choices = ['340-440_Angstrom_Exponent',
                            '440-675_Angstrom_Exponent',
                            '440-870_Angstrom_Exponent']
        angstrom_cols = np.select(condlist, angstrom_choices)

        wlen_choices = [440, 675, 870]
        wlen_from = np.select(condlist, wlen_choices)
        AOD_Extinction_from_cols = [f'AOD_Extinction-Total[{w}nm]' for w in wlen_from]
        AOD_Scattering_from_cols = [f'Scattering_AOD[{w}nm]' for w in wlen_from]

        g_from_choices = [440, 440, 675]
        g_to_choices = [675, 675, 870]
        gwl_from = np.select(condlist, g_from_choices)
        gwl_to = np.select(condlist, g_to_choices)

        # Further refine g_from/g_to for wl > 870
        wl_mask = (effective_wavelengths > 870)
        gwl_from[wl_mask], gwl_to[wl_mask] = 870, 1020

        g_from_cols = [f'Asymmetry_Factor-Total[{g}nm]' for g in gwl_from]
        g_to_cols = [f'Asymmetry_Factor-Total[{g}nm]' for g in gwl_to]

        # Interpolate/Extrapolate AODs and g
        AOD_Extinction_cols = [f'AOD_Extinction_{f}band' for f in filter_names]
        AOD_Scattering_cols = [f'AOD_Scattering_{f}band' for f in filter_names]
        Asymmetry_Factor_cols = [f'Asymmetry_Factor_{f}band' for f in filter_names]

        lambda_over_lambda0_alpha = ((effective_wavelengths.values / wlen_from).reshape(1, -1)) \
                                    ** (-df_Mie[angstrom_cols].values)

        df_Mie[AOD_Extinction_cols] = df_Mie[AOD_Extinction_from_cols].values * lambda_over_lambda0_alpha
        df_Mie[AOD_Scattering_cols] = df_Mie[AOD_Scattering_from_cols].values * lambda_over_lambda0_alpha
        df_Mie[Asymmetry_Factor_cols] = linear_interpolation(effective_wavelengths.values,
                                                             gwl_from, df_Mie[g_from_cols].values,
                                                             gwl_to, df_Mie[g_to_cols].values)

        # Now, pick closest AERONET Site and month and transport to df_final!
        MieColnames = AOD_Extinction_cols + AOD_Scattering_cols + Asymmetry_Factor_cols
        def getMieParameters(row):
            obsMonth = row['time_UTC'].month
            df_Mie_closest = df_Mie[df_Mie.Month == obsMonth].iloc[0]
            MieParameters = df_Mie_closest[MieColnames].values
            return MieParameters

        df_final[MieColnames] = df_final.apply(lambda row: getMieParameters(row), axis=1, result_type='expand')

        # Save to dataframe as lookup table, and we are done!
        df_final.to_csv(LUT_FILEPATH_MOON, index=False)

        print(f'>> AERONET Aerosol parameters computed. LUT exported to {LUT_FILEPATH_MOON} ({time.time()-start:.2f} s)')
    else:
        print('[MP_Wrapper_Datacenter] Moon background metric LUT already exists.')