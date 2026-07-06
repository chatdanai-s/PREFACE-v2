# Wrapper script for Part 7 of the pipeline, to set up multiprocessing jobs.
# The transit predictions + plots are made by MultiprocessingProcess, which this script imports.
from astropy.utils import data
from astropy.utils.iers import IERS_A_URL_MIRROR
import os
import time
import urllib.error
import numpy as np
import datetime as dt
import pandas as pd
import ephem

# Multiprocessing
from joblib import Parallel, delayed, parallel_config
from tqdm import tqdm  # For multiprocessing progress

# Other pipeline imports
from preface import P2_MultiprocessingProcess
from preface.P1_ModCheck import creation_date

# Timezone handling
import astropy.units as u
from astropy.coordinates import EarthLocation, CIRS
from timezonefinder import TimezoneFinder

# Moonlight background handling modules are lazy imported.
# SkyCoord lookup table handling
from astropy.time import Time, TimeDelta
import astropy.coordinates as ac
from astropy.coordinates import AltAz, get_sun, get_body, SkyCoord
from scipy.interpolate import interp1d


# Get parameters associated with UBVRI/ugriz filters
def getFilterParams(CSV_core_folder, filter_name):
    df_filters = pd.read_csv(CSV_core_folder / 'filter_information.csv')
    row_used = df_filters.loc[df_filters['filter'] == filter_name].iloc[0]
    effective_wavelength = float(row_used['effective_wavelength'])      # Retrieved from Bessell+98 and Fukugita+96
    zeropoint_irradiance = float(row_used['f_lamdba_zeropoint_e-11'])   # Retrieved from Bessell+98 and Fukugita+96
    allen_full_moon_mag  = float(row_used['allen_fullmag'])  # Sampled from LIME_tbx results generated during Oct2025-May2026 for TNT ULTRASPEC

    return effective_wavelength, zeropoint_irradiance, allen_full_moon_mag

# Calculates apparent magnitude of the Moon based on Allen 1976, page 44
def calcMoonMag(full_moon_mag, phase_angle_deg):
    allen_mag = full_moon_mag + 0.026 * np.abs(phase_angle_deg) + 4e-9 * phase_angle_deg**4
    return allen_mag

# Input in degrees
def haversine_distances(lat1, lon1, lat2_arr, lon2_arr):
    R = 6371.0  # Earth radius in km
    lat1, lon1 = np.radians(lat1), np.radians(lon1)
    lat2_arr, lon2_arr = np.radians(lat2_arr), np.radians(lon2_arr)

    dlat = lat2_arr - lat1
    dlon = lon2_arr - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2_arr) * np.sin(dlon/2)**2
    distances = 2 * R * np.arcsin(np.sqrt(a))
    return distances

# Interpolation for asymmetry factor
def linear_interpolation(x, x1, y1, x2, y2):
    x = x.reshape(1, -1)
    y = y1 + (x - x1) * (y2 - y1)/(x2 - x1)
    return y

# Cublic spline interpolation for altaz from obs_times to dense_times
# Converts to cartesian coordinates first to safely handle:
# Az wraparound (0/360) and Alt near zenith (~90 deg) singularity
def interpolate_altaz(obs_times, altazs, dense_times):
    alt_rad = altazs.alt.rad
    az_rad  = altazs.az.rad

    # AltAz -> Cartesian unit vectors
    x = np.cos(alt_rad) * np.cos(az_rad)
    y = np.cos(alt_rad) * np.sin(az_rad)
    z = np.sin(alt_rad)

    jd_coarse = obs_times.jd
    jd_dense  = dense_times.jd

    # Interpolate each Cartesian component independently
    x_i = interp1d(jd_coarse, x, kind='cubic', fill_value='extrapolate')(jd_dense)
    y_i = interp1d(jd_coarse, y, kind='cubic', fill_value='extrapolate')(jd_dense)
    z_i = interp1d(jd_coarse, z, kind='cubic', fill_value='extrapolate')(jd_dense)

    # Re-normalise (cubic spline doesn't preserve unit length)
    norm = np.sqrt(x_i**2 + y_i**2 + z_i**2)
    x_i /= norm
    y_i /= norm
    z_i /= norm

    # Cartesian -> AltAz
    alt_interp = np.rad2deg(np.arcsin(np.clip(z_i, -1, 1)))
    az_interp  = np.rad2deg(np.arctan2(y_i, x_i)) % 360

    return alt_interp, az_interp


# Phase Two needs a cached IERS table to perform co-ordinate transformations, which auto-updates every 7 days.
# A manual download on a shorter timescale is needed to prevent multiple checks/downloads when multiprocessing.
# The hash used to denote the file also updates every 7 days, so this must also be tracked and updated in a .txt file.
def download_IERS(CSV_core_folder):
    from astropy.utils import iers
    iers.conf.auto_max_age = None
    iers.conf.auto_download = False

    # IERS_A_URL = 'https://datacenter.iers.org/data/9/finals2000A.all' # Website does not respond as of 30 June 2026
    IERS_A_URL = 'https://maia.usno.navy.mil/ser7/finals2000A.all'
    IERS_PATH  = CSV_core_folder / 'IERSPath.txt'

    if data.is_url_in_cache(IERS_A_URL):
        cached_path = data.download_file(IERS_A_URL, cache=True, show_progress=True)  # Get the name of the cached file.
        with open(IERS_PATH, 'w') as IERS_File:    # Write location and new hash to .txt file for future ref.
            c = IERS_File.write(str(cached_path))

    else:
        print('[P2_MultiprocessingWrapper] Retrieving IERS table from mirror...')
        try:    # Try/except blocks are needed in case of bad internet connection. #JustEduroamThings
            r = urllib.request.Request(IERS_A_URL_MIRROR, headers={'User-Agent': 'astropy/iers', 'Accept': '*/*'})
            u_= urllib.request.urlopen(r, timeout=30).read().decode('utf-8')
            mirror_path = CSV_core_folder / 'finals2000A.all'
            with open(IERS_PATH, 'w') as IERS_File:
                c = IERS_File.write(str(u_))
            with open(IERS_PATH, 'w') as IERS_File:
                c = IERS_File.write(str(mirror_path))
            print('>> Retrieval complete.')

        except urllib.error.URLError:
            print('>> The request timed out; trying again...')
            cached_path = data.download_file(IERS_A_URL, cache=True, show_progress=True)  # Get the new file.
            with open(IERS_PATH, 'w') as IERS_File:
                c = IERS_File.write(str(cached_path))
            print('>> Retrieval complete!')

    # Now that we definitely have a cached copy, we need to check its age.
    with open(IERS_PATH, 'r') as IERS_File:
       c = IERS_File.read().replace('\n','')        # Reads name of cached IERS table.

    IERS_creation = creation_date(IERS_PATH) 
    IERS_creation = dt.datetime.fromtimestamp(IERS_creation)  # Converts date to human-readable timestamp.

    if (not os.path.isfile(c)) or (dt.datetime.now() - IERS_creation >= dt.timedelta(days=7)):  # If old, we update.
        def retrieveIERS():
            data.clear_download_cache(hashorurl=IERS_A_URL)                 # Otherwise the next line just returns the name and location of cached file!
            cached_path = data.download_file(IERS_A_URL, cache=True, show_progress=True)
            with open(IERS_PATH, 'w') as IERS_File:                         # Write location and new hash to .txt file for future ref.
                c = IERS_File.write(str(cached_path))
            print('>> Retrieval complete!')

        print('[P2_MultiprocessingWrapper] Retrieving updated IERS table: please stand by...')
        try:
            retrieveIERS()
        except urllib.error.URLError:
            print('>> The request timed out; trying again...')
            retrieveIERS()
    else:
        print('[P2_MultiprocessingWrapper] IERS table is up to date.')

    # Explicitly load the cached table so astropy uses this exact file instead of
    # lazily resolving IERS_Auto on first transform
    with open(IERS_PATH, 'r') as IERS_File:
        c = IERS_File.read().replace('\n','')
    table = iers.IERS_A.open(c)
    iers.earth_orientation_table.set(table)


# Make directory for output graphs, if it doesn't exist already.
def make_OutputGraphDir(output_folder, instrument, filter_name, obs_start, obs_end):
    folder_name = rf'{instrument}_{filter_name}-band_{obs_start.strftime("%b-%d-%Y")}_to_{obs_end.strftime("%b-%d-%Y")}'
    graph_folder_path = (output_folder / "phase_2" / "graph" / folder_name)
    graph_folder_path.mkdir(parents=True, exist_ok=True) 
    return graph_folder_path

# Find telescope location and therefore local time zone for local time handling
# As much as I would like to just store the utc offset in hours, daylight savings is unfortunately a thing.
def get_LocAndTimezoneStr(scope_df, scope_idx):
    Loc = EarthLocation(lat=scope_df['Lat'][scope_idx] * u.deg,
                        lon=scope_df['Long'][scope_idx] * u.deg,
                        height=scope_df['Alt'][scope_idx] * u.m)  # Location info of telescope
    tf = TimezoneFinder()
    latvalue, lonvalue, altvalue = Loc.lat.value, Loc.lon.value, Loc.height.value
    timezone_str = tf.timezone_at(lat=latvalue, lng=lonvalue)

    return Loc, timezone_str


# Creates lookup tables (LUT) containing local AltAz positions of the Sun and Moon
# First create a table with 10-min precision, then interpolate for 1-min precision.
def make_LUT_AltAzs(CSV_intermediate_folder, instrument, obs_start: dt.datetime, obs_end: dt.datetime, Loc):
    LUT_FOLDER = CSV_intermediate_folder / "lut_altaz"
    LUT_FOLDER.mkdir(parents=True, exist_ok=True)
    LUT_FILENAME = f'sun_moon_altazs_at_{instrument}_TDB{obs_start.strftime("%b-%d-%Y")}_to_TDB{obs_end.strftime("%b-%d-%Y")}.parquet'
    LUT_FILEPATH_ALTAZS = LUT_FOLDER / LUT_FILENAME

    if not os.path.exists(LUT_FILEPATH_ALTAZS):
        print('[P2_MultiprocessingWrapper] Creating lookup table for Sun and Moon position in local observer frame. This may take a while!')

        obs_start_LUT = Time(obs_start - dt.timedelta(days=2), format='datetime', scale='tdb')
        obs_end_LUT = Time(obs_end + dt.timedelta(days=2), format='datetime', scale='tdb')

        timestep_10m = TimeDelta(60*10, format='sec')  # 10-minute sampling
        n_steps_10m = int(((obs_end_LUT - obs_start_LUT) / timestep_10m).value) + 1
        obstimes_10m = obs_start_LUT + timestep_10m * np.arange(n_steps_10m)

        timestep_1m = TimeDelta(60, format='sec')    # 1-minute sampling
        n_steps_1m = int(((obs_end_LUT - obs_start_LUT) / timestep_1m).value) + 1
        obstimes_1m = obs_start_LUT + timestep_1m * np.arange(n_steps_1m)

        frame = AltAz(obstime=obstimes_10m, location=Loc)  # Converts from sky frame (RA/Dec) to local observer frame (Alt/Az)

        # Where is the sun at obstimes?
        start = time.time()

        sun = get_sun(obstimes_10m)
        sun_SkyCoord = SkyCoord(ra=sun.ra, dec=sun.dec, unit="deg", frame='icrs')
        sunaltazs = sun_SkyCoord.transform_to(frame)

        time_taken = time.time() - start
        print(f'>> AltAz pairs created for Sun. ({time_taken:.2f} s)')

        # Where is the Moon at obstimes?
        start = time.time()

        moon = get_body("moon", obstimes_10m)
        moon_SkyCoord = SkyCoord(ra=moon.ra, dec=moon.dec, unit="deg", frame='icrs')
        moonaltazs = moon_SkyCoord.transform_to(frame)   # AltAz pairs created for moon

        time_taken = time.time() - start
        print(f'>> AltAz pairs created for Moon. ({time_taken:.2f} s)')

        # Interpolate altazs
        sun_alt, sun_az = interpolate_altaz(obstimes_10m, sunaltazs, obstimes_1m)
        moon_alt, moon_az = interpolate_altaz(obstimes_10m, moonaltazs, obstimes_1m)


        # Build dataframe with these altazs and export
        start = time.time()

        df_altazs = pd.DataFrame({
            "obstime": obstimes_1m.to_datetime(),
            "sun_alt": np.round(sun_alt, 3).astype(np.float32),
            "sun_az": np.round(sun_az, 3).astype(np.float32),
            "moon_alt": np.round(moon_alt, 3).astype(np.float32),
            "moon_az": np.round(moon_az, 3).astype(np.float32),
        })
        df_altazs.to_parquet(LUT_FILEPATH_ALTAZS, engine="pyarrow", compression="snappy", index=False)

        time_taken = time.time() - start
        print(f'>> AltAz LUT exported to {LUT_FILEPATH_ALTAZS}. ({time_taken:.2f} s)')

    else:
        print('[P2_MultiprocessingWrapper] Appropriate Sun-Moon AltAz LUT already exists.')

    return LUT_FILEPATH_ALTAZS


# Creates a lookup table (LUT) containing hourly precomputed lunar brightness in mags and Mie scattering parameters
# as a function of effective band wavelength (that is pulled from another lookup table).
def make_LUT_Moon(CSV_core_folder, CSV_intermediate_folder, instrument, obs_start, obs_end, toggle_moonlight_noise):
    LUT_FOLDER = CSV_intermediate_folder / "lut_moonmags"
    LUT_FOLDER.mkdir(parents=True, exist_ok=True)
    LUT_FILENAME = f'lunar_TOA_Mags_at_{instrument}_UTC{obs_start.strftime("%b-%d-%Y")}_to_UTC{obs_end.strftime("%b-%d-%Y")}.parquet'
    LUT_FILEPATH_MOON = LUT_FOLDER / LUT_FILENAME

    if toggle_moonlight_noise and not os.path.exists(LUT_FILEPATH_MOON):
        print('[P2_MultiprocessingWrapper] Creating lookup table for moon background metric. This may take a while!')

        start = time.time()

        # Sample all hours from obs_start to obs_end
        obs_start_LUT = obs_start - dt.timedelta(days=2)
        obs_start_LUT.replace(hour=0, minute=0, second=0)

        obs_end_LUT = obs_end + dt.timedelta(days=2)
        obs_end_LUT.replace(hour=0, minute=0, second=0)

        n_hours = int((obs_end_LUT - obs_start_LUT).total_seconds() // 3600) + 1

        hourly_times = np.array([obs_start_LUT + dt.timedelta(hours=i)
                                 for i in range(n_hours)])
        
        # Get phase angles
        phase_angles = np.zeros_like(hourly_times)
        moon = ephem.Moon()
        for i, t in enumerate(hourly_times):
            moon.compute(t)  # Assume TDB ~ UTC, which ~10 min error is sufficient for hourly precision
            fraction_of_illumination = moon.moon_phase

            phase_angle_deg = np.degrees(np.arccos(2*fraction_of_illumination - 1))
            phase_angles[i] = phase_angle_deg

        # Get filter information
        df_filters = pd.read_csv(CSV_core_folder / 'filter_information.csv')
        filter_names = df_filters['filter']
        allen_full_moon_mags = df_filters['allen_fullmag']

        # Create dataframe and get lunar magnitudes per phase angles
        df_moon = pd.DataFrame({
            "time_UTC": hourly_times,
            "phase_angle_deg": phase_angles.astype(np.float32),
        })
        for filter_name, full_mag in zip(filter_names, allen_full_moon_mags):
            df_moon[f"{filter_name}mag"] = calcMoonMag(full_mag, phase_angles).astype(np.float32)


        # Save to dataframe as lookup table, and we are done!
        df_moon.to_parquet(LUT_FILEPATH_MOON, engine="pyarrow", compression="snappy", index=False)

        time_taken = time.time()-start
        print(f'>> Moonlight mag. LUT exported to {LUT_FILEPATH_MOON} ({time_taken:.2f} s)')

    elif toggle_moonlight_noise == False:
        pass
    else:
        print('[P2_MultiprocessingWrapper] Moon background metric LUT already exists.')

    return LUT_FILEPATH_MOON


# Queues and executes Phase Two transit predictor for all planets.
def P2Wrap(CSV_core_folder, CSV_intermediate_folder, output_folder,
           scope_df, scope_idx, instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode, viable_cumulative_cut,
           toggle_moonlight_noise, scattering_aod, absorption_aod, asymmetry_factor, moonlight_amplification_factor,
           toggle_graph_outputs, event_weight_graph_threshold,
           obs_start, obs_end,
           toggle_multiprocessing, cores_used):

    print('[P2_MultiprocessingWrapper] Pre-multiprocessing procedures initiating...')

    download_IERS(CSV_core_folder)
    graph_folder_path = make_OutputGraphDir(output_folder, instrument, filter_name, obs_start, obs_end)
    Loc, timezone_str = get_LocAndTimezoneStr(scope_df, scope_idx)
    LUT_FILEPATH_ALTAZS = make_LUT_AltAzs(CSV_intermediate_folder, instrument, obs_start, obs_end, Loc)
    LUT_FILEPATH_MOON = make_LUT_Moon(CSV_core_folder, CSV_intermediate_folder, instrument, obs_start, obs_end, toggle_moonlight_noise) 


    # Multiprocessing start
    if toggle_multiprocessing == True:
        cores_actually_used = cores_used
        print(f'\n+++ Phase Two: Multiprocessing engaged for {cores_actually_used} cores. +++\n')
    elif toggle_multiprocessing == False:
        cores_actually_used = 1
        print(f'\n+++ Phase Two: Engaged for 1 core. +++\n')
    

    # Directory of all Phase Two Inputs (dont touch!)
    sky_noise_text = 'Y-SkyNoise' if toggle_sky_noise else 'N-SkyNoise'
    defocus_text = 'Y-Defocus' if toggle_defocus else 'N-Defocus'
    config_str = f'{instrument}_{filter_name}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}'

    filename_pattern = (
        rf"{config_str}_{metric_mode}-Mode_{viable_cumulative_cut*100}%_Cut_*.csv"
    )
    jobs = list((CSV_intermediate_folder / "phase_2_inputs").glob(filename_pattern))

    # Start multiprocessing and wrap with tqdm to get progress bar
    with parallel_config(backend='loky', prefer='processes'):
        job = Parallel(
            n_jobs=cores_actually_used,
            pre_dispatch='all',
            return_as='generator_unordered'
        )(
            delayed(P2_MultiprocessingProcess.P2Predictor)(
                CSV_core_folder, csv_initiated, output_folder,
                obs_start, obs_end, scope_df, scope_idx,
                instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode,
                toggle_moonlight_noise, scattering_aod, absorption_aod, asymmetry_factor, moonlight_amplification_factor,
                toggle_graph_outputs, event_weight_graph_threshold,
                Loc, timezone_str, graph_folder_path,
                LUT_FILEPATH_ALTAZS, LUT_FILEPATH_MOON
            ) for csv_initiated in jobs
        )

        for _ in tqdm(job, total=len(jobs), desc=">> CSVs completed"):
            pass

    print(f'\n+++ Multiprocessing complete! +++\n')

    return filename_pattern, cores_actually_used
    