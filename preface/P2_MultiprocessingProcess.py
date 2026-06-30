# Part 7 of overall pipeline.
# Also known as Phase Two, transit predictor.
# Take in .csv for each planet and find viable events in a given observing window.
# Now telescope, filter, location and time-dependent.

# Over-riding astropy_mpl_style makes things pretty!
import sys
import os
import datetime as dt
from zoneinfo import ZoneInfo

import pandas as pd
import numpy as np
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.style as mplstyle
from matplotlib.colors import LinearSegmentedColormap, Normalize
import ephem

# Switch matplotlib to non-interactive environment and other performance commands
plt.switch_backend('Agg')
mplstyle.use('fast')
rcParams['path.simplify'] = True
rcParams['path.simplify_threshold'] = 1.0

from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
plt.rcParams.update({'font.size': 18})  # Font size over-ride
from astropy.time import Time
import astropy.units as u
import astropy.coordinates as ac
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from dateutil.relativedelta import relativedelta

# Forces astropy to use IERS cache
from astropy.utils import iers
iers.conf.auto_max_age = None
iers.conf.auto_download = False

# Cache all LUTs outside main function, as to not reopen everytime
_lut_cache = {}
def _get_shared_luts(LUT_FILEPATH_ALTAZS, LUT_FILEPATH_MOON, toggle_moonlight_noise, Loc):
    cache_key = (str(LUT_FILEPATH_ALTAZS),
                 str(LUT_FILEPATH_MOON))

    if cache_key not in _lut_cache:
        # This block only runs the FIRST time this worker sees this path.
        altaz_LUT = pd.read_parquet(LUT_FILEPATH_ALTAZS, engine="pyarrow")
        obstime_LUT = Time(altaz_LUT["obstime"].to_numpy(dtype=str), format="isot")
        sunaltazs_LUT = ac.SkyCoord(alt=altaz_LUT['sun_alt'], az=altaz_LUT['sun_az'], unit='deg',
                                    obstime=obstime_LUT, location=Loc, frame='altaz')
        moonaltazs_LUT = ac.SkyCoord(alt=altaz_LUT['moon_alt'], az=altaz_LUT['moon_az'], unit='deg',
                                    obstime=obstime_LUT, location=Loc, frame='altaz')

        moon_LUT = pd.read_parquet(LUT_FILEPATH_MOON, engine="pyarrow") if toggle_moonlight_noise else None

        _lut_cache[cache_key] = (sunaltazs_LUT, moonaltazs_LUT, moon_LUT)

    # Every subsequent call (2nd, 3rd, ... Nth CSV handled by this worker) hits this line directly.
    return _lut_cache[cache_key]


# Altitude constants needed
minAlt   = 30        # Min obs. altitude for air mass = 2.
maxAlt   = 88        # Max obs. altitude for the Liverpool instrument
nightAlt = -18       # Marker for night start/end altitude

# Blank arrays for use with time cut/observation ranking
transit_start_times = []
transit_start_times_err = []
x_hold_fall = []
x_hold_rise = []
minAlt_cross_hold = []

# Output plot configurations
PLOT_DECIMATION = 2   # Display-only thinning of target/moon path points
# Define custom colormap
colors = ['blue', 'aqua']
cmap = LinearSegmentedColormap.from_list('target_cmap', colors, N=180)
norm = Normalize(vmin=0, vmax=360)


# Specify sky co-ordinates of target appropriate units.
def AltChecker(row): 
    LockOn = ac.SkyCoord(
        row['RA:HH'] + (row['RA:MM']/60) + (row['RA:SS']/3600),
        row['Dec:DD'] + np.sign(row['Dec:DD']) * (row['Dec:MM']/60 + row['Dec:SS']/3600),
        unit=(u.hourangle, u.deg),
        frame='icrs'
    )
    return LockOn.ra.hour, LockOn.dec.degree

# Conversion from BJD_TDB to UTC
def BJD2UTC(time, LockOn, Loc):
    time_BJD = Time(time, scale='tdb', location=Loc)   # BJD_TDB (implicit)
    ltt_bary = time_BJD.light_travel_time(LockOn)      # JD to BJD conversion factor (light travel time)
    JD_TDB = time_BJD - ltt_bary                       # JD_TDB
    JD_UTC = JD_TDB.utc                                # JD_UTC

    return JD_UTC

# Given T1 (datetime), what is the closest midnight (Time)?
def closest_midnight(T1):
    if T1.hour >= 12:
        midnight = T1.replace(hour=0, minute=0, second=0, microsecond=0) + relativedelta(days=1)
    else:
        midnight = T1.replace(hour=0, minute=0, second=0, microsecond=0)
    midnight = Time(midnight)   # Astropy object conversion for compatibility.
    return midnight


# Master function for Phase Two - calculation of event times/metrics.
def P2Predictor(
    CSV_core_folder, csv_initiated, output_folder,
    obs_start, obs_end, scope_df, scope_idx,
    instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode,
    toggle_moonlight_noise, scattering_aod, absorption_aod, asymmetry_factor, moonlight_amplification_factor,
    toggle_graph_outputs, event_weight_graph_threshold,
    Loc, timezone_str, graph_folder_path,
    LUT_FILEPATH_ALTAZS, LUT_FILEPATH_MOON
    ):

    # Filename of this process
    csv_name = os.path.basename(csv_initiated)
    # Get relevant (cached) LUTs
    sunaltazs_LUT, moonaltazs_LUT, moon_LUT = _get_shared_luts(
        LUT_FILEPATH_ALTAZS, LUT_FILEPATH_MOON, toggle_moonlight_noise, Loc
    )


    # Internal conversion to recover midpoint of first observation in BJD_TDB.    
    def T1_Calc(row):
        
        period_d = row['P (day)']
        period_sec = (period_d * u.d).to(u.s).value
        eph_obj = Time(row['T0 (HJD or BJD)'], format='jd').datetime    # Scale is not JD, but the format (245XXXX.XXX...) is!
        T0_sec_err = (row['T0_err'] * u.d).to(u.s)      # Conversion of ephemeris error.
        P_sec_err = (row['P_err'] * u.d).to(u.s)        # Conversion of error in planet's determined period.

        # Recovers transit start time from published midpoint, returning a timedelta object.
        transit_start = eph_obj - dt.timedelta(0.5 * row['T14'])
        
        # Builds range of epoch numbers relevant for target/window combination.
        # epoch_start finds final transit before our window begins, so you gotta add +1!
        epoch_start = (obs_start - transit_start).total_seconds() // period_sec + 1 
        epoch_end = (obs_end - transit_start).total_seconds() // period_sec
        epoch_range = np.arange(epoch_start, epoch_end+1)

        # Increase the number of orbits n to find future transit start times in our window.
        for epoch in epoch_range:
            # Calculation of individual transit start times
            transit_start_time = transit_start + (epoch * dt.timedelta(period_d))
            transit_start_times.append(transit_start_time)
            
            # Calculation of total error at epoch n.
            # Uncertainty should scale overall as no. of orbits n, so this needs to be inside the squared term!
            transit_start_time_err = np.sqrt(T0_sec_err**2 + (epoch * P_sec_err)**2).value              
            transit_start_times_err.append(transit_start_time_err)
    

    # Calculation of key times during transit.
    def Times_Calc(row):
        # Locked to index 0 except T1, because the other values are NaN at index 1 onwards
        period_d = df['P (day)'][0]
        Rs = df['R*'][0] * u.solRad
        Rp = (df['Rp'][0] * u.R_jup).to(u.solRad)
        T1 = row['T1']
        T14 = df['T14'][0]
        b = df['Impact Parameter'][0]
        a_Calc = df['a_Calc'][0]
        
        # Calculates midpoint (T0), end of transit (T4) and length of mid-transit event (T23).
        T0 = T1 + dt.timedelta(0.5 * T14)                
        T4 = T1 + dt.timedelta(T14)   

        try:
            with np.errstate(invalid='raise'):        
                T23 = (period_d/np.pi) * np.arcsin(
                      np.sqrt((Rs - Rp)**2 - (b * Rs)**2) / (a_Calc * u.AU).to(u.solRad)
                )
                # Recovers length of ingress/egress. abs is a necessary guess if T23 > T14 (unphysical soln.)     
                T12 = abs(T14 - T23.value) / 2
                T2 = T1 + dt.timedelta(T12)     # Time at end of ingress
                T3 = T4 - dt.timedelta(T12)     # Time at start of egress
                
                # TESS targets have dodgy parameters - this calculation will catch them.
                # Forces a FloatingPointError so the transit is treated as grazing, passed to except block.
                np.sqrt((T3 - T2).total_seconds())
        

        # Grazing transits have no T23, so a different approach is needed.
        except FloatingPointError:
            T23 = 0      
            T2 = T0  # Time at end of ingress
            T3 = T0  # Time at start of egress

        except ValueError:
            sys.exit(f"Time error for {csv_name} ({planet_name}) transiting at {T1} - flag me!")

        # A solid baseline is important for a good transit curve!
        baseline_left = T1 - dt.timedelta(hours=1)
        baseline_right = T4 + dt.timedelta(hours=1) 

        return T2, T0, T3, T4, baseline_left, baseline_right


    # Adds moon noise metric which is essentially SNR with yes moon / SNR with no moon, at start of observation.
    # This works because the metrics are directly proportional to SNR.
    def MoonNoiseMetric(lunar_altaz, target_altaz,
                        moon_mag, AOD_scatter, AOD_absorption, asymmetry_factor,
                        dif_BS, delta_midnight, dif_BE, filter_name):
        from preface.P1_RankMaker import findSkyB
        from preface.P2_MultiprocessingWrapper import getFilterParams

        # Acquire effective wavelength
        effective_wavelength, _, _ = getFilterParams(CSV_core_folder, filter_name)

        # Calculate moon background in mag/arcsec^2, then to mag at aperture much like the same way as skyB.
        # Based on Winkler (2022) model.
        lunar_zenith = lunar_altaz.zen.rad
        target_zenith = target_altaz.zen.rad
        coord_targ = ac.SkyCoord(target_altaz.az, target_altaz.alt, unit="deg")
        coord_moon = ac.SkyCoord(lunar_altaz.az, lunar_altaz.alt, unit="deg")
        theta = coord_targ.separation(coord_moon).to(u.rad).value

        # timemask removes out-of-observation values AND zmask removes unphysical values when zeniths go beyond 90 deg
        timemask = (delta_midnight.value >= dif_BS) & (delta_midnight.value <= dif_BE)
        max_zenith = np.pi/2 - 0.002  # 0.002 prevents gradation underflow to zero
        zmask = (lunar_zenith <= max_zenith) & (target_zenith <= max_zenith)
        mask = timemask & zmask

        # If zmask excludes all lunar_zenith, then it basically means the moon and target is
        # never present over the horizon together during observation day at all -- Consider neglible SNR change.
        if zmask.any() == False:
            return 1
        
        # If zmask AND timemask excludes all times, then consider masking only for zenith
        # This will give a lower bound for moon noise metric score instead.
        elif mask.any() == False:
            mask = zmask

        lunar_zenith = lunar_zenith[mask]
        target_zenith = target_zenith[mask]
        theta = theta[mask]

        # Scattering parameters
        scaleHeight = 8500       # Of Earth, in m
        g = asymmetry_factor     # Asymmetry factor for in Mie scattering
        tau_R = np.exp(-Loc.height.value/scaleHeight) * (1.229e+10) * effective_wavelength**(-4.05)  # Eqn 13
        tau_M = AOD_scatter
        tau_sum = tau_R + tau_M  # Total scattering optical depth
        tau_tot = tau_sum + AOD_absorption  # Total optical depth

        # Secant of zenith angles and other precomputed constants
        sec_zeta = 1/np.cos(target_zenith)
        sec_z = 1/np.cos(lunar_zenith)

        gradation = sec_zeta * (np.exp(-tau_tot*sec_zeta) - np.exp(-tau_tot*sec_z)) / (sec_z - sec_zeta)

        # Rayleigh and Mie scattering phase functions calculation, eqns 9 and 10.
        phase_R = 0.0571 * (1.06 + np.cos(theta)**2)  # Rayleigh, take chi = 0.0148 and precalculated (9)
        phase_M = (1 - g**2) / (1 + g**2 - 2*g*np.cos(theta))**(1.5) / (4 * np.pi)  # Mie
        phase = (tau_R * phase_R + tau_M * phase_M) / tau_sum   # Total phase function (1/sr)
        phase = (phase * 1/u.sr).to(1/u.arcsec**2).value

        # Calculate moonlight intensity decrease at target from ONE-time scattering, has unit 1/arcsec^2
        moon_intensity_decrease = phase * gradation * tau_sum/tau_tot

        # With a little bit of math, you can disregard calling of zeropoint magnitude completely
        # when converting from mag to flux, then flux/arcsec^2 back to mag/arcsec^2.
        # Returns moon contribution to sky background in mag/arcsec^2.
        with np.errstate(divide='ignore', invalid='ignore'):
            moon_intensity_mag = moon_mag - 2.5 * np.log10(moon_intensity_decrease)
        moon_intensity_mag[~np.isfinite(moon_intensity_mag)] = np.nan
        
        # Calculate noise metric
        moon_flux = 10**(-0.4 * (
            findSkyB(moon_intensity_mag, scope_df, scope_idx, toggle_defocus, run_mode) - moonlight_amplification_factor
            )
        )
        if toggle_sky_noise == True:
            sky_flux = 10**(-0.4 * findSkyB(scope_df[f'msky_{filter_name}'].iloc[scope_idx],
                                            scope_df, scope_idx, toggle_defocus, run_mode))
        else:
            sky_flux = 0
        target_flux = 10**(-0.4 * df[f'{filter_name}mag'][0])
        
        MoonNoiseMetric = 1 / np.sqrt(1 + moon_flux/(sky_flux + target_flux))
        return np.nanmin(MoonNoiseMetric)  # Take min value (signifies highest SNR reduction)


    ##### WE GET INTO CALCULATION AND PLOTTING HERE #####
    def EventMetric(row):
        from preface.P2_MultiprocessingWrapper import interpolate_altaz
 
        # Chooses appropriate midnight for each useful time (for graphing purposes).
        # Use relativedelta to avoid end of months/years breaking the code.
        T1 = row['T1']
        T14 = df['T14'][0]  # It's NAN in second row onwards
        targ_samples = 1801

        # Set closest midnight (BJD) associated with start of transit T1
        midnight = closest_midnight(T1)

        # Retrieve relevent delta_midnights and observation times in range (midday to evening)
        delta_midnight_targ = np.linspace(-12, 18, targ_samples) * u.hour  # Time sampled 30h from midnight for target (1 min)
        delta_midnight_sep  = np.linspace(-12, 18, 21) * u.hour            # Time sampled 30h from midnight for moon sep (90 min)

        targ_times = midnight + delta_midnight_targ   # Times for Moon/target calc
        sep_times  = midnight + delta_midnight_sep    # Times for Moon/target separation calc


        # Where is the sun and moon from midday to evening? 
        # LUTs are used here because it makes the pipeline way faster
        obstime_num = sunaltazs_LUT.obstime.mjd  # Convert to array for indexing mask
        targ_times_num = targ_times.mjd
        sep_times_num  = sep_times.mjd

        targ_timemask = (obstime_num >= targ_times_num[0]) & (obstime_num <= targ_times_num[-1])
        sep_timemask  = np.isin(obstime_num, sep_times_num)

        sun_altazs = sunaltazs_LUT[targ_timemask]      # AltAz pairs created for sun (1 min)
        moon_altazs = moonaltazs_LUT[targ_timemask]    # AltAz pairs created for moon (1 min)
        moonsep_altazs = moonaltazs_LUT[sep_timemask]  # AltAz pairs created for moon (90 min)


        # Where is the target from midday to evening?
        # Retrieve altazs for 10-minute precision then interpolate much like Sun and Moon
        targ_times_10m = targ_times_10m_list[row.name]
        targ_altazs_10m = targ_altazs_10m_by_row[row.name]
        targ_alt_1m, targ_az_1m = interpolate_altaz(targ_times_10m, targ_altazs_10m, targ_times)
        targ_altazs = ac.SkyCoord(alt=targ_alt_1m, az=targ_az_1m, unit="deg",
                                  obstime=targ_times, location=Loc, frame='altaz')
        
        sep_times_in_targ_times = np.isin(targ_times_num, sep_times_num)
        targsep_altazs = targ_altazs[sep_times_in_targ_times]

        # Calculates separation between target and Moon at 21 intervals.
        coord_targ = ac.SkyCoord(targsep_altazs.az, targsep_altazs.alt, unit="deg")
        coord_moon = ac.SkyCoord(moonsep_altazs.az, moonsep_altazs.alt, unit="deg")

        seps = coord_targ.separation(coord_moon).value
        seps_int = seps.astype(int)


        ### LUNAR ILLUMINATION CALC (Returns moon phase in %)
        m = ephem.Moon()
        m.compute(T1)
        phase = m.moon_phase * 100


        ### CLEARANCE PROTOCOLS FOR MULTIPLE-EVENT RUNNING
        x_hold_fall.clear()
        x_hold_rise.clear()
        minAlt_cross_hold.clear()

        # When will the transit, baseline and night begin and end -- relative to midnight? 
        # Finds differences either side of 00:00:00.
        # Note: Take only the value when plotting! You can't plot a unit.
        DIF_PRECISION = 6
        midnight_mjd = midnight.mjd  # float, compute once

        def hours_relative_to_midnight(time):
            dif = (Time(time).mjd - midnight_mjd) * 24  # MJD is in days, *24 = hours
            return round(dif, DIF_PRECISION)

        dif_BS = hours_relative_to_midnight(row['Baseline_Start'])
        dif_T1 = hours_relative_to_midnight(T1)
        dif_T2 = hours_relative_to_midnight(row['T2'])
        dif_T0 = hours_relative_to_midnight(row['T0'])
        dif_T3 = hours_relative_to_midnight(row['T3'])
        dif_T4 = hours_relative_to_midnight(row['T4'])
        dif_BE = hours_relative_to_midnight(row['Baseline_End'])


        # When does each event night start and end? 
        # Checks when sun path crosses twilight marker at y=-18, appends to x_hold_fall or _rise (depending on gradient).
        # x_holds_fall and _rise indicate delta_midnight in minutes where sunrise-sunfall crossover occurs
        nightAlt_cross_check = sun_altazs.alt.value - nightAlt
        nightAlt_cross_check_idxs = np.arange(len(nightAlt_cross_check) - 1)

        for i in nightAlt_cross_check_idxs:
            if nightAlt_cross_check[i] * nightAlt_cross_check[i+1] <= 0:       # If the moment sun crosses over y=-18 at i to i+1
                x_ = delta_midnight_targ.to('min').value[i]  # Time in delta_midnight_targ (minutes) where the crossing happens 

                if sun_altazs.alt.value[i+1] < sun_altazs.alt.value[i]:      # If sunfall
                    x_hold_fall.append(x_)
                elif sun_altazs.alt.value[i+1] >= sun_altazs.alt.value[i]:   # If sunrise
                    x_hold_rise.append(x_)

        # For "two-night events", with a bit at the start and end!
        # It is the second condition in particular that catches most of these.
        first_sunfall_cut = (x_hold_rise[0] < x_hold_fall[0])
        nfalls, nrises = len(x_hold_fall), len(x_hold_rise)

        T0_before_UTC_midnight = (dif_T0 < 0)
        T0_after_UTC_midnight = (dif_T0 > 0)

        # All condition keys are coded such that (T)ransit, sun(R)ise, sun(F)all.
        # (Fall) |---(Transit)--- (Rise) ---day--- (Fall) ---|
        conditions_TRF = [nfalls == nrises, first_sunfall_cut, T0_before_UTC_midnight]
        # |--- (Rise) ---day--- (Fall) ---(Transit)---| (Rise)
        conditions_RFT = [nfalls == nrises, first_sunfall_cut, T0_after_UTC_midnight]
        # |--- (Fall) --- (Transit) --- (Rise) ---day--- (Fall) ---| (Rise)
        conditions_FTRF = [nfalls > nrises, T0_before_UTC_midnight]
        
        if all(conditions_TRF):         # One fall, one rise, transit in first night.
            x_hold_fall[0] = -13*60
        elif all(conditions_RFT):       # One fall, one rise, transit in second night.
            x_hold_rise[0] = +21*60
        elif all(conditions_FTRF):      # Two falls, one rise, transit in first night.
            x_hold_rise.append(+21*60)
        else:
            pass
        
        # Recompute
        nfalls, nrises = len(x_hold_fall), len(x_hold_rise)

        if nrises > nfalls:
            night_length_mins = abs(x_hold_rise[1] - x_hold_fall[0])
        elif nfalls > nrises:
            night_length_mins = abs(x_hold_rise[0] - x_hold_fall[0])  
        elif nfalls == nrises:  # Normal 
            night_length_mins = abs(x_hold_rise[0] - x_hold_fall[0])


        # |--- (Fall) ------ (Rise) ---day--- (Fall) ---(Transit)---|
        conditions_FRFT = [nfalls > nrises, T0_after_UTC_midnight]
        # |--- (Rise) ---day--- (Fall) ---(Transit)--- (Rise) ---|
        conditions_RFTR = [nrises > nfalls, T0_after_UTC_midnight]
        # |---(Transit)--- (Rise) ---day--- (Fall) ------ (Rise) ---|
        conditions_TRFR = [nrises > nfalls, T0_before_UTC_midnight]

        if all(conditions_FRFT):              # Two falls, one rise, transit in second night.
            x_hold_fall[0] = x_hold_fall[1]
            x_hold_rise[0] = x_hold_fall[1] + night_length_mins
        elif all(conditions_RFTR):            # One fall, two rises, transit in second night.
            x_hold_rise[0] = x_hold_rise[1]    
        elif all(conditions_TRFR):            # One fall, two rises, transit in first night.
            sunfall_one = x_hold_rise[0] - night_length_mins
            x_hold_fall.insert(0, sunfall_one)
        else:
            pass

        night_midpoint = x_hold_rise[0] - night_length_mins/2      # Night midpoint
        sunfall_min, sunrise_min = x_hold_fall[0], x_hold_rise[0]
        # Sunfall and sunrise time relative to midnight
        dif_sunfall = round(sunfall_min/60, DIF_PRECISION)
        dif_sunrise = round(sunrise_min/60, DIF_PRECISION)


        ### HEIGHT AND NIGHT SPOT CHECKS
        sampling_interval_min = 1
        T14_min = T14 * 60 * 24

        # Transit timespace check
        transit_dur = np.arange(0, T14_min + 120, sampling_interval_min) * u.min    # Generates transit time-space
        transit_times = Time(row['Baseline_Start']) + transit_dur                   # Date and Time where transit occurs
        transit_times_num = transit_times.mjd
        transit_times_mask = (targ_times_num >= transit_times_num[0]) & (targ_times_num <= transit_times_num[-1])
        targ_altazs_T = targ_altazs[transit_times_mask]   # AltAz pairs during transit created
        
        # Night timespace check
        night_dur = np.arange(sunfall_min, sunrise_min, sampling_interval_min) * u.min    # Generates night time-space
        night_times = midnight + night_dur                                   # Date and Time where relevant night occurs
        night_times_num = night_times.mjd
        night_times_mask = (targ_times_num >= night_times_num[0]) & (targ_times_num <= night_times_num[-1])
        targ_altazs_N = targ_altazs[night_times_mask]     # AltAz pairs during relevant night created

        # Transit-Night overlap (TNO) timespace check
        TNO_times_mask = (transit_times_mask & night_times_mask)
        targ_altazs_TNO = targ_altazs[TNO_times_mask]

        # Is the target at observable height anywhere in the run?
        target_reaches_min_alt = np.any(targ_altazs_T.alt.value >= minAlt) or np.any(targ_altazs_N.alt.value >= minAlt)
        # Is the target at observable height throughout transit+baseline?      
        target_above_min_alt_entire_transit = np.all(targ_altazs_T.alt.value >= minAlt) 
        # Is the target at observable height at some point in the night?
        target_reaches_min_alt_at_night  = np.any(targ_altazs_N.alt.value >= minAlt)
        # Is the target observable at height, at night?
        target_reaches_min_alt_during_transit_at_night = np.any(targ_altazs_TNO.alt.value >= minAlt)


        # IMPLEMENT INTERNAL_RANKS   
        # Bank of generator expressions (F denotes "strict", for which target must always be sufficiently high)
        strict_altitude_conditions = [target_above_min_alt_entire_transit, target_reaches_min_alt_at_night]
        lax_altitude_conditions = [target_reaches_min_alt, target_reaches_min_alt_during_transit_at_night]

        transit_midnight_distance = abs(dif_T0*60 - night_midpoint)
        sufficient_baseline_margin = ( (night_length_mins - T14_min)/2 - transit_midnight_distance >= 60 )
        sufficient_transit_margin = ( (night_length_mins - T14_min)/2 - transit_midnight_distance >= 0 )
        transit_at_night_exists = ( transit_midnight_distance <= (T14_min + night_length_mins)/2 )

        cond_03_F = [sufficient_baseline_margin, *strict_altitude_conditions]
        cond_03_P = [sufficient_baseline_margin, *lax_altitude_conditions]

        cond_02_F = [sufficient_transit_margin, *strict_altitude_conditions]
        cond_02_P = [sufficient_transit_margin, *lax_altitude_conditions]

        cond_01_F = [transit_at_night_exists, *strict_altitude_conditions]
        cond_01_P = [transit_at_night_exists, *lax_altitude_conditions]

        # Assigns internal rank. If one condition is untrue, an all statement will return 'False'.
        # Events which fail an all statement are passed down to the next statement until 'True' is returned and they stop.
        if all(cond_03_F):
            Internal_Rank = '03_F'
        elif all(cond_03_P):
            Internal_Rank = '03_P'
        elif all(cond_02_F):
            Internal_Rank = '02_F'
        elif all(cond_02_P):
            Internal_Rank = '02_P'
        elif all(cond_01_F):
            Internal_Rank = '01_F'
        elif all(cond_01_P):
            Internal_Rank = '01_P'
        else:
            Internal_Rank = 'X'   


        # Choose airmass integral limits L1, L2 and return key times.
        # We have one integration to do, for transit + baselines for air mass marker.
        
        # If the target is always above 30 degrees, there are no crossing points to worry about!
        if 'F' in Internal_Rank:
            global L1, L2
            try:
                # FOR TRANSIT + BASELINE
                L1 = max(dif_BS, dif_sunfall)
                L2 = min(dif_BE, dif_sunrise)

                # Get-out for super-marginal events (eg HAT-P-66).
                if L2 < L1:
                    Internal_Rank = 'X'
                else:
                    pass

            except NameError:
                sys.exit(f"[MP_Process] {csv_name} ({planet_name}) has thrown an exception for a full event at {T1} - look at me!")
                
        # If there are crossing points (lax altitude condition), passed to bank of generator expressions.
        # any() will stop iterating as soon as it returns a "True" case.
        # For lower integral limit:
        elif 'P' in Internal_Rank:
            try:
                targ_alt_T_start = targ_altazs_T.alt.value[0]
                targ_alt_T_end = targ_altazs_T.alt.value[-1]

                minAlt_crossing_exists = (targ_alt_T_start < minAlt) or (targ_alt_T_end < minAlt)
                if minAlt_crossing_exists:
                    minAlt_cross_check = targ_altazs_T.alt.value - minAlt
                    minAlt_cross_check_idxs = np.arange(len(minAlt_cross_check) - 1)

                    for i in minAlt_cross_check_idxs:
                        if (minAlt_cross_check[i] == 0) or (minAlt_cross_check[i] * minAlt_cross_check[i+1] < 0):
                            minAlt_cross_hr_relative_to_transit = transit_dur.to('hr').value[i] - 1   # relative to start of transit (-1 because transit_dur includes baseline)
                            dif_minAlt_cross = dif_T1 + minAlt_cross_hr_relative_to_transit           # convert to be minAlt cross time relative to midnight

                            dif_minAlt_cross = round(dif_minAlt_cross, DIF_PRECISION)
                            minAlt_cross_hold.append(dif_minAlt_cross)
         

                # For targets that are barely visible, we get two crossing points!
                # These observing limits are easier to calculate.
                if len(minAlt_cross_hold) == 2:
                    L1 = max(minAlt_cross_hold[0], dif_sunfall)
                    L2 = min(minAlt_cross_hold[1], dif_sunrise)


                # One crossing point, calculate observing limits:
                elif len(minAlt_cross_hold) == 1:
                    # Stores every value to diagnose issue if integration limit determination goes wrong
                    events = [
                        ("dif_BS", dif_BS),
                        ("dif_sunfall", dif_sunfall),
                        ("dif_T1", dif_T1),
                        ("dif_minAlt_cross", dif_minAlt_cross),
                        ("dif_T4", dif_T4),
                        ("dif_sunrise", dif_sunrise),
                        ("dif_BE", dif_BE),
                    ]
                    events.sort(key=lambda x: x[1])
                    debugString = (
                        f"\n>> planet_name={planet_name}\n"
                        + "\n".join(f">> {name}={value}" for name, value in events)
                        + f"\n>> targ_alt_T_start={targ_alt_T_start}"
                        + f"\n>> targ_alt_T_end={targ_alt_T_end}"
                        + f"\n>> minAlt={minAlt}"
                    )
                    
                    choose_L1_sunfall_1 = (dif_BS <= dif_sunfall) and (targ_alt_T_end < minAlt)
                    choose_L1_sunfall_2 = (dif_T1 <= dif_sunfall <= dif_minAlt_cross <= dif_T4 <= dif_sunrise) and (targ_alt_T_end < minAlt)
                    choose_L1_sunfall_3 = (dif_T1 <= dif_sunfall <= dif_minAlt_cross <= dif_sunrise <= dif_T4) and (targ_alt_T_end < minAlt)
                    choose_L1_sunfall_4 = (dif_minAlt_cross <= dif_sunfall)
                    choose_L1_sunfall = choose_L1_sunfall_1 or choose_L1_sunfall_2 or choose_L1_sunfall_3 or choose_L1_sunfall_4
                    
                    choose_L1_minAlt_cross_1 = (dif_sunfall <= dif_minAlt_cross <= dif_BE) and (dif_BS < dif_sunfall)
                    choose_L1_minAlt_cross_2 = (dif_sunfall <= dif_minAlt_cross <= dif_BE) and (dif_BS < dif_sunrise) and (targ_alt_T_end > minAlt)
                    choose_L1_minAlt_cross_3 = (dif_BS <= dif_sunfall <= dif_minAlt_cross <= dif_BE) and (targ_alt_T_end > minAlt)
                    choose_L1_minAlt_cross = choose_L1_minAlt_cross_1 or choose_L1_minAlt_cross_2 or choose_L1_minAlt_cross_3                

                    choose_L1_baseline_start = (dif_BS <= dif_sunrise <= dif_minAlt_cross) or (dif_BS <= dif_minAlt_cross <= dif_sunrise)

                    # Now we get L1
                    if choose_L1_sunfall:
                        L1 = dif_sunfall
                    elif choose_L1_minAlt_cross:
                        L1 = dif_minAlt_cross
                    elif choose_L1_baseline_start:
                        L1 = dif_BS   
                    else:
                        print(debugString)
                        sys.exit(f'[P2_MultiprocessingProcess] Bad lower integration limit for {csv_name} ({planet_name}) transting at {T1} - flag me!')
                    
                    #############################################################################
    
                    choose_L2_baseline_end = (dif_BS <= dif_minAlt_cross <= dif_BE <= dif_sunrise) and (targ_alt_T_end > minAlt)
                    
                    choose_L2_minAlt_cross = (dif_minAlt_cross <= dif_sunrise <= dif_BE) or (dif_minAlt_cross < dif_BE)  

                    choose_L2_sunrise_1 = (dif_BS <= dif_minAlt_cross <= dif_sunrise <= dif_BE) and (targ_alt_T_end > minAlt)
                    choose_L2_sunrise_2 = (dif_sunfall <= dif_T1 <= dif_minAlt_cross <= dif_sunrise <= dif_T4) and (targ_alt_T_end > minAlt)
                    choose_L2_sunrise_3 = (dif_sunrise <= dif_minAlt_cross)
                    choose_L2_sunrise = choose_L2_sunrise_1 or choose_L2_sunrise_2 or choose_L2_sunrise_3

                    # Now we get L2
                    if choose_L2_baseline_end: 
                        L2 = dif_BE 
                    elif choose_L2_sunrise:
                        L2 = dif_sunrise
                    elif choose_L2_minAlt_cross:
                        L2 = dif_minAlt_cross
                    else:
                        print(debugString)
                        sys.exit(f'[P2_MultiprocessingProcess] Bad upper integration limit for {csv_name} ({planet_name}) transting at {T1} - flag me!')

                # Get-out or super-marginal events (eg HAT-P-66).
                if L2 <= L1:
                    Internal_Rank = 'X'
                else:
                    pass

            except NameError:
                sys.exit(f'\n[P2_MultiprocessingProcess] {csv_name} ({planet_name}) has thrown an exception for a partial event at {T1} - look at me!')
                
        else:
            pass  # This will throw an UnboundLocalError if an event has been incorrectly flagged.     
       

        # Choose all segment limits for event visibility metric.
        if Internal_Rank != 'X':
            def findVisibilityLimits(dif_from, dif_to):
                not_visible_during_transit = (dif_from <= L1 and dif_to <= L1) or (dif_from >= L2 and dif_to >= L2)
                if not_visible_during_transit == False:  # If at least part is INSIDE of visibility limit
                    L_from, L_to = max(dif_from, L1), min(dif_to, L2)
                elif not_visible_during_transit == True:
                    L_from, L_to = 0, 0
                
                return L_from, L_to
            
            # First baseline visibility limits
            L_BS1, L_BS2 = findVisibilityLimits(dif_BS, dif_T1) 
            # Ingress visibility limits
            L_In1, L_In2 = findVisibilityLimits(dif_T1, dif_T2) 
            # Full transit limits
            L_Tr1, L_Tr2 = findVisibilityLimits(dif_T2, dif_T3) 
            # Egress limits
            L_Eg1, L_Eg2 = findVisibilityLimits(dif_T3, dif_T4)
            # Second baseline limits
            L_BE1, L_BE2 = findVisibilityLimits(dif_T4, dif_BE)


            # What times should you put on your observing proposal? For this, UTC is convenient.
            # All previous calculations use BJD_TDB - these are best for publications!
            # This step was previously handled by P2_Timesplitter, but has been shrunk to go here.
            ProperStart = midnight.datetime + dt.timedelta(0,0,0,0,0,L1)
            ProperEnd = midnight.datetime + dt.timedelta(0,0,0,0,0,L2)
            ObsStart_JD, ObsEnd_JD, midnight_UTC = BJD2UTC([ProperStart, ProperEnd, midnight.datetime],
                                                           LockOn, Loc)
            
            # Recovery of base for air mass integration.
            # If it's a two-night event, a slightly extended integration base is needed. Try/except catches this.

            def Air_Mass_Marker_Calc(L1, L2, INT_val):
                Air_Mass_Marker = max(0,
                        3 * ((L2-L1) / INT_val - 2/3)
                        ) # Returns an air mass weighting between 0 and 1, where 30 deg = 0.   
                return Air_Mass_Marker
            
            try:
                epsrel = 1e-5
                INT_y = abs(np.sin(np.deg2rad(targ_altazs.alt.degree)))**(-0.6)

                with np.errstate(invalid='raise'):
                    f = interp1d(delta_midnight_targ, INT_y, kind='cubic')
                    INT_val, INT_err = integrate.quad(f, L1, L2, epsrel=epsrel)
                    Air_Mass_Marker = Air_Mass_Marker_Calc(L1, L2, INT_val)
                     
            except ValueError:
                extra_delta_midnight = np.linspace(-13, 19, targ_samples) * u.hour
                f = interp1d(extra_delta_midnight, INT_y, kind='cubic')
                try:
                    with np.errstate(invalid='raise'):
                        INT_val, INT_err = integrate.quad(f, L1, L2, epsrel=epsrel)
                        Air_Mass_Marker = Air_Mass_Marker_Calc(L1, L2, INT_val)

                except ValueError:  # Second exception for very long transiters (eg. Kepler-432)
                    INT_val, INT_err = integrate.quad(f, L1-24, L2-24, epsrel=epsrel)
                    Air_Mass_Marker = Air_Mass_Marker_Calc(L1-24, L2-24, INT_val)


            # Weight by percentage of event + baseline captured.
            EVENT_WEIGHT_PRECISION = DIF_PRECISION - 1
            try:
                with np.errstate(invalid='raise'):
                    WBase = ((L_BS2-L_BS1) + (L_BE2-L_BE1)) / 2            # Baseline
                    WTrans = (L_Tr2-L_Tr1) / (dif_T3-dif_T2)               # Full depth
                    WInOut = ((L_Eg2-L_Eg1) + (L_In2-L_In1)) / (2 * (dif_T2-dif_T1))    # Ingress-egress

                    event_weight = np.around(WBase * WTrans * WInOut, EVENT_WEIGHT_PRECISION)
                   
            except FloatingPointError: # Case grazing transits -- No need to consider WTrans!
                WBase = ((L_BS2-L_BS1) + (L_BE2-L_BE1)) / 2                 # Baseline
                WTrans = np.nan                                             # Full depth (never achieved for grazing events!)
                WInOut = ((L_Eg2-L_Eg1) + (L_In2-L_In1)) / (2 * (dif_T2-dif_T1))        # Ingress-egress

                event_weight = np.around(WBase * WInOut, EVENT_WEIGHT_PRECISION)

            except ZeroDivisionError:   
                sys.exit(f'\n[P2_MultiprocessingProcess] ZeroDivisionError for {csv_name} ({planet_name}) transiting at {T1} - flag me!')

            WBase, WTrans, WInOut = round(WBase, EVENT_WEIGHT_PRECISION), round(WTrans, EVENT_WEIGHT_PRECISION), round(WInOut, EVENT_WEIGHT_PRECISION)

            # Acquire moon noise metric
            if toggle_moonlight_noise == True:
                LUT_times = pd.to_datetime(moon_LUT['time_UTC'])
                ObsStart_JD_dt = ObsStart_JD.to_datetime()
                idx_ObsStart = (LUT_times - ObsStart_JD_dt).abs().idxmin() 

                moon_mag = moon_LUT[f'{filter_name}mag'].iloc[idx_ObsStart]

                moon_noise_metric = MoonNoiseMetric(moon_altazs, targ_altazs,
                                                    moon_mag, scattering_aod, absorption_aod, asymmetry_factor,
                                                    dif_BS, delta_midnight_targ, dif_BE, filter_name)    
            elif toggle_moonlight_noise == False:
                moon_noise_metric = 1  # SNR doesnt change if no noise source added

        else:
            ObsStart_JD = np.nan
            ObsEnd_JD = np.nan
            Air_Mass_Marker = 0  # Unobservable 
            WBase = 0
            WTrans = 0
            WInOut = 0
            event_weight = 0
            moon_noise_metric = 0 # No need to calculate as it is not used


        # Create plot only if event_weight >= event_weight_graph_threshold to reduce unneccessary images and CPU load
        if (toggle_graph_outputs == True) and (event_weight >= event_weight_graph_threshold):

            # Shift time axis from BJD to local time for that telescope (for graphing purposes)
            # First, find timedelta where midnight_UTC - midnight_BJD in hours
            delta_h_UTC_BJD = (midnight_UTC - midnight).to('hour').value
            # Next, find utc offset where midnight_local - midnight_UTC in hours
            utc_offset = midnight_UTC.to_datetime(timezone=ZoneInfo(timezone_str)).utcoffset() / dt.timedelta(hours=1)
            # Thus, we get offset = local - BJD
            offset = delta_h_UTC_BJD + utc_offset

            # Shift delta_midnight and Moon_step to local time since midnight: local = BJD + offset
            delta_midnight_targ_local = delta_midnight_targ + (offset * u.hour)
            delta_midnight_sep_local = delta_midnight_sep + (offset * u.hour)
            xlim = (delta_midnight_targ_local.min().value, delta_midnight_targ_local.max().value)

            # Altitude mask to only plot whatever is in plot
            moonaltmask = (moon_altazs.alt.value >= -1)
            targaltmask = (targ_altazs.alt.value >= -1)

            # Plot target path and transit viability for all events.
            # Scatter c argument denotes colour, lw=linewidth, default s=20, zorder stops fill from mucking up the colour of your other lines.
            # Use fig.close(fig) to stop all the figures from eating your RAM.
            fig, ax = plt.subplots(figsize=(18, 10))

            ax.plot(delta_midnight_targ_local[moonaltmask][::PLOT_DECIMATION],
                    moon_altazs.alt[moonaltmask][::PLOT_DECIMATION],
                    color='darkorange', zorder=2, linestyle='-.', linewidth=2.5,
                    label='Lunar Path')                        # Moon path
            
            sc = ax.scatter(delta_midnight_targ_local[targaltmask][::PLOT_DECIMATION],
                            targ_altazs.alt[targaltmask][::PLOT_DECIMATION],
                            c=targ_altazs.az.value[targaltmask][::PLOT_DECIMATION],
                            cmap=cmap, norm=norm, zorder=3, lw=0, s=15,
                            antialiased=False)   # Target path
            ax.plot([], [], color='blue', label='Target Path') # Dummy plot to put on legend

            ax.scatter(delta_midnight_sep_local, targsep_altazs.alt,
                       facecolors='white', edgecolors='blue', zorder=4, s=100, linewidths=2, marker='o',
                       label='Lunar Separation (deg)')         # Separation scatterplot
            
            # Attaches separations to markers along target path.
            for w, txt in enumerate(seps_int):
                ax.annotate(seps_int[w], (delta_midnight_sep_local.value[w], targsep_altazs.alt.value[w]),
                            color='violet', xytext=(-7.5, -27), textcoords='offset pixels')
                
            ax.fill_between(delta_midnight_targ_local.to('hr').value, 0, 90, (sun_altazs.alt <= nightAlt*u.deg),
                            facecolor='black', alpha=1, zorder=0, label='Night')   # Bounded by astronomical twilight
            
            # Red fill denotes transit in progress, superimposed on other zones. alpha denotes transparency.
            ax.axvspan(dif_BS+offset, dif_BE+offset, 0,90, alpha=0.30, facecolor='r', zorder=1,
                       label='Observation')
            ax.axvspan(dif_T1+offset, dif_T4+offset, 0,90, alpha=0.30, facecolor='r', edgecolor='r', hatch='X', zorder=1,
                       label='Transit Active')
            
            ax.axhspan(0, minAlt, xlim[0], xlim[1], facecolor='grey', alpha=0.2, zorder=5,
                       label='Unobservable altitude')
            ax.axhspan(maxAlt ,90, xlim[0], xlim[1], facecolor='grey', alpha=0.2, zorder=5)

            fig.colorbar(sc, ax=ax).set_label('Azimuth [deg]')
            ax.legend(loc='upper right')

            textpos = xlim[1] + 5
            ax.text(textpos, 85, f'Moon Illumination:\n{phase:.0f}%')
            ax.text(textpos, 75, f'Closest Lunar\nSeparation: {min(seps_int):.0f} deg')

            # Custom xtick formatting
            ax.set_xlim(xlim)
            tick_start = 2 * np.floor(np.rint(xlim[0]) / 2) # Force ticks to be always even
            tick_end = 2 * np.ceil(np.rint(xlim[1]) / 2)
            ax.set_xticks(np.arange(tick_start+2, tick_end, 2))
            
            def format_tick(val, pos):
                hour = int(val % 24)  # Hour of date
                if hour == 0:
                    # Add val hours to midnight (UTC -> local) to get actual date on tick
                    date_midnight = midnight_UTC + dt.timedelta(hours=(utc_offset+val))
                    return date_midnight.strftime('00\n%d %b')
                else:
                    return f'{hour:02d}' # Return only hour
            ax.xaxis.set_major_formatter(FuncFormatter(format_tick))

            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_ylim(0, 90)
            ax.set_xlabel(f'Local time in hours (UTC+{utc_offset:.1f})')
            ax.set_ylabel('Altitude [deg]')

            ax.set_title(f'{planet_name} Transit Observability for {instrument} from\n' +\
                         f'{(ObsStart_JD + dt.timedelta(hours=utc_offset)).strftime("%b %d %Y %H:%M")} to ' +\
                         f'{(ObsEnd_JD   + dt.timedelta(hours=utc_offset)).strftime("%b %d %Y %H:%M")} (UTC+{utc_offset:.1f})')

            fig.savefig(graph_folder_path / f'{event_weight:.3f}_{Internal_Rank}_{planet_name}_{instrument}_{T1.strftime("%b-%d-%Y-%H-%M-%S")}.jpg',
                        dpi=70, transparent=False, facecolor='white', edgecolor='black', format='jpg',
                        pil_kwargs={'quality': 75, 'optimize': False, 'progressive': False})

            plt.close(fig)  # This is still needed to release system resources after plot is saved. Must be explicit!
        
        return phase, min(seps), Internal_Rank, ObsStart_JD, ObsEnd_JD, Air_Mass_Marker, WBase, WTrans, WInOut, event_weight, moon_noise_metric
    

    def FinalMetric(row):
        final_metric = df[metric_mode][0] * row['Air_Mass_Metric'] * row['Event_Weight'] * row['Moon_Noise_Metric']
        return final_metric


    ##### DEFINITIONS END HERE #####  
    # csv file containing one planet.
    df_cols = ['Planet', 'R*', 'Rp', 'RA:HH', 'RA:MM', 'RA:SS', 'Dec:DD', 'Dec:MM', 'Dec:SS', f'{filter_name}mag',
        'T14', 'Depth', 'T0 (HJD or BJD)', 'T0_err', 'P (day)', 'P_err',
        'Impact Parameter', 'a_Calc', 'Previous Study Flag',
        'TSM', 'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']
    df = pd.read_csv(csv_initiated, usecols=df_cols, skipinitialspace=True)

    # Altitude spot-check - will the target ever be high enough to be observed?
    df[['RA:HMS', 'Dec:Deg']] = df.apply(lambda row: AltChecker(row), axis=1, result_type='expand')
    planet_name = df['Planet'][0]
    
    if minAlt >= 90 - abs(scope_df['Lat'].iloc[scope_idx] - df['Dec:Deg'][0]):  # Case physically unobservable target
        pass
    else:
        # Fill in transit_start_times and transit_start_times_err lists
        T1_Calc(df.iloc[0])   # There is only one row anyways

        if len(transit_start_times) == 0:  # Case no transits in observation window:           
            pass

        else:
            # Create T1 and T1_err columns according to Figure 3.6
            # Reindex original di to create all these NaN rows
            df = df.reindex(np.arange(len(transit_start_times)))
            # Concatenate T1 and T1_err
            df = pd.concat([df,
                            pd.DataFrame({'T1': transit_start_times, 'T1_err': transit_start_times_err})],
                            axis=1)

            # Optimization routine: Run LockOn.transform_to(frame) ONCE per CSV initiated (instead of once per row)
            # While this makes code harder to read, the performance gain is well worth it
            # A big bottleneck is actually the overhead for each .transform_to call

            # Target coord in RA/Dec
            LockOn = ac.SkyCoord(df['RA:HMS'][0], df['Dec:Deg'][0], unit=(u.hourangle, u.deg), frame='icrs')

            # List of target times to transform
            targ_samples = 1801
            min_interval = 10
            targ_times_10m_list = [
                (closest_midnight(T1) + np.linspace(-12, 18, targ_samples) * u.hour)[::min_interval]
                for T1 in df['T1']
            ]
            concat_targ_times_10m = Time(   # Transform to vectorizable form
                np.concatenate([t.jd for t in targ_times_10m_list]),
                format='jd', scale='utc'
            )  
            frame = ac.AltAz(obstime=concat_targ_times_10m, location=Loc)
            targ_altazs_10m_all = LockOn.transform_to(frame)

            n5m = targ_samples//min_interval + 1
            targ_altazs_10m_by_row = [
                targ_altazs_10m_all[i*n5m:(i+1)*n5m]
                for i in range(len(df))
            ]  # Call targ_altazs_10m at that row by targ_altazs_5m_by_row[row.name]

            
            # Writes additional key times and results
            df[['T2', 'T0', 'T3', 'T4', 'Baseline_Start', 'Baseline_End']] = df.apply(lambda row: Times_Calc(row),
                                                                                      axis=1, result_type='expand')
            df[['Lunar_Illumination', 'Closest_Lunar_Approach', 'Internal Rank', 'Observation_Start_(UTC)', 'Observation_End_(UTC)',
                'Air_Mass_Metric', 'Baseline_Weight', 'Transit_Curve_Weight', 'Ingress-Egress_Weight', 'Event_Weight',
                'Moon_Noise_Metric']] = df.apply(lambda row: EventMetric(row),
                                                 axis=1, result_type='expand')

            # Final metrics, then export
            df[f'Final_{metric_mode}'] = df.apply(lambda row: FinalMetric(row), axis=1)
            df = df.drop(columns=f'{filter_name}mag')

            df.to_csv(output_folder / 'phase_2' / 'individual_planets'/ csv_name, index=False)
    
    # CLEARANCE PROTOCOLS FOR MULTIPLE-TARGET RUNNING   
    transit_start_times.clear()
    transit_start_times_err.clear()
