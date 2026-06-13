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
from astropy.table import QTable
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from dateutil.relativedelta import relativedelta

# Master function for Phase Two - calculation of event times/metrics.
def P2Predictor(csvcorepath, csvinpath, csvoutpath,
                ObsStart, ObsEnd, ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
                Moon_Noise, Scattering_AOD, Absorption_AOD, Asymmetry_Factor, Amplification_Factor, EventPercent_minThreshold,
                Loc, timezone_str):
    # Filename of this process
    Name = os.path.basename(csvinpath)

    # Constants needed
    PT8 = 30        # Min obs. altitude for air mass = 2.
    LIV = 88        # Max obs. altitude for the Liverpool instrument
    Cross = -18     # Marker for night start/end altitude

    # Blank arrays for use with time cut/observation ranking
    TargetTimes = []
    ErrTimes = []
    x_hold_fall = []
    x_hold_rise = []
    Cross_hold = []

    # Lookup table for moon noise calculations
    if Moon_Noise == 'Y_Moon':
        moonlightParams_name = rf'{csvcorepath}/Lookup_tables/lunar_TOA_Mags_at_{Inst}_{ObsStart.strftime("%b %d %Y")}_to_{ObsEnd.strftime("%b %d %Y")}.csv'
        moonlightParams_LUT = pd.read_csv(moonlightParams_name)

    # Lookup table for sun/moon AltAz skycoords
    AltAz_name = rf'{csvcorepath}/Lookup_tables/sun_moon_altazs_at_{Inst}_{ObsStart.strftime("%b %d %Y")}_to_{ObsEnd.strftime("%b %d %Y")}.ecsv'
    AltAz_LUT = QTable.read(AltAz_name, format='ascii.ecsv')

    sunaltazs_LUT = ac.SkyCoord(alt=AltAz_LUT['sun_alt'], az=AltAz_LUT['sun_az'], unit='deg',
                                obstime=AltAz_LUT['obstime'], location=Loc, frame='altaz')
    moonaltazs_LUT = ac.SkyCoord(alt=AltAz_LUT['moon_alt'], az=AltAz_LUT['moon_az'], unit='deg',
                                 obstime=AltAz_LUT['obstime'], location=Loc, frame='altaz')

    def AltChecker(row): 
        # Specify sky co-ordinates of target.
        LockOn = ac.SkyCoord(row['RA:HH'] + (row['RA:MM']/60) + (row['RA:SS']/3600),
                             row['Dec:DD'] + np.sign(row['Dec:DD']) * (row['Dec:MM']/60 + row['Dec:SS']/3600),
                             unit=(u.hourangle, u.deg), frame='icrs')
        
        return LockOn.ra.hour, LockOn.dec.degree

    def T1Calc(row):
        # Internal conversion to recover midpoint of first observation in BJD_TDB.    
        period_d = row['P (day)'] 
        eph_obj = Time(row['T0 (HJD or BJD)'], format='jd').datetime    # Scale is not JD, but the format (245XXXX.XXX...) is!
        ESecErr = (row['T0_err'] * u.d).to(u.s)       # Conversion of ephemeris error.
        PSecErr = (row['P_err'] * u.d).to(u.s)        # Conversion of error in planet's determined period.

        # Recovers transit start time from published midpoint, returning a timedelta object.
        StartTime = eph_obj - dt.timedelta(0.5 * row['T14'])
        
        # Builds range of n (epoch number) relevant for target/window combination.
        # nstart finds final transit before our window begins, so you gotta add +1!
        nStart = (ObsStart - StartTime).total_seconds() // (period_d * u.d).to(u.s).value + 1 
        nEnd = (ObsEnd - StartTime).total_seconds() // (period_d * u.d).to(u.s).value
        Range = np.arange(nStart, nEnd + 1)

        # Increase the number of orbits n to find future transit start times in our window.
        for n in Range:
            # Calculation of individual transit start times
            NextTransitTime = StartTime + (n * dt.timedelta(period_d))
            TargetTimes.append(NextTransitTime)
            
            # Calculation of total error at epoch n.
            # Uncertainty should scale overall as no. of orbits n, so this needs to be inside the squared term!
            ErrTot = np.sqrt(ESecErr**2 + (n * PSecErr)**2).value              
            ErrTimes.append(ErrTot)
    

    # Calculation of key times.
    def TimesCalc(row):
        # Locked to index 0 except T1 because the other values are NaN at index 1 onwards
        period_d = di['P (day)'][0]
        Rs = di['R*'][0] * u.solRad
        Rp = (di['Rp'][0] * u.R_jup).to(u.solRad)
        T1 = row['T1']
        T14 = di['T14'][0]
        b = di['Impact Parameter'][0]
        a_Calc = di['a_Calc'][0]
        
        # Calculates midpoint (T0), end of transit (T4) and length of transit event.
        TransitMid = T1 + dt.timedelta(0.5 * T14)                
        TransitEnd = T1 + dt.timedelta(T14)   

        try:
            with np.errstate(invalid='raise'):        
                T23 = (period_d/np.pi) * np.arcsin(
                      np.sqrt((Rs - Rp)**2 - (b * Rs)**2) / (a_Calc * u.AU).to(u.solRad)
                )
                # Recovers length of ingress/egress. abs is a necessary guess if T23 > T14 (unphysical soln.)     
                T12 = abs(T14 - T23.value) / 2
                T2 = T1 + dt.timedelta(T12)             # Time at end of ingress
                T3 = TransitEnd - dt.timedelta(T12)     # Time at start of egress
                
                # TESS targets have dodgy parameters - this calculation will catch them.
                # Forces a FloatingPointError so the transit is treated as grazing, passed to except block.
                np.sqrt((T3 - T2).total_seconds())
        
        # Grazing transits have no T23, so a different approach is needed.
        except FloatingPointError:
            T23 = 0      
            T2 = TransitMid # Time at end of ingress
            T3 = TransitMid # Time at start of egress

        except ValueError:
            sys.exit(f"Time error for {Name} ({PlanetName}) transiting at {T1} - flag me!")

        # A solid baseline is important for a good transit curve!
        BaseOne = T1 - dt.timedelta(hours=1)
        BaseTwo = TransitEnd + dt.timedelta(hours=1) 

        return T2, TransitMid, T3, TransitEnd, BaseOne, BaseTwo

    # Conversion from BJD_TDB to UTC
    def BJD2UTC(time, LockOn):
        time_BJD = Time(time, scale='tdb', location=Loc)   # BJD_TDB (implicit)
        ltt_bary = time_BJD.light_travel_time(LockOn)      # JD to BJD conversion factor (light travel time)
        JD_TDB = time_BJD - ltt_bary                       # JD_TDB
        JD_UTC = JD_TDB.utc                                # JD_UTC

        return JD_UTC

    # Adds moon noise metric which is essentially SNR of yes moon / SNR of no moon, at start of observation.
    # This works because the metrics are directly proportional to SNR.
    def MoonNoiseMetric(lunar_altaz, target_altaz,
                        moon_mag, AOD_scatter, AOD_absorption, asymmetry_factor,
                        dif_BS, delta_midnight, dif_BE, filter):
        from .P1_RankMaker import findSkyB
        from .P2_MultiprocessingProcess import getFilterParams

        # Acquire effective wavelength
        effective_wavelength, _, _ = getFilterParams(csvcorepath, filter)

        # Calculate moon background in mag/arcsec^2, then to mag at aperture much like the same way as skyB.
        # Based on Winkler (2022) model.
        # TOA moon flux is computed using LIME model with Allen extrapolation.
        lunar_zenith = lunar_altaz.zen.rad
        target_zenith = target_altaz.zen.rad
        coord_targ = ac.SkyCoord(target_altaz.az, target_altaz.alt, unit="deg")
        coord_moon = ac.SkyCoord(lunar_altaz.az, lunar_altaz.alt, unit="deg")
        theta = coord_targ.separation(coord_moon).to(u.rad).value

        # Remove out-of-observation values AND unphysical values when zeniths go beyond 90 deg
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
        tau_R = np.exp(-Loc.height.value/scaleHeight) * (1.229e10) * effective_wavelength**(-4.05)  # Eqn 13
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
        moon_intensity_mag = moon_mag - 2.5 * np.log10(moon_intensity_decrease)
        
        # Calculate noise metric
        moon_flux = 10**(-0.4 * (findSkyB(moon_intensity_mag, ds, S, Defocus, Run_Mode) - Amplification_Factor))
        if Add_Noise == 'Y_Noise':
            sky_flux = 10**(-0.4 * findSkyB(ds[f'msky_{Filter}'].iloc[S], ds, S, Defocus, Run_Mode))
        else:
            sky_flux = 0
        target_flux = 10**(-0.4 * di[f'{Filter}mag'][0])
        
        MoonNoiseMetric = (1 + moon_flux/(sky_flux + target_flux)) ** (-0.5)
        return np.min(MoonNoiseMetric)  # Take min value (signifies highest SNR reduction)


    ##### WE GET INTO CALCULATION AND PLOTTING HERE #####
    def EventMetric(row):
        from P2MP_Wrapper_Datacenter import interpolate_altaz

        # Specify co-ordinates of target in sky coords, as previously converted.
        LockOn = ac.SkyCoord(di['RA:HMS'][0], di['Dec:Deg'][0], unit=(u.hourangle, u.deg), frame='icrs')
        
        # Chooses appropriate midnight for each useful time (for graphing purposes).
        # Use relativedelta to avoid end of months/years breaking the code.
        T1 = row['T1']
        T14 = di['T14'][0]  # It's NAN in second row onwards
        t_samples = 1801

        # Set closest midnight (BJD) associated with start of transit T1
        if T1.hour >= 12:
            MidnightSet = T1.replace(hour=0, minute=0, second=0, microsecond=0) + relativedelta(days=1)
        else:
            MidnightSet = T1.replace(hour=0, minute=0, second=0, microsecond=0)
        midnight = Time(MidnightSet)                            # Astropy object conversion for compatibility.

        delta_midnight = np.linspace(-12, 18, t_samples) * u.hour  # Time sampled 30h from midnight for target (1 min)
        Moon_step = np.linspace(-12, 18, 21) * u.hour              # Time sampled 30h from midnight for moon sep (90 min)

        times = midnight + delta_midnight          # Times for Moon/target calc
        Moon_times = midnight + Moon_step          # Times for Moon/target separation calc

        # Where is the sun and moon from midday to evening? 
        # LUTs are used here because it makes the pipeline way faster (2-3x faster, seriously)
        obstime_num = sunaltazs_LUT.obstime.mjd  # Convert to array for indexing mask
        times_num = times.mjd
        moon_times_num = Moon_times.mjd

        timemask = (obstime_num >= times_num[0]) & (obstime_num <= times_num[-1])
        septimemask = np.isin(obstime_num, moon_times_num)

        sunaltazs = sunaltazs_LUT[timemask]           # AltAz pairs created for sun
        moonaltazs = moonaltazs_LUT[timemask]         # AltAz pairs created for moon (1 min)
        moonsepaltazs = moonaltazs_LUT[septimemask]   # AltAz pairs created for moon (90 min)

        # Where is the target from midday to evening?
        # For very, very significant speed optimization, use LockOn.transform_to() exactly once (instead of 5)
        # Calculate for 5 minute precision then interpolate much like Sun and Moon
        times_5m = times[::5]
        frame = ac.AltAz(obstime=times_5m, location=Loc)
        Targetaltazs = LockOn.transform_to(frame)        # AltAz pairs created for target
        Target_alt, Target_az = interpolate_altaz(times_5m, Targetaltazs, times)
        Targetaltazs = ac.SkyCoord(alt=Target_alt, az=Target_az, unit="deg",
                                   obstime=times, location=Loc, frame='altaz')
        
        moon_times_in_times = np.isin(times_num, moon_times_num)
        Targetsepaltazs = Targetaltazs[moon_times_in_times]

        # Calculates separation between target and Moon at 21 intervals.
        coord_targ = ac.SkyCoord(Targetsepaltazs.az, Targetsepaltazs.alt, unit="deg")
        coord_moon = ac.SkyCoord(moonsepaltazs.az, moonsepaltazs.alt, unit="deg")

        Sep = coord_targ.separation(coord_moon).value
        FormSep = Sep.astype(int)

        ### LUNAR ILLUMINATION CALC (Returns moon phase in %)
        m = ephem.Moon()
        m.compute(T1)
        Phase = m.moon_phase * 100

        ### CLEARANCE PROTOCOLS FOR MULTIPLE-EVENT RUNNING
        x_hold_fall.clear()
        x_hold_rise.clear()
        Cross_hold.clear()
        

        # When will the transit, baseline and night begin and end -- relative to midnight? 
        # Finds differences either side of 00:00:00.
        # NOTE: Take only the value when plotting! You can't plot a unit.
        def hours_relative_to_midnight(time):
            return (Time(time) - midnight).to('hr').value

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
        CrossCheck = sunaltazs.alt.value - Cross

        for i in np.arange(len(CrossCheck) - 1):
            if CrossCheck[i] * CrossCheck[i+1] <= 0:    # If the moment sun crosses over y=-18 at i
                x_ = delta_midnight.to('min').value[i]  # Time in delta_midnight where the crossing happens 

                if sunaltazs.alt.value[i+1] < sunaltazs.alt.value[i]:      # If sunfall
                    x_hold_fall.append(x_)
                elif sunaltazs.alt.value[i+1] >= sunaltazs.alt.value[i]:   # If sunrise
                    x_hold_rise.append(x_)

        # For "two-night events", with a bit at the start and end!
        # It is the second condition in particular, which catches these.
        NGen1 = [len(x_hold_fall) == len(x_hold_rise), x_hold_fall[0] > x_hold_rise[0], dif_T0 < 0]
        NGen2 = [len(x_hold_fall) == len(x_hold_rise), x_hold_fall[0] > x_hold_rise[0], dif_T0 > 0]
        NGen3 = [len(x_hold_fall) > len(x_hold_rise), dif_T0 < 0]  # Condition for two falls, one rise, transit in first night.
        
        if all(NGen1):
            x_hold_fall[0] = -780
        elif all(NGen2):
            x_hold_rise[0] = 1160
        elif all(NGen3):
            x_hold_rise.append(1160)
        else:
            pass

        if len(x_hold_rise) > len(x_hold_fall):
            TNight = abs(x_hold_rise[1] - x_hold_fall[0])
        elif len(x_hold_rise) < len(x_hold_fall):
            TNight = abs(x_hold_rise[0] - x_hold_fall[1])  
        else:
            TNight = abs(x_hold_rise[0] - x_hold_fall[0])  # Night length
        
        NGen4 = [len(x_hold_fall) > len(x_hold_rise), dif_T0 > 0]
        NGen6 = [len(x_hold_fall) < len(x_hold_rise), dif_T0 > 0]
        NGen5 =  len(x_hold_fall) < len(x_hold_rise)

        if all(NGen4):              # Two falls, one rise, transit during second night.
            x_hold_fall[0] = x_hold_fall[1]
            x_hold_rise[0] = x_hold_fall[1] + TNight
        elif all(NGen6):            # Exception for one fall, two rises, transit in second night.
            x_hold_rise[0] = x_hold_rise[1]    
        elif NGen5 == True:         # Exception for one fall, two rises, transit in first night.
            FallOne = x_hold_rise[0] - TNight
            x_hold_fall.insert(0, FallOne)
        else:
            pass

        TNight0 = x_hold_rise[0] - TNight/2 # Night midpoint


        ### HEIGHT AND NIGHT SPOT CHECKS
        sampling_interval_min = 1

        # Transit timespace check
        TDur = np.arange(0, (T14*1440) + 120, sampling_interval_min) * u.min    # Generates transit time-space
        Ttimes = Time(row['Baseline_Start']) + TDur              # Date and Time where transit occurs
        Ttimes_num = Ttimes.mjd
        Ttimes_mask = (times_num >= Ttimes_num[0]) & (times_num <= Ttimes_num[-1])
        SC = Targetaltazs[Ttimes_mask]                           # AltAz pairs during transit created
        
        # Night timespace check
        Nt = np.arange(x_hold_fall[0], x_hold_rise[0], sampling_interval_min) * u.min     # Generates night time-space
        Ntimes = midnight + Nt                                   # Date and Time where relevant night occurs
        Ntimes_num = Ntimes.mjd
        Ntimes_mask = (times_num >= Ntimes_num[0]) & (times_num <= Ntimes_num[-1])
        NC = Targetaltazs[Ntimes_mask]                           # AltAz pairs during relevant night created

        HSC = np.any(SC.alt.value >= PT8) or np.any(NC.alt.value >= PT8)   # Is the target at observable height anywhere in the run?        
        HFC = np.all(SC.alt.value >= PT8)  # Is the target at observable height throughout the run?
        NSC = np.any(NC.alt.value >= PT8)  # Is the target at observable height at some point in the night?

        # Transit-night overlap timespace check
        ingress_sunfall = [dif_T1, x_hold_fall[0]/60]
        egress_sunrise = [dif_T4, x_hold_rise[0]/60]   
        transit_night_overlap = midnight + np.arange(max(ingress_sunfall), min(egress_sunrise), sampling_interval_min/60) * u.hr
        # Sometimes, the transit happens just at day and transit_night_overlap doesn't exist.
        if transit_night_overlap.size > 0:
            TNO_num = transit_night_overlap.mjd
            TNO_mask = (times_num >= TNO_num[0]) & (times_num <= TNO_num[-1])
            KC = Targetaltazs[TNO_mask]        # AltAz pairs during relevant transit-night overlap
            KSC = np.any(KC.alt.value >= PT8)  # Is the target observable at height, at night? (K2-99 exception.)
        else:
            KSC = False

        # IMPLEMENT CONDITIONS 1-3        
        # Bank of generator expressions (S denotes "strict", for which target must always be sufficiently high)
        C1Gen = [abs(dif_T0*60 - TNight0) <= (T14*1440 + TNight)/2, HSC, KSC]
        C1GenS = [abs(dif_T0*60 - TNight0) <= (T14*1440 + TNight)/2, HFC, NSC]
        C2Gen = [T14*1440 + 2*abs(dif_T0*60 - TNight0) <= TNight, HSC, KSC]
        C2GenS = [T14*1440 + 2*abs(dif_T0*60 - TNight0) <= TNight, HFC, NSC]
        C3Gen = [(TNight - T14*1440)/2 - abs(dif_T0*60 - TNight0) >= 60, HSC, KSC]
        C3GenS = [(TNight - T14*1440)/2 - abs(dif_T0*60 - TNight0) >= 60, HFC, NSC]

        # Assigns internal rank. If one condition is untrue, an all statement will return 'False'.
        # Events which fail an all statement are passed down to the next statement until 'True' is returned and they stop.
        if all(C3GenS):
            Internal_Rank = '03_F'
        elif all(C3Gen):
            Internal_Rank = '03_P'
        elif all(C2GenS):
            Internal_Rank = '02_F'
        elif all(C2Gen):
            Internal_Rank = '02_P'
        elif all(C1GenS):
            Internal_Rank = '01_F'
        elif all(C1Gen):
            Internal_Rank = '01_P'
        else:
            Internal_Rank = 'X'   

        # Choose integral limits L1, L2 and return key times.
        # We have one integration to do, for transit + baselines for air mass marker.
        
        # If the target is always above 30 degrees, there are no crossing points to worry about!
        if 'F' in Internal_Rank:
            global L1, L2
            try:
                # FOR TRANSIT + BASELINE
                L1 = max(dif_BS, x_hold_fall[0]/60)
                L2 = min(dif_BE, x_hold_rise[0]/60)

                # Get-out for super-marginal events (eg HAT-P-66).
                if L2 < L1:
                    Internal_Rank = 'X'
                else:
                    pass

            except NameError:
                sys.exit(f"[MP_Process] {Name} ({PlanetName}) has thrown an exception for a full event at {T1} - look at me!")
                
        # If there are crossing points, passed to bank of generator expressions.
        # any() will stop iterating as soon as it returns a "True" case.
        # For lower integral limit:
        elif 'P' in Internal_Rank:
            try:     
                PGen = [SC.alt.value[0] < PT8, SC.alt.value[-1] < PT8]
                if any(PGen):
                    PT8CC = SC.alt.value - PT8
                    for i in np.arange(len(PT8CC) - 1):
                        if PT8CC[i] == 0 or PT8CC[i] * PT8CC[i+1] < 0:
                            LT1 = TDur.to('hr').value[i] - 1    # relative to start of transit
                            LTC1 = dif_T1 + LT1                 # convert to be relative to midnight.

                            Cross_hold.append(LTC1)
         
                # For targets that are barely visible, we get two crossing points!
                # These observing limits are easier to calculate.
                if len(Cross_hold) == 2:
                    L1 = max(Cross_hold[0], x_hold_fall[0] / 60)
                    L2 = min(Cross_hold[1], x_hold_rise[0] / 60)

                # One crossing point, calculate observing limits:
                elif len(Cross_hold) == 1:
                    PGen2 = [dif_BS < x_hold_fall[0]/60, SC.alt.value[-1] < PT8]

                    # For transits that are too long on fall side for one night.
                    PGen12 = [dif_T1 < x_hold_fall[0]/60 < LTC1 < dif_T4 < x_hold_rise[0]/60,
                              SC.alt.value[-1] < PT8]

                    # This condition is for transits that are too long for 1 night!
                    PGen15 = [dif_T1 < x_hold_fall[0]/60 < LTC1 < x_hold_rise[0]/60 < dif_T4,
                              SC.alt.value[-1] < PT8]
 
                    PGen3 = [all(PGen2), all(PGen12), all(PGen15),
                             x_hold_fall[0]/60 > LTC1]
                    
                    PGen4 = [x_hold_fall[0]/60 < LTC1 < dif_BE,
                             x_hold_fall[0]/60 > dif_BS]
                    
                    PGen5 = [dif_BS < x_hold_rise[0]/60 < LTC1,
                             dif_BS < LTC1 < x_hold_rise[0]/60]
                    
                    PGen6 = [x_hold_fall[0]/60 < LTC1 < dif_BE,
                             x_hold_rise[0]/60 > dif_BS,
                             SC.alt.value[-1] > PT8]
                    
                    PGen7 = [dif_BS < x_hold_fall[0]/60 < LTC1 < dif_BE,
                             SC.alt.value[-1] > PT8]
                    
                    PGen14= [all(PGen4), all(PGen6), all(PGen7)]

                    # Now we get L1 and L2
                    # any () & all() only accept one argument, so it is necessary to nest generator expressions.
                    # Pass denotes that no action should be taken. != means 'is not equal to'.
                    if any(PGen3):
                        L1 = x_hold_fall[0]/60
                    elif any(PGen14):
                        L1 = LTC1
                    elif any(PGen5):
                        L1 = dif_BS   
                    else:
                        sys.exit(f'[MP_Process] Bad lower integration limit for {Name} ({PlanetName}) transting at {T1} - flag me!')
                    
                    #############################################################################
    
                    PGen8 = [round(dif_BS, 6) <= round(LTC1, 6) < dif_BE < x_hold_rise[0]/60,
                             SC.alt.value[-1] > PT8]    
                    PGen9 = [LTC1 < x_hold_rise[0]/60 < dif_BE,
                             LTC1 < dif_BE]   
                    PGen10= [round(dif_BS,6) <= round(LTC1, 6) < x_hold_rise[0]/60 < dif_BE,
                             SC.alt.value[-1] > PT8]
                    PGen13= [x_hold_fall[0]/60 < dif_T1 < LTC1 < x_hold_rise[0]/60 < dif_T4,
                             SC.alt.value[-1] > PT8]  # For transits that are too long on rise side for one night.
                    PGen11= [all(PGen10), all(PGen13),
                             LTC1 > x_hold_rise[0]/60]
    
                    if all(PGen8): 
                        L2 = dif_BE 
                    elif any(PGen11):
                        L2 = x_hold_rise[0]/60.
                    elif any(PGen9):
                        L2 = LTC1
                    else:
                        sys.exit(f'[MP_Process] Bad upper integration limit for {Name} ({PlanetName}) transting at {T1} - flag me!')

                # Get-out for super-marginal events (eg HAT-P-66).
                if L2 < L1:
                    Internal_Rank = 'X'
                else:
                    pass

            except NameError:
                sys.exit(f'[MP_Process] {Name} ({PlanetName}) has thrown an exception for a partial event at {T1} - look at me!')
                
        else:
            pass  # This will throw an UnboundLocalError if an event has been incorrectly flagged.     
       
        # Choose all segment limits for event visibility metric.
        if Internal_Rank != 'X':
            def findVisibilityLimits(dif_from, dif_to):
                PGen = (dif_from <= L1 and dif_to <= L1) or (dif_from >= L2 and dif_to >= L2)
                if PGen == False:  # If at least part is INSIDE of visibility limit
                    L_from, L_to = max(dif_from, L1), min(dif_to, L2)
                elif PGen == True:
                    L_from, L_to = 0, 0
                
                return L_from, L_to
            
            # First baseline visibility limits
            L3, L4 = findVisibilityLimits(dif_BS, dif_T1)
            # Ingress visibility limits
            L5, L6 = findVisibilityLimits(dif_T1, dif_T2)
            # Full transit limits
            L7, L8 = findVisibilityLimits(dif_T2, dif_T3)
            # Egress limits
            L9, L10 = findVisibilityLimits(dif_T3, dif_T4)
            # Second baseline limits
            L11, L12 = findVisibilityLimits(dif_T4, dif_BE)


            # What times should you put on your observing proposal? For this, UTC is convenient.
            # All previous calculations use BJD_TDB - these are best for publications!
            # This step was previously handled by P2_Timesplitter, but has been shrunk to go here.
            PropStart = midnight.datetime + dt.timedelta(0,0,0,0,0,L1)
            ObsStart_JD = BJD2UTC(PropStart, LockOn)
            
            PropEnd = midnight.datetime + dt.timedelta(0,0,0,0,0,L2)
            ObsEnd_JD = BJD2UTC(PropEnd, LockOn)
            
            # Recovery of base for air mass integration.
            # If it's a two-night event, a slightly extended integration base is needed. Try/except catches this.
            try:
                with np.errstate(invalid='raise'):
                    f3 = interp1d(delta_midnight,
                                  abs(np.sin(np.deg2rad(Targetaltazs.alt.degree)))**(-0.6),
                                  kind='cubic')
                    Sval, Serr = integrate.quad(f3, L1, L2, epsrel=1e-5)
                    Air_Mass_Marker = max(0, 3 * ((L2-L1) / Sval - 2/3)) # Returns an air mass weighting between 0 and 1, where 30 deg = 0.   
                     
            except ValueError:
                ext_delta_midnight = np.linspace(-13, 19, t_samples) * u.hour
                f3 = interp1d(ext_delta_midnight,
                              abs(np.sin(np.deg2rad(Targetaltazs.alt.degree)))**(-0.6),
                              kind='cubic')
                try:
                    with np.errstate(invalid='raise'):
                        Sval, Serr = integrate.quad(f3, L1, L2, epsrel=1e-5)
                        Air_Mass_Marker = max(0, 3 * ((L2-L1) / Sval - 2/3)) # Returns an air mass weighting between 0 and 1, where 30 deg = 0.
                except ValueError: # Second exception for very long transiters (eg. Kepler-432)
                    Sval, Serr = integrate.quad(f3, L1-24, L2-24, epsrel=1e-5)
                    Air_Mass_Marker = max(0, 3 * (((L2-24) - (L1-24)) / Sval - 2/3))

            # Weight by percentage of event + baseline captured.    
            try:
                with np.errstate(invalid='raise'):
                    WBase = ((L4-L3) + (L12-L11)) / 2                  # Baseline
                    WTrans = (L8-L7) / (dif_T3-dif_T2)                     # Full depth
                    WInOut = ((L10-L9) + (L6-L5)) / (2 * (dif_T2-dif_T1))  # Ingress-egress

                    EventPercent = np.around(WBase * WTrans * WInOut, 10)
                   
            except FloatingPointError: # Case grazing transits -- No need to consider WTrans!
                WBase = ((L4-L3) + (L12-L11)) / 2                  # Baseline
                WTrans = np.nan                                    # Full depth (never achieved for grazing events!)
                WInOut = ((L10-L9) + (L6-L5)) / (2 * (dif_T2-dif_T1))  # Ingress-egress

                EventPercent = np.around(WBase * WInOut, 10)

            except ZeroDivisionError:   
                sys.exit(f'[MP_Process] ZeroDivisionError for {Name} ({PlanetName}) transiting at {T1} - flag me!')

            # Acquire moon noise metric
            if Moon_Noise == 'Y_Moon':
                LUT_times = pd.to_datetime(moonlightParams_LUT['time_UTC'])
                ObsStart_JD_dt = ObsStart_JD.to_datetime()
                idx_ObsStart = (LUT_times - ObsStart_JD_dt).abs().idxmin() 

                moon_mag = moonlightParams_LUT[f'{Filter}mag'].iloc[idx_ObsStart]
                Scattering_AOD_AERONET = moonlightParams_LUT[f'AOD_Scattering_{Filter}band'].iloc[idx_ObsStart]
                Extinction_AOD_AERONET = moonlightParams_LUT[f'AOD_Extinction_{Filter}band'].iloc[idx_ObsStart]

                if Scattering_AOD == 'Default':
                    AOD_scatter = Scattering_AOD_AERONET
                else:
                    AOD_scatter = Scattering_AOD

                if Absorption_AOD == 'Default':
                    AOD_absorption = Extinction_AOD_AERONET - Scattering_AOD_AERONET
                else:
                    AOD_absorption = Absorption_AOD

                if Asymmetry_Factor == 'Default':
                    AsymmFactor = moonlightParams_LUT[f'Asymmetry_Factor_{Filter}band'].iloc[idx_ObsStart]
                else:
                    AsymmFactor = Asymmetry_Factor
            
                Moon_Noise_Metric = MoonNoiseMetric(moonaltazs, Targetaltazs,
                                                    moon_mag, AOD_scatter, AOD_absorption, AsymmFactor,
                                                    dif_BS, delta_midnight, dif_BE, Filter)    
            elif Moon_Noise == 'N_Moon':
                Moon_Noise_Metric = 1  # SNR doesnt change if no noise source added

        else:
            ObsStart_JD = np.nan
            ObsEnd_JD = np.nan
            Air_Mass_Marker = 0  # Unobservable 
            WBase = 0
            WTrans = 0
            WInOut = 0
            EventPercent = 0
            Moon_Noise_Metric = 0 # No need to calculate as it is not used
        

        # Create plot only if EventPercent >= minThreshold to reduce unnneccessary images and CPU load
        if EventPercent >= EventPercent_minThreshold:
            # Shift time axis from BJD to local time for that telescope (for graphing purposes)
            # First, find timedelta where midnight_UTC - midnight_BJD in hours
            midnight_UTC = BJD2UTC(midnight, LockOn)
            delta_h_UTC_BJD = (midnight_UTC - midnight).to('hour').value
            # Next, find utc offset where midnight_local - midnight_UTC in hours
            utc_offset = midnight_UTC.to_datetime(timezone=ZoneInfo(timezone_str)).utcoffset() / dt.timedelta(hours=1)
            # Thus, we get offset = local - BJD
            offset = delta_h_UTC_BJD + utc_offset

            # Shift delta_midnight and Moon_step to local time since midnight: local = BJD + offset
            delta_midnight_local = delta_midnight + offset * u.hour
            Moon_step_local = Moon_step + offset * u.hour
            xlim = (delta_midnight_local.min().value, delta_midnight_local.max().value)

            # Define custom colormap
            colors = ['blue', 'aqua']
            cmap = LinearSegmentedColormap.from_list('target_cmap', colors, N=180)
            norm = Normalize(vmin=0, vmax=360)

            # Altitude mask to only plot whatever is in plot
            moonaltmask = (moonaltazs.alt.value >= -1)
            targaltmask = (Targetaltazs.alt.value >= -1)

            # Plot target path and transit viability for all events.
            # Scatter c argument denotes colour, lw=linewidth, default s=20, zorder stops fill from mucking up the colour of your other lines.
            # Use fig.close(fig) to stop all the figures from eating your RAM.
            plt.ioff()
            fig, ax = plt.subplots(figsize=(18, 10))

            ax.plot(delta_midnight_local[moonaltmask], moonaltazs.alt[moonaltmask], markevery=20,
                    color='darkorange', zorder=2, linestyle='-.', linewidth=2.5,
                    label='Lunar Path')                        # Moon path
            sc = ax.scatter(delta_midnight_local[targaltmask], Targetaltazs.alt[targaltmask], c=Targetaltazs.az.value[targaltmask], 
                            cmap=cmap, norm=norm, zorder=3, lw=0, s=15)   # Target path
            ax.plot([], [], color='blue', label='Target Path') # Dummy plot to put on legend
            ax.scatter(Moon_step_local, Targetsepaltazs.alt,
                       facecolors='white', edgecolors='blue', zorder=4, s=100, linewidths=2, marker='o',
                       label='Lunar Separation (deg)')      # Separation scatterplot
            
            # Attaches separations to markers along target path.
            for w, txt in enumerate(FormSep):
                ax.annotate(FormSep[w], (Moon_step_local.value[w], Targetsepaltazs.alt.value[w]),
                            color='violet', xytext=(-7.5, -27), textcoords='offset pixels')
                
            ax.fill_between(delta_midnight_local.to('hr').value, 0, 90, sunaltazs.alt <= -18*u.deg,
                            facecolor='black', alpha=0.95, zorder=0, label='Night')   # Bounded by astronomical twilight
            
            # Red fill denotes transit in progress, superimposed on other zones. alpha denotes transparency.
            ax.axvspan(dif_BS+offset, dif_BE+offset, 0,90, alpha=0.3, facecolor='r', zorder=1,
                       label='Observation')
            ax.axvspan(dif_T1+offset, dif_T4+offset, 0,90, alpha=0.3, facecolor='r', edgecolor='r', hatch='X', zorder=1,
                       label='Transit Active')
            
            ax.axhspan(0, PT8, xlim[0], xlim[1], facecolor='grey', alpha=0.2, zorder=5,
                       label='Unobservable altitude')
            ax.axhspan(LIV,90, xlim[0], xlim[1], facecolor='grey', alpha=0.2, zorder=5)

            fig.colorbar(sc, ax=ax).set_label('Azimuth [deg]')
            ax.legend(loc='upper right')

            textpos = xlim[1] + 5
            ax.text(textpos, 85, f'Moon Illumination:\n{Phase:.0f}%')
            ax.text(textpos, 75, f'Closest Lunar\nSeparation: {min(FormSep):.0f} deg')

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

            ax.set_title(f'{PlanetName} Transit Observability for {Inst} from\n' +\
                         f'{(ObsStart_JD + dt.timedelta(hours=utc_offset)).strftime("%b %d %Y %H:%M")} to ' +\
                         f'{(ObsEnd_JD + dt.timedelta(hours=utc_offset)).strftime("%b %d %Y %H:%M")} (UTC+{utc_offset:.1f})')

            fig.savefig(rf'{PlotPath}/{EventPercent:.3f}_{Internal_Rank}_{Name}_{Inst}_{T1.strftime("%b_%d_%Y_%H_%M_%S")}.jpg',
                        dpi=100, transparent=False, facecolor='white', edgecolor='black', format='jpg')

            plt.close(fig)  # This is still needed to release system resources after plot is saved. Must be explicit!
        
        return Phase, min(Sep), Internal_Rank, ObsStart_JD, ObsEnd_JD, Air_Mass_Marker, WBase, WTrans, WInOut, EventPercent, Moon_Noise_Metric
        # Note the conversion to UTC for our observing start/end times!
    
    def FinalMetric(row):
        FinMet = di[Metric_Mode][0] * row['Air_Mass_Metric'] * row['Event_Weight'] * row['Moon_Noise_Metric']
        return FinMet


    ##### DEFINITIONS END HERE #####
    # Graphs will live here.
    PlotPath = rf'{csvoutpath}/{Inst}_{Filter}-band_{ObsStart.strftime("%b %d %Y")}_to_{ObsEnd.strftime("%b %d %Y")}'
    
    # csv file containing one planet.
    Cols=['Planet', 'R*', 'Rp', 'RA:HH', 'RA:MM', 'RA:SS', 'Dec:DD', 'Dec:MM', 'Dec:SS', f'{Filter}mag',
          'T14', 'Depth', 'T0 (HJD or BJD)', 'T0_err', 'P (day)', 'P_err',
          'Impact Parameter', 'a_Calc', 'Previous Study Flag',
          'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']
    di = pd.read_csv(csvinpath, usecols=Cols, skipinitialspace=True)

    # Altitude spot-check - will the target ever be high enough to be observed?
    di[['RA:HMS', 'Dec:Deg']] = di.apply(lambda row: AltChecker(row), axis=1, result_type='expand')
    PlanetName = di['Planet'][0]
    
    if PT8 >= 90 - abs(ds['Lat'].iloc[S] - di['Dec:Deg'][0]):  # Case physically unobservable target
        pass
    else:
        # Create TargetTimes and ErrTimes list
        T1Calc(di.iloc[0])  # There is only one row anyways

        if len(TargetTimes) == 0:  # Case no transits in observation window:           
            pass

        else:
            # Create T1 and T1_err columns according to Figure 3.6
            # Reindex original di to create all these NaN rows
            di = di.reindex(np.arange(len(TargetTimes)))
            # Concatenate T1 and T1_err
            di = pd.concat([di,
                            pd.DataFrame({'T1': TargetTimes, 'T1_err': ErrTimes})],
                            axis=1)
                            
            # Writes additional key times and results
            di[['T2', 'T0', 'T3', 'T4', 'Baseline_Start', 'Baseline_End']] = di.apply(lambda row: TimesCalc(row),
                                                                                      axis=1, result_type='expand')
            di[['Lunar_Illumination', 'Closest_Lunar_Approach', 'Internal Rank', 'Observation_Start_(UTC)', 'Observation_End_(UTC)',
                'Air_Mass_Metric', 'Baseline_Weight', 'Transit_Curve_Weight', 'Ingress-Egress_Weight', 'Event_Weight',
                'Moon_Noise_Metric']] = di.apply(lambda row: EventMetric(row),
                                                 axis=1, result_type='expand')  
            
            # Final metrics, then export
            di[f'Final_{Metric_Mode}'] = di.apply(lambda row: FinalMetric(row), axis=1)
            di = di.drop(columns=f'{Filter}mag')

            di.to_csv(rf'{csvoutpath}/{Name}', index=False)
    
    # CLEARANCE PROTOCOLS FOR MULTIPLE-TARGET RUNNING   
    TargetTimes.clear()
    ErrTimes.clear()
