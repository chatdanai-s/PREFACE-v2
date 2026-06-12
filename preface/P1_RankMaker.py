# Part 4 of overall pipeline
# Takes in WorkingTEPSet, calculates exposure times, cadences and ranks for a chosen telescope/filter combination.
# Outputs a full ranked .csv for further processing.

import os
import sys
import pandas as pd
import numpy as np
import astropy.units as u
import datetime as dt

from preface.P1_ModCheck import creation_date

# Calculation of sky background. Spectrographs use a square aperture instead (as part of the slit, 10" width).
# Defocusing means we take in more of the sky background!
# This is identical for all ranking parameter functions.
# Put here because the function is later reused in Phase Two!
def findSkyB(mu_sky, scope_df, scope_idx, toggle_defocus, run_mode):
    try:
        if toggle_defocus == True:
            area = np.pi / 4 * scope_df['Theta_DF'].iloc[scope_idx]**2
        elif toggle_defocus != False:
            sys.exit('[RankMaker] Invalid parameter input. Please check toggle_defocus.')

        elif run_mode in ['Spectral_Half_Well', 'IR_Half_Well']:
            area = 10 * 4 * scope_df['Theta_see'].iloc[scope_idx]
        elif run_mode == 'Half_Well':
            area = np.pi * 4 * scope_df['Theta_see'].iloc[scope_idx]**2
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check run_mode.')

        SkyBackground = -2.5 * np.log10(area) + mu_sky

    except IndexError:
        sys.exit('[RankMaker] Telescope name not recognised; a full list is given in Scope.csv.')

    return SkyBackground


# The central rank-making function, containing several sub-routines (listed first).
# csvpath is the the path to the core folder, as usual
def RankMaker(CSV_core_folder, CSV_intermediate_folder,
              scope_df, scope_idx, instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus):
    
    sky_noise_text = 'Y-SkyNoise' if toggle_sky_noise else 'N-SkyNoise'
    defocus_text = 'Y-Defocus' if toggle_defocus else 'N-Defocus'

    fulltepset_csv_path = CSV_core_folder / 'FullTEPSet.csv'
    workingtep_csv_path = CSV_core_folder / 'WorkingTEPSet.csv'
    rankedtep_csv_path = CSV_intermediate_folder / 'ranked_tep_sets' / f'RankedTEPSet_{instrument}_{filter_name}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}.csv'
    (CSV_intermediate_folder / 'ranked_tep_sets').mkdir(parents=True, exist_ok=True)

    df_tep = pd.read_csv(workingtep_csv_path, skipinitialspace=True)

    # Calculation for exposure times, SNR and cadence to fill half of a CCD well.
    # If noise is on, we must include sky backgrounds - they will non-negligibly suppress faint targets.
    # Add these with star magnitudes to get an effective composite magnitude, m_comp.
    def TexpCalc_HalfWell(row):
        mtarget = row[f'{filter_name}mag']
        mzp = scope_df[f'mzp_{filter_name}'].iloc[scope_idx]

        # m_comp evaluation based on Add_Noise
        if toggle_sky_noise == True:
            m_comp = -2.5 * np.log10(10**(-0.4 * mtarget) + 10**(-0.4 * SkyB))
        elif toggle_sky_noise == False:
            m_comp = mtarget
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check toggle_sky_noise.')
        
        # theta_see and texp evaluation based on Defocus
        Omega_pix = scope_df['Omega_pix'].iloc[scope_idx]
        Half_Well = scope_df['Half_Well'].iloc[scope_idx]
        Gain = scope_df['Gain'].iloc[scope_idx]

        if toggle_defocus == True:
            theta = scope_df['Theta_DF'].iloc[scope_idx]
            Texp = (np.pi * (theta/2)**2 / Omega_pix) * (Half_Well / Gain) * 10**(-0.4 * (mzp - m_comp))
        elif toggle_defocus == False:
            theta = scope_df['Theta_see'].iloc[scope_idx]
            Texp = (np.pi * (theta*2)**2 / Omega_pix) / (4 * np.log(2)) * (Half_Well / Gain) * 10**(-0.4 * (mzp - m_comp))            
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check toggle_defocus.')

        # SNR and Data count from exposure time
        SNR = np.sqrt(Texp) * 10**(-0.2 * (mtarget - mzp))
        Overhead = scope_df['Overhead'].iloc[scope_idx]
        Npoints = (row['T14'] * u.d).to(u.s) / ((Texp + Overhead) * u.s)
        
        return Texp, SNR, Npoints
    
    # Calculation for exposure times, SNR and cadence to fill half of a CCD well. (SPECTRO CASE)
    # Our flux is dispersed across many spectroscopic elements, which can each be modelled as a 2D Gaussian.
    # LT Sprat's zero-points are already per element, so this doesn't need to be modelled twice.
    def TexpCalc_HalfWellSpectro(row):
        mtarget = row[f'{filter_name}mag']

        if toggle_sky_noise == True:
            m_comp = -2.5 * np.log10(10**(-0.4 * mtarget) + 10**(-0.4 * SkyB))
        elif toggle_sky_noise == False:
            m_comp = mtarget

        resolution_element_exists = (scope_df['Element?'].iloc[scope_idx] == 'Y')
        if resolution_element_exists:
            Elements = scope_df['Lambda_Range'].iloc[scope_idx] / scope_df['Dispersion'].iloc[scope_idx]
        else:
            Elements = 1

        mzp = scope_df[f'mzp_{filter_name}'].iloc[scope_idx]
        Theta_see = scope_df['Theta_see'].iloc[scope_idx]
        Half_Well = scope_df['Half_Well'].iloc[scope_idx]
        Gain = scope_df['Gain'].iloc[scope_idx]

        m_5Sig = scope_df['m_5Sig'].iloc[scope_idx]
        Lambda_Cent = scope_df['Lambda_Cent'].iloc[scope_idx]
        Res = scope_df['Res'].iloc[scope_idx]
        
        F = 10**(0.4 * (mzp - m_comp)) / (Elements * np.pi) * np.log(2) / Theta_see**2
        
        Texp = (Half_Well / Gain) / F
        # Texp = (ds['Half_Well'].iloc[scope_idx] / ds['Gain'].iloc[scope_idx]) / 25 * \
        #        3600 * 10**(-0.4 * (ds['m_5Sig'].iloc[scope_idx] - mtarget))
        
        # SNR and Data count
        SNR = np.sqrt(25 * (Texp/3600) * 10**(0.4*(m_5Sig - mtarget)) * Lambda_Cent / Res)
        Overhead = scope_df['Res'].iloc[scope_idx]
        Npoints = (row['T14'] * u.d).to(u.s) / ((Texp + Overhead) * u.s)

        return Texp, SNR, Npoints
    

    # Creation and writing of two ranking parameters, for single events with and without T-dependence respectively.
    # Habitable zone studies are only interested in a narrow range of Teq - becomes constant.
    # Now for known and unknown masses - if a sensible mass is present, it is always used.
    # Otherwise, a two-track mass-radius relation is used. n=2 for smaller than Saturn, n=0 for larger.
    def Priority(row, habitable_study=False):
        # Mass-radius variables
        Mp = row['Mp']
        Rp = row['Rp']
        Merr = max(row['Mp_err_up'], row['Mp_err_dn']) / row['Mp']

        # Misc variables
        texp = row[f'Exposure Time for {run_mode}']
        texp_min = scope_df['texp_min'].iloc[scope_idx]
        tover = scope_df['Overhead'].iloc[scope_idx]
        Depth = row['Depth']
        T14_d = row['T14']
        T14_sec = ((T14_d*u.d).to(u.s)).value
        Lambda_range = scope_df['Lambda_Range'].iloc[scope_idx]
        scaleH_earth = 8000

        if habitable_study == True:
            Teq_Calc = 1
        elif habitable_study == False:
            Teq_Calc = row['Teq_Calc']

        # Magnitude variables
        mzp = scope_df[f'mzp_{filter_name}'].iloc[scope_idx]
        mtarget = row[f'{filter_name}mag']
        
        # Are we using the sky backgrounds and out-of-transit noise?
        if toggle_sky_noise == True:
            C_y = scope_df['C_y'].iloc[scope_idx]
            Aperture = scope_df['Aperture'].iloc[scope_idx]
            Scope_Altitude = scope_df['Alt'].iloc[scope_idx]
            Half_Well = scope_df['Half_Well'].iloc[scope_idx]
            Gain = scope_df['Gain'].iloc[scope_idx]

            Noise_Scint = 0.0032 * C_y * Aperture**(-2/3) * texp**(-0.5) * np.exp(-Scope_Altitude / scaleH_earth) * (Half_Well / Gain)
            Noise_R = 10**(0.4 * (mzp-mtarget)) / 10000
            Noise_Exp = np.sqrt(texp * \
                                (10**(0.4 * (mzp-mtarget)) + 10**(0.4 * (mzp-SkyB))) * \
                                (1 + 12 * T14_d) + \
                                Noise_Scint**2 + Noise_R**2
                               )

            Signal = (Lambda_range / 200)**(-0.5) * \
                     10**(0.4 * (mzp-mtarget)) * texp * Teq_Calc * (Depth/100) * \
                     np.sqrt(T14_sec / (texp + tover))
            
            # ULTRASPEC has a minimum exposure time for CCD windowing - some targets will subsequently not be practical.
            if texp <= texp_min:
                Rank = 0.00
            elif ~np.isnan(Mp) and (Merr <= 0.4):
                Rank = Signal * (Rp/Mp) / Noise_Exp
            elif np.isnan(Mp) or np.isnan(Merr) or (Merr > 0.4):
                Rank = Signal * max(Rp/0.8, 0.8/Rp) / Noise_Exp
            
            else:
                sys.exit(f'[RankMaker] Planet fallen through generators at {row.Planet} - flag me!')

            return Rank
        
        elif toggle_sky_noise == False:
            C_T = (Lambda_range/200)**(-0.5) * 10**(0.2 * mzp) * np.sqrt(texp / (texp + tover))
            Flux = 10**(-0.2 * mtarget)

            if texp <= texp_min:
                Rank = 0.00
            elif ~np.isnan(Mp) and (Merr <= 0.4):
                Rank = C_T * Flux * T14_sec * Teq_Calc * (Depth/100) * Rp / Mp
            elif np.isnan(Mp) or np.isnan(Merr) or (Merr > 0.4):
                if Rp >= 0.8:
                    n = 0
                elif Rp < 0.8:
                    n = 2

                Rank = C_T * Flux * T14_sec * Teq_Calc * (Depth/100) * Rp**(1-n) / 0.8

            else:
                sys.exit(f'[RankMaker] Planet fallen through generators at {row.Planet} - flag me!')
            
            return Rank
    
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check toggle_sky_noise.') 
    

    # Ranking parameters for multiple-transit studies in chosen run mode.
    # Masses (or lack thereof) are now in at this point.
    def MultiTransit_Priority(row):
        MultiTransit_Rank = row['Rank'] / np.sqrt(row['P (day)'])
        return MultiTransit_Rank
        
    def MultiTransitHabitable_Priority(row):
        MultiTransitHabitable_Rank = row['Habitable_Rank'] / np.sqrt(row['P (day)'])
        return MultiTransitHabitable_Rank

    # Transmission Spectroscopy Metric (TSM) as defined by Kempton+18
    def TSM(row):
        # Initialize variables
        Rp = row['Rp']
        Rs = row['R*']
        Mp = row['Mp']  # Mp_Calc will be implemented in WorkingTEPSetBuilder soon but not now
        Teq = row['Teq'] if pd.notna(row['Teq']) else row['Teq_Calc']
        mJ = row['Jmag']
        ScaleFactor = 1

        # Convert to appropriate units
        Rp = (Rp * u.jupiterRad).to(u.earthRad).value
        Mp = (Mp * u.jupiterMass).to(u.earthMass).value

        TSM = ScaleFactor * (Rp**3 * Teq) / (Mp * Rs**2) * 10**(-0.2 * mJ)
        return TSM

    ########################################################################################################
    ## END OF DEFINITIONS ##

    print("[RankMaker] Creating the ranked TEPSet file...")
    config_str = f'{instrument}, {filter_name}-band for {run_mode}, {sky_noise_text}, {defocus_text}'
    
    fulltep_creation = creation_date(fulltepset_csv_path)
    fulltep_creation = dt.datetime.fromtimestamp(fulltep_creation)  # Converts date to human-readable timestamp.    
    
    # Set to fire if WorkingTEPSet was recently modified, so will always fire when downloaded .csvs are updated.
    # Alternatively, will also fire when a new valid telescope/filter combo is requested.
    forty_seconds_since_creation = ( dt.datetime.now() - fulltep_creation >= dt.timedelta(0, 40.) )

    if forty_seconds_since_creation or not os.path.exists(rankedtep_csv_path):
        # Get Sky background        
        mu_sky = scope_df[f'msky_{filter_name}'].iloc[scope_idx]
        SkyB = findSkyB(mu_sky, scope_df, scope_idx, toggle_defocus, run_mode)
     
        # Applies functions to filtered data set and creates new columns.
        # For functions returning multiple variables, count over them as a Series to write to individual columns.         
        if run_mode in ['Spectral_Half_Well', 'IR_Half_Well']:
            colnames = ['Exposure Time for Spectral_Half_Well', 'SNR (Spectral_Half_Well)', 'Cadence (Spectral_Half_Well)']
            df_tep[colnames] = df_tep.apply(lambda row: TexpCalc_HalfWellSpectro(row),
                                            axis=1, result_type='expand')
            
        elif run_mode == 'Half_Well':
            colnames = ['Exposure Time for Half_Well', 'SNR (Half_Well)', 'Cadence (Half_Well)']
            df_tep[colnames] = df_tep.apply(lambda row: TexpCalc_HalfWell(row),
                                            axis=1, result_type='expand')    
                      
        else:
            sys.exit("[RankMaker] Please select either 'Half_Well' for photometry, 'Spectral_Half_Well' for optical spectroscopy or 'IR_Half_Well' for NIR spectroscopy.")

        # Apply our functions.
        df_tep['Rank'] = df_tep.apply(lambda row: Priority(row, habitable_study=False), axis=1)
        df_tep['Habitable_Rank'] = df_tep.apply(lambda row: Priority(row, habitable_study=True), axis=1)
        df_tep['Multi_Transit_Rank'] = df_tep.apply(lambda row: MultiTransit_Priority(row), axis=1)
        df_tep['Multi_Transit_Habitable_Rank'] = df_tep.apply(lambda row: MultiTransitHabitable_Priority(row), axis=1)
        df_tep['TSM'] = df_tep.apply(lambda row: TSM(row), axis=1)

        # Sort by chosen column (signal strength parameter), high to low.
        df_tep = df_tep.sort_values(by='Rank', ascending=False)
        
        # Kill-switch to fire if an invalid telescope-filter combination is selected.
        if pd.isnull(max(df_tep['Rank'])):
            sys.exit('[RankMaker] Invalid telescope-filter combo selected; please check!')
        else:
            pass
        
        # Writes dataframe to csvs.
        df_tep.to_csv(rankedtep_csv_path, index=False)

        print(rf'[RankMaker] Ranked TEP list for {config_str} generated.')

    else:
        print(rf'[RankMaker] Ranked TEP list for {config_str} is up to date.')
    