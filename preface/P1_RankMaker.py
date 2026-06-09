# Part 4 of overall pipeline
# Takes in WorkingTEPSet, calculates exposure times, cadences and ranks for a chosen telescope/filter combination.
# Outputs a full ranked .csv for further processing.

import os
import sys
import pandas as pd
import numpy as np
import astropy.units as u
import datetime as dt

from ModCheck import creation_date

# Calculation of sky background. Spectrographs use a square aperture instead (as part of the slit, 10" width).
# Defocusing means we take in more of the sky background!
# This is identical for all ranking parameter functions.
# Put here because the function is later reused in Phase Two!
def findSkyB(mu_sky, ds, S, Defocus, Run_Mode):
    try:
        if Defocus == 'Y':
            area = np.pi / 4 * ds['Theta_DF'].iloc[S]**2
        elif Defocus != 'N':
            sys.exit('[RankMaker] Invalid parameter input. Please check Defocus.')

        elif Run_Mode in ['Spectral_Half_Well', 'IR_Half_Well']:
            area = 10 * 4 * ds['Theta_see'].iloc[S]
        elif Run_Mode == 'Half_Well':
            area = np.pi * 4 * ds['Theta_see'].iloc[S]**2
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check Run_Mode.')

        SkyB = -2.5 * np.log10(area) + mu_sky

    except IndexError:
        sys.exit('[RankMaker] Telescope name not recognised; a full list is given in Scope.csv.')

    return SkyB

# The central rank-making function, containing several sub-routines (listed first).
# csvpath is the the path to the core folder, as usual
def RankMaker(csvpath, csvrankpath, ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus):
    di = pd.read_csv(rf'{csvpath}/WorkingTEPSet.csv', skipinitialspace=True)

    # Calculation for exposure times, SNR and cadence to fill half of a CCD well.
    # If noise is on, we must include sky backgrounds - they will non-negligibly suppress faint targets.
    # Add these with star magnitudes to get an effective composite magnitude, m_comp.
    def TexpCalcHalfWell(row):
        mtarget = row[f'{Filter}mag']
        mzp = ds[f'mzp_{Filter}'].iloc[S]

        # m_comp evaluation based on Add_Noise
        if Add_Noise == 'Y_Noise':
            m_comp = -2.5 * np.log10(10**(-0.4 * mtarget) + 10**(-0.4 * SkyB))
        elif Add_Noise == 'N_Noise':
            m_comp = mtarget
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check Add_Noise.')
        
        # theta_see and texp evaluation based on Defocus
        if Defocus == 'Y':
            theta = ds['Theta_DF'].iloc[S]
            Texp = (np.pi * (theta/2)**2 / ds['Omega_pix'].iloc[S]) * \
                   (ds['Half_Well'].iloc[S] / ds['Gain'].iloc[S]) * \
                   10**(-0.4 * (mzp - m_comp))

        elif Defocus == 'N':
            theta = ds['Theta_see'].iloc[S]
            Texp = (np.pi * (theta*2)**2 / ds['Omega_pix'].iloc[S]) / (4 * np.log(2)) *\
                   (ds['Half_Well'].iloc[S] / ds['Gain'].iloc[S]) * \
                   10**(-0.4 * (mzp - m_comp))            
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check Denoise.')

        # SNR and Data count from exposure time
        S_N = np.sqrt(Texp) * 10**(-0.2 * (mtarget - mzp))
        Npts = (row['T14'] * u.d).to(u.s) / ((Texp + ds['Overhead'].iloc[S]) * u.s)
        
        return Texp, S_N, Npts
    
    # Calculation for exposure times, SNR and cadence to fill half of a CCD well. (SPECTRO CASE)
    # Our flux is dispersed across many spectroscopic elements, which can each be modelled as a 2D Gaussian.
    # LT Sprat's zero-points are already per element, so this doesn't need to be modelled twice.
    def TexpCalcHalfWellSpectro(row):
        mtarget = row[f'{Filter}mag']

        if Add_Noise == 'Y_Noise':
            m_comp = -2.5 * np.log10(10**(-0.4 * mtarget) + 10**(-0.4 * SkyB))
        elif Add_Noise == 'N_Noise':
            m_comp = mtarget

        if ds['Element?'].iloc[S] == 'Y':
            Elements = ds['Lambda_Range'].iloc[S] / ds['Dispersion'].iloc[S]
        else:
            Elements = 1

        F = 10**(0.4 * (ds[f'mzp_{Filter}'].iloc[S] - m_comp)) / (Elements * np.pi) * \
            np.log(2) / ds['Theta_see'].iloc[S]**2
        
        Texp = (ds['Half_Well'].iloc[S] / ds['Gain'].iloc[S]) / F
        # Texp = (ds['Half_Well'].iloc[S] / ds['Gain'].iloc[S]) / 25 * \
        #        3600 * 10**(-0.4 * (ds['m_5Sig'].iloc[S] - mtarget))
        
        S_N = np.sqrt(25 * (Texp/3600) * 10**(0.4*(ds['m_5Sig'].iloc[S] - mtarget)) * \
                      ds['Lambda_Cent'].iloc[S] / ds['Res'].iloc[S]
                      )
        
        Npts = (row['T14'] * u.d).to(u.s) / ((Texp + ds['Overhead'].iloc[S]) * u.s)

        return Texp, S_N, Npts
    

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
        texp = row[f'Exposure Time for {Run_Mode}']
        texp_min = ds['texp_min'].iloc[S]
        tover = ds['Overhead'].iloc[S]
        Depth = row['Depth']
        T14 = row['T14']
        Lambda_range = ds['Lambda_Range'].iloc[S]
        scaleH_earth = 8000

        if habitable_study == True:
            Teq_Calc = 1
        elif habitable_study == False:
            Teq_Calc = row['Teq_Calc']

        # Magnitude variables
        mzp = ds[f'mzp_{Filter}'].iloc[S]
        mtarget = row[f'{Filter}mag']
        
        # Are we using the sky backgrounds and out-of-transit noise?
        if Add_Noise == 'Y_Noise':
            Noise_Scint = 0.0032 * ds['C_y'].iloc[S] * ds['Aperture'].iloc[S]**(-2/3) * \
                          texp**(-0.5) * np.exp(-ds['Alt'].iloc[S] / scaleH_earth) * \
                          ds['Half_Well'].iloc[S] / ds['Gain'].iloc[S]

            Noise_R = 10**(0.4 * (mzp-mtarget)) / 10000

            Noise_Exp = np.sqrt(texp * \
                                (10**(0.4 * (mzp-mtarget)) + 10**(0.4 * (mzp-SkyB))) * \
                                (1 + 12 * T14) + \
                                Noise_Scint**2 + Noise_R**2
                               )

            Signal = (Lambda_range / 200)**(-0.5) * \
                     10**(0.4 * (mzp-mtarget)) * texp * Teq_Calc * (Depth/100) * \
                     np.sqrt(((T14 * u.d).to(u.s)).value / (texp + tover))
            
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
        
        elif Add_Noise == 'N_Noise':
            C_T = (Lambda_range/200)**(-0.5) * 10**(0.2 * mzp) * np.sqrt(texp / (texp + tover))
            Flux = 10**(-0.2 * mtarget)

            if texp <= texp_min:
                Rank = 0.00
            elif ~np.isnan(Mp) and (Merr <= 0.4):
                Rank = C_T * Flux * ((T14*u.d).to(u.s)).value * Teq_Calc * (Depth/100) * Rp / Mp
            elif np.isnan(Mp) or np.isnan(Merr) or (Merr > 0.4):
                if Rp >= 0.8:
                    n = 0
                elif Rp < 0.8:
                    n = 2

                Rank = C_T * Flux * ((T14*u.d).to(u.s)).value * Teq_Calc * (Depth/100) * Rp**(1-n) / 0.8

            else:
                sys.exit(f'[RankMaker] Planet fallen through generators at {row.Planet} - flag me!')
            return Rank
    
        else:
            sys.exit('[RankMaker] Invalid parameter input. Please check Add_Noise.') 
    

    # Ranking parameters for multiple-transit studies in chosen run mode.
    # Masses (or lack thereof) are now in at this point.
    def MultiTPriority(row):
        MTRank = row['Rank'] / np.sqrt(row['P (day)'])
        return MTRank
        
    def MultiTHPriority(row):
        MTHRank = row['Habitable_Rank'] / np.sqrt(row['P (day)'])
        return MTHRank
            

    ########################################################################################################
    ## END OF DEFINITIONS ##

    print("[RankMaker] Creating the ranked TEPSet file...")
    
    D = creation_date(rf'{csvpath}/FullTEPSet.csv')
    LegD = dt.datetime.fromtimestamp(D)  # Converts date to human-readable timestamp.    
    
    # Set to fire if WorkingTEPSet was recently modifed, so will always fire when downloaded .csvs are updated.
    # Alternatively, will also fire when a new valid telescope/filter combo is requested.
    Gen1 = [dt.datetime.now() - LegD >= dt.timedelta(0, 40.), 
            not os.path.exists(rf'{csvrankpath}/RankedTEPSet_{Inst}_{Filter}-band_for_{Run_Mode}_{Add_Noise}_{Defocus}.csv')
           ]
    # Gen2 = pd.notnull(ds['m_5Sig'].iloc[S])

    # if True:
    if any(Gen1):
        # Get Sky background        
        SkyB = findSkyB(ds[f'msky_{Filter}'].iloc[S], ds, S, Defocus, Run_Mode)
     
        # Applies functions to filtered data set and creates new columns.
        # For functions returning multiple variables, count over them as a Series to write to individual columns.         
        if Run_Mode in ['Spectral_Half_Well', 'IR_Half_Well']:
            colnames = ['Exposure Time for Spectral_Half_Well', 'SNR (Spectral_Half_Well)','Cadence (Spectral_Half_Well)']
            di[colnames] = di.apply(lambda row: TexpCalcHalfWellSpectro(row),
                                    axis=1, result_type='expand')
            
        elif Run_Mode == 'Half_Well':
            colnames = ['Exposure Time for Half_Well', 'SNR (Half_Well)','Cadence (Half_Well)']
            di[colnames] = di.apply(lambda row: TexpCalcHalfWell(row),
                                    axis=1, result_type='expand')    
                      
        else:
            sys.exit('[RankMaker] Please select either Half_Well for photometry, Spectral_Half_Well for optical spectroscopy or IR_Half_Well for NIR spectroscopy.')

        # Apply our functions.
        di['Rank'] = di.apply(lambda row: Priority(row, habitable_study=False), axis=1)    
        di['Habitable_Rank'] = di.apply(lambda row: Priority(row, habitable_study=True), axis=1)    
        di['Multi_Transit_Rank'] = di.apply(lambda row: MultiTPriority(row), axis=1)    
        di['Multi_Transit_Habitable_Rank'] = di.apply(lambda row: MultiTHPriority(row), axis=1)    
        
        # Sort by chosen column (signal strength parameter), high to low.
        di = di.sort_values(by='Rank', ascending=False)
        
        # Kill-switch to fire if an invalid telescope-filter combination is selected.
        if pd.isnull(max(di['Rank'])):
            sys.exit('[RankMaker] Invalid telescope-filter combo selected; please check!')
        else:
            pass
        
        # Writes dataframe to csvs.
        di.to_csv(rf'{csvrankpath}/RankedTEPSet_{Inst}_{Filter}-band_for_{Run_Mode}_{Add_Noise}_{Defocus}.csv',
                  index=False)

        print(rf'[RankMaker] Ranked list for {Inst}, {Filter}-band for {Run_Mode}, {Add_Noise}, {Defocus}-defocus generated.')

    else:
        print(rf'[RankMaker] Ranked list for {Inst}, {Filter}-band for {Run_Mode}, {Add_Noise}, {Defocus}-defocus is up to date.')
    