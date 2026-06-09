# Part 0 of overall pipeline.
# Check if all inputs are valid.
import os
import sys
from datetime import timedelta
import numpy as np
import pandas as pd
import joblib


def isNumber(var):
    return isinstance(var, (int, float)) and not np.isnan(var)

def InputCheck(csvcorepath, csvrankpath, csvcutpath, csvinpath, csvoutpath, csvendpath, csvscorepath,
               Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
               Moon_Noise, Scattering_AOD, Absorption_AOD, Asymmetry_Factor, Amplification_Factor,
               EventWeight_GraphThreshold, ObsStart, ObsEnd, ncores):
    InputErrorFlag = False

    # Check if packages are properly installed (except lime_tbx)
    for pkg in ['astropy', 'matplotlib', 'scipy', 'tqdm', 'ephem', 'joblib']:
        try:
            __import__(pkg)
        except ImportError:
            InputErrorFlag = True
            print(f"[InputCheck] Missing package detected: {pkg}")


    # Check if file paths exist
    for path in [csvcorepath, csvrankpath, csvcutpath, csvinpath, csvoutpath, csvendpath, csvscorepath]:
        if not os.path.exists(path):
            InputErrorFlag = True
            print(f"[InputCheck] Invalid folder path detected. {path} does not exist!")

    # Creates lookup table folder for csvcorepath
    if os.path.exists(csvcorepath):
        LUT_folder = rf'{csvcorepath}/Lookup_tables'
        if not os.path.exists(LUT_folder):
            os.mkdir(LUT_folder)


    # Check instrument configurations
    try:
        ds = pd.read_csv(rf'{csvcorepath}/Scope.csv')
        S = np.where(ds['Telescope'] == Inst)[0][0]

        mzp = ds[f'mzp_{Filter}'].iloc[S]
        msky = ds[f'msky_{Filter}'].iloc[S]
        if not (isNumber(mzp) and isNumber(msky)):
            InputErrorFlag = True
            print('[InputCheck] Invalid Filter detected. A full list is given in Scope.csv.')
    
    except IndexError:
        InputErrorFlag = True
        print('[InputCheck] Invalid Telescope name detected. A full list is given in Scope.csv.')
    except KeyError:
        InputErrorFlag = True
        print('[InputCheck] Invalid Filter detected. A full list is given in Scope.csv.')    

    if Run_Mode not in ['Half_Well', 'Spectral_Half_Well', 'IR_Half_Well']:
        InputErrorFlag = True
        print("[InputCheck] Invalid Run_Mode detected. Is it 'Half_Well', 'Spectral_Half_Well', or 'IR_Half_Well'?")

    try:
        Lambda_Cent = ds['Lambda_Cent'].iloc[S]
        Dispersion = ds['Dispersion'].iloc[S]
        if (Run_Mode in ['Spectral_Half_Well', 'IR_Half_Well']) and not (isNumber(Lambda_Cent) and isNumber(Dispersion)):
            InputErrorFlag = True
            print(f"[InputCheck] Run_Mode is not available for this telescope! Set Run_Mode = 'Half_Well'.")
    except UnboundLocalError:
        if Run_Mode in ['Spectral_Half_Well', 'IR_Half_Well']:
            print("[InputCheck] Run_Mode validity cannot be checked due to invalid telescope index.")

    if Add_Noise not in ['Y_Noise', 'N_Noise']:
        InputErrorFlag = True
        print("[InputCheck] Invalid Add_Noise detected. Is it 'Y_Noise' or 'N_Noise'?")
    
    if Defocus not in ['Y', 'N']:
        InputErrorFlag = True
        print("[InputCheck] Invalid Defocus detected. Is it 'Y' or 'N'?")

    try:
        theta_DF = ds['Theta_DF'].iloc[S]
        if (Defocus == 'Y') and not isNumber(theta_DF):
            InputErrorFlag = True
            print("[InputCheck] Defocus is not available for this telescope! Set Defocus = 'N'.")
    except UnboundLocalError:
        if Defocus == 'Y':
            print("[InputCheck] Defocus validity cannot be checked due to invalid telescope index.")

    if Metric_Mode not in ['Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']:
        InputErrorFlag = True
        print("[InputCheck] Invalid Metric_Mode detected. Is it 'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', or 'Multi_Transit_Habitable_Rank'?") 

    if (not isNumber(ViableCut)) or (not (0 < ViableCut <= 1)):
        InputErrorFlag = True
        print("[InputCheck] Invalid ViableCut detected. Is it a number between 0 and 1?") 


    # Check moon background configurations
    if Moon_Noise not in ['Y_Moon', 'N_Moon']:
        InputErrorFlag = True
        print("[InputCheck] Invalid Moon_Noise detected. Is it 'Y_Moon' or 'N_Moon'?")

    if Moon_Noise == 'Y_Moon':
        # Check validity of lime_tbx setup
        try:
            import lime_tbx
            from lime_tbx.common.logger import get_logger
            from lime_tbx.persistence.local_storage.appdata import get_appdata_folder
            logger = get_logger()
            appdata_folder = get_appdata_folder(logger)

            main_eocfi_path = rf'{appdata_folder}\eocfi_data'
            main_kernels_path = rf'{appdata_folder}\kernels'

            for path in [main_eocfi_path, main_kernels_path]:
                if not os.path.exists(path):
                    InputErrorFlag = True
                    print(f"[InputCheck] Invalid folder paths for lime_tbx modeling detected. {path} does not exist!")
                    print(' '*13 + 'Please refer to README_lime_tbx_installation_guide.txt on installing lime_tbx.')

        except ImportError as e:
            InputErrorFlag = True
            print(f"[InputCheck] Missing package detected: {e.name}")
            print(' '*13 + 'Please refer to README_lime_tbx_installation_guide.txt on installing lime_tbx.')

        # Check other variables associated with moon background
        available_filters = ['U','B','V','R','I','u','g','r','i','z']
        if Filter not in available_filters:
            InputErrorFlag = True
            print("[InputCheck] Invalid Filter for Moon_noise detected. Is it UBVRI or ugriz?")

        if (Scattering_AOD != 'Default') and (not isNumber(Scattering_AOD) or Scattering_AOD < 0):
            InputErrorFlag = True
            print("[InputCheck] Invalid Scattering_AOD detected. Is it 'Default' or at least 0?")

        if (Absorption_AOD != 'Default') and (not isNumber(Absorption_AOD) or Absorption_AOD < 0):
            InputErrorFlag = True
            print("[InputCheck] Invalid Absorption_AOD detected. Is it 'Default' or at least 0?")

        if (Asymmetry_Factor != 'Default') and (not isNumber(Asymmetry_Factor) or not (-1 < Asymmetry_Factor < 1)):
            InputErrorFlag = True
            print("[InputCheck] Invalid Asymmetry_Factor detected. Is it 'Default' or between -1 and +1? (0.5-0.8 recommended)")

        if not isNumber(Amplification_Factor):
            InputErrorFlag = True
            print("[InputCheck] Invalid Amplification_Factor detected. Is it a number? (5 is default)")


    # Check EventWeight_GraphThreshold
    if (not isNumber(EventWeight_GraphThreshold)) or (not (0 <= EventWeight_GraphThreshold)):
        InputErrorFlag = True
        print('[InputCheck] Invalid EventWeight_GraphThreshold detected. Is it a number at least 0?')

    # Check observation time and multiprocessing configurations
    if ObsEnd <= ObsStart:
        InputErrorFlag = True
        print('[InputCheck] Invalid Observation times. ObsEnd is earlier than ObsStart!')
    elif ObsEnd - ObsStart < timedelta(hours=1):
        InputErrorFlag = True
        print('[InputCheck] Invalid Observation times. Observation period must be at least one hour!')

    cpu_count = joblib.cpu_count()
    if not isNumber(ncores) or (int(ncores) != ncores) or (ncores <= 0):
        InputErrorFlag = True
        print('[InputCheck] Invalid ncores detected. Is it a positive integer?')
    elif ncores > cpu_count:
        InputErrorFlag = True
        print('[InputCheck] Invalid ncores detected. More CPU cores than available specified.')
    elif ncores == cpu_count:
        print('[InputCheck] WARNING: All CPU cores will be called during multiprocessing. Your machine will likely lock up!')
        continueFlag = input(' '*13 + 'Are you sure you want to continue? [Y/N]: ')
        if continueFlag == 'Y':
            pass
        else:
            sys.exit('[InputCheck] Pipeline terminated.')

    # Check InputErrorFlag
    if InputErrorFlag == True:
        sys.exit('[InputCheck] Invalid input(s) found -- Pipeline terminated.')
        