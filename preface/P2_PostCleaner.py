# Part 8 of the pipeline - concatenates all events for all viable planets into final data products.

# Import from future to access this!
# Note: Imports from future must always go first.

import os
import glob

import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

from .P1_ImpactMerger import BWriter


# With both phases fired, all outputs must be brought together.
# Concatenation function
def get_merged_csv(filelist, ncores):
    dfs = Parallel(n_jobs=ncores, backend='loky')(delayed(pd.read_csv)(f)
                                                  for f in tqdm(filelist, desc='CSVs merged'))
    df = pd.concat(dfs, ignore_index=True)
    return df


def Cleaner(csvrankpath, csvoutpath, csvendpath, csvscorepath,
            ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
            ObsStart, ObsEnd, ncores):
    print("[P2_PostCleaner] Merging all CSV outputs...")

    # COMMON ROUTE: First, combine all output CSVs to one big dataframe
    # P2 input finders (don't touch these!)
    filename_pattern = rf"{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{Metric_Mode}_{ViableCut}_Cut_*.csv"
    # Create filename pattern to glob that includes the path of all our output files.
    fmask = os.path.normpath(os.path.join(csvoutpath, filename_pattern))
    
    bigdf = get_merged_csv(glob.glob(fmask), ncores)  # Get all .csvs that match our pattern and cut together.
    # glob is not fussy about the order in which it reads files in - we have to do a sort by final metric!
    Cols = ['Planet', 'R*', 'Rp', 'RA:HH', 'RA:MM', 'RA:SS', 'Dec:DD', 'Dec:MM', 'Dec:SS',
            'T14', 'Depth', 'T0 (HJD or BJD)', 'T0_err', 'P (day)', 'P_err',
            'Impact Parameter', 'a_Calc', 'Previous Study Flag',
            'TSM', 'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank',
            'RA:HMS', 'Dec:Deg']        # Cols to Forward Fill NA
    bigdf[Cols] = bigdf[Cols].ffill()   # Forward filling NA
    
    bigdf = bigdf[bigdf['Internal Rank'] != 'X'] # Discard unobservable events.
    bigdf = bigdf[bigdf[f'Final_{Metric_Mode}'] != 0]
    bigdf.sort_values(by=[Metric_Mode, f'Final_{Metric_Mode}'], ascending=[False, False], inplace=True)  # Sort bigdf appropriately


    # Then, grab a few extra columns which are handy for observers to have -- dq
    dq = pd.read_csv(rf'{csvrankpath}/RankedTEPSet_{Inst}_{Filter}-band_for_{Run_Mode}_{Add_Noise}_{Defocus}.csv',
                     usecols=['Planet','Teff','[Fe/H]','M*','logg*','Mp','gp','Teq', f'{Filter}mag', f'Exposure Time for {Run_Mode}'])
    dq['Table Position'] = dq.index + 1


    # ROUTE 1: Full list of ranked events sorted by final metric -- dw
    dw = dq.sort_values(by='Planet')
    
    # Recovers list of target names appearing in both data sets then format
    dw = pd.merge(bigdf, dw, how='outer', indicator=True)
    dw = dw.drop(['_merge'], axis=1)
    dw = dw.dropna(subset=['R*'])
    
    dw['FormattedPlanet'] = dw.apply(lambda row: BWriter(row), axis=1)
    dw['FormattedPlanet'] = dw['FormattedPlanet'].replace({'-00':'-', '-0':'-',
                                                           '_00':'_', '_0':'_'},
                                                  regex=True)    
    dw['FormattedPlanet'] = dw['FormattedPlanet'].str.lower()    # exoplanet.eu uses lower-case.
    
    
    # These planets orbit one component of a binary - must be specified individually.
    dw['FormattedPlanet'] = dw['FormattedPlanet'].str.replace('wasp-77', 'wasp-77_a')
    
    # Sort the order of our columns.
    dw = dw[['Planet','FormattedPlanet','Teff','[Fe/H]','M*','R*','logg*', f'{Filter}mag','Mp','Rp','gp','Teq',
             'RA:HH','RA:MM','RA:SS','RA:HMS','Dec:DD','Dec:MM','Dec:SS','Dec:Deg',
             'T14','Depth','T0 (HJD or BJD)','T0_err','P (day)','P_err','Impact Parameter','a_Calc',
             'Previous Study Flag','Rank','Habitable_Rank','Multi_Transit_Rank','Multi_Transit_Habitable_Rank','Internal Rank',
             'Baseline_Start','T1','T1_err','T2','T0','T3','T4','Baseline_End','Lunar_Illumination','Closest_Lunar_Approach',
             'Observation_Start_(UTC)','Observation_End_(UTC)','Air_Mass_Metric','Baseline_Weight','Transit_Curve_Weight',
             'Ingress-Egress_Weight','Event_Weight','Moon_Noise_Metric',f'Final_{Metric_Mode}', f'Exposure Time for {Run_Mode}']]
    
    dw.sort_values(by=f'Final_{Metric_Mode}', ascending=False, inplace=True)    # Sort by our final metric.    

    # Export the final DataFrame to a csv file.
    dw.to_csv(rf'{csvendpath}/Full_Event_List_for_{Inst}_{Filter}-band,{Run_Mode}_{Metric_Mode}-modes,{Add_Noise},{Defocus},{ViableCut}_Cut, {ObsStart.strftime("%b_%d_%Y")}_to_{ObsEnd.strftime("%b_%d_%Y")}.csv', 
              index=False)
    print("[P2_PostCleaner] Full list of ranked events created.")


    # ROUTE 2: Full list of cumulative scores -- dc
    # Cumulative observability calculation, also carrying forward previous study flag.
    du = bigdf.filter(['Planet', 'Previous Study Flag'], axis=1)
    du.drop_duplicates(subset='Planet', keep='first', inplace=True)
    du.reset_index(inplace=True)  # du is a df of planets and its study flag score
    del du['index']
    
    B = []
    for i in np.arange(len(du)):
        bigdg = bigdf[bigdf['Planet'].str.contains(du['Planet'][i])]   # Recover all events for that planet specifically
        bigdh = bigdg[bigdg['Event_Weight'] >= 0.5]                    # These are events for which at least all transit + one baseline is visible.
        CumObsScore = bigdh[f'Final_{Metric_Mode}'].sum()              # Sum scores of all the 'good' events from previous.
        B.append([du['Planet'][i], du['Previous Study Flag'][i], CumObsScore, bigdg.shape[0], bigdh.shape[0]])  # Append to a list.
    
    COS = pd.DataFrame(B)  # Convert to DataFrame called Cumulative Obseravability Score.
    COS.columns = ['Planet', 'Previous Study Flag', f'CumObsScore_{Inst}', f'{Inst}_Event_No', f'{Inst}_Event_No (W>=0.5)']
    COS = COS[COS[f'{Inst}_Event_No (W>=0.5)'] != 0]                        # Drop all planets with no good events.
    COS.sort_values(f'CumObsScore_{Inst}', ascending=False, inplace=True)   # Sort and tidy!
    
    # Recover ranked table positions for cumulative score .csv.
    dc = pd.merge(COS, dq, how='inner', indicator=True)
   
    # Rearrange columns
    dc = dc[['Planet', 'Previous Study Flag', f'CumObsScore_{Inst}', f'{Inst}_Event_No', f'{Inst}_Event_No (W>=0.5)', 
             'Table Position']]
              
    dc.to_csv(rf'{csvscorepath}/CumulativeObsScore_for_{Inst}_{Filter}-band,{Run_Mode}_{Metric_Mode}-modes,{Add_Noise},{Defocus},{ViableCut}_Cut, {ObsStart.strftime("%b_%d_%Y")}_to_{ObsEnd.strftime("%b_%d_%Y")}.csv',
              index=False)
    
    print('[P2_PostCleaner] Cumulative observability scores calculated.')
