# Part 8 of the pipeline - concatenates all events for all viable planets into final data products.

# Import from future to access this!
# Note: Imports from future must always go first.
import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

from preface.P1_ImpactMerger import BWriter


# With both phases fired, all outputs must be brought together.
# Concatenation function
def get_merged_csv(filelist, cores_used):
    dfs = Parallel(n_jobs=cores_used, backend='loky')(
        delayed(pd.read_csv)(f) for f in tqdm(filelist, desc='>> CSVs merged')
    )
    return pd.concat(dfs, ignore_index=True)


def Cleaner(CSV_intermediate_folder, filename_pattern, output_folder,
            instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode, viable_cumulative_cut,
            obs_start, obs_end, cores_actually_used):
    
    print("[P2_PostCleaner] Merging all CSV outputs...")

    # COMMON ROUTE: First, combine all output CSVs to one big dataframe
    filelist = list(
        (output_folder / 'phase_2' / 'individual_planets').glob(filename_pattern)
    )
    bigdf = get_merged_csv(filelist, cores_actually_used)  # Get all .csvs that match our pattern and cut together.


    # glob is not fussy about the order in which it reads files in - we have to do a sort by final metric!
    fillNA_Cols = ['Planet', 'R*', 'Rp', 'RA:HH', 'RA:MM', 'RA:SS', 'Dec:DD', 'Dec:MM', 'Dec:SS',
        'T14', 'Depth', 'T0 (HJD or BJD)', 'T0_err', 'P (day)', 'P_err',
        'Impact Parameter', 'a_Calc', 'Previous Study Flag',
        'TSM', 'Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank',
        'RA:HMS', 'Dec:Deg']                          # Cols to Forward Fill NA
    bigdf[fillNA_Cols] = bigdf[fillNA_Cols].ffill()   # Forward filling NA

    bigdf = bigdf[bigdf['Internal Rank'] != 'X'] # Discard unobservable events.
    bigdf = bigdf[bigdf[f'Final_{metric_mode}'] != 0]
    bigdf.sort_values(by=[metric_mode, f'Final_{metric_mode}'], ascending=[False, False], inplace=True)  # Sort bigdf appropriately


    # Then, grab a few extra columns which are handy for observers to have -- ranked_df
    sky_noise_text = 'Y-SkyNoise' if toggle_sky_noise else 'N-SkyNoise'
    defocus_text = 'Y-Defocus' if toggle_defocus else 'N-Defocus'
    rankedtep_csv_path = CSV_intermediate_folder / 'ranked_tep_sets' / f'RankedTEPSet_{instrument}_{filter_name}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}.csv'

    ranked_df = pd.read_csv(rankedtep_csv_path,
        usecols=['Planet','Teff','[Fe/H]','M*','logg*','Mp','gp','Teq','Teq_Calc','Mp_Calc', f'{filter_name}mag', f'Exposure Time for {run_mode}'])
    ranked_df['Table Position'] = ranked_df.index + 1


    # ROUTE 1: Full list of ranked events sorted by final metric -- full_df
    full_df = ranked_df.sort_values(by='Planet')
    
    # Recovers list of target names appearing in both data sets then format
    full_df = pd.merge(bigdf, full_df, how='outer', indicator=True)
    full_df = full_df.drop(['_merge'], axis=1)
    full_df = full_df.dropna(subset=['R*'])
    
    full_df['FormattedPlanet'] = full_df.apply(lambda row: BWriter(row), axis=1)
    full_df['FormattedPlanet'] = full_df['FormattedPlanet'].replace({'-00':'-', '-0':'-',
                                                                    '_00':'_', '_0':'_'},
                                                                    regex=True)    
    full_df['FormattedPlanet'] = full_df['FormattedPlanet'].str.lower()    # exoplanet.eu uses lower-case.
    
    
    # These planets orbit one component of a binary - must be specified individually.
    full_df['FormattedPlanet'] = full_df['FormattedPlanet'].str.replace('wasp-77', 'wasp-77_a')
    
    # Sort the order of our columns.
    full_df = full_df[
        ['Planet','FormattedPlanet','Teff','[Fe/H]','M*','R*','logg*',f'{filter_name}mag',
        'Mp','Mp_Calc','Rp','gp','Teq','Teq_Calc',
        'RA:HH','RA:MM','RA:SS','RA:HMS','Dec:DD','Dec:MM','Dec:SS','Dec:Deg',
        'T14','Depth','T0 (HJD or BJD)','T0_err','P (day)','P_err','Impact Parameter','a_Calc',
        'Previous Study Flag','TSM','Rank','Habitable_Rank','Multi_Transit_Rank','Multi_Transit_Habitable_Rank','Internal Rank',
        'Baseline_Start','T1','T1_err','T2','T0','T3','T4','Baseline_End','Lunar_Illumination','Closest_Lunar_Approach',
        'Observation_Start_(UTC)','Observation_End_(UTC)','Air_Mass_Metric','Baseline_Weight','Transit_Curve_Weight',
        'Ingress-Egress_Weight','Event_Weight','Moon_Noise_Metric',f'Final_{metric_mode}',f'Exposure Time for {run_mode}']
    ]
    
    full_df.sort_values(by=f'Final_{metric_mode}', ascending=False, inplace=True)    # Sort by our final metric.    

    # Export the final DataFrame to a csv file.
    fullranked_folder = (output_folder / "phase_2" / "full_ranked_event_list")
    fullranked_folder.mkdir(parents=True, exist_ok=True)
    
    fullranked_csv_name = (
        f"Full_Event_List_for_{instrument}_"
        f"{filter_name}-band_"
        f"for_{run_mode}_{metric_mode}-modes_"
        f"{sky_noise_text}_{defocus_text}_"
        f"{viable_cumulative_cut*100}%_cut_"
        f"from_{obs_start.strftime('%b-%d-%Y')}_"
        f"to_{obs_end.strftime('%b-%d-%Y')}.csv"
    )
    full_df.to_csv(fullranked_folder / fullranked_csv_name, index=False)

    print("[P2_PostCleaner] Full list of ranked events created.")


    # ROUTE 2: Full list of cumulative observability scores (COS) -- cumu_df
    # Cumulative observability calculation, also carrying forward previous study flag.

    # flag_df is a df of planets and its study flag score
    flag_df = bigdf.filter(['Planet', 'Previous Study Flag'], axis=1)
    flag_df.drop_duplicates(subset='Planet', keep='first', inplace=True)
    flag_df.reset_index(inplace=True)    
    del flag_df['index']
    
    COS_data = []
    flag_df_idxs = np.arange(len(flag_df))

    for i in flag_df_idxs:
        planet_name = flag_df['Planet'][i]
        previous_study_flag = flag_df['Previous Study Flag'][i]

        bigdg = bigdf[bigdf['Planet'].str.contains(planet_name)]   # Recover all events for that planet specifically
        all_event_numbers = bigdg.shape[0]

        bigdh = bigdg[bigdg['Event_Weight'] >= 0.5]                # These are events for which at least all transit + one baseline is visible.
        good_event_numbers = bigdh.shape[0]

        CumObsScore = bigdh[f'Final_{metric_mode}'].sum()          # Sum scores of all the 'good' events from previous.

        COS_data.append([planet_name, previous_study_flag, CumObsScore, all_event_numbers, good_event_numbers])  # Append to a list.
    
    COS = pd.DataFrame(COS_data)  # Convert to DataFrame called Cumulative Obseravability Score.
    COS.columns = ['Planet', 'Previous Study Flag', f'CumObsScore_{instrument}', f'{instrument}_Event_No', f'{instrument}_Event_No (W>=0.5)']
    COS = COS[COS[f'{instrument}_Event_No (W>=0.5)'] != 0]                        # Drop all planets with no good events.
    COS.sort_values(f'CumObsScore_{instrument}', ascending=False, inplace=True)   # Sort and tidy!
    

    # Recover ranked table positions for cumulative score .csv.
    cumu_df = pd.merge(COS, ranked_df, how='inner', indicator=True)
   
    # Rearrange columns
    cumu_df = cumu_df[
        ['Planet', 'Previous Study Flag', f'CumObsScore_{instrument}', f'{instrument}_Event_No', f'{instrument}_Event_No (W>=0.5)', 'Table Position']
    ]
    
    cumulative_folder = (output_folder / "phase_2" / "cumulative_observability_scores")
    cumulative_folder.mkdir(parents=True, exist_ok=True)

    csv_config_suffix = fullranked_csv_name.removeprefix("Full_Event_List_for_")
    cumulative_csv_name = f'CumulativeObsScore_for_{csv_config_suffix}.csv'

    cumu_df.to_csv(cumulative_folder / cumulative_csv_name, index=False)
    
    print('[P2_PostCleaner] Cumulative observability scores calculated.')
