# Part 4.1 of overall pipeline.
# Takes in WorkingTEPSet, calculates exposure times, cadences and ranks for a chosen calibration telescope/filter combination.
# This then provides the cut-off to be subsequently used in the wider pipeline.
# This cut-off must be chosen by the user - at present, a physically-motivated cut is not possible.
# Outputs a small .csv of D_min values for wider use.
import pandas as pd
import numpy as np

#The central rank-making function.
def RankMaker(CSV_core_folder, CSV_intermediate_folder,
              scope_df, instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode, viable_cumulative_cut):

    # Parameters for RankMaker for calibrating instruments
    if run_mode == 'Spectral_Half_Well':
        Inst_cal = 'VLT FORS2 (200kHz) 600RI+19'  # Calibrating instrument
        Filter_cal = '600RI+19'
        SkyNoise_cal = toggle_sky_noise
        Defocus_cal = toggle_defocus

    elif run_mode == 'IR_Half_Well':
        Inst_cal = 'VLT KMOS_HK'
        Filter_cal = 'HK'
        SkyNoise_cal = True
        Defocus_cal = toggle_defocus

    else:  # Half_Well Run_Mode
        Inst_cal = instrument
        Filter_cal = filter_name
        SkyNoise_cal = toggle_sky_noise
        Defocus_cal = toggle_defocus

    cal_idx = np.where(scope_df['Telescope'] == Inst_cal)[0][0]  # Index of calibrating Inst

    sky_noise_text = 'Y-SkyNoise' if SkyNoise_cal else 'N-SkyNoise'
    defocus_text = 'Y-Defocus' if Defocus_cal else 'N-Defocus'
    calibrator_csv_path = (
        CSV_intermediate_folder / "ranked_tep_sets"
        / f"RankedTEPSet_{Inst_cal}_{Filter_cal}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}.csv"
    )
    rankcut_csv_path = (
        CSV_intermediate_folder / "minranks_cut_sets"
        / f"minRanks_{Inst_cal}_{Filter_cal}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}_{viable_cumulative_cut*100}%_cut.csv"
    )
    (CSV_intermediate_folder / "minranks_cut_sets").mkdir(parents=True, exist_ok=True)

    # Run RankMaker on Calibrating instrument if configuration is different
    different_calibrating_configuration = (Inst_cal != instrument) or (Filter_cal != filter_name) \
                                       or (SkyNoise_cal != toggle_sky_noise) or (Defocus_cal != toggle_defocus)
    
    if different_calibrating_configuration == True:
        import preface.P1_RankMaker as P1_RankMaker
        print("[Cutter] Running RankMaker on calibrating instrument...")
        P1_RankMaker.RankMaker(CSV_core_folder, CSV_intermediate_folder,
                               scope_df, cal_idx, Inst_cal, Filter_cal, run_mode, SkyNoise_cal, Defocus_cal)

    df_cal = pd.read_csv(calibrator_csv_path, skipinitialspace=True)


    # 'Index' column is only there to make the CSV look like the legacy ver as much as possible.
    # This is the list to store cut thresholds.
    RMin_cols = ['Index', 'RMin_Rank', 'RMin_Habitable_Rank',
                 'RMin_Multi_Transit_Rank', 'RMin_Multi_Transit_Habitable_Rank']
    RMin_vals = [1]

    def CumCutter(metric_col):
        # Sort by chosen column (signal strength parameter), high to low.
        df_cal_temp = df_cal.sort_values(by=metric_col, ascending=False)[[metric_col]]  # [[col]] keeps dj a dataframe
        df_cal_temp.reset_index(inplace=True)
        df_cal_temp['CumSum'] = df_cal_temp[metric_col].cumsum() # Cumulative sum calculation
        df_cal_temp['CumSumFrac'] = df_cal_temp['CumSum'] / df_cal_temp['CumSum'].max()

        df_cal_temp = df_cal_temp[df_cal_temp['CumSumFrac'] <= viable_cumulative_cut]
        RMin = df_cal_temp[metric_col].min()
        RMin_vals.append(RMin)
    
    # For loop on all metrics
    for metric in ['Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']:
        CumCutter(metric)

    
    # Brings together into DataFrame and provides cutoff value
    df = pd.DataFrame([RMin_vals], columns=RMin_cols)
    Rmin = df[f'RMin_{metric_mode}'][0]

    # Write to CSV.
    df.to_csv(rankcut_csv_path, index=False)     
    print(rf'[Cutter] {viable_cumulative_cut*100}% cuts for {Inst_cal}, {Filter_cal}-band for {run_mode}, {sky_noise_text}, {defocus_text} generated.')

    return Rmin
