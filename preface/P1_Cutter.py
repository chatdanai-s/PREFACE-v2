# Part 4.1 of overall pipeline.
# Takes in WorkingTEPSet, calculates exposure times, cadences and ranks for a chosen calibration telescope/filter combination.
# This then provides the cut-off to be subsequently used in the wider pipeline.
# This cut-off must be chosen by the user - at present, a physically-motivated cut is not possible.
# Outputs a small .csv of D_min values for wider use.
import pandas as pd
import numpy as np

#The central rank-making function.
def RankMaker(csvcorepath, csvrankpath, csvcutpath,
              ds, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut):
    
    # Parameters for RankMaker for calibrating instruments
    ds = pd.read_csv(rf'{csvcorepath}/Scope.csv') 

    if Run_Mode == 'Spectral_Half_Well':
        Inst_cal = 'VLT FORS2 (200kHz) 600RI+19'  # Calibrating instrument
        Filter_cal = '600RI+19'
        Add_Noise_cal = Add_Noise
        Defocus_cal = Defocus

    elif Run_Mode == 'IR_Half_Well':
        Inst_cal = 'VLT KMOS_HK'
        Filter_cal = 'HK'
        Add_Noise_cal = 'Y_Noise'
        Defocus_cal = Defocus

    else:  # Half_Well Run_Mode
        Inst_cal = Inst
        Filter_cal = Filter
        Add_Noise_cal = Add_Noise
        Defocus_cal = Defocus

    S_cal = np.where(ds['Telescope'] == Inst_cal)[0][0]  # Index of calibrating Inst
    
    # Run RankMaker on Calibrating instrument if configuration is different
    if (Inst_cal != Inst) or (Filter_cal != Filter) or (Add_Noise_cal != Add_Noise) or (Defocus_cal != Defocus):
        import P1_RankMaker
        print("[Cutter] Running RankMaker on calibrating instrument...")
        P1_RankMaker.RankMaker(csvcorepath, csvrankpath, ds,
                               S_cal, Inst_cal, Filter_cal, Run_Mode, Add_Noise_cal, Defocus_cal)

    di = pd.read_csv(rf'{csvrankpath}/RankedTEPSet_{Inst_cal}_{Filter_cal}-band_for_{Run_Mode}_{Add_Noise_cal}_{Defocus_cal}.csv',
                     skipinitialspace=True)

    # 'Index' column is only there to make the CSV look like the legacy ver as much as possible.
    # This is the list to store cut thresholds.
    RMin_cols = ['Index', 'RMin_Rank', 'RMin_Habitable_Rank',
                 'RMin_Multi_Transit_Rank', 'RMin_Multi_Transit_Habitable_Rank']
    RMin_vals = [1]

    def CumCutter(col):
        # Sort by chosen column (signal strength parameter), high to low.
        dj = di.sort_values(by=col, ascending=False)[[col]]  # [[col]] keeps dj a dataframe
        dj.reset_index(inplace=True)
        dj['CumSum'] = dj[col].cumsum() # Cumulative sum calculation
        dj['SumFrac'] = dj['CumSum'] / dj['CumSum'].max()

        dj = dj[dj['SumFrac'] <= ViableCut]
        RMin = dj[col].min()
        RMin_vals.append(RMin)
    
    # For loop on all metrics
    for metric in ['Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']:
        CumCutter(metric)

    # Brings together into DataFrame and provides cutoff value
    df = pd.DataFrame([RMin_vals], columns=RMin_cols)
    RMin = df[f'RMin_{Metric_Mode}'][0]

    # Write to CSV.
    df.to_csv(rf'{csvcutpath}/Rmin_{Inst_cal}_{Filter_cal}-band_for_{Run_Mode},{Add_Noise_cal},{Defocus_cal},{ViableCut*100}%_cut.csv',
              index=False)       
    print(rf'[Cutter] {ViableCut*100}% cuts for {Inst_cal}, {Filter_cal}-band for {Run_Mode}, {Add_Noise_cal}, {Defocus_cal}-defocus generated.')

    return RMin
