# Part 4.2 of overall pipeline.
# Takes in WorkingTEPSet, calculates exposure times, cadences and ranks for a chosen calibration telescope/filter combination.
# This then provides the cut-off to be subsequently used in the wider pipeline.
# This cut-off must be chosen by the user - at present, a physically-motivated cut is not possible.
# Outputs a small .csv of D_min values for wider use.

# This module should be used for spectroscopic observations - it doesn't make much sense to compare FORS2 spectroscopy to
# photometry from a small telescope!
import pandas as pd

# The central rank-making function.
def RankMaker(csvpath, csvrankpath, csvcutpath,
              ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, ViableCut):
    
    # Run RankMaker on Calibrating instrument
    import P1_RankMaker
    print("Running RankMaker on calibrating instrument...")
    P1_RankMaker.RankMaker(csvpath, csvrankpath,
                           ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus)

    di = pd.read_csv(rf'{csvrankpath}/RankedTEPSet_{Inst}_{Filter}-band_for_{Run_Mode}_{Add_Noise}_{Defocus}.csv',
                     skipinitialspace=True)
    
    # 'Index' column is only there to make the CSV look like the legacy ver as much as possible.
    # These are the list to store cut thresholds.
    RMin_cols = ['Index', 'RMin_Rank', 'RMin_Habitable_Rank',
                 'RMin_Multi_Transit_Rank', 'RMin_Multi_Transit_Habitable_Rank']
    RMin_vals = [1]

    def CumCutter(col):
        # Sort by chosen column (signal strength parameter), high to low.
        dj = di[col].sort_values(by=col, ascending=False)
        dj.reset_index(inplace=True)
        dj['CumSum'] = dj[col].cumsum() # Cumulative sum calculation
        dj['SumFrac'] = dj['CumSum'] / dj['CumSum'][-1]

        RMin = dj[dj['SumFrac'] <= ViableCut]
        RMin_vals.append(RMin)
    
    # For loop on all metrics
    for metric in ['Rank', 'Habitable_Rank', 'Multi_Transit_Rank', 'Multi_Transit_Habitable_Rank']:
        CumCutter(metric)

    # Brings together into DataFrame.
    df = pd.DataFrame([RMin_vals], columns=RMin_cols)

    # Write to CSV.
    df.to_csv(rf'{csvcutpath}/Rmin_{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{ViableCut*100}%_cut.csv',
              index=False)       
    print(rf'{ViableCut*100}% cuts for {Inst}, {Filter}-band for {Run_Mode}, {Add_Noise}, {Defocus}-defocus generated.')
