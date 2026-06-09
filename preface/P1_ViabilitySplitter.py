# Part 5 of overall pipeline.
# Small module to take the top planets of our ranked set and split into separate files.
# Two files made - viable planets, and not-viable planets.
import os
import fnmatch
import pandas as pd

csvtoppath = '../CSV_Bank/Top_Sets'
csvlowpath = '../CSV_Bank/Unusable_Sets'


# Main function
def Splitter(csvrankpath, csvinpath, csvoutpath,
             Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, RMin):
    df = pd.read_csv(rf'{csvrankpath}/RankedTEPSet_{Inst}_{Filter}-band_for_{Run_Mode}_{Add_Noise}_{Defocus}.csv')
    df = df.sort_values(by=Metric_Mode, ascending=False)
    df['CumSum'] = df[Metric_Mode].cumsum()
    df['SumFrac'] = df['CumSum'] / df['CumSum'].max()

    df_top = df[df[Metric_Mode] >= RMin]    # Takes all planets ranked above the cumulative cut.
    df_low = df[df[Metric_Mode] < RMin]     # Takes all planets ranked below the cumulative cut.

    print("[ViabilitySplitter] Outputting rank threshold-separated TEPSets to .csv now..." )
    df_top.to_csv(rf'{csvtoppath}/TopTEPSet_{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{Metric_Mode}_Mode,{ViableCut}_Cut.csv',
                  index=False)
    df_low.to_csv(rf'{csvlowpath}/UnusableTargets_{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{Metric_Mode}_Mode,{ViableCut}_Cut.csv',
                  index=False)
    
    print(rf'[ViabilitySplitter] New viability division ({ViableCut} fraction cut) for {Metric_Mode}-mode successful.')
        
    # Split top set into CSVs for individual targets, to be used by Phase Two.
    
    # Clear the Input and Output Parts folders if a previous run for this instrument/filter combination has been executed.
    # This should prevent any conflict as TEPCat is periodically updated and planet ranks are revised.
    def clear_configparts(configpath):
        configfiles = [os.path.join(path, f)
                       for path, dirs, files in os.walk(configpath)
                       for f in fnmatch.filter(files,
                                               f'{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{Metric_Mode}_{ViableCut}_Cut_*.csv')
        ]
        for f in configfiles:
            os.remove(f)

    clear_configparts(csvinpath)
    clear_configparts(csvoutpath)
        
    print('[ViabilitySplitter] Previous input/output wipe successful.')
    
    # Now we can safely write our new parts!
    df_top.reset_index()
    rows = df_top.shape[0]
    for row in range(rows):
        df.iloc[[row]].to_csv(rf'{csvinpath}/{Inst}_{Filter}-band_for_{Run_Mode},{Add_Noise},{Defocus},{Metric_Mode}_{ViableCut}_Cut_{row}.csv',
                              index=False)
        
    print('Phase One Complete.\n')
