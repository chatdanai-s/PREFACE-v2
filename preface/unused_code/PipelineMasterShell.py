# Master script for pipeline for Python 3.11
import numpy as np
import pandas as pd
import datetime as dt
import joblib  # For parallel computing

# Please note the other nonstandard packages used in this pipeline include:
# ephem, joblib, tqdm, lime_tbx
# Refer to ./README_lime_tbx_installation_guide.txt on installing lime_tbx.

if __name__ == '__main__':  # Needed in case multiprocessing runs into an infinite loop
    ### Paths to .csv files
    csvcorepath = '../CSV_Bank/Core_Files'      # Used to collect, store, and process essential parameters
    csvrankpath = '../CSV_Bank/Ranked_Sets'     # Used to store ranked results
    csvcutpath = '../CSV_Bank/Place_Cutters'    # Used to store threashold-cut results
    csvinpath = '../CSV_Bank/P2_CSV_InputParts'        # Output of Phase 1; Input of Phase 2
    csvoutpath = '../CSV_Bank/P2_CSV_OutputParts'      # Output for individual planet CSVs
    csvendpath = '../CSV_Bank/Final_Outputs'           # Final output of the pipeline -- For individual transits
    csvscorepath = '../CSV_Bank/Cumulative_ObsScores'  # Final output of the pipeline -- For cumulative obseravability scores

    ### Choose your instrument, filter and desired running modes.
    Inst = 'TNT ULTRASPEC'  # A wrong/missing name here will give an IndexError.
    Filter = 'r'            # BVRI usually, ugriz for TNT, RISE/SPRAT_Red/SPRAT_Blue for LT. Check in Scope.csv
    Run_Mode = 'Half_Well'  # Half_Well, Spectral_Half_Well or IR_Half_Well (latter for VLT/KMOS observations)
    Add_Noise = 'Y_Noise'   # Toggles additional noise terms (sky backgrounds and out-of-transit); Accepts Y_Noise or N_Noise
    Defocus = 'Y'           # Y or N to toggle defocus. 
    Metric_Mode = 'Rank'    # Rank, Habitable_Rank, Multi_Transit_Rank, Multi_Transit_Habitable_Rank
    ViableCut = 0.97        # Take this top fraction of the cumulative distribution of planet metrics for the calibration instrument.

    ### Choose your configurations for moon background considerations. (Only supports UBVRI and ugriz bands)
    # 'Default' argument is available for Scattering_AOD, Extinction_AOD, and Asymmetry_Factor, which picks the
    # median aggregated by month of year at closest AERONET site, interpolated at wavelength.
    Moon_Noise = 'Y_Moon'   # Toggles moon background considerations; Accepts Y_Moon or N_Moon.
    Scattering_AOD = 'Default'      # Aerosol optical depth for scattering. Set as 0 if no aerosol scattering.
    Absorption_AOD = 'Default'      # Aerosol optical depth for absorption. Set as 0 if no aerosol absorption.
    Asymmetry_Factor = 'Default'    # Usually 0.5 - 0.8 for atmospheric aerosols.
    Amplification_Factor = 5        # The higher, the lower the score for brighter Moons. 5 is default. 0 implies no amplification.

    ### For Phase Two image outputs, what is the minimum EventWeight required for image generation?
    # 0 for all, 0.5 recommended, 1 for only complete transit signals, >1 if none
    EventWeight_GraphThreshold = 0.5

    ### Specify your desired time (in BJD) of observation! Can be future or past.
    # Format: Y,M,D,H,M,s
    ObsStart = dt.datetime(2025,10,1,0,0,0)
    ObsEnd = dt.datetime(2026,5,31,0,0,0)

    ### How many cores do you wish to call on?
    # PROTIP: Don't use them all, or your machine will likely lock up. Leave one spare!
    ncores = joblib.cpu_count() - 2
 

    ### Check input validity
    import InputCheck
    InputCheck.InputCheck(csvcorepath, csvrankpath, csvcutpath, csvinpath, csvoutpath, csvendpath, csvscorepath,
                          Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
                          Moon_Noise, Scattering_AOD, Absorption_AOD, Asymmetry_Factor, Amplification_Factor,
                          EventWeight_GraphThreshold, ObsStart, ObsEnd, ncores)
    
    # Retrive scope.csv and index for telescope parameters
    ds = pd.read_csv(rf'{csvcorepath}/Scope.csv')
    S = np.where(ds['Telescope'] == Inst)[0][0]

    ### Fire Phase One of the pipeline.
    ### Call each script as a module, then fire the main function associated with that script.
    import ModCheck
    ModCheck.Check(csvcorepath)
    import ImpactMerger
    ImpactMerger.ExoOrgImpacts(csvcorepath)
    import ExoplanetseuImpactMerger
    ExoplanetseuImpactMerger.ExoeuImpacts(csvcorepath)
    import P1_WorkingTEPSetBuilder
    P1_WorkingTEPSetBuilder.WorkBuilder(csvcorepath)

    import P1_RankMaker
    P1_RankMaker.RankMaker(csvcorepath, csvrankpath,
                           ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus)

    import P1_Cutter   # Extracts cut-off metric value RMin.
    RMin = P1_Cutter.RankMaker(csvcorepath, csvrankpath, csvcutpath,
                               ds, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut)

    import P1ViabilitySplitter
    P1ViabilitySplitter.Splitter(csvrankpath, csvinpath, csvoutpath,
                                 Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, RMin)


    ### Fire Phase Two of the pipeline.
    import P2MP_Wrapper_Datacenter
    P2MP_Wrapper_Datacenter.P2Wrap(csvcorepath, csvinpath, csvoutpath,
                                   ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
                                   Moon_Noise, Scattering_AOD, Absorption_AOD, Asymmetry_Factor, Amplification_Factor,
                                   EventWeight_GraphThreshold, ObsStart, ObsEnd, ncores)

    import P2PostCleaner
    P2PostCleaner.Cleaner(csvrankpath, csvoutpath, csvendpath, csvscorepath,
                          ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut,
                          ObsStart, ObsEnd, ncores)
