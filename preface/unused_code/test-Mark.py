import numpy as np
import pandas as pd
import datetime as dt
#import multiprocessing as mp
import multiprocess as mp

import ModCheck
import ImpactMerger
import ExoplanetseuImpactMerger
import P1_WorkingTEPSetBuilder
import P1_RankMaker
import P1_SpecCutter
import P1_Cutter
import P1ViabilitySplitter
import P2MP_Wrapper_Datacenter_NoMP #This one explicitly points to the mirror!
import P2PostCleaner

csvcorepath='../CSV_Bank/Core_Files'
csvcutpath = '../CSV_Bank/Place_Cutters'#'/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/Place_Cutters'


csvpath='../CSV_Bank/Core_Files' #'/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/Core_Files'
#Choose your instrument, filter and desired running modes.
Inst='Sierra Remote' #A wrong/missing name here will give an IndexError.
Filter='V' #BVRI usually, ugriz for TNT, RISE/SPRAT_Red/SPRAT_Blue for LT
Run_Mode='Half_Well' #Half_Well, Spectral_Half_Well or IR_Half_Well (latter for VLT/KMOS observations)
Add_Noise='Y_Noise' #Toggles additional noise terms (sky backgrounds and out-of-transit)
Defocus='N' #Y or N to toggle defocus. 
Metric_Mode='Rank' #Rank, Habitable_Rank, Multi_Transit_Rank, 'Multi_Transit_Habitable_Rank'
ViableCut=0.985 #Take this top fraction of the cumulative distribution of planet metrics for the calibration instrument.

ObsStart=dt.datetime(2022,12,1,0,0,0) #Specify your desired time of observation! Can be future or past.
ObsEnd=dt.datetime(2022,12,10,0,0,0) #Format Y,M,D,H,M,s

#Read in telescope.csv and selects index of chosen instrument. 
ds=pd.read_csv(r'%s/Scope.csv' % csvcorepath)
S = np.where(ds['Telescope'] == Inst)[0]

ncores=(mp.cpu_count()-1)



#ModCheck.Check() 
#ModCheck is Not needed in this version

ImpactMerger.ExoOrgImpacts()
ExoplanetseuImpactMerger.ExoeuImpacts()
P1_WorkingTEPSetBuilder.WorkBuilder()
P1_RankMaker.RankMaker(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus)

if Run_Mode == 'Spectral_Half_Well':
    P1_SpecCutter.RankMaker(ds, np.where(ds['Telescope'] == 'VLT FORS2 (200kHz) 600RI+19')[0], 'VLT FORS2 (200kHz) 600RI+19', 
      '600RI+19', 'Spectral_Half_Well', Add_Noise, Defocus, ViableCut)
    dd=pd.read_csv((r'%s/Rmin_%s_%s-band_for_%s,%s,%s,%s%%_cut.csv') % (csvcutpath, 'VLT FORS2 (200kHz) 600RI+19', '600RI+19', 
    'Spectral_Half_Well', Add_Noise, Defocus, ViableCut))
    RMin=dd['RMin_%s' % Metric_Mode][0]
else:
    P1_Cutter.RankMaker(ds, np.where(ds['Telescope'] == 'MuSCAT2_Dark')[0], 'MuSCAT2_Dark', 
      'r', 'Half_Well', Add_Noise, 'Y', ViableCut)
    dd=pd.read_csv((r'%s/Rmin_%s_%s-band_for_%s,%s,%s,%s%%_cut.csv') % (csvcutpath, 'MuSCAT2_Dark', 'r', 
      'Half_Well', Add_Noise, 'Y', ViableCut))
    RMin=dd['RMin_%s' % Metric_Mode][0]

P1ViabilitySplitter.Splitter(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, RMin)
#'''


#----- PHASE 2 -----
# P2MP_Wrapper_Datacenter takes around 90 minutes to run.

if __name__ == "__main__":
	P2MP_Wrapper_Datacenter_NoMP.P2Wrap(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, ObsStart, ObsEnd, ncores)


P2PostCleaner.Cleaner(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, ObsStart, ObsEnd)
