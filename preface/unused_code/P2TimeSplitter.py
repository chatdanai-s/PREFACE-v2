# Part 6 of overall pipeline.
# Takes in subset of viable targets, recovers T0 in JD for each.
# Now dependent on instrument, filter and observing location.
# Then splits into many .csvs, one for each viable planet, to be handled by P2.
import os
import platform
#import glob
import fnmatch

import numpy as np
import pandas as pd
import datetime as dt
import astropy.coordinates as ac
import astropy.units as u

from functools import wraps
from astropy import time
from bisect import bisect_left

csvinpath='/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_InputParts'
csvoutpath='/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_OutputParts'
csvtoppath='/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/Top_Sets'

BJD_List=[]

#Print behaviour stored as a variable.
p = print

def get_print():
    return print

#Cross-platform check for file creation/modification date
def creation_date(path_to_file):
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file) 
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            return stat.st_mtime
        
#Counting function that executes every time the BJD_JD converter is called.
def counter(func):
    @wraps(func)
    def wrapper(row,ds,S):
        wrapper.count += 1    # executed every time the wrapped function is called
        return func(row,ds,S)
    wrapper.count = 0         # executed only once in decorator definition time
    return wrapper

#Which generated BJD is closest to the true value?
def takeClosest(HJD_List,q):
    
    pos=bisect_left(HJD_List,q)
    if pos == 0:
        return HJD_List[0]
    if pos == len(HJD_List):
        return HJD_List[-1]
    before = HJD_List[pos-1]
    after = HJD_List[pos]
    if after - q < q - before:
        return after
    else:
        return before
    
def TimeSplit(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut):
    
    U = creation_date(r'%s/TopTEPSet_%s_%s-band_for_%s,%s,%s,%s_Mode,%s_Cut.csv' % 
                      (csvtoppath, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut))
    LegU = dt.datetime.fromtimestamp(U) #Converts date to human-readable timestamp.    
    #p(LegU)
    
    #Set to fire when TopTEPSet was recently modifed, so will always fire when downloaded .csvs are updated.

    if dt.datetime.now() - LegU <= dt.timedelta(0, 10.):
    
        db=pd.read_csv(r'%s/TopTEPSet_%s_%s-band_for_%s,%s,%s,%s_Mode,%s_Cut.csv' % 
                       (csvtoppath, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut))
        
        #Clear the Input and Output Parts folders if a previous run for this instrument/filter combination has been executed.
        #This should prevent any conflict as TEPCat is periodically updated and planet ranks are revised.
        configfilesin= [os.path.join(path, f) for path, csvnames, files in os.walk(csvinpath) 
            for f in fnmatch.filter(files, r'%s_%s-band_for_%s,%s,%s,%s_%s_Cut_*.csv' % 
                                    (Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut))]
        for q in configfilesin:    
            os.remove(q)

        configfilesout= [os.path.join(path, f) for path, csvnames, files in os.walk(csvoutpath) 
            for f in fnmatch.filter(files, r'%s_%s-band_for_%s,%s,%s,%s_%s_Cut_*.csv' % 
                                    (Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut))]
        for q in configfilesout:    
            os.remove(q)
            
        p('Previous input/output wipe successful.')
        
        @counter #Count for this function!
        def BJD_JD_Machine(row,ds,S):
               
            #Specify co-ordinates of target.
            HD = ac.SkyCoord((row['RA:HH']+(row['RA:MM']/60)+(row['RA:SS']/3600)), (row['Dec:DD']+(row['Dec:MM']/60)+(row['Dec:SS']/3600)),
                             unit=(u.hourangle, u.deg), frame='icrs')
            T_Targ=row['T0 (HJD or BJD)'] #Target BJD from .csv.
            JD_Diff=0.006944 #10 minutes in Julian time scale - the JD will be somewhere in this range.
            Loc=ac.EarthLocation(lat=(ds['Lat'][S])*u.deg, lon=(ds['Long'][S])*u.deg, height=(ds['Alt'][S])*u.m) #Location info
            
            #Create Julian time range to look over.
            T=np.linspace((T_Targ-JD_Diff),(T_Targ+JD_Diff),1200.)
            #print T
            
            #Iterates over Julian dates to get set of test BJDs.
            for i in T:
                times=time.Time([i], format='jd',scale='tdb', location=Loc)
                ltt_bary = times.light_travel_time(HD)
                BJD_Entry = time.Time([i], format='jd', scale='tdb')+ltt_bary
                BJD_List.append(BJD_Entry)
            
            #Returns closest BJD value and its index within the list.    
            #print 'Closest BJD is %s' % takeClosest(BJD_List, time.Time(T_Targ, format='jd', scale='tdb'))
            Y=BJD_List.index(takeClosest(BJD_List, time.Time(T_Targ, format='jd', scale='tdb')))
            #print S
            JD_True=time.Time(T[Y], format='jd', scale='tdb') #Finds corresponding JD value.
            p(r'(%s - %s of %s)' % (JD_True.utc, BJD_JD_Machine.count, db.shape[0])) 
            #Note the conversion to UTC - this is what the astroutils online calc does.
            
            del BJD_List[:]
            
            return JD_True
        
        p('Recovering JD ephemerides for all viable targets - this may take a while!')
        db['T0 (JD)']=db.apply(lambda row: BJD_JD_Machine(row,ds,S), axis=1)
        p('All ephemerides recovered!')
        
        for y in range(db.shape[0]):
            db[db.index.isin([y])].to_csv((r'%s/%s_%s-band_for_%s,%s,%s,%s_%s_Cut_%s.csv') % 
              (csvinpath, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, y))
            
    else:
        p('We have times.')
            
    return get_print()