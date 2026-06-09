#Wrapper script for Part 7 of the pipeline, to set up multiprocessing jobs.
#The transit predictions + plots are made by P2MP_Process, which this script imports.

#Print is a keyword in python 2.X, but function-ness was added in 3.X.
#Import from future to access this!
#NOTE: Imports from future must always go first.
from __future__ import print_function
from astropy.utils import data

import os
import platform
import time
import glob
import urllib.error
import P2MP_Process

import numpy as np
import datetime as dt
import multiprocessing as mp

#Paths to .csv files
csvcorepath = '/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/Core_Files'
csvinpath = '/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_InputParts'
csvoutpath = '/mirror/scratch/jmorgan/Codes/Pipeline/CSV_Bank/P2_CSV_OutputParts'

#Our multiprocesses will live here.
jobs = []

#Print behaviour stored as a variable.
p = print

def get_print():
    return print

#Cross-platform check for file creation/modification date
def creation_date(path_to_file):
    if platform.system() == 'Windows':
        return os.path.getctime(path_to_file) #Windows is easy.
    else:
        stat = os.stat(path_to_file)
        try:
            return stat.st_birthtime
        except AttributeError:
            #The Linux kernel does not offer an easy way to get file creation dates.
            #Instead, go for when the file was last modified.
            #But~! If our .csv is never modified by us, this will just return the download date!
            return stat.st_mtime
        
#Sets up subprocess to queue planets for each core, otherwise they all fire at once! This kills the kernel. :-(
def _subprocess(ObsStart,ObsEnd,ds,S, Inst, Filter, Run_Mode, Add_Noise, Defocus,
                                       Metric_Mode, ViableCut, JobChunk, Core_No):    
    for Planet_No, csv_in in enumerate(JobChunk[Core_No]): 
        Dir, Name = csv_in.split('Parts/')
        P2MP_Process.P2Predictor(csv_in, Name, ObsStart, ObsEnd, ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus,
                                   Metric_Mode, ViableCut, JobChunk, Core_No, Planet_No)
            
#Queues and executes Phase Two transit predictor for all planets.
def P2Wrap(ds, S, Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut, ObsStart, ObsEnd, ncores):
    #P2 input finders (don't touch these!)
    stub = '%s_%s-band_for_%s,%s,%s,%s_%s_Cut_' % (Inst, Filter, Run_Mode, Add_Noise, Defocus, Metric_Mode, ViableCut)
    fnpattern = stub + '*.csv'
    inpattern_wpath = os.path.normpath(os.path.join(csvinpath, fnpattern))
    
    #Make directory for output graphs, if it doesn't exist already.
    if not os.path.exists('/mirror/scratch/jmorgan/Codes/Pipeline/P2_Outputs/%s_%s-band_%s_to_%s' % (Inst, Filter, ObsStart.strftime('%b %d %Y'), ObsEnd.strftime('%b %d %Y'))):
        os.mkdir('/mirror/scratch/jmorgan/Codes/Pipeline/P2_Outputs/%s_%s-band_%s_to_%s' % (Inst, Filter, ObsStart.strftime('%b %d %Y'), ObsEnd.strftime('%b %d %Y')))
    
    #Phase Two needs a cached IERS table to perform co-ordinate transformations, which auto-updates every 7 days.
    #A manual download on a shorter timescale is needed to prevent multiple checks/downloads when multiprocessing.
    #The hash used to denote the file also updates every 7 days, so this must also be tracked and updated in a .txt file.
    if data.is_url_in_cache('http://maia.usno.navy.mil/ser7/finals2000A.all'):
        data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True) #Get the name of the cached file.
        with open((r'%s/IERSPath.txt' % csvcorepath),'w' ) as IERSFile: #Write location and new hash to .txt file for future ref.
            c=IERSFile.write(str(data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True)))
    else:
        p('Hold up - no cached IERS table is present! Retrieving...')
        try: #Try/except blocks are needed in case of bad internet connection. #JustEduroamThings
            data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True) #Get the new file.
            with open((r'%s/IERSPath.txt' % csvcorepath),'w' ) as IERSFile:
                c=IERSFile.write(str(data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True)))
            p('Retrieval complete.')
        except urllib.error.URLError:
            p('The request timed out; trying again...')
            data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True) #Get the new file.
            with open((r'%s/IERSPath.txt' % csvcorepath),'w' ) as IERSFile:
                c=IERSFile.write(str(data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True)))
            p('Retrieval complete.')
            
    #Now that we definitely have a cached copy, we need to check its age.
    with open((r'%s/IERSPath.txt' % csvcorepath),'r' ) as IERSFile:
        c=IERSFile.read().replace('\n','') #Reads name of cached IERS table.
    A = creation_date(c) 
    LegA = dt.datetime.fromtimestamp(A) #Converts date to human-readable timestamp.
    p(LegA)
    
    if dt.datetime.now() - LegA >= dt.timedelta(1.): #If old, we update.
        p('Retrieving updated IERS table: please stand by...')
        try:
            data.clear_download_cache(hashorurl='http://maia.usno.navy.mil/ser7/finals2000A.all') #Otherwise the next line just returns the name and location of cached file!
            data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True) #Get the new file.
            with open((r'%s/IERSPath.txt' % csvcorepath),'w' ) as IERSFile: #Write location and new hash to .txt file for future ref.
                c=IERSFile.write(str(data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True)))
            p('Retrieval complete!')
        except urllib.error.URLError:
            p('The request timed out; trying again...')
            data.clear_download_cache(hashorurl='http://maia.usno.navy.mil/ser7/finals2000A.all') #Otherwise the next line just returns the name and location of cached file!
            data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True) #Get the new file.
            with open((r'%s/IERSPath.txt' % csvcorepath),'w' ) as IERSFile: #Write location and new hash to .txt file for future ref.
                c=IERSFile.write(str(data.download_file('http://maia.usno.navy.mil/ser7/finals2000A.all', cache=True, show_progress=True)))
            p('Retrieval complete!')
    else:
        p('IERS table is up to date.')
    
    #Test to try and make sure no further IERS downloads are triggered.
    time.sleep(3)
    
    JobChunk=(np.array_split(glob.glob(inpattern_wpath),ncores)) #Splits target list into even chunks, one for each core.
    p('+++Multiprocessing engaged for %s cores.+++' % ncores)
    
    #Targeting our subprocess, iterate over available cores.
    for Core_No in range(ncores):    
        process = mp.Process(target=_subprocess, args=(ObsStart,ObsEnd,ds,S, Inst, Filter, Run_Mode, Add_Noise, Defocus,
                                   Metric_Mode, ViableCut, JobChunk, Core_No))
        os.system("taskset -p -c %d %d" % (Core_No % ncores, os.getpid())) #Manually assigns subset of targets to each core. Linux-only?
        jobs.append(process) #Each planet becomes one child process, appended to jobs list.
        process.start() #Executes process. 
        
    for job in jobs:
        job.join() #Block subsequent code from running until all processes (planets) are completed.
        
    return get_print() #Return all printed statements as part of main line.