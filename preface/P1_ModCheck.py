# Part 1 of overall pipeline.
# Checks modification date of three TEPCat .csvs, downloads new copies and builds if necessary.
import os
import platform
import urllib.request

import numpy as np
import pandas as pd
import datetime as dt


# urls for TEPCat files
target_url_1 = 'http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'
target_url_2 = 'http://www.astro.keele.ac.uk/jkt/tepcat/kepplanets-csv.csv'
target_url_3 = 'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv'

# Headers to be used in master csv. Star properties first, then planet.
Head = ['Planet',
        'Teff', 'Teff_err_up', 'Teff_err_dn',
        '[Fe/H]', '[Fe/H]_err_up', '[Fe/H]_err_dn',
        'M*', 'M*_err_up', 'M*_err_dn',
        'R*', 'R*_err_up', 'R*_err_dn',
        'logg*', 'logg*_err_up', 'logg*_err_dn',
        'rho*', 'rho*_err_up', 'rho*_err_dn',
        'Period',
        'e', 'e_err_up', 'e_err_dn',
        'a (AU)', 'a_err_up', 'a_err_dn',
        'Mp', 'Mp_err_up', 'Mp_err_dn',
        'Rp','Rp_err_up', 'Rp_err_dn',
        'gp', 'gp_err_up', 'gp_err_dn',
        'rhop', 'rhop_err_up', 'rhop_err_dn',
        'Teq', 'Teq_err_up', 'Teq_err_dn',
        'Discovery_Reference',
        'Recent_Reference'
        ]
H1 = np.arange(0, len(Head), 1)     # Indices of Head
H2 = H1[np.arange(len(H1)) != 19]   # H1, but excludes 'Period' header (Obs P (day) will be used)

Obs = ['Planet', 'Type',
       'RA:HH', 'RA:MM', 'RA:SS',
       'Dec:DD', 'Dec:MM', 'Dec:SS',
       'Vmag', 'Kmag',
       'T14', 'Depth',
       'T0 (HJD or BJD)', 'T0_err',
       'P (day)', 'P_err',
       'Ephemeris_Reference'
       ]
O1 = np.arange(1, len(Obs), 1)      # Indices of Obs with 'Planet' excluded


# Cross-platform check for file creation/modification date, shamelessly cribbed from StackOverflow.
# Code is modified to get earliest of file creation and modification date instead,
# because first run can cause issues with file creation metadata
def creation_date(path_to_file):   
    if platform.system() == 'Windows':  # For Windows
        last_time = os.path.getmtime(path_to_file)
        return last_time
    else:
        stat = os.stat(path_to_file)
        try:
            last_time = min(stat.st_birthtime, stat.st_mtime)  # For macOS
            return last_time 
        except AttributeError:
            # The Linux kernel does not offer an easy way to get file creation dates.
            # Instead, go for when the file was last modified.
            # But~! If our .csv is never modified by us, this will just return the download date!
            return stat.st_mtime


# Keep the whole thing in a function so the parent script can call it.
# csvpath argument MUST be csvcorepath from the master shell.
def Check(csvpath):
    # When were the files created/modified?
    A = creation_date(rf'{csvpath}/allplanets-csv.csv')
    K = creation_date(rf'{csvpath}/kepplanets-csv.csv')
    O = creation_date(rf'{csvpath}/observables.csv')
    LegA = dt.datetime.fromtimestamp(A)  # Converts date to human-readable timestamp.
    LegK = dt.datetime.fromtimestamp(K)
    LegO = dt.datetime.fromtimestamp(O)
    
    # Retrieve new files if a week has elapsed between this code running and last download.
    # Flip all instances of >= to <= for bug-hunting.
    now = dt.datetime.now()

    if now - LegA >= dt.timedelta(7.):
        print('[ModCheck] Retrieving updated allplanets-csv.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_1, rf'{csvpath}/allplanets-csv.csv')
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] allplanets-csv.csv is up to date.')
        
    if now - LegK >= dt.timedelta(7.):
        print('[ModCheck] Retrieving updated kepplanets-csv.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_2, rf'{csvpath}/kepplanets-csv.csv')    
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] kepplanets-csv.csv is up to date.')
        
    if now - LegO >= dt.timedelta(7.):
        print('[ModCheck] Retrieving updated observables.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_3, rf'{csvpath}/observables.csv')
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] observables.csv is up to date.')
    

    # DAN: Code should check if FullTEPSet.csv exists. <-- Already accounted for
    # If there are issues with creating the file, change the condition to "if True:" to guarantee file creation.
    # if True:
    if (now - LegA >= dt.timedelta(7.)) or (not os.path.exists(rf'{csvpath}/FullTEPSet.csv')):
        # Reads in three TEPCat csvs. 
        # I reject the headers (with header=None) and substitute my own!
        # In planet .csvs, -1 is NaN -- pandas must be told this.        
        # Some columns are read in as int64 -- they must be forced to be floats.
        make_float64 = ['Teff_err_up', 'Teff_err_dn',
                         'Teq', 'Teq_err_up', 'Teq_err_dn']

        df = pd.read_csv(rf'{csvpath}/allplanets-csv.csv',
                         skipinitialspace=True,
                         names=Head, usecols=H2, header=None, skiprows=1,
                         dtype={key: np.float64 for key in make_float64},
                         na_values=['-1'])

        dg = pd.read_csv(rf'{csvpath}/kepplanets-csv.csv',
                         skipinitialspace=True,
                         names=Head, usecols=H2, header=None,
                         dtype={key: np.float64 for key in make_float64},
                         na_values=['-1'])
        
        dg1 = pd.read_csv(rf'{csvpath}/observables.csv',
                          skipinitialspace=True,
                          names=Obs, usecols=O1, header=None, skiprows=1)

        # Cuts together, cleans and reads out master .csv file.
        # gp_err_dn contains a bad value (2..9), and so must be handled here rather than at read-in.
        # NaN k-mags must be specified here to preserve -1 co-ordinate values elsewhere.
        dh = pd.concat([df,dg], ignore_index=True)
        dh1 = pd.concat([dh,dg1], axis=1)                   # Set axis=1 to combine column-wise.        

        # Hardcoded fixes (Due to unpredictability of error formatting)
        # If there are more hardcoding needed, put them here
        dh1['gp_err_dn'] = dh1['gp_err_dn'].apply(pd.to_numeric, errors='coerce')   # Makes 2..9 go away!
        dh1['rhop_err_dn'] = dh1['rhop_err_dn'].replace('0.3-', 0.3)                # Fix for TOI-561d

        # Discard unusable entries
        dh1 = dh1[(dh1['T14'] != 0.) | (dh1['Depth'] != 0.)]      # Discard entries with zero T14 OR depth
        dh1 = dh1[~dh1['Type'].str.contains('BD')]                # Filter out brown dwarves
        dh1 = dh1[~dh1['Planet'].str.contains('EPIC_201702477')]  # This object is also a brown dwarf     

        # Set NaN appropriately
        dh1['gp_err_dn'] = dh1['gp_err_dn'].replace(-1., np.nan)
        dh1['Kmag'] = dh1['Kmag'].replace(-1., np.nan)
        dh1['Mp'] = dh1['Mp'].replace(0, np.nan)        # Correctly sets missing masses to avoid dividing by 0

        # Whitespace elimination
        for colname in ['Planet', 'Discovery_Reference', 'Recent_Reference', 'Ephemeris_Reference']:
            dh1[colname] = dh1[colname].str.strip()

        # Export
        dh1.to_csv(rf'{csvpath}/FullTEPSet.csv', index=False)     
        print('[ModCheck] New version of FullTEPSet.csv constructed.')
    
    else:
        print('[ModCheck] Current version of FullTEPSet.csv is up to date.')
