# Part 1 of overall pipeline.
# Checks modification date of three TEPCat .csvs, downloads new copies and builds if necessary.
import os
import platform
import urllib.request

import numpy as np
import pandas as pd
import datetime as dt

# urls for TEPCat files
target_url_a = 'http://www.astro.keele.ac.uk/jkt/tepcat/allplanets-csv.csv'
target_url_k = 'http://www.astro.keele.ac.uk/jkt/tepcat/kepplanets-csv.csv'
target_url_o = 'http://www.astro.keele.ac.uk/jkt/tepcat/observables.csv'

# Headers to be used in master csv. Star properties first, then planet.
planets_cols = ['Planet',
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
planets_cols_idx = np.arange(0, len(planets_cols), 1)     # Indices of planets_cols
planets_cols_idx_to_use = planets_cols_idx[np.arange(len(planets_cols_idx)) != 19]   # Excludes 'Period' header (Obs P (day) will be used)

obs_cols = ['Planet', 'Type',
    'RA:HH', 'RA:MM', 'RA:SS',
    'Dec:DD', 'Dec:MM', 'Dec:SS',
    'Vmag', 'Kmag',
    'T14', 'Depth',
    'T0 (HJD or BJD)', 'T0_err',
    'P (day)', 'P_err',
    'Ephemeris_Reference'
]
obs_cols_idx_to_use = np.arange(1, len(obs_cols), 1)      # Indices of obs_cols with 'Planet' excluded


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
# CSV_core_folder argument MUST be CSV_core_folder from the master shell.
def Check(CSV_core_folder):
    # Locate files
    allplanets_csv_path = CSV_core_folder / 'allplanets-csv.csv'
    kepplanets_csv_path = CSV_core_folder / 'kepplanets-csv.csv'
    observable_csv_path = CSV_core_folder / 'observables.csv'
    fulltepset_csv_path = CSV_core_folder / 'FullTEPSet.csv'

    # When were the files created/modified?
    creation_date_A = creation_date(allplanets_csv_path)
    creation_date_K = creation_date(kepplanets_csv_path)
    creation_date_O = creation_date(observable_csv_path)

    creation_date_A = dt.datetime.fromtimestamp(creation_date_A)  # Converts date to human-readable timestamp.
    creation_date_K = dt.datetime.fromtimestamp(creation_date_K)
    creation_date_O = dt.datetime.fromtimestamp(creation_date_O)
    
    # Retrieve new files if a week has elapsed between this code running and last download.
    # Flip all instances of >= to <= for bug-hunting.
    now = dt.datetime.now()
    week_elapsed_since_creation = dt.timedelta(7.)

    if now - creation_date_A >= week_elapsed_since_creation:
        print('[ModCheck] Retrieving updated allplanets-csv.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_a, allplanets_csv_path)
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] allplanets-csv.csv is up to date.')
        
    if now - creation_date_K >= week_elapsed_since_creation:
        print('[ModCheck] Retrieving updated kepplanets-csv.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_k, kepplanets_csv_path)    
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] kepplanets-csv.csv is up to date.')
        
    if now - creation_date_O >= week_elapsed_since_creation:
        print('[ModCheck] Retrieving updated observables.csv file: please stand by...')
        urllib.request.urlretrieve(target_url_o, observable_csv_path)
        print('>> Retrieval complete!')
    else:
        print('[ModCheck] observables.csv is up to date.')
    

    # DAN: Code should check if FullTEPSet.csv exists. <-- Already accounted for
    # If there are issues with creating the file, change the condition to "if True:" to guarantee file creation.
    # if True:
    if (now - creation_date_A >= week_elapsed_since_creation) or (not os.path.exists(fulltepset_csv_path)):
        # Reads in three TEPCat csvs. 
        # I reject the headers (with header=None) and substitute my own!
        # In planet .csvs, -1 is NaN -- pandas must be told this.
        # Some columns are read in as int64 -- they must be forced to be floats.
        make_float64 = ['Teff_err_up', 'Teff_err_dn',
                        'Teq', 'Teq_err_up', 'Teq_err_dn']

        # df dg dg1
        df_a = pd.read_csv(allplanets_csv_path,
                           skipinitialspace=True,
                           names=planets_cols, usecols=planets_cols_idx_to_use, header=None, skiprows=1,
                           dtype={key: np.float64 for key in make_float64},
                           na_values=['-1'])

        df_k = pd.read_csv(kepplanets_csv_path,
                           skipinitialspace=True,
                           names=planets_cols, usecols=planets_cols_idx_to_use, header=None, skiprows=0,
                           dtype={key: np.float64 for key in make_float64},
                           na_values=['-1'])
        
        df_o = pd.read_csv(observable_csv_path,
                           skipinitialspace=True,
                           names=obs_cols, usecols=obs_cols_idx_to_use, header=None, skiprows=1)

        # Cuts together, cleans and reads out master .csv file.
        # gp_err_dn contains a bad value (2..9), and so must be handled here rather than at read-in.
        # NaN k-mags must be specified here to preserve -1 co-ordinate values elsewhere.
        df_fulltep = pd.concat([df_a, df_k], ignore_index=True)
        df_fulltep = pd.concat([df_fulltep, df_o], axis=1)                   # Set axis=1 to combine column-wise.        

        # Hardcoded fixes (Due to unpredictability of error formatting)
        # If there are more hardcoding needed, put them here
        df_fulltep['gp_err_dn'] = df_fulltep['gp_err_dn'].apply(pd.to_numeric, errors='coerce')   # Makes 2..9 go away!
        df_fulltep['rhop_err_dn'] = df_fulltep['rhop_err_dn'].replace('0.3-', 0.3)                # Fix for TOI-561d

        # Discard unusable entries
        usable_entries_mask = (
            ((df_fulltep['T14'] != 0.0) | (df_fulltep['Depth'] != 0.0))     # Discard entries with zero T14 OR depth
            & ~df_fulltep['Type'].str.contains('BD')                        # Filter out brown dwarves
            & ~df_fulltep['Planet'].str.contains('EPIC_201702477')          # This object is also a brown dwarf
        )
        df_fulltep = df_fulltep[usable_entries_mask]      

        # Set NaN appropriately
        df_fulltep['gp_err_dn'] = df_fulltep['gp_err_dn'].replace(-1., np.nan)
        df_fulltep['Kmag'] = df_fulltep['Kmag'].replace(-1., np.nan)
        df_fulltep['Mp'] = df_fulltep['Mp'].replace(0, np.nan)          # Correctly sets missing masses to avoid dividing by 0

        # Whitespace elimination
        for colname in ['Planet', 'Discovery_Reference', 'Recent_Reference', 'Ephemeris_Reference']:
            df_fulltep[colname] = df_fulltep[colname].str.strip()

        # Export
        df_fulltep.to_csv(fulltepset_csv_path, index=False)     
        print('[ModCheck] New version of FullTEPSet.csv constructed.')
    
    else:
        print('[ModCheck] Current version of FullTEPSet.csv is up to date.')
