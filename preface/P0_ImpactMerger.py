# Part 2.1 of overall pipeline.
# TEPCat does not store impact parameters, so they must be corralled from elsewhere.

import pandas as pd

# Small function to massage names such that both catalogues use same naming convention.
def BWriter(row):
    if row['Planet'][-1].isdigit() == True:
        row['Planet'] = row['Planet'] + 'b'
        
    return row['Planet']

# Main function
# csvpath argument MUST be csvcorepath from the master shell.
def ExoOrgImpacts(csvpath):
    # Takes in all names in FullTEPSet and processes them.
    da = pd.read_csv(rf'{csvpath}/FullTEPSet.csv')
    da.insert(1, 'Planet_temp', da['Planet'])   # Temporary column

    da['Planet'] = da['Planet'].replace({'-00':'-', '-0':'-',
                                         '_00':'_', '_0':'_'},
                                        regex=True
                                        )
    da['Planet'] = da.apply(lambda row: BWriter(row), axis=1)

    # Takes in table of transiting planets from exoplanets.org
    # exoplanets.org is a catalogue built by hand, so no auto-updater here.
    # Also, exoplanets.org is not updated with planets published after June 2018
    colnames = ['NAME', 'B']
    db = pd.read_csv(rf'{csvpath}/ExoOrgList.csv', usecols=colnames)
    db.rename(columns={'NAME': 'Planet', 'B': 'Impact Parameter'}, inplace=True)

    # Process planet names before merging
    db['Planet'] = db['Planet'].replace({'55 Cnc e':'55_Cnc_e', 'OGLE2':'OGLE',
                                         'BD +20 594': 'K2-56b', 
                                         r'([A-Z]) ': r'\1_',
                                         'A b':'b',
                                         ' ':''
                                         },
                                        regex=True
                                        )
    
    # Append impact parameter to original dataframe
    dc = pd.merge(da, db, on='Planet', how='left')
    dc = dc.drop('Planet', axis=1)                      # Recover old planet names
    dc = dc.rename(columns={'Planet_temp': 'Planet'})  

    # Export
    dc.to_csv(rf'{csvpath}/FullTEPSetWithExoOrgImpacts.csv', index=False)

    print('[ImpactMerger] Impact parameters from exoplanets.org recovered.')
