# Part 2.2 of overall pipeline.
# Further b-values can be retrieved from exoplanets.eu, particularly for new objects.

import pandas as pd
from ImpactMerger import BWriter

# Main function
# csvpath argument MUST be csvcorepath from the master shell.
def ExoeuImpacts(csvpath):
    # Takes in FullTEPSet, processes names to find those with missing b.
    da = pd.read_csv(rf'{csvpath}/FullTEPSetWithExoOrgImpacts.csv')
    da.insert(1, 'Planet_temp', da['Planet'])   # Temporary column

    da['Planet'] = da['Planet'].replace({'-00':'-', '-0':'-',
                                         '_00':'_', '_0':'_'},
                                        regex=True)
    da['Planet'] = da.apply(lambda row: BWriter(row), axis=1)
    
    # Takes in table of transiting planets from exoplanets.eu, processes names to find b values
    colnames = ['name', 'impact_parameter']
    db = pd.read_csv(rf'{csvpath}/exoplaneteu_catalog.csv', usecols=colnames)
    db.rename(columns={'name': 'Planet',
                       'impact_parameter': 'Impact Parameter_2'},
              inplace=True)

    # This scheme is not perfect yet, needs polishing
    db['Planet'] = db['Planet'].replace({'AU Mic ':'AU_Mic_',
                                         'A b':'b', ' A':'', ' B':'', ' (AB)':'',
                                         r'([A-Z]) ': r'\1_',
                                         ' ':'',
                                         'EPIC_22881391b':'EPIC_228813918b'
                                         },
                                        regex=True)

    # Append impact parameter to missing values in original dataframe then export
    dc = pd.merge(da, db, on='Planet', how='left')
    dc['Impact Parameter'] = dc['Impact Parameter'].fillna(dc['Impact Parameter_2'])
    dc = dc.drop(columns=['Planet', 'Impact Parameter_2'])
    dc = dc.rename(columns={'Planet_temp': 'Planet'})   # Recover old planet names  

    # Reads back out with:
    dc.to_csv(rf'{csvpath}/FullTEPSetWithAllImpacts.csv', index=False)

    print('[ExoplanetseuImpactMerger] Impact parameters from exoplanets.eu recovered.')
