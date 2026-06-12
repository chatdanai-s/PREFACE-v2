# Part 2.2 of overall pipeline.
# Further b-values can be retrieved from exoplanets.eu, particularly for new objects.

import pandas as pd
from preface.P1_ImpactMerger import BWriter

# Main function

def ExoeuImpacts(CSV_core_folder):
    # Locate files
    tep_with_exob_csv_path = CSV_core_folder / 'FullTEPSetWithExoOrgImpacts.csv'
    exoplanets_eu_csv_path = CSV_core_folder / 'exoplaneteu_catalog.csv'
    tep_with_allb_csv_path = CSV_core_folder / 'FullTEPSetWithAllImpacts.csv'


    # Takes in FullTEPSet, processes names to find those with missing b.
    df_tep_exob = pd.read_csv(tep_with_exob_csv_path)
    df_tep_exob.insert(1, 'Planet_temp', df_tep_exob['Planet'])   # Temporary column

    df_tep_exob['Planet'] = df_tep_exob['Planet'].replace({'-00':'-', '-0':'-',
                                                           '_00':'_', '_0':'_'},
                                                           regex=True)
    df_tep_exob['Planet'] = df_tep_exob.apply(lambda row: BWriter(row), axis=1)
    

    # Takes in table of transiting planets from exoplanets.eu, processes names to find b values
    colnames = ['name', 'impact_parameter']
    df_eu = pd.read_csv(exoplanets_eu_csv_path, usecols=colnames)
    df_eu.rename(columns={'name': 'Planet',
                          'impact_parameter': 'Impact Parameter_from_eu'},
                 inplace=True)

    # This scheme is not perfect yet, needs polishing
    df_eu['Planet'] = df_eu['Planet'].replace({'AU Mic ':'AU_Mic_',
                                             'A b':'b', ' A':'', ' B':'', ' (AB)':'',
                                             r'([A-Z]) ': r'\1_',
                                             ' ':'',
                                             'EPIC_22881391b':'EPIC_228813918b'
                                             },
                                             regex=True)


    # Append impact parameter to missing values in original dataframe then export
    df_tep_allb = pd.merge(df_tep_exob, df_eu, on='Planet', how='left')
    df_tep_allb['Impact Parameter'] = df_tep_allb['Impact Parameter'].fillna(df_tep_allb['Impact Parameter_from_eu'])
    df_tep_allb = df_tep_allb.drop(columns=['Planet', 'Impact Parameter_from_eu'])
    df_tep_allb = df_tep_allb.rename(columns={'Planet_temp': 'Planet'})   # Recover old planet names  

    # Reads back out with:
    df_tep_allb.to_csv(tep_with_allb_csv_path, index=False)

    print('[ExoplanetseuImpactMerger] Impact parameters from exoplanets.eu recovered.')
