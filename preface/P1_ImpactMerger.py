# Part 2.1 of overall pipeline.
# TEPCat does not store impact parameters, so they must be corralled from elsewhere.

import pandas as pd

# Small function to massage names such that both catalogues use same naming convention.
def BWriter(row):
    if row['Planet'][-1].isdigit() == True:
        row['Planet'] = row['Planet'] + 'b'
        
    return row['Planet']

# Main function
def ExoOrgImpacts(CSV_core_folder):
    # Locate files
    fulltepset_csv_path = CSV_core_folder / 'FullTEPSet.csv'
    exoorglist_csv_path = CSV_core_folder / 'ExoOrgList.csv'
    tep_with_exoorg_csv_path = CSV_core_folder / 'FullTEPSetWithExoOrgImpacts.csv'


    # Takes in all names in FullTEPSet and processes them.
    df_fulltep = pd.read_csv(fulltepset_csv_path)
    df_fulltep.insert(1, 'Planet_temp', df_fulltep['Planet'])   # Temporary column to be used in final planet name

    # Setup planet column as key
    df_fulltep['Planet'] = df_fulltep['Planet'].replace({'-00':'-', '-0':'-',
                                                         '_00':'_', '_0':'_'},
                                                         regex=True
                                                         )
    df_fulltep['Planet'] = df_fulltep.apply(lambda row: BWriter(row), axis=1)


    # Takes in table of transiting planets from exoplanets.org
    # exoplanets.org is a catalogue built by hand, so no auto-updater here.
    # Also, exoplanets.org is not updated with planets published after June 2018
    colnames = ['NAME', 'B']
    df_exoorg = pd.read_csv(exoorglist_csv_path, usecols=colnames)
    df_exoorg.rename(columns={'NAME': 'Planet', 'B': 'Impact Parameter'}, inplace=True)

    # Process planet names before merging
    df_exoorg['Planet'] = df_exoorg['Planet'].replace({'55 Cnc e':'55_Cnc_e', 'OGLE2':'OGLE',
                                                       'BD +20 594': 'K2-56b', 
                                                       r'([A-Z]) ': r'\1_',
                                                       'A b':'b',
                                                       ' ':''
                                                       },
                                                       regex=True
                                                       )
    

    # Append impact parameter to original dataframe
    df_combined = pd.merge(df_fulltep, df_exoorg, on='Planet', how='left')
    df_combined = df_combined.drop('Planet', axis=1)                      # Recover old planet names
    df_combined = df_combined.rename(columns={'Planet_temp': 'Planet'})  

    # Export
    df_combined.to_csv(tep_with_exoorg_csv_path, index=False)

    print('[ImpactMerger] Impact parameters from exoplanets.org recovered.')
