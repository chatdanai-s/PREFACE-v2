# Part 5 of overall pipeline.
# Small module to take the top planets of our ranked set and split into separate files.
# Two files made - viable planets, and not-viable planets.
import os
import fnmatch
import pandas as pd

# Main function
def Splitter(CSV_core_folder, CSV_intermediate_folder, output_folder,
             instrument, filter_name, run_mode, toggle_sky_noise, toggle_defocus, metric_mode, viable_cumulative_cut, RMin):
    
    sky_noise_text = 'Y-SkyNoise' if toggle_sky_noise else 'N-SkyNoise'
    defocus_text = 'Y-Defocus' if toggle_defocus else 'N-Defocus'
    config_str = f'{instrument}_{filter_name}-band_for_{run_mode}_{sky_noise_text}_{defocus_text}'

    rankedtep_csv_path = CSV_intermediate_folder / 'ranked_tep_sets' / f'RankedTEPSet_{config_str}.csv'

    # Create output folders if doesnt exist
    viable_folder = (output_folder / "phase_1" / "viable_target_list")
    nonviable_folder = (output_folder / "phase_1" / "nonviable_target_list")

    viable_folder.mkdir(parents=True, exist_ok=True)  
    nonviable_folder.mkdir(parents=True, exist_ok=True)

    # Import ranked TEP file and split into viable and nonviable files
    df_ranked = pd.read_csv(rankedtep_csv_path)
    df_ranked = df_ranked.sort_values(by=metric_mode, ascending=False)
    df_ranked['CumSum'] = df_ranked[metric_mode].cumsum()
    df_ranked['SumFrac'] = df_ranked['CumSum'] / df_ranked['CumSum'].max()

    df_top = df_ranked[df_ranked[metric_mode] >= RMin]    # Takes all planets ranked above the cumulative cut.
    df_low = df_ranked[df_ranked[metric_mode] < RMin]     # Takes all planets ranked below the cumulative cut.

    # Export phase 1 output
    print("[ViabilitySplitter] Outputting rank threshold-separated TEPSets to .csv now..." )
    df_top.to_csv(viable_folder / f'TopTEPSet_{config_str}_{metric_mode}-Mode_{viable_cumulative_cut*100}%_Cut.csv',
                  index=False)
    df_low.to_csv(nonviable_folder / f'UnusableTargets_{config_str}_{metric_mode}-Mode_{viable_cumulative_cut*100}%_Cut.csv',
                  index=False)
    
    print(rf'[ViabilitySplitter] New viability division ({viable_cumulative_cut} fraction cut) for {metric_mode}-mode successful.')
        
    
    # Split top set into CSVs for individual targets, to be used by Phase Two.

    # Clear the Input and Output Parts folders if a previous run for this instrument/filter combination has been executed.
    # This should prevent any conflict as TEPCat is periodically updated and planet ranks are revised.
    individual_csv_initial_filename = f'{config_str}_{metric_mode}-Mode_{viable_cumulative_cut*100}%_Cut'

    phase_2_input_folder = (CSV_intermediate_folder / 'phase_2_inputs')
    phase_2_output_individual_folder = (output_folder / 'phase_2' / 'individual_planets')
    phase_2_output_individual_folder.mkdir(parents=True, exist_ok=True)  # Create directory if doesnt exist

    def clear_configparts(configpath):
        configfiles = [os.path.join(path, f)
                       for path, dirs, files in os.walk(configpath)
                       for f in fnmatch.filter(files, f'{individual_csv_initial_filename}_*.csv')
        ]
        for f in configfiles:
            os.remove(f)

    clear_configparts(phase_2_input_folder)
    clear_configparts(phase_2_output_individual_folder)
        
    print('[ViabilitySplitter] Previous input/output wipe successful.')
    

    # Now we can safely write our new parts!
    df_top.reset_index()
    rows = df_top.shape[0]
    for row in range(rows):
        df_top.iloc[[row]].to_csv(phase_2_input_folder / f'{individual_csv_initial_filename}_{row}.csv', index=False)

    print('Phase One of PREFACE Complete.\n')
