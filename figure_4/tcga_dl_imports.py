import os
import re
import numpy as np
import pandas as pd


def read_dl_results(output_folder):
    files = os.listdir(output_folder)
    unique_training_ids = [re.search('ID_([a-z0-9]*)_', f) for f in files]
    unique_training_ids = [e.group(0).replace('ID_', '').replace('_', '') 
                           for e in unique_training_ids if e is not None]
    unique_training_ids = np.unique(unique_training_ids).astype(str)
    print(len(unique_training_ids))

    # Load results
    global_rank_df = pd.DataFrame()
    all_rank_files = {training_id:[e for e in files 
                                    if (training_id in e) and (test in e) and ('ustat.txt' in e)][0]
                      for training_id in unique_training_ids}
    all_rank_scores = {
        training_id: pd.read_csv(output_folder + f, 
                                 header=None if data_type == 'TCGA' else 0,
                                 index_col=0)
        for training_id, f in all_rank_files.items()
    }
    all_rank_scores = {
        training_id: df[1] if data_type == 'TCGA' else df['PR-PD']
        for training_id, df in all_rank_scores.items()
    }
    rank_scores_df = pd.DataFrame(all_rank_scores).T
    rank_scores_df['AUC'] = rank_scores_df['ustat'] / rank_scores_df['product_samples']
    
    median_classifier_ID = rank_scores_df[rank_scores_df['AUC'] <= np.median(rank_scores_df['AUC'])].sort_values('AUC').index[-1]
    
    if type_agg == 'median':
        baseline_prediction = 'baseline_predicted_AUC_ID_%s_%sMann-Whitney-ls_prediction.csv'%(
            median_classifier_ID,
            '%s_'%(GDSC_drug_id) if GDSC_drug_id is not None else ''
        )
        baseline_prediction = pd.read_csv(output_folder + baseline_prediction)
    elif type_agg == 'aggregated':
        baseline_prediction_df = []
        for x in os.listdir(output_folder):
            if ('baseline_predicted_AUC_ID' not in x) or ('_prediction.csv' not in x):
                continue

            baseline_prediction_df.append(pd.read_csv(output_folder + x))

        baseline_prediction = pd.concat(baseline_prediction_df, axis=1)
        baseline_prediction['aggregated'] = np.mean(baseline_prediction['predicted'], axis=1)
        baseline_prediction = baseline_prediction.loc[:,~baseline_prediction.columns.duplicated()]
        
    baseline_prediction.columns = ['bcr_patient_barcode', 'measure_of_response', 'predicted']
    return baseline_prediction, rank_scores_df