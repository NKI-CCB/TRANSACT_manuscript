drop_intermediate = False

def compute_ROC_curve(binary_df):
    fpr, tpr, tr_thresholds = roc_curve(
        binary_df['measure_of_response'].replace('Non Responder', 1).replace('Responder', 0),
        binary_df['predicted'].values,
        drop_intermediate=drop_intermediate
    )
    AUC = roc_auc_score(
        binary_df['measure_of_response'].replace('Non Responder', 1).replace('Responder', 0),
        binary_df['predicted'].values
    )
    
    return fpr, tpr, tr_thresholds, AUC


def compute_significance(df):
    two_sided_test = scipy.stats.mannwhitneyu(
        df[df['measure_of_response'] == 'Non Responder']['predicted'].values,
        df[df['measure_of_response'] == 'Responder']['predicted'].values,
        alternative='two-sided'
    )
    one_sided_test = scipy.stats.mannwhitneyu(
        df[df['measure_of_response'] == 'Non Responder']['predicted'].values,
        df[df['measure_of_response'] == 'Responder']['predicted'].values,
        alternative='greater'
    )
    
    product = df[df['measure_of_response'] == 'Non Responder'].shape[0]
    product = product * df[df['measure_of_response'] == 'Responder'].shape[0]
    
    return {
        'one-sided':{'ustat': one_sided_test[0], 
                     'p-val': one_sided_test[1],
                     'label': '%.3f [%0.2f]'%(one_sided_test[1], one_sided_test[0]/product)},
        'two-sided':{'ustat': two_sided_test[0], 
                     'p-val': two_sided_test[1],
                     'label': '%.3f [%0.2f]'%(two_sided_test[1], two_sided_test[0]/product)}
    }