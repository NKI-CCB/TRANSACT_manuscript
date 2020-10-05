# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

Computing proportion of each component for one value of gamma.
"""

import sys
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from precise.pv_computation import PVComputation
from precise.interpolation import Interpolation
from precise.matrix_operations import _center_kernel, _right_center_kernel, _left_center_kernel
from precise.kernel_computer import KernelComputer
from precise.PRECISE import PRECISE

sys.path.insert(0, '../src/')
from taylor_expansion import *

def compute_proportion(gamma,
                    n_pc,
                    n_pv,
                    normalized_data_df,
                    source_data_key,
                    target_data_key,
                    n_interpolation=100,
                    n_jobs=30):
    
    # Fit domain adaptation classifier
    print('FIT CLASSIFIER')
    PRECISE_clf = PRECISE(kernel='rbf',
                        kernel_params={'gamma':gamma},
                        n_components=n_pc,
                        n_jobs=n_jobs,
                        verbose=1)
    PRECISE_clf.fit(normalized_data_df[source_data_key],
                    normalized_data_df[target_data_key],
                    n_pv=n_pv,
                    step=n_interpolation,
                    with_interpolation=True)
    
    # Angular coefficients
    optimal_time = PRECISE_clf.optimal_time
    source_angular = PRECISE_clf.interpolation_._gamma_interpolations(optimal_time)
    target_angular = PRECISE_clf.interpolation_._xi_interpolations(optimal_time)

    source_angular = np.diag(source_angular)
    target_angular = np.diag(target_angular)
    
    # Gaussian depth, i.e. offset
    print('COMPUTE GAUSSIAN DEPTH')
    source_gaussian_depth, target_gaussian_depth = compute_gaussian_depth(gamma, normalized_data_df, source_data_key, target_data_key)
    
    sigma_offset_source = PRECISE_clf.principal_vectors_.gamma_coef['source'].dot(source_gaussian_depth)
    sigma_offset_target = PRECISE_clf.principal_vectors_.gamma_coef['target'].dot(target_gaussian_depth)
    sigma_offset = source_angular.dot(sigma_offset_source) + target_angular.dot(sigma_offset_target)
    
    offset_contribution = {
        'source': np.square(sigma_offset_source),
        'target': np.square(sigma_offset_target),
        'consensus': np.square(sigma_offset_source)
    }
    
    # Linear term
    print('COMPUTE LINEAR CONTRIBUTION')
    p = normalized_data_df[source_data_key].shape[1]
    source_basis_eval = Parallel(n_jobs=n_jobs, verbose=1)(delayed(basis_function)(normalized_data_df[source_data_key],
                                                                                   i,
                                                                                   gamma)
                                                       for i in range(p))
    target_basis_eval = Parallel(n_jobs=n_jobs, verbose=1)(delayed(basis_function)(normalized_data_df[target_data_key],
                                                                                   i,
                                                                                   gamma)
                                                       for i in range(p))
    target_basis_eval = np.array(target_basis_eval).T
    source_basis_eval = np.array(source_basis_eval).T
    
    sigma_linear_source = PRECISE_clf.principal_vectors_.gamma_coef['source'].dot(source_basis_eval)
    sigma_linear_target = PRECISE_clf.principal_vectors_.gamma_coef['target'].dot(target_basis_eval)

    sigma_linear = source_angular.dot(sigma_linear_source) + target_angular.dot(sigma_linear_target)
    linear_contribution = {
        'source': np.square(np.linalg.norm(sigma_linear_source, axis=1)),
        'target': np.square(np.linalg.norm(sigma_linear_target, axis=1)),
        'consensus': np.square(np.linalg.norm(sigma_linear, axis=1))
    }
    
    # Interaction term
    print('COMPUTE INTERACTION CONTRIBUTION')
    loadings = Parallel(n_jobs=n_jobs, verbose=1)(delayed(interaction_loadings_genes)(i, gamma,
                                                                                    PRECISE_clf,
                                                                                    normalized_data_df,
                                                                                    source_data_key,
                                                                                    target_data_key)
                                               for i in range(p))

    loadings = np.array([np.array(e) for e in loadings])
    source_interaction_contribution = [np.sum(np.square(e[:,0,:]), axis=0) for e in loadings]
    source_interaction_contribution = np.array(source_interaction_contribution)
    source_interaction_contribution = np.sum(source_interaction_contribution, axis=0)
    
    target_interaction_contribution = [np.sum(np.square(e[:,1,:]), axis=0) for e in loadings]
    target_interaction_contribution = np.array(target_interaction_contribution)
    target_interaction_contribution = np.sum(target_interaction_contribution, axis=0)
    
    consensus_interaction_contribution = [np.sum(np.square(e[:,2,:]), axis=0) for e in loadings]
    consensus_interaction_contribution = np.array(consensus_interaction_contribution)
    consensus_interaction_contribution = np.sum(consensus_interaction_contribution, axis=0)
    
    interaction_contribution = {
        'source': source_interaction_contribution,
        'target': target_interaction_contribution,
        'consensus': consensus_interaction_contribution
    }
    
    return PRECISE_clf, {
        'offset': offset_contribution,
        'linear': linear_contribution, 
        'interaction': interaction_contribution
    }
        