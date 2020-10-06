# -*- coding: utf-8 -*-
"""
@author: Soufiane Mourragui

Functions supporting Taylor expansion of Gaussian kernel.
"""

import sys
import pandas as pd
import numpy as np


def compute_gaussian_depth(gamma, normalized_data_df, source_data_key, target_data_key):
    source_gauss_depth = np.square(np.linalg.norm(normalized_data_df[source_data_key], axis=1))
    source_gauss_depth = np.exp(- gamma * source_gauss_depth)

    target_gauss_depth = np.square(np.linalg.norm(normalized_data_df[target_data_key], axis=1))
    target_gauss_depth = np.exp(- gamma * target_gauss_depth)
    
    return source_gauss_depth, target_gauss_depth


def basis_function(x, i, rbf_gamma):
    
    norm = np.square(np.linalg.norm(x, axis=1))
    constant = np.sqrt(2*rbf_gamma)
    y = np.multiply(constant * x[:,i], np.exp(- rbf_gamma * norm))
    
    return y - np.mean(y)


def interaction_loading(i, j, rbf_gamma, norm_source, norm_target, clf, normalized_data_df, source_data_key, target_data_key):
    constant = 2*rbf_gamma
    if i == j:
        constant /= 2
    
    optimal_time = clf.optimal_time
    source_angular = clf.interpolation_._gamma_interpolations(optimal_time)
    target_angular = clf.interpolation_._xi_interpolations(optimal_time)
    source_angular = np.diag(source_angular)
    target_angular = np.diag(target_angular)
    
    X_source = normalized_data_df[source_data_key]
    X_source = np.multiply(X_source[:,i], X_source[:,j])
    X_source = np.multiply(X_source, norm_source)
    X_source = clf.principal_vectors_.gamma_coef['source'].dot(X_source)
    
    X_target = normalized_data_df[target_data_key]
    X_target = np.multiply(X_target[:,i], X_target[:,j])
    X_target = np.multiply(X_target, norm_target)
    X_target = clf.principal_vectors_.gamma_coef['target'].dot(X_target)
    
    return constant* np.array([
        X_source,
        X_target,
        source_angular.dot(X_source) + target_angular.dot(X_target)
    ])


def interaction_loadings_genes(i, rbf_gamma, clf, normalized_data_df, source_data_key, target_data_key):
    loadings = []
    
    norm_source = np.exp(-rbf_gamma*np.square(np.linalg.norm(normalized_data_df[source_data_key], axis=1)))
    norm_target = np.exp(-rbf_gamma*np.square(np.linalg.norm(normalized_data_df[target_data_key], axis=1)))
    for j in range(i+1):
        x = interaction_loading(i, j, rbf_gamma, norm_source, norm_target, clf, normalized_data_df, source_data_key, target_data_key)
        loadings.append(x)
    
    return loadings