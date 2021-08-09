kernel_surnames = [
    'linear_centered_standardized',
    'rbf_gamma_1_centered_standardized',
    'rbf_gamma_2_centered_standardized',
    'rbf_gamma_3_centered_standardized',
    'rbf_gamma_4_centered_standardized',
    'rbf_gamma_5_centered_standardized',
    'rbf_gamma_6_centered_standardized',
    'rbf_gamma_7_centered_standardized',
    'rbf_gamma_star'
]


kernel_names = ['linear', 'rbf', 'rbf', 'rbf', 'rbf', 'rbf', 'rbf', 'rbf', 'rbf']
kernel_param = [
    {},
    {'gamma': 10**(-5)},
    {'gamma': 10**(-4.5)},
    {'gamma': 10**(-4)},
    {'gamma': 10**(-3.5)},
    {'gamma': 10**(-3)},
    {'gamma': 10**(-2.5)},
    {'gamma': 10**(-2)},
    {'gamma': 5.*10**(-4)}
]

kernel_param = {k:p for k,p in zip(kernel_surnames, kernel_param)}

number_pc = {
    'source': 70,
    'target': 50
}

n_iter_dl = 50

n_pv = [20, 20, 20, 20, 20, 20, 20, 20, 20]
n_pv = {k:p for k,p in zip(kernel_surnames, n_pv)}
n_interpolation = 100

order = [
    'uncorrected_EN',
    '1',
    'linear_centered_standardized',
    '2',
    'uncorrected_network',
    'combat_network',
    '3',
    'KRR_rbf_gamma_1_centered_standardized',
    'KRR_rbf_gamma_2_centered_standardized',
    'KRR_rbf_gamma_3_centered_standardized',
    'KRR_rbf_gamma_4_centered_standardized',
    'KRR_rbf_gamma_5_centered_standardized',
    'KRR_rbf_gamma_6_centered_standardized',
    'KRR_rbf_gamma_7_centered_standardized',
    '4',
    'rbf_gamma_1_centered_standardized',
    'rbf_gamma_2_centered_standardized',
    'rbf_gamma_3_centered_standardized',
    'rbf_gamma_4_centered_standardized',
    'rbf_gamma_5_centered_standardized',
    'rbf_gamma_6_centered_standardized',
    'rbf_gamma_7_centered_standardized'
]

labels = [
    'Elastic Net',
    '',
    'PRECISE',
    '',
    'DL',
    'ComBat + DL',
    '',
    r'$\gamma$=1$\times$$10^{-5}$',
    r'$\gamma$=3$\times$$10^{-4}$',
    r'$\gamma$=1$\times$$10^{-4}$',
    r'$\gamma$=3$\times$$10^{-3}$',
    r'$\gamma$=1$\times$$10^{-3}$',
    r'$\gamma$=3$\times$$10^{-2}$',
    r'$\gamma$=1$\times$$10^{-2}$',
    '',
    r'$\gamma$=1$\times$$10^{-5}$',
    r'$\gamma$=3$\times$$10^{-4}$',
    r'$\gamma$=1$\times$$10^{-4}$',
    r'$\gamma$=3$\times$$10^{-3}$',
    r'$\gamma$=1$\times$$10^{-3}$',
    r'$\gamma$=3$\times$$10^{-2}$',
    r'$\gamma$=1$\times$$10^{-2}$',
]