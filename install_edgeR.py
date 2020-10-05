import rpy2.robjects as robjects

robjects.r('''
    if (!requireNamespace("BiocManager"))
        install.packages("BiocManager")
    BiocManager::install(version="3.8")''')

robjects.r('''
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("edgeR", version = "3.8")''')