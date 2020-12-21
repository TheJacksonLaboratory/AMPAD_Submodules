#!/usr/bin/env python3

# Import library
import numpy as np
import pandas as pd

# Global definitions for names in DLPFC
CASE_CONTROL = ["Diagnosis", "Diagnosis_MCI"]
SUBTYPES = ["A", "B"]
MODULES = ["DLPFCturquoise", "DLPFCblue", "DLPFCyellow", "DLPFCbrown"]
SUBMODULES = [
    "DLPFCturquoise_1", "DLPFCturquoise_2", "DLPFCblue_1", "DLPFCblue_2", "DLPFCblue_3", 
    "DLPFCblue_4", "DLPFCyellow_1", "DLPFCyellow_2", "DLPFCyellow_3", "DLPFCbrown_1", "DLPFCbrown_2"
]

# Define a common function for loading in association results
def read_assoc_files(col, assoc_names):
    assocs = None
    for name in assoc_names:
        print("Reading %s.pheno.output.ps" % name)
        df = pd.read_csv(
            "clean_data/5_gwas/5a_gwas_results/" + name + ".pheno.output.ps",
            sep="\t", header=None, names=["SNP", "CHR", "Beta", "P"], index_col=0
        )
        if type(assocs) == pd.DataFrame:
            assocs[name] = df[col]
        else:
            assocs = pd.DataFrame(df[col])
            assocs.columns = [name]
            
    return assocs

# Read in association files' p-values as matrix
assoc_mtx = read_assoc_files("P", CASE_CONTROL + SUBTYPES + MODULES + SUBMODULES)
assoc_mtx.to_csv("clean_data/11_genetic_arch/assoc_mtx.csv")

# Read in association files' beta values as matrix
assoc_mtx = read_assoc_files("Beta", CASE_CONTROL + SUBTYPES + MODULES + SUBMODULES)
assoc_mtx.to_csv("clean_data/11_genetic_arch/assoc_mtx_beta.csv")
