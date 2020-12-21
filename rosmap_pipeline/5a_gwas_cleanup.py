#!/usr/bin/env python

# Notes
# * Scripts integrated from work by Cai John
# * Data taken from the Amp-AD Repositories

# Goal
# * Take in .ps files from EMMAX and reduce to most significant SNPs
# * Save file back for analysis in R

# Usage:
# python3 -W ignore 5a_gwas_cleanup.py

# Import libraries
import os

# Static variables
SOURCE_DIRECTORY = "clean_data/5_gwas/5a_gwas_results"
RESULTS_DIRECTORY = "clean_data/5_gwas/5a_cleaning"

# Retrieve all .ps files in 16_gwas/
gwas_files = {}
for file in os.listdir(SOURCE_DIRECTORY):
    if file.endswith(".ps"):
        gwas_files[file] = os.path.join(SOURCE_DIRECTORY, file)

# Open each file and only choose lines that 
# have a p<0.05 significance
for file in gwas_files:

    print("Creating modified file for %s" % file)

    # Check if file already exists and skip
    if os.path.isfile(os.path.join(RESULTS_DIRECTORY, (file + ".sig"))):
        print("%s already exists. Delete file to recreate it" % os.path.join(RESULTS_DIRECTORY, (file + ".sig")))
        continue

    # Keep track of line number
    line_number = 0

    # Open data file to read and modified file to write
    with open(gwas_files[file], "r") as read, open(os.path.join(RESULTS_DIRECTORY, (file + ".sig")), "w") as write:
        
        # Iterate through reading file
        for line in read:
            
            line_number += 1
            print("Line %s: SNP %s" % (line_number, line.split("\t")[0]), end="\r")

            # Retreive p-value and write if less than 0.05
            p_val = float(line.split("\t")[3])
            if p_val < 0.05:
                write.write(line)
    
    print("Written to new file %s" % (file + ".sig"))
    print()

