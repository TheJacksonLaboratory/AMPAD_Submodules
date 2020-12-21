#!/usr/bin/env python

# Notes
# * Scripts integrated from work by Cai John
# * Data taken from the Amp-AD Repositories

# Goal
# * Read in data from NHGRI-EBI GWAS Catalog
# * Locate closest alleles for all suggestive alleles
# * Append information for closest alleles

# Import libraries
import os
import pandas
import pdb

# Static variables
NHGRI_EBI = "raw_data/NHGRI_EBI_data.tsv"
SOURCE_DIRECTORY = "clean_data/5_gwas/5d_gwas_variants/"
RESULTS_DIRECTORY = "clean_data/5_gwas/5e_gwas_lookup/"

# Retrieve all suggestive SNPs .csv files
snp_files = {}
for file in os.listdir(SOURCE_DIRECTORY):
    if file.endswith(".csv"):
        snp_files[file] = os.path.join(SOURCE_DIRECTORY, file)

print("*** Determining Closest Known Loci for each SNP ***")

# Open each file and compare each SNP to pre-recorded SNPs
for file in snp_files:

    print("Analyzing file %s" % snp_files[file])

    # Read in file
    snps = pandas.read_csv(snp_files[file])

    # Create data frame to store output
    snps_out = snps.copy()
    snps_out["Closest.SNP.BP"] = pandas.Series([-1] * len(snps_out), index=snps_out.index)
    snps_out["Closest.Distance"] = pandas.Series([-1] * len(snps_out), index=snps_out.index)
    snps_out["Reported.Genes"] = pandas.Series([""] * len(snps_out), index=snps_out.index)
    snps_out["Mapped.Genes"] = pandas.Series([""] * len(snps_out), index=snps_out.index)
    snps_out["PubMed.ID"] = pandas.Series([""] * len(snps_out), index=snps_out.index)
    snps_out["NHGRI.EBI.Line.Number"] = pandas.Series([-1] * len(snps_out), index=snps_out.index)

    # Iterate over each SNP
    for snp_index in range(len(snps)): 

        print("\tWorking on SNP: %s" % snps["SNP"][snp_index])

        # Keep track of line with smallest SNP distance
        smallest_data = []
        smallest_dist = 10**9
        line_number = -1

        # Open connection to NHGRI-EBI data
        with open(NHGRI_EBI, "r") as nhgri:
            
            line_counter = 0
            # Iterate through every line
            for line in nhgri:
                line_data = line.strip().split("\t")
                line_counter += 1
                # Make sure that chromosomes match (at least)
                if str(line_data[11]) == str(snps["CHR"][snp_index]):
                    if abs(int(line_data[12]) - int(snps["BP"][snp_index])) < smallest_dist:
                        smallest_dist = abs(int(line_data[12]) - int(snps["BP"][snp_index]))
                        smallest_data = line_data
                        line_number = line_counter

        # Add data into new data frame
        snps_out["Closest.SNP.BP"][snp_index] = int(smallest_data[12])
        snps_out["Closest.Distance"][snp_index] = smallest_dist
        snps_out["Reported.Genes"][snp_index] = smallest_data[13]
        snps_out["Mapped.Genes"][snp_index] = smallest_data[14]
        snps_out["PubMed.ID"][snp_index] = smallest_data[1]
        snps_out["NHGRI.EBI.Line.Number"][snp_index] = line_number

    # Write out file
    file_name = "_".join(file.split(".")[0].split("_")[0:3]) + ".csv"
    print("Creating file %s" % file_name)

    snps_out.to_csv(os.path.join(RESULTS_DIRECTORY, file_name), header=True, index=False)
