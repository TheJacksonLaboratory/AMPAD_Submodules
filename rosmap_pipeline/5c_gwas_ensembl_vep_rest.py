#!/usr/bin/env python

# Notes
# * Scripts integrated from work by Cai John
# * Data taken from the Amp-AD Repositories

# Goal
# * Read in significant SNPs and their RefSNP IDs
# * Query REST API from ENSEMBL's Variant Effects Predictor 
# * Tabulate results and write into CSV file

# Import libraries
import os
import pandas
import requests
import sys
import json

# Static variables
SOURCE_DIRECTORY = "clean_data/5_gwas/5b_gwas_setup/"
RESULTS_DIRECTORY = "clean_data/5_gwas/5c_gwas_vep/"

# Retrieve all suggestive SNPs .csv files
snp_files = {}
for file in os.listdir(SOURCE_DIRECTORY):
    if file.endswith("_suggestive_SNPs.csv"):
        snp_files[file] = os.path.join(SOURCE_DIRECTORY, file)

# Open each file and query REST API for associated effects
for file in snp_files:

    print("Analyzing file %s" % snp_files[file])
    
    # Check if file already exists and skip
    if os.path.isfile(os.path.join(RESULTS_DIRECTORY, snp_files[file])):
        print("%s already exists. Delete file to recreate it" % os.path.join(RESULTS_DIRECTORY, snp_files[file]))
        continue

    # Read in file
    snps = pandas.read_csv(snp_files[file])
    snps = snps.dropna()

    # Define REST API parameters
    server = "http://rest.ensembl.org"
    ext = "/vep/human/id"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    # The ENSEMBL REST API only accepts up to 200 variants per request
    # Split variants into sets of 200
    resp_json = []
    for start_index in range(0, len(snps["RefSNP.ID"]), 200):

        start = start_index
        end = min(start_index + 200, len(snps["RefSNP.ID"]))
        # Compile IDs into JSON string
        post_json = {
            "ids": (snps["RefSNP.ID"][start:end]).tolist()
        }
        post_json_str = json.dumps(post_json)

        # Post request
        req = requests.post(server + ext, headers=headers, data=post_json_str)

        # Decode JSON from response
        resp_json += req.json()

    # Populate data frame with resulting data
    # Response JSON has the following format:
    # [
    #   {
    #       'input': 'rs23498103',
    #       'transcript_consequences': {
    #           'gene_symbol': 'JSDF',
    #           'biotype': 'lincRNA',
    #           'consequence_terms': [ 'intron_variant', 'non_coding_transcript_variant' ]
    #       }
    #   },
    #   ...
    # ]

    # Create new data frame to store data
    cons = pandas.DataFrame(columns=["RefSNP.ID", "Gene.Symbol", "Biotype", "Consequence.Terms"])
    cons_count = 0

    for snp_dict in resp_json:

        # Check to make sure SNP has transcript consequences
        if "transcript_consequences" in snp_dict:

            # Loop through all transcript consequences
            for consequence in snp_dict["transcript_consequences"]:

                print("\t SNP %s" % snp_dict["input"])
                
                # Create a new row in data frame with information
                # Make sure data that we are looking for is available in the transcript consequence
                cons.loc[cons_count] = [
                    snp_dict["input"],
                    consequence["gene_symbol"] if "gene_symbol" in consequence else "",
                    consequence["biotype"] if "biotype" in consequence else "",
                    ",".join(consequence["consequence_terms"]) if "consequence_terms" in consequence else ""
                ]
                cons_count += 1

    out_file = file.split("_suggestive_SNPs.csv")[0]
    print("Creating file %s" % (out_file + "_VEP.csv"))

    # Write consequences to CSV file
    cons.to_csv(os.path.join(RESULTS_DIRECTORY, (out_file + "_VEP.csv")), header=True, index=False)
