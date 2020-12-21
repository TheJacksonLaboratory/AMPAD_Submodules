# Run iterativeWGCNA for Mayo
mkdir ../results/1_cleaning/mayo/cohort/
qsub -F "../../../../clean_data/1_cleaning/mayo_iterativeWGCNA.tsv ../results/1_cleaning/mayo/cohort/" 1d_iterative_WGCNA_module_individual.sh

# Run iterativeWGCNA for Mayo
mkdir ../results/1_cleaning/mssm/cohort/
qsub -F "../../../../clean_data/1_cleaning/mssm_iterativeWGCNA.tsv ../results/1_cleaning/mssm/cohort/" 1d_iterative_WGCNA_module_individual.sh

# Run iterativeWGCNA for Mayo
mkdir ../results/1_cleaning/rosmap/cohort/
qsub -F "../../../../clean_data/1_cleaning/rosmap_iterativeWGCNA.tsv ../results/1_cleaning/rosmap/cohort/" 1d_iterative_WGCNA_module_individual.sh
