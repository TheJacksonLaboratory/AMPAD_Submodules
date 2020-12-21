# Iterate through all module files in Mayo
for filename in ../clean_data/1_cleaning/iterativeWGCNA_mayo/*.tsv; do
  
  # Create directory to store results
  mkdir ../results/1_cleaning/mayo/${filename##*/}
  
  # Run iterativeWGCNA
  qsub -F "../../../../clean_data/1_cleaning/iterativeWGCNA_mayo/${filename##*/} ../clean_data/1_cleaning/mayo/${filename##*/}" 1d_iterative_WGCNA_module_individual.sh
  
done

# Iterate through all module files in MSSM
for filename in ../clean_data/1_cleaning/iterativeWGCNA_mssm/*.tsv; do
  
  # Create directory to store results
  mkdir ../results/1_cleaning/mssm/${filename##*/}
  
  # Run iterativeWGCNA
  qsub -F "../../../../clean_data/1_cleaning/iterativeWGCNA_mssm/${filename##*/} ../clean_data/1_cleaning/mssm/${filename##*/}" 1d_iterative_WGCNA_module_individual.sh
  
done

# Iterate through all module files in ROSMAP
for filename in ../clean_data/1_cleaning/iterativeWGCNA_rosmap/*.tsv; do
  
  # Create directory to store results
  mkdir ../results/1_cleaning/rosmap/${filename##*/}
  
  # Run iterativeWGCNA
  qsub -F "../../../../clean_data/1_cleaning/iterativeWGCNA_rosmap/${filename##*/} ../clean_data/1_cleaning/rosmap/${filename##*/}" 1d_iterative_WGCNA_module_individual.sh
  
done
