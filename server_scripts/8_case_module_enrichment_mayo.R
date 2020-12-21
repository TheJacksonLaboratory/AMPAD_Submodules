# Notes
# * Scripts integrated from work by Cai John
# * Data taken from the Amp-AD Repositories

# Goal
# * Determine similarity between cases and modules
# * Assign modules of enrichment to each case

#------------------------------------------------
# Static
#------------------------------------------------

# Clear environment
rm(list = ls())
setwd("../")

# Number of permutations for permutation testing
N.PERMS <- 100

#------------------------------------------------
# Libraries
#------------------------------------------------

library(dplyr)
library(parallel)

#------------------------------------------------
# Load Data
#------------------------------------------------

# Load eigengenes for each cohort
mayo.eigensamples <- readRDS("clean_data/2_build_svd_sets/mayo_eigensamples.RDS")

# Load cleaned data for each cohort
mayo.cleaned <- readRDS("clean_data/1_cleaning/mayo_cleaned.RDS")

# Retrieve number of cores
num.cores <- detectCores()

#------------------------------------------------
# Define Projection
#------------------------------------------------

# Takes two vectors (sample and eigensample) of equal
# length. Calculates the magnitude of the projection of
# the sample onto the eigensample, maintaining the sign
# of the dot product as an indication of directionality
# when compared to the eigensample.
projection <- function(sample, eigensample) {
  
  # Calculate dot product between eigensample and sample
  dot <- (sample * eigensample) %>% 
    sum()
  
  # Calculate magnitude of eigensample
  eigensample.mag <- (eigensample * eigensample) %>%
    sum() %>%
    sqrt()
  
  # Calculate projection magnitude (with direction)
  proj <- dot / eigensample.mag
  
  # Normalize projection by length of eigensample
  proj <- proj / eigensample.mag
  
  return(proj)
}

#------------------------------------------------
# Mayo Projections
#------------------------------------------------

# Create matrix to hold value of projections
mayo.samples <- matrix(nrow=nrow(mayo.cleaned), ncol=ncol(mayo.eigensamples))
rownames(mayo.samples) <- mayo.cleaned$Patient
colnames(mayo.samples) <- colnames(mayo.eigensamples)

# Iterate through each patient and each submodule
cat("Calculating Projections for the Mayo Cohort")
for (patient.index in 1:nrow(mayo.cleaned)) {
  
  cat(paste0("Calculating Projections for Patient ", patient.index, "\n"))
  
  for (submodule.index in 1:ncol(mayo.eigensamples)) {
    
    cat(paste0("Submodule ", colnames(mayo.eigensamples)[submodule.index]))
    
    # Retrieve eigensample for the submodule
    # Retrieve sample for the given patient
    mayo.samples[patient.index, submodule.index] <- projection(
      sample=mayo.cleaned[patient.index,-(1:19)] %>%
        as.numeric(),
      eigensample=mayo.eigensamples[,submodule.index]
    )
    
    # Produce permutations to calculate projection null distribution
    perms.projs <- mclapply(
      1:N.PERMS,
      FUN=function(i) projection(
        sample=mayo.cleaned[patient.index,-(1:19)] %>%
          as.numeric() %>%
          sample(., nrow(mayo.eigensamples)),
        eigensample=mayo.eigensamples[,submodule.index]
      ),
      mc.cores=num.cores
    )
    perms.projs <- unlist(perms.projs) %>%
      as.numeric()
    
    perms.projs <- perms.projs[order(perms.projs)]
    proj.upper <- perms.projs[round(N.PERMS * 0.975)]
    proj.lower <- perms.projs[round(N.PERMS * 0.025)]
    
    # If projection is p>0.05, consider it not significant
    mayo.samples[patient.index, submodule.index] <- mayo.samples[patient.index, submodule.index] / abs(mayo.samples[patient.index, submodule.index])
    if (mayo.samples[patient.index, submodule.index] > proj.lower && mayo.samples[patient.index, submodule.index] < proj.upper) {
      mayo.samples[patient.index, submodule.index] <- 0
    }
  }
  cat("\n")
}

# Save Mayo Samples
saveRDS(mayo.samples, file="results/8_case_module_enrichment/mayo_sample_module_enrichment.RDS")
