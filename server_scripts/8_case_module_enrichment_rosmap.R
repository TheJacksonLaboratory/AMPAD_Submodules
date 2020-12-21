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
rosmap.eigensamples <- readRDS("clean_data/2_build_svd_sets/rosmap_eigensamples.RDS")

# Load cleaned data for each cohort
rosmap.cleaned <- readRDS("clean_data/1_cleaning/rosmap_cleaned.RDS")

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
# ROSMAP Projections
#------------------------------------------------

# Create matrix to hold value of projections
rosmap.samples <- matrix(nrow=nrow(rosmap.cleaned), ncol=ncol(rosmap.eigensamples))
rownames(rosmap.samples) <- rosmap.cleaned$Patient
colnames(rosmap.samples) <- colnames(rosmap.eigensamples)

# Iterate through each patient and each submodule
cat("Calculating Projections for the ROSMAP Cohort")
for (patient.index in 1:nrow(rosmap.cleaned)) {
  
  cat(paste0("Calculating Projections for Patient ", patient.index, "\n"))
  
  for (submodule.index in 1:ncol(rosmap.eigensamples)) {
    
    cat(paste0("Submodule ", colnames(rosmap.eigensamples)[submodule.index]))
    
    # Retrieve eigensample for the submodule
    # Retrieve sample for the given patient
    rosmap.samples[patient.index, submodule.index] <- projection(
      sample=rosmap.cleaned[patient.index,-(1:19)] %>%
        as.numeric(),
      eigensample=rosmap.eigensamples[,submodule.index]
    )
    
    # Produce permutations to calculate projection null distribution
    perms.projs <- mclapply(
      1:N.PERMS,
      FUN=function(i) projection(
        sample=rosmap.cleaned[patient.index,-(1:19)] %>%
          as.numeric() %>%
          sample(., nrow(rosmap.eigensamples)),
        eigensample=rosmap.eigensamples[,submodule.index]
      ),
      mc.cores=num.cores
    )
    perms.projs <- unlist(perms.projs) %>%
      as.numeric()
    
    perms.projs <- perms.projs[order(perms.projs)]
    proj.upper <- perms.projs[round(N.PERMS * 0.975)]
    proj.lower <- perms.projs[round(N.PERMS * 0.025)]
    
    # If projection is p>0.05, consider it not significant
    if (rosmap.samples[patient.index, submodule.index] > proj.lower && rosmap.samples[patient.index, submodule.index] < proj.upper) {
      rosmap.samples[patient.index, submodule.index] <- 0
    }
    rosmap.samples[patient.index, submodule.index] <- rosmap.samples[patient.index, submodule.index] / abs(rosmap.samples[patient.index, submodule.index])
  }
  cat("\n")
}

# Save ROSMAP Samples
saveRDS(rosmap.samples, file="results/8_case_module_enrichment/rosmap_sample_module_enrichment.RDS")

