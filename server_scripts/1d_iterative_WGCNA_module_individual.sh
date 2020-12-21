#PBS -l nodes=1:ppn=1,walltime=48:00:00
#PBS -q batch

# Change to the working directory of the script
cd $PBS_O_WORKDIR

# Load modules for iterativeWGCNA
module load R/3.5.1
module load python/3.6.5

# Run iterativeWGCNA
iterativeWGCNA -i $1 -o $2 --verbose --wgcnaParameters power=6,minKMEtoStay=0.6,minCoreKME=0.6,minCoreKMESize=100