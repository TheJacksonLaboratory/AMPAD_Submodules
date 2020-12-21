#PBS -l nodes=1:ppn=20,walltime=48:00:00
#PBS -m bea
#PBS -M nikhil.milind@jax.org
#PBS -q batch

cd $PBS_O_WORKDIR

module load R/3.4.4

R --no-restore --no-save --quiet < 8_case_module_enrichment_mayo.R