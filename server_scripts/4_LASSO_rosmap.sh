#PBS -l nodes=1:ppn=1,walltime=48:00:00
#PBS -m bea
#PBS -M nikhil.milind@jax.org
#PBS -q batch

cd $PBS_O_WORKDIR

module load R/3.4.4

R --no-restore --no-save --quiet < 4_LASSO_rosmap.R