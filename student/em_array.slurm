#!/bin/bash
# Example Slurm array job: one energy minimisation per array task.
# 1) make a list of prepared systems:   ls -d systems/*/ > systems.txt
# 2) submit:                             sbatch --array=1-$(wc -l < systems.txt)%50 em_array.slurm
# Adapt partition, time and modules to your cluster (this example follows the Snellius setup used for IspG).
#SBATCH --job-name=fes_em
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=rome
#SBATCH --output=logs/%x.%A_%a.out

set -euo pipefail
module load 2024 GROMACS/2024.3-foss-2024a

SYS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" systems.txt)
MDP=$(pwd)/minim.mdp
cd "$SYS"
gmx_mpi grompp -f "$MDP" -c system.gro -p system.top -o em.tpr -maxwarn 0 > grompp.log 2>&1
srun gmx_mpi mdrun -deffnm em -ntomp "$SLURM_CPUS_PER_TASK" > mdrun.log 2>&1
