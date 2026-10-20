#!/bin/sh
#
#SBATCH --job-name=AF3_job_10_YDIT_ECOLI
#SBATCH --partition=gpu-a100
#SBATCH --account=research-as-bn
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=21
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --mail-type=ALL

module load 2025
module load cuda/12.5

singularity exec \
    --nv \
    --bind /scratch/lucasverhage/alphafold3_db:/root/public_databases \
    --bind $PWD/data:/root/af_input \
    --bind $PWD:/root/af_output_job_4_YDIT_ECOLI \
    --bind $HOME/AF3_MODEL_PARAMETERS_DIR:/root/models \
    /projects/alphafold/alphafold3.0.1.sif \
    python /app/alphafold/run_alphafold.py \
    --json_path=/root/af_input/job_10.json \
    --model_dir=/root/models \
    --db_dir=/root/public_databases \
    --output_dir=/root/af_output_job_4_YDIT_ECOLI/output_job_4_YDIT_ECOLI