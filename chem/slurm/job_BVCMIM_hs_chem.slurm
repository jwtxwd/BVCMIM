#!/bin/bash
# --------------------------------------------------------------
### PART 1: Requests resources to run your job.
# --------------------------------------------------------------
###Set the partition where the job will run
#SBATCH --partition=public

###Set your qos that you want to run under
###SBATCH --qos=proj-546

### Optional. Set the job name
#SBATCH --job-name=jia_job

### Optional. Set the output filename.
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=slurm_out/%x-%j.out
#SBATCH --error=slurm_out/%x-%j.err
### REQUIRED. Specify the PI group for this job
### SBATCH --account=manager
### Optional. Request email when job begins and ends
### Specify high priority jobs
### SBATCH --qos=user_qos_manager
#SBATCH --mail-type=ALL
### Optional. Specify email address to use for notification
#SBATCH --mail-user=jiawi@mail.uc.edu
### REQUIRED. Set the partition for your job.
### SBATCH -p general

### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=1
#SBATCH --nodes=1
###SBATCH --cpus-per-task=40
###SBATCH --gres=gpu
###SBATCH --gres=gpu:a100:1
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=100gb
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=24:00:00

# --------------------------------------------------------------
### PART 2: Executes bash commands to run your job
# --------------------------------------------------------------
### SLURM Inherits your environment. cd $SLURM_SUBMIT_DIR not needed
pwd; hostname; 

date
### Load required modules/libraries if needed
module load R/4.2.1
#source /home/zhouw6/.conda/envs/pinn
#cd /home/zhouw6/.conda/envs/pinn

#source activate /home/zhouw6/.conda/envs/torch
#cd /N/lustre/project/proj-262/PINN_DI/Codes

#python DG_PINN_NS_loss_weighting.py

#python DG_PINN_KdV_3.py

#python DG_PINN_1d_wave_equation_loss_weighting.py

#python DG_PINN_beam_with_analytic_solution_losswe.py

#python DG_PINN_beam_adaptive_weighting2.py
#python DG_PINN_1d_wave_equation_hard_constrain.py
#python DG_PINN_1d_wave_equation_case3_hard_constrain.py
#python DG_PINN_1d_wave_equation_case3_soft_constrain.py
#python DG_PINN_1d_wave_equation_case3_ill_constrain.py

#python DG_PINN_1d_wave_equation_soft_constrain_NTK_L2.py
R CMD BATCH BVCMIM/chem/R/BVCMIM_hs_chem.R

date
~