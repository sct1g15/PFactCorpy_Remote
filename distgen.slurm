#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=10
#SBATCH --time=24:00:00

source /local/software/conda/miniconda-py3-new/bin/activate /scratch/sct1g15/mydocuments/conda/pymol
module load mpich/3.2.1/intel

mpirun -np 400 python distgen2.py --trajin /mainfs/scratch/sct1g15/mydocuments/RPC/3PWL_A0201_free_310K/500ns_con/md_0_1_protein_molurcenter_uratom_urfit_time.xtc --topolin /mainfs/scratch/sct1g15/mydocuments/RPC/3PWL_A0201_free_310K/500ns_con/md_0_1_protein.pdb
