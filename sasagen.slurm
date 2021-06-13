#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=01:00:00

source /local/software/conda/miniconda-py3-new/bin/activate /scratch/sct1g15/mydocuments/conda/pymol

python sasagen.py --trajin /mainfs/scratch/sct1g15/mydocuments/RPC/3PWL_A0201_free_310K/500ns_con/md_0_1_protein_molurcenter_uratom_urfit_time.xtc --topolin /mainfs/scratch/sct1g15/mydocuments/RPC/3PWL_A0201_free_310K/500ns_con/md_0_1_protein.pdb
