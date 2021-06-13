#!/bin/bash

#SBATCH --ntasks-per-node=40
#SBATCH --nodes=2
#SBATCH --time=00:10:00
#SBATCH --partition=scavenger

source /local/software/conda/miniconda-py3-new/bin/activate /scratch/sct1g15/mydocuments/conda/pymol
module load mpich/3.2.1/intel

mpirun -np 80 python HIST_M8_modelfitting_MPI.py -x '20180705 A0201 4 state.xlsx' -t 298 -p 7 -st 'UV EXP' 
