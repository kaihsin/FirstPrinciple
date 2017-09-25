#!/bin/sh
#BSUB -a mvapich
#BSUB -n 2 
#BSUB -J vasp
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q  test              # depend on -n
#BSUB -R 'span[ptile=48]' # uncomment when multi-nodes
#BSUB -R "ipathavail==0"  # uncomment when multi-nodes

source /pkg/mpi/intel/12/mvapich2-1.8/bin/mpivars.sh

VASP=/home/$USER/vasp_ntu

mpiexec -np ${LSB_DJOB_NUMPROC} ${VASP}
