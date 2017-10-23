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

VASP=/home/$USER/bin/vasp_ntu

for a in 2.7 2.75 2.8 2.85 2.9 2.95 3.0
do
echo "a= $a"

cat >POSCAR <<!
Fe bulk bcc
$a
   -0.5  0.5  0.5
    0.5 -0.5  0.5
    0.5  0.5 -0.5
Si
 1
Direct
 0.00000000 0.00000000 0.00000000
!
mpiexec -np ${LSB_DJOB_NUMPROC} ${VASP}

#E=`grep F OSZICAR| awk '{printf $5}'`


E=`grep "energy without entropy" OUTCAR |tail -1|awk '{printf "%12.6f\n",$5}'`
echo $a $E >> ev.dat
done
