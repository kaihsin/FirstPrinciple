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

for a in 3.924654726 3.965115084 4.005575442 4.0460358 4.086496158 4.126956516 4.167416874 4.207877232

do
echo "a= $a"

cat >POSCAR <<!
Na
$a
1.0000000000         0.0000000000         0.0000000000
0.0000000000         1.0000000000         0.0000000000
0.0000000000         0.0000000000         1.0000000000
Na
2
Direct
0.000000000         0.000000000         0.000000000
0.500000000         0.500000000         0.500000000
!
mpiexec -np ${LSB_DJOB_NUMPROC} ${VASP}

#E=`grep F OSZICAR| awk '{printf $5}'`


E=`grep "energy without entropy" OUTCAR |tail -1|awk '{printf "%12.6f\n",$5}'`
echo $a $E >> ev.dat
done
