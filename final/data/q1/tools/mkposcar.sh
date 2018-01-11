cp POSCAR.tmp POSCAR 
sed -i -e "s/LATTICE_CONST/$1/g" POSCAR
