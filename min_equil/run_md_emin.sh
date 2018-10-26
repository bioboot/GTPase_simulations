nproc=8
if test $# -ge 1; then
   nproc=$1
fi
mpirun -np $nproc $AMBERHOME/bin/sander.MPI -O -i min_1.sander -p sys_gtp_box.prmtop -c sys_gtp_box.inpcrd -ref sys_gtp_box.inpcrd -r min_1.rst -o min_1.out -x min_1.mdcrd

mpirun -np $nproc $AMBERHOME/bin/sander.MPI -O -i min_2.sander -p sys_gtp_box.prmtop -c min_1.rst -ref min_1.rst -r min_2.rst -o min_2.out -x min_2.mdcrd

mpirun -np $nproc $AMBERHOME/bin/sander.MPI -O -i min_3.sander -p sys_gtp_box.prmtop -c min_2.rst -ref min_2.rst -r min_3.rst -o min_3.out -x min_3.mdcrd

mpirun -np $nproc $AMBERHOME/bin/sander.MPI -O -i min_4.sander -p sys_gtp_box.prmtop -c min_3.rst -r min_4.rst -o min_4.out -x min_4.mdcrd

mpirun -np $nproc $AMBERHOME/bin/sander.MPI -O -i dyna_heat.sander -p sys_gtp_box.prmtop -c min_4.rst -r dyna_heat.rst -o dyna_heat.out -x dyna_heat.mdcrd -e dyna_heat.mden

mpirun -np $nproc $AMBERHOME/bin/pmemd.MPI -O -i dyna_equil.sander -p sys_gtp_box.prmtop -c dyna_heat.rst -r dyna_equil.rst -o dyna_equil.out -x dyna_equil.mdcrd -e dyna_equil.mden

