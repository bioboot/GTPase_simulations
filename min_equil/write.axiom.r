write.axiom <- 
  function(prefix = "r", emin = FALSE, iter=1:10, 
           cluster="axiom", runtime="72:00:00",
           memsize="4000mb", nnode=1, ncpu=1, ngpu=1, new.vel = FALSE) {
  if (nchar(prefix) > 8) {
    warning(paste("Your filename 'prefix' is over 8 characters in length,\n\t",
                  "consider shortning so iteration digits are visible with qstat"))
  }

  if(ncpu < ngpu) ncpu = ngpu
  post.pmemd=""
  mpi=""
  if(ncpu>1) {
     post.pmemd=".MPI"
     mpi="$MPIRUN -np $NPROCS "
  }

  if(cluster=="local") {
    if(emin) {
       cat("#!/bin/bash
./run_md_emin.sh ", ncpu, sep="", file = paste("submit_", prefix, ".sh", sep=""))
    } else {
       ## Run pmemd on local machine 
       commfiles <- paste(prefix,".",sprintf("%02.0f", iter),".sh",sep="")
       submitfile <- paste("submit_",prefix,".sh",sep="")
       
       for(i in 1:length(iter)) {
         cat(paste("#!/bin/bash
prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"
NPROCS=",ncpu,"
AMBER=$AMBERHOME/bin/pmemd.cuda",post.pmemd,"
MPIRUN=$MPI_HOME/bin/mpirun
AMBER_ARGS=\"-O -i ", ifelse(i==1 && new.vel, "dyna_prod0.sander", "dyna_prod.sander"), " -p sys_gtp_box.prmtop -c dyna.$prv.rst -r dyna.$cur.rst -o dyna.$cur.out -x dyna.$cur.traj.nc -e dyna.$cur.mden -inf dyna.$cur.inf -amd dyna.$cur.amdlog\"
echo Running cmd: ",mpi,"$AMBER $AMBER_ARGS
", mpi, "$AMBER $AMBER_ARGS
echo Finished at time: `date`\n", sep=""),
           file=commfiles[i])
       }
   
       ##-- Write a master submission shell script for dependent jobs
       head <- "#!/bin/bash\n"
       middle <- paste("./", commfiles, "\n", collapse="", sep="")
       cat(paste(head, middle, sep=""), file=submitfile)
    }
    system("chmod a+x *.sh")
  }
  else if(cluster=="axiom") {
    if(emin) {
       commfile = paste(prefix, ".pbs", sep="")
       cat("#PBS -S /bin/sh
#PBS -M hyangl@umich.edu
#PBS -A grant_lab
#PBS -q first
#PBS -j oe
#PBS -o out.$PBS_JOBNAME.log
#PBS -V
#PBS -l nodes=",nnode,":ppn=",ncpu,",mem=",memsize,",walltime=",runtime,"

echo Running job name $PBS_JOBNAME with ID $PBS_JOBID on host $PBS_O_HOST
echo Nodes for run:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`

./run_md_emin.sh $NPROCS
", sep = "", file = commfile)
       cat("#!/bin/bash
qsub ", commfile, "\n", sep = "", file = paste("submit_pbs_", prefix, ".sh", sep=""))
    } else {
       ## Run pmemd on AXIOM
       pbsfiles <- paste(prefix,".",sprintf("%02.0f", iter),".pbs",sep="")
       submitfile <- paste("submit_pbs_",prefix,".sh",sep="")
   
       for(i in 1:length(iter)) {
   
         cat(paste("#PBS -S /bin/sh
#PBS -M hyangl@umich.edu
###PBS -A bjgrant_flux
#PBS -A grant_lab
###PBS -N test
#PBS -q first
###PBS -m abe
#PBS -j oe
#PBS -o out.$PBS_JOBNAME.log
###PBS -e err.$PBS_JOBNAME.log
#PBS -V
#PBS -l nodes=",nnode,":ppn=",ncpu,":gpus=",ngpu,",mem=",memsize,",walltime=",runtime,"

echo Running job name $PBS_JOBNAME with ID $PBS_JOBID on host $PBS_O_HOST
echo Nodes for run:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS cores
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Time is `date`
echo Directory is `pwd`

prv=",sprintf("%02.0f", (iter[i]-1)),"
cur=",sprintf("%02.0f",  iter[i]),"

AMBER=$AMBERHOME/bin/pmemd.cuda",post.pmemd,"
MPIRUN=$MPI_HOME/bin/mpirun
AMBER_ARGS=\"-O -i ", ifelse(i==1 && new.vel, "dyna_prod0.sander", "dyna_prod.sander"), " -p sys_gtp_box.prmtop -c dyna.$prv.rst -r dyna.$cur.rst -o dyna.$cur.out -x dyna.$cur.traj.nc -e dyna.$cur.mden -inf dyna.$cur.inf -amd dyna.$cur.amdlog\"

echo Running cmd: ",mpi,"$AMBER $AMBER_ARGS
",mpi,"$AMBER $AMBER_ARGS
echo Finished at time: `date`\n", sep=""),
           file=pbsfiles[i])
       }
   
       ##-- Write a master submission shell script for dependent jobs
       ## pbsfiles <- paste("r",".",sprintf("%02.0f", iter),".pbs",sep="")
       pbsids <- paste("r", iter, sep="")
       head <- paste("#!/bin/bash \n\n",
                     pbsids[1],"=`qsub ",pbsfiles[1],"`\necho $r1\n", sep="")
   
       middle <- paste(paste(pbsids[-1],"=`qsub -W depend=afterok:$",
                       pbsids[-length(pbsids)]," ",pbsfiles[-1],"`\necho $",
                       pbsids[-1], "\n", sep=""), collapse="", sep="")
   
       cat(head, middle, file=submitfile, sep="")
     }
     system("chmod a+x submit_*.sh")
  }


#######

}
