pdbfile=/Users/hyangl/project/ras/results/2015/1127_cMD_Ras/pdb/wt/5P21.pdb

bprod=FALSE
bhead=FALSE
bnowat=FALSE
bmhc=FALSE
bnopka=FALSE

dirn=md_`basename $pdbfile .pdb`
if test ! -r $dirn; then mkdir $dirn; fi
mdpath=`pwd`/$dirn
if test $bprod == "FALSE"; then cp $pdbfile $mdpath; fi
pdbfile=`basename $pdbfile`
cd $mdpath

cp $MD_SCRIPTS/scripts/setup_01.r ./
cp $MD_SCRIPTS/scripts/write.pdb.R ./

R -e "source(\"setup_01.r\"); source(\"write.pdb.R\"); setup_01(\"$pdbfile\")" &> pka.log

nline=`grep -n -i -e warning  pka.log`
if grep -i -e error pka.log; then
   exit 1;
fi
awk '/^\| HIS \|/&&$7!="NA"{print $4, $7}' pka.log

# check files and define variables
protfile=(`ls ????_prot.pdb 2> /dev/null`)
if test ${#protfile[@]} -eq 0 || test ${#protfile[@]} -gt 1; then
   echo Ambiguous protein files: ????_prot.pdb
   exit 1
fi

if test $bnopka == "TRUE"; then
   if ! cp *_4pka.pdb $protfile; then
      echo Error: failed copy file
      exit 1
   fi
fi

head=`echo $protfile | sed 's/^\(....\).*/\1/'`
ligfiles=(`ls "$head"_gdp.pdb "$head"_gtp.pdb "$head"_1kx.pdb "$head"_mg.pdb 2> /dev/null`)
lig=(`echo ${ligfiles[@]} | sed -e "s/$head//g" -e 's/.pdb//g'`)

# 2. check HIS state; can be skipped
echo 'Checking HIS state; Can be skipped...'
hisresid=( `awk '/^\| HIS \|/&&$7!="NA"{printf "%d ", $4}' pka.log` )
# Check HIS in the output structure (Default: out)
# Meanwhile, check the output PKA of HIS in setup_01.r run
# Decide HIS state, and rename them in the PDB file
rm -f check_his.pml
echo "reinitialize" >> check_his.pml
echo "load pka/out" >> check_his.pml
j=1
for i in ${hisresid[@]}; do
   echo "select his$j, i. $i" >> check_his.pml
   echo "select pp$j, byres ((his$j and n. nd1+ne2) around 5)" >> check_his.pml
   let j=j+1
done
echo "orient his1" >> check_his.pml
echo "hide (all)" >> check_his.pml
echo "as sticks, his*" >> check_his.pml
echo "as sticks, pp*" >> check_his.pml

#3. download force field files
echo Download force field files...
if test ! -r ff; then
   mkdir ff
#wget http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/cof/frcmod.phos
#wget http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/cof/GTP.prep
#wget http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/cof/GDP.prep
   cp $MD_SCRIPTS/data/frcmod.phos ff/
   cp $MD_SCRIPTS/data/GTP.prep ff/
   cp $MD_SCRIPTS/data/GDP.prep ff/
   #sed 's/ gtp  INT/ GTP  INT/' GTP.prep > ff/GTP.prep
   #sed 's/ gdp  INT/ GDP  INT/' GDP.prep > ff/GDP.prep

   if test $bmhc == TRUE; then
      # abacavir
      cp $MD_SCRIPTS/data/aba.off ff/
      cp $MD_SCRIPTS/data/aba.frcmod ff/
   fi
fi


echo Parepare input files for tleap...
cp $MD_SCRIPTS/scripts/run_gen_leapfile.r ./

bsolvate=TRUE
if test $bnowat == TRUE; then bsolvate=FALSE; fi
if test $bmhc == TRUE; then
   # For MHC
   R -e "source(\"run_gen_leapfile.r\"); \
     gen_var_complex_mhc(\"$protfile\", \
        unlist(strsplit(\"${ligfiles[@]}\", \" \")), solvate=$bsolvate)"
else
   R -e "source(\"run_gen_leapfile.r\"); \
     write.leap(\"$protfile\", unlist(strsplit(\"${ligfiles[@]}\", \" \")), \
        file=\"setup_02.leap\", ext=\"\", solvate=$bsolvate, ccbond=FALSE)"
fi

rm -f leap.log

# 4.5 generate pdb files
echo Generating PDB files...
cp $MD_SCRIPTS/scripts/run_gen_pdb.sh ./

#######

# up
rsync -av ../md_5P21/ hyangl@axiom.ccmb.med.umich.edu:/home/hyangl/ras/2015/1127_cMD_Ras/md_5P21/

tleap -s -f setup_02.leap
sh run_gen_pdb.sh

# down
rsync -av hyangl@axiom.ccmb.med.umich.edu:/home/hyangl/ras/2015/1127_cMD_Ras/md_5P21/ ../md_5P21/


######

#5. EM, heating, and equilibration
# If you want to add more (e.g. control) simulations, be sure to set lig=_gtp here !!

echo Running EM, heating, and equilibration...
if test ! -r min_equil; then mkdir min_equil; fi
cd min_equil

# Here if you have multiple runs for one state, start from equilibration
mkdir run01
cd run01
## for run02-05, just copy files in run01 and modify the random seed in dyna_heat.sander
## echo $RANDOM

cp $MD_SCRIPTS/scripts/min_1.sander ./
cp $MD_SCRIPTS/scripts/min_2.sander ./
cp $MD_SCRIPTS/scripts/min_3.sander ./
cp $MD_SCRIPTS/scripts/min_4.sander ./
# cp $MD_SCRIPTS/scripts/dyna_heat.sander ./
## Here I use $RANDOM to assign the random seed! :)
sed -e "s/71277/$RANDOM/" $MD_SCRIPTS/scripts/dyna_heat.sander > dyna_heat.sander
cp $MD_SCRIPTS/scripts/dyna_equil.sander ./
cp ../../sys"$lig"_box.prmtop ./
cp ../../sys"$lig"_box.inpcrd ./
cp ../../sys"$lig"_nowat.prmtop ./
cp ../../sys"$lig"_nowat.pdb ./
sed -e "s/sys_box/sys${lig}_box/g" $MD_SCRIPTS/scripts/write.axiom.r > write.axiom.r
sed -e "s/sys_box/sys${lig}_box/g" $MD_SCRIPTS/scripts/run_md_emin.sh > run_md_emin.sh
chmod a+x run_md_emin.sh

lig2=gtp
if test ${#lig[@]} -ge 1; then lig2=`echo $lig | sed 's/^_//'`; fi
R -e 'source("write.axiom.r")' -e "write.axiom(\"$lig2\", emin=TRUE, cluster=\"axiom\", memsize=\"4000mb\", nnode=1, ncpu=8)"

cd ../../

# up
rsync -av ../md_5P21/ hyangl@axiom.ccmb.med.umich.edu:/home/hyangl/ras/2015/1127_cMD_Ras/md_5P21/




########### This is part II ##########
## after equilibration, set up production!

# down
rsync -av hyangl@axiom.ccmb.med.umich.edu:/home/hyangl/ras/2015/1127_cMD_Ras/md_5P21/ ../md_5P21/


echo Preparing files for production...
if test ! -r production; then mkdir production; fi
cd production

if test ! -r 00Readme || ! grep -q -e "^run01:" 00Readme; then
   echo "run01: 100ns cMD simulation" >> 00Readme
fi

if test ! -r run01; then mkdir run01; fi
cd run01

lig=_gtp
cp ../../sys"$lig"_nowat.pdb ./
cp ../../sys"$lig"_box.inpcrd ./
cp ../../sys"$lig"_box.prmtop ./
cp ../../sys"$lig"_nowat.prmtop ./
cp ../../min_equil/run01/dyna_equil.rst ./dyna.00.rst

pdb=./sys"$lig"_nowat.pdb
natom=`sed -n '2p' ./sys"$lig"_box.inpcrd`
eptot=`sed -n '/A V E R A G E S/,/R M S  F/p' ../../min_equil/run01/dyna_equil.out |\
       sed -n '/EPtot/p' | sed 's/.*EPtot[[:space:]]*=[[:space:]]*\(.*\)/\1/'`
dihed=`sed -n '/A V E R A G E S/,/R M S  F/p' ../../min_equil/run01/dyna_equil.out |\
       sed -n '/DIHED/p' | sed 's/.*DIHED[[:space:]]*=[[:space:]]*\(.*\)/\1/'`
sed -e "s%pdb[[:space:]]*<-[[:space:]]*read.pdb.*%pdb <- read.pdb(\"$pdb\")%" \
    -e "s/NATOM[[:space:]]*=.*/NATOM = $natom/" \
    -e "s/EPtot[[:space:]]*=.*/EPtot = $eptot/" \
    -e "s/DIHED[[:space:]]*=.*/DIHED = $dihed/" \
       $MD_SCRIPTS/scripts/estimatePara.r > estimatePara.r
R -e 'source("estimatePara.r")'
cp $MD_SCRIPTS/scripts/catPara.r ./
para=`Rscript catPara.r`

# if cMD, uncomment following line
para=`echo $para | sed 's/.*ntwprt/ntwprt/g'`

sed '/iamd=3/,$'d $MD_SCRIPTS/scripts/dyna_prod.sander > ./dyna_prod.sander
echo " $para" >> ./dyna_prod.sander
echo " &end" >> ./dyna_prod.sander
sed '/iamd=3/,$'d $MD_SCRIPTS/scripts/dyna_prod0.sander > ./dyna_prod0.sander
echo " $para" >> ./dyna_prod0.sander
echo " &end" >> ./dyna_prod0.sander

sed -e "s/sys_box/sys${lig}_box/g" $MD_SCRIPTS/scripts/write.axiom.r > write.axiom.r

# !!!!! The origanal one doesn't work!(md_group_02: .01.pbs instead apo.01.pbs!) 
# let lig2=apo!
lig2=gtp
# if test ${#lig[@]} -ge 1; then lig2=`echo $lig | sed 's/^_//'`; fi
R -e 'source("write.axiom.r")' -e "write.axiom(\"$lig2\", iter=1:16, cluster=\"axiom\", memsize=\"4000mb\", nnode=1, ncpu=1, ngpu=1)"
cd ../../

# up
rsync -av ../md_5P21/ hyangl@axiom.ccmb.med.umich.edu:/home/hyangl/ras/2015/1127_cMD_Ras/md_5P21/

## Above is for run01, for run02-05, just copy the directory run01 and change the corresponding dyna.00.rst
## cp ../../min_equil/run01/dyna_equil.rst ./dyna.00.rst




############  check #############

## check equilibration!
mkdir analysis
cd analysis
cp $MD_SCRIPTS/scripts/run_analysis_basic.sh ./
cp $MD_SCRIPTS/scripts/process_mdout.perl ./
cp $MD_SCRIPTS/scripts/plot_analysis_basic.r ./
# check trj of equilibration!
# vmd load sys_gtp_box.prmtop (AMBER7 Parm)
# vmd load dyna_equil.mdcrd (AMBER Coordinates with Periodic Box


## unwrap coordinates of trajectories
cp $MD_SCRIPTS/scripts/run_oktraj.sh production/run01/
# modify the number of chunks (20->16)
./run_oktraj.sh


## some modification for mac
## date: 10/26/2015

# 1. sed space
# sed uses \s to represent space in linux
# but it doesn't work on mac!
# [[:space:]] can universally represent space!
# use the following commend to modify 5P21_E232A_setup_cMD.sh !:P
# :%s/\\s/[[:space:]]/gc
# e.g.
sed -e "s%pdb[[:space:]]*<-[[:space:]]*read.pdb.*%pdb <- read.pdb(\"$pdb\")%" \
    -e "s/NATOM[[:space:]]*=.*/NATOM = $natom/" \
    -e "s/EPtot[[:space:]]*=.*/EPtot = $eptot/" \
    -e "s/DIHED[[:space:]]*=.*/DIHED = $dihed/" \
$MD_SCRIPTS/scripts/estimatePara.r > estimatePara.r

# 2. sed newlines
# it seems \n cannot be used by mac to present enter
# Xinqiu helped me modify the part generating $para

# 3. pdb2pqr
# pdb2pqr v1.8 (v2.0 doesn't generate the out.propka file..)
# modify setup_01.r: change pdb2pqr to pdb2pqr.py




