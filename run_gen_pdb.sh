# generate pdb file 
files=`ls *nowat*.prmtop 2>/dev/null`
for i in ${files[@]}; do
   head=`echo $i | sed 's/\.prmtop//'`
   $AMBERHOME/bin/ambpdb -p $head.prmtop < $head.inpcrd > $head.pdb
done
