#!/bin/bash

## script to unwrap coordinates of trajectories

if [ ! -r origin.traj ]; then
   mkdir origin.traj
   mv dyna.??.traj.nc origin.traj
fi
# Change the number of chunks (20) if necessary
for((i=1; i<=16; i++)); do
  j=`echo $i | awk '{printf("%02d",$1)}'`
  echo trajin origin.traj/dyna.$j.traj.nc > ptraj.in
  echo 'center !:WAT,GDP,GTP,MG2,Na+ origin mass' >> ptraj.in
  echo image origin center >> ptraj.in
  echo trajout dyna.$j.ok.traj.nc netcdf  >> ptraj.in
  # Change the pdb file name if necessary
  cpptraj ./sys_gtp_nowat.pdb < ptraj.in
  ln -s ./dyna.$j.ok.traj.nc dyna.$j.traj.nc
done

