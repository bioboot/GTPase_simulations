#!/bin/bash 

r1=`qsub gtp.01.pbs`
echo $r1
r2=`qsub -W depend=afterok:$r1 gtp.02.pbs`
echo $r2
r3=`qsub -W depend=afterok:$r2 gtp.03.pbs`
echo $r3
r4=`qsub -W depend=afterok:$r3 gtp.04.pbs`
echo $r4
r5=`qsub -W depend=afterok:$r4 gtp.05.pbs`
echo $r5
r6=`qsub -W depend=afterok:$r5 gtp.06.pbs`
echo $r6
r7=`qsub -W depend=afterok:$r6 gtp.07.pbs`
echo $r7
r8=`qsub -W depend=afterok:$r7 gtp.08.pbs`
echo $r8
r9=`qsub -W depend=afterok:$r8 gtp.09.pbs`
echo $r9
r10=`qsub -W depend=afterok:$r9 gtp.10.pbs`
echo $r10
r11=`qsub -W depend=afterok:$r10 gtp.11.pbs`
echo $r11
r12=`qsub -W depend=afterok:$r11 gtp.12.pbs`
echo $r12
r13=`qsub -W depend=afterok:$r12 gtp.13.pbs`
echo $r13
r14=`qsub -W depend=afterok:$r13 gtp.14.pbs`
echo $r14
r15=`qsub -W depend=afterok:$r14 gtp.15.pbs`
echo $r15
r16=`qsub -W depend=afterok:$r15 gtp.16.pbs`
echo $r16
