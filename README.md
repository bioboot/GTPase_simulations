# Comparative structural dynamic analysis of GTPases

This repository contains molecular simulation input files for the project **Comparative structural dynamic analysis of GTPases**
by Hongyang Li, Xin-Qiu Yao and Barry J. Grant

The results of this work are currently *in press* at PLoS Computational Biology. A pre-print of this work is available in the [BioRxiv](https://www.biorxiv.org/content/early/2018/07/16/370197).

## Dependencies  

* [AMBER12](http://ambermd.org/) molecular simulation package.  
* AMBER force field [ff99SB](http://ambermd.org/AmberModels.php).
* [PROPKA](https://github.com/jensengroup/propka-3.1) software for predicting the pKa values of ionizable groups in proteins.  

## Description of major file   

* End-to-end bash file: `setup_cMD.sh`  
* Directory for energy minimization: `min_equil`  
* Directory for final MD simulation: `production`  
* Additional system input files as AMBER format `inpcrd`, `prmtop` and `pdb`   

For questions about individual files please contact [Hongyang](hyangl@umich.edu) and for project questions contact [Barry](bjgrant@ucsd.edu). 
