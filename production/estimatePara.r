##

##- Determine initial aMD simulation paramaters (Amber12)

##   The input values for EPtot, DIHED and NATOM

##   come from cMD out/log file. NRES can come from protein only PDB

## 

 

library(bio3d)

pdb <- read.pdb("./sys_gtp_nowat.pdb")

NRES = sum(pdb$calpha) ## Number of protein residues

NATOM_PRO  =  nrow(pdb$atom) ## Number of total atoms
NATOM =  31168

#EPtot = -99044.8987
#DIHED = 1760.9184
EPtot = -99044.8987
DIHED = 1760.9184


##-- Total potential values, use a factor of 0.2

##   (for lower boost use value btwn 0.15-0.19)

alphaP = NATOM*(1/5)

EthreshP = EPtot + alphaP

 

##-- Dihedral potential values

##   (use factor of 3.5 kcal/mol/residue => Approximate energy contribution per d.f.)

alphaD = (1/5)*(3.5*NRES)

EthreshD = (3.5*NRES) + DIHED

 

##-- For a higher acceleration, add to EthreshD multiples of alphaD.

EthreshD.2 = EthreshD + alphaD

EthreshD.3 = EthreshD.2 + alphaD

EthreshD.4 = EthreshD.3 + alphaD

EthreshD.4 = EthreshD.4 + alphaD 

 

 

cmd <- paste(" iamd=3,\n   EthreshD = ",

             c(EthreshD,  EthreshD.2, EthreshD.3, EthreshD.4),

             ", alphaD = ", alphaD,

             ",\n   EthreshP = ", EthreshP, 

             ", alphaP = ", alphaP, ",\n  ntwprt = ", NATOM_PRO, ",\n", sep="")

 

out <- list(EthreshP = EthreshP, alphaP = alphaP,

           EthreshD = EthreshD, alphaD = alphaD,

           EthreshD.2 = EthreshD.2,

           EthreshD.3 = EthreshD.3,

           EthreshD.4 = EthreshD.4, cmd=cmd)

 

cat(cmd)
