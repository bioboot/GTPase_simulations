"wiki.tbl" <-
function(mat) {
  cat(paste("| **",
            paste(colnames(mat), collapse = "** | **"),
            "** |",sep=""), sep="\n")
  for(i in 1:nrow(mat)) {
    cat(paste("|",paste(mat[i,], collapse = " | "),"|"),sep="\n")
  }
}

read.propka <- function(file) {

  ##-- examine out.propka
  ## ~/software/pdb2pqr/pdb2pqr.py --with-ph=7 --chain --salt --ffout=amber --ff=amber 1BE9.pdb out
  ## pka <- read.propka("out.propka")
  ## wiki.tbl(pka$his)

  ## modified for PDB2PQR V3.0  - yao


  ## lineformat
  first <- c(1, 4, 10, 17, 25, 50, 68, 86)
  last  <- c(3, 9, 16, 24, 49, 67, 85, 103)

  split.string <- function(x) {
    x <- substring(x, first, last)
    x[nchar(x) == 0] <- as.character(NA)
    x
  }
  trim <- function (s) {
    s <- sub("^ +", "", s)
    s <- sub(" +$", "", s)
    s[(s=="")]<-NA
    s
  }

  raw.lines <- readLines(file)

  type <- substring(raw.lines,1,6)
  blank <- type==""
  type <- type[!blank]; raw.lines <- raw.lines[!blank]
  col <- which(type==" RESID")
  start <- col+2
  end <- which(type=="SUMMAR")-2
  ##trim(split.string(raw.lines[start]))

  mat <- matrix(trim(sapply(raw.lines[start:end], split.string)),byrow=TRUE, ncol=length(first))
  colnames(mat) <- c("resid","resno","pka","location","desolv","hside","hback","coulombic")
  return(list(pka=mat, his=mat[mat[,"resid"]=="HIS",]))
}

write.ligand <- function(pdb, lig="GDP", file) {
  # Supporting ligands: GDP, GSP, GNP, GTP, MG, 1KX

##   Remove water and hydrogen; don't renumber residues, otherwise resno for 
##   ligands incorrect (a bug?)
  npdb <- convert.pdb(pdb, type="amber", renumber=FALSE, rm.wat=TRUE, rm.h=TRUE)
  hpdb <- NULL
  if(lig %in% npdb$atom[, "resid"]) {
    hetatom <- npdb$atom[npdb$atom[, "resid"]==lig,,drop = FALSE ]
   #? chain=A is not necessary?
   #hetatom <- npdb$atom[npdb$atom[,"resid"]=="1KX" & npdb$het[,"chain"]=="A",]
     # NOTE: check if one atom or not

    switch(lig, 
      GSP = { 
         hetatom[, "resid"] <- "GTP"
         hetatom[hetatom[,"elety"]=="S1G","elety"] <- "O1G"
         hetatom[,"elety"] <- sub("(.*)'$","\\1*", hetatom[,"elety"])
         outn <- "gtp"
      },
      GNP = {
         hetatom[, "resid"] <- "GTP"
         hetatom[hetatom[,"elety"]=="N3B","elety"] <- "O3B"
         hetatom[,"elety"] <- sub("(.*)'$","\\1*", hetatom[,"elety"])
         outn <- "gtp"
      },
      GTP = {
         hetatom[,"elety"] <- sub("(.*)'$","\\1*", hetatom[,"elety"])
         outn <- "gtp"
      },
      GDP = {
         hetatom[,"elety"] <- sub("(.*)'$","\\1*", hetatom[,"elety"])
         outn <- "gdp"
      },
      MG = {
         hetatom[, "resid"] <- "MG2"
         outn <- "mg"
      },
      "1KX" = {
         outn <- "1kx"
      },
      stop("Unknow ligand")
    )
    hetxyz  <- as.numeric(t(hetatom[, c("x","y","z")]))

    hpdb$atom <- hetatom
    hpdb$xyz  <- hetxyz
     
    write.pdb(hpdb, resno=vec2resno(1:sum(!duplicated(hetatom[,"resno"])), hetatom[,"resno"]), 
              eleno=(1:nrow(hetatom)), file=paste(file,"_", outn,".pdb",sep=""))
  }
}

##-- Read input structure
##   (Normally this is a single 'core-fitted' chain)
setup_01 <- function(infile) {

  require(bio3d)
  
  id <- tolower( substr(basename(infile), 1, 4) )
  pdb <- read.pdb( infile )
  
  
  ##-- Format protein for pKa determination
  pkfile <- paste(id,"_4pka.pdb",sep="")
  ##   Remove water and hydrogen; don't renumber residues, otherwise resno for 
  ##   ligands incorrect (a bug?)
  opdb <- convert.pdb(pdb, type="pdb", renumber=FALSE, rm.wat=TRUE, rm.h=TRUE)
  # remove ligands
  spdb <- trim.pdb(opdb, inds=atom.select(opdb, "protein"))
  write.pdb(spdb, file=pkfile)
  
  
  ##-- Format Ligand for Leap setup
  # It doesn't matter if ligand does NOT exist
  # GDP
  write.ligand(pdb, "GDP", file=id)
  
  # GTP
  write.ligand(pdb, "GTP", file=id)
  write.ligand(pdb, "GNP", file=id)
  write.ligand(pdb, "GSP", file=id)
  
  # MG
  write.ligand(pdb, "MG", file=id)
  
  # 1KX or abacavir
  write.ligand(pdb, "1KX", file=id)
  
  
  ##-- Run pdb2pqr & read output with read.propka()
  system( "if test ! -r pka; then mkdir pka; fi" )
  # --chain, keep the chain IDs
  system( paste("pdb2pqr.py --chain --with-ph=7 --ff=amber", pkfile,"pka/out") )
  pka <- read.propka("pka/out.propka")
  wiki.tbl(pka$his)
  
  ##-- Note.
  ## For each Histidine residue examine its pKa value in outfile.propka.
  ## If the pH < pKa then you may consider protonating both ND1 and NE2
  ## (i.e. making the residue HIP).
  ## Otherwise you should chose between HID (delta), HIE (epsilon).
  ##
  ## Choosing between HID and HIE often entails analyzing what is around
  ## each ring nitrogen in each histidine.
  ## If there is an electron donor (O) nearby this may indicate that the
  ## nitrogen should be protonated.
  ## If there is an electron acceptor (H) nearby perhaps it should not
  ## be protonated.
  ##
  ## Once you have come to a justifiable decision, replace all HIS
  ## residue names in your pdbfile with one of the HID, HIE or HIP options.
  ## [N.B. Remember to keep a record of your chosen protonation states!]
  ##
  
  ##-- Change HIS residue type and remove hydrogens
  nels <- read.pqr("pka/out")
  # now renumber is fine because only peptide considered
  # must remain hydrogen for checking HIS protonation state
  new <- convert.pdb(nels, type="amber", renumber=TRUE, rm.h=FALSE)
  his.inds <- atom.select(new, resid="HID")
  if(length(his.inds$atom)>0) {
     his.resno <- as.numeric(new$atom[his.inds$atom, "resno"])
     ii <- bounds(his.resno, dup.inds=TRUE)
     his.type <- apply(ii, 1, function(x) {
        bhd1 = FALSE
        bhe2 = FALSE
        for(j in x[2]:x[3]) {
           if(new$atom[his.inds$atom[j], "elety"] == "HD1") bhd1=TRUE
           if(new$atom[his.inds$atom[j], "elety"] == "HE2") bhe2=TRUE
        }
        if(bhd1 && bhe2) {
           return ("HIP")
        } else if(bhd1) {
           return ("HID") 
        } else if(bhe2) {
           return ("HIE")
        } else {
           stop("No proton in either ND1 or NE2")
        }
     } )
     new$atom[his.inds$atom, "resid"] <- vec2resno(his.type, his.resno)
  }
  
  # remove hydrogen
  new <- convert.pdb(new, type="amber", renumber=TRUE, rm.h=TRUE)
  write.pdb(new, file=paste(id,"_prot.pdb",sep=""))
} # end funtion setup_01


##-- Run Leap
## see "setup_02.leap"


##-- Run sander min1-4, heat and equil

##-- Run production dynamics on flux, group cluster, or gpu
