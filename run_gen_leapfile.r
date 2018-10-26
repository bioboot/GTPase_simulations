write.leap <- function(pdbfile, ligfiles=NULL, file="out.leap", 
   ext=paste("_", sub("\\.pdb.*", "", basename(pdbfile)), sep=""), 
   solvate=TRUE, ccbond=FALSE, cc.refpdb=NULL, cc.refresno=NULL) {

   head <- sub("\\.pdb.*", "", basename(pdbfile))
   ligs <- sub(".*_(.*).pdb", "\\1", ligfiles)
   if(length(ligs)==0 || nchar(ligs)==0) ligs <- NULL
   cat("logFile leap", ext, ".log\n", file=file, sep="")
   cat("source leaprc.ff99SB\n", file=file, append=TRUE)
   cat("source leaprc.gaff\n", file=file, append=TRUE)
   cat("addAtomTypes{{\"O3\" \"O\" \"sp2\"}}\n", file=file, append=TRUE)
   cat("\n", file=file, append=TRUE)
   if(file.exists("ff/aba.frcmod")) {
      cat("loadamberparams ff/aba.frcmod\n", file=file, append=TRUE)
      cat("loadoff ff/aba.off\n", file=file, append=TRUE)
   }
   if(file.exists("ff/GDP.prep"))
      cat("loadamberprep ff/GDP.prep\n", file=file, append=TRUE)
   if(file.exists("ff/GTP.prep"))
      cat("loadamberprep ff/GTP.prep\n", file=file, append=TRUE)
   if(file.exists("ff/frcmod.phos"))
      cat("loadamberparams ff/frcmod.phos\n", file=file, append=TRUE)
   cat("\n", file=file, append=TRUE)
   if(ccbond) {
      require(bio3d)
      # get the equivalent CYS-CYS bonds
      pdb <- read.pdb(pdbfile)
      pdb2 <- read.pdb(cc.refpdb)
      aln <- seqaln.pair(seqbind(pdbseq(pdb), pdbseq(pdb2)), id=c(pdbfile, cc.refpdb))
      pdbs <- read.fasta.pdb(aln)
      inds <- which(pdbs$resno[2,] %in% cc.refresno)
      # check
#      pdbs$ali[2, inds]
      bcys <- pdbs$resno[1, inds]
      cat("CYS-CYS bonds correspondence\n")
      cat("  Reference: ", paste(cc.refresno, collapse=", "), "\n", sep="")
      cat("  Current:   ", paste(bcys, collapse=", "), "\n", sep="")
      if(any(is.na(bcys)))
         stop("in finding equivalent S-S bonds")
      
      # Modify pdb file: CYS -> CYX
      npdb <- pdb
      cys.inds <- atom.select(npdb, resno=bcys)
      npdb$atom[cys.inds$atom, "resid"] <- "CYX"
      pdbfile <- paste(head, "_cyx.pdb", sep="")
      write.pdb(npdb, file=pdbfile)
   }
   cat("prot=loadpdb ", pdbfile, "\n", sep="", file=file, append=TRUE)
   if(!is.null(ligs)) {
      cat(paste(ligs, "=loadpdb ", ligfiles, sep="", collapse="\n"), 
            "\n", sep="", file=file, append=TRUE) 
   }
   cat("\n", file=file, append=TRUE)
   if(ccbond) {
      bcys <- matrix(bcys, nr=2) 
      cat(paste("bond prot.",bcys[1,],".SG prot.", bcys[2,], ".SG", 
          sep="", collapse="\n"), "\n", sep="",  file=file, append=TRUE)
   }
   cat("\n", file=file, append=TRUE)
   cat("##-- make SYSTEM one\n", file=file, append=TRUE)   
   cat("saveamberparm prot sys_nowat", ext, ".prmtop sys_nowat", ext, ".inpcrd\n", sep="", file=file, append=TRUE)
   cat("\n", file=file, append=TRUE)
   if(!is.null(ligs))  {
      cat("## System 2\n", file=file, append=TRUE)
      cat("sys=combine { prot", ligs, "}\n", file=file, append=TRUE)
      cat("saveamberparm sys sys_", ligs[1], "_nowat", ext, ".prmtop sys_", 
          ligs[1], "_nowat", ext, ".inpcrd\n", sep="", file=file, append=TRUE)
      cat("\n", file=file, append=TRUE)
   }
   if(solvate) {
      cat("##-- IONS\n", file=file, append=TRUE)
      cat("addions prot Na+ 0\n", file=file, append=TRUE)
      cat("#addions prot Cl- 0\n", file=file, append=TRUE)
      cat("\n", file=file, append=TRUE)
      cat("##-- SOLVATE\n", file=file, append=TRUE)
      cat("solvatebox prot TIP3PBOX 12\n", file=file, append=TRUE)
      cat("saveamberparm prot sys_box", ext, ".prmtop sys_box", ext, ".inpcrd\n", sep="", file=file, append=TRUE)
      cat("\n", file=file, append=TRUE)
      if(!is.null(ligs)) {
         cat("addions sys Na+ 0\n", file=file, append=TRUE)
         cat("#addions sys Cl- 0\n", file=file, append=TRUE)
         cat("solvatebox sys TIP3PBOX 12\n", file=file, append=TRUE)
         cat("saveamberparm sys sys_", ligs[1], "_box", ext, ".prmtop sys_", 
               ligs[1], "_box", ext, ".inpcrd\n", sep="", file=file, append=TRUE)
         cat("\n", file=file, append=TRUE)
      } 
   }
   cat("quit\n", file=file, append=TRUE)
}
gen_var_complex_mhc <- function(pdbfile, ligfiles=NULL, solvate=FALSE) {
   require(bio3d)
   head <- sub("\\.pdb.*", "", basename(pdbfile))
   fnm.nopep <- paste(head, "_nopep.pdb", sep="")
   fnm.polyG <- paste(head, "_polyG.pdb", sep="")
   fnm.polyA <- paste(head, "_polyA.pdb", sep="")
   pdb <- read.pdb(pdbfile)
   lens <- table(pdb$atom[pdb$calpha, "chain"])
   cid <- NULL
   for(i in 1:length(lens)) {
      if(lens[i] <= 40 ) {
         cid <- names(lens[i])
         break
      }
   }
   if(is.null(cid)) {
      warning("No peptide found")
      write.pdb(pdb, file=fnm.nopep)
      write.leap(pdbfile, ligfiles=ligfiles, file="setup_02.leap", 
        ext="", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
      write.leap(fnm.nopep, ligfiles=ligfiles, file="setup_02_nopep.leap", 
        ext="_nopep", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
   } else {    
      cid.hla <- unique(pdb$atom[!(pdb$atom[,"chain"] %in% cid), "chain"])
      spdb <- trim.pdb(pdb, inds=atom.select(pdb, chain=cid.hla))
      write.pdb(spdb, file=fnm.nopep)
    
      inds.nopep <- atom.select(pdb, chain=cid.hla)
      inds.g <- atom.select(pdb, chain=cid, elety=c("N", "CA", "C", "O", "OXT"))
      inds.a <- atom.select(pdb, chain=cid, elety=c("N", "CA", "CB", "C", "O", "OXT"))
      inds.polyG <- combine.sel(inds.nopep, inds.g, "+")
      spdb <- trim.pdb(pdb, inds=inds.polyG)
      spdb$atom[spdb$atom[, "chain"]==cid, "resid"] <- "GLY"
      write.pdb(spdb, file=fnm.polyG)
    
      inds.polyA <- combine.sel(inds.nopep, inds.a, "+")
      spdb <- trim.pdb(pdb, inds=inds.polyA)
#      inds <- atom.select(spdb, chain=cid)
#      inds <- combine.sel(inds, atom.select(spdb, chain=cid, resid="GLY"), "-")
#      spdb$atom[inds$atom, "resid"] <- "ALA"
      spdb$atom[spdb$atom[, "chain"]==cid, "resid"] <- "ALA"
      write.pdb(spdb, file=fnm.polyA)

      write.leap(pdbfile, ligfiles=ligfiles, file="setup_02.leap", 
        ext="", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
      write.leap(fnm.nopep, ligfiles=ligfiles, file="setup_02_nopep.leap", 
        ext="_nopep", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
      write.leap(fnm.polyG, ligfiles=ligfiles, file="setup_02_polyG.leap", 
        ext="_polyG", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
      write.leap(fnm.polyA, ligfiles=ligfiles, file="setup_02_polyA.leap", 
        ext="_polyA", solvate=solvate, ccbond=TRUE, 
        cc.refpdb="3upr_prot.pdb", cc.refresno=c(100,163,202,258,300,355))
   }
}
