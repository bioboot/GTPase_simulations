invisible(capture.output(source("estimatePara.r")))
cmd <- paste("iamd=3,\n   EthreshD = ", EthreshD, ", alphaD = ", alphaD, 
             ",\n   EthreshP = ", EthreshP, ", alphaP = ", alphaP, 
             ",\n  ntwprt = ", NATOM_PRO, ",\n", sep="")
cat(cmd)
