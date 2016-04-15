library("rpx")
library("mzID")
library("MSnID")
library("MSnbase")

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args)!=4) {
  stop("Areguments ", call.=FALSE)
} else {
  score.treshold  = as.numeric(args[1])
  error.treshold  = as.numeric(args[2])
  fdr             = as.double(args[3])
  iteration       = as.numeric(args[4])
}

print(paste("score_treshold:", score.treshold))
print(paste("error_treshold:", error.treshold))
print(paste("fdr:", fdr))
print(paste("iteration:", iteration))

setwd("/root/data/")

mz.files <- list.files(path = ".", pattern ="mzid", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

####################################### FILTERING ###################################################

prj               <- MSnID(".")
mzid              <- read_mzIDs(object = prj, mzids = mz.files)
mzid              <- assess_termini(mzid)
mzid              <- assess_missed_cleavages(mzid)
show(mzid)
#print(paste("Before Filtering:",length(accessions(prj))))
prj               <- apply_filter(mzid, "numIrregCleavages == 0")
prj               <- apply_filter(prj, "numMissCleavages < 2")
prj$msmsScore     <- -log10(prj$`MS-GF:SpecEValue`)
prj               <- correct_peak_selection(prj)
prj$massError     <- abs(mass_measurement_error(prj)) # ppm
fObj              <- MSnIDFilter(prj)
fObj$msmsScore    <- list(comparison=">", threshold=score.treshold)
fObj$massError    <- list(comparison="<", threshold=error.treshold)
fObj.grid         <- optimize_filter(fObj, prj, fdr.max=fdr,
                                     method="Grid", level="peptide", n.iter=iteration)
set.seed(0)
fObj.sann         <- optimize_filter(fObj.grid, prj, fdr.max=fdr,
                                     method="SANN", level="peptide", n.iter=iteration)
prj               <- apply_filter(prj, fObj.sann)

prj <- apply_filter(prj, "!grepl('Contaminant', accession)")
prj <- apply_filter(prj, "!grepl('XXX_', accession)")
#print(paste("After Filtering:",length(accessions(prj))))

####################################### QUANTIFICATION ###################################################
msnset  <- as(prj, "MSnSet")
msnset  <- combineFeatures(msnset, fData(msnset)$accession, redundancy.handler="unique", fun="sum",cv=FALSE)
exprs.table     <- exprs(msnset)
exprs.table     <- cbind(Protein=row.names(exprs.table), as.data.frame(exprs.table))
          
write.table(exprs.table, file="SpectrumCount.txt", sep="\t", row.names = F, col.names = T)
save(msnset, file = "msnset.rda")
save.image(file="scquant_results.RData")
