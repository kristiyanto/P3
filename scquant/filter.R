library("rpx")
library("mzID")
library("MSnID")
library("MSnbase")

if (length(args)!=4) {
  stop("Areguments ", call.=FALSE)
} else if (length(args)==1) {
  score.treshold  = args[1]
  error.treshold  = args[2]
  fdr             = args[3]
  iteration       = args[4]
}

setwd("/root/data/")

mz.files <- list.files(path = ".", pattern ="mzid", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

####################################### FILTERING ###################################################

prj               <- MSnID(".")
prj               <- read_mzIDs(object = prj, mzids = mz.files)
prj               <- assess_termini(prj)
prj               <- assess_missed_cleavages(prj)
show(prj)
#print(paste("Before Filtering:",length(accessions(prj))))
prj               <- apply_filter(prj, "numIrregCleavages == 0")
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
print(paste("After Filtering:",length(accessions(prj))))

####################################### QUANTIFICATION ###################################################
msnset  <- combineFeatures(as(prj, "MSnSet"), fData(out)$accession, redundancy.handler="unique", fun="sum",cv=FALSE)
rm(prj)
gc()
tmp     <- exprs(msnset)
tmp     <- tmp[tmp>0]
save.table(tmp, file="SpectrumCount.txt", sep="\t")
save(msnset, file = "OUT.RData")
