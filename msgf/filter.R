library("rpx")
library("mzID")
library("MSnID")
library(MSnbase)

if (length(args)!=3) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
  mz.file     = args[1]
  mzml.file   = args[2]
  output.file = args[3]
}

score.treshold  = 7.0
error.treshold  = 20
fdr             = 0.01
iteration       = 5000

setwd("/root/data/")

####################################### FILTERING ###################################################

prj               <- MSnID()
prj               <- read_mzIDs(object = prj, mzids = mz.file)
prj               <- assess_termini(prj)
prj               <- assess_missed_cleavages(prj)
show(prj)
print(paste("Before Filtering:",length(accessions(prj))))
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

####################################### IDENTIFICATION ###################################################
msexp.raw         <- readMSData(mzml.file)
msexp             <- addIdentificationData(msexp.raw, id = prj)
idSummary(msexp)
msexp             <- removeNoId(msexp)
msexp             <- removeMultipleAssignment(msexp)

save(prj, file = "OUT.RData")
