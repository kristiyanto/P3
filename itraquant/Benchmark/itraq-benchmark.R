## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Benchmark for iTRAQ4 
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)

library(MSnbase)
library(mzID)
library(stringr)
library(dplyr)
start.time = Sys.time()

evalue_treshold = 75
pNA = 0
quant_method = trap
combine_by = mean

####################################### READ FILE FROM PNNL ###################################################

setwd("~/Documents/GITHUB/P3/itraquant/Benchmark/itraq-benchmark/PNNL")
mzml.sp        <- list.files(path = ".", pattern ="txt", all.files = F, 
                             full.names = F, recursive = F, ignore.case = T, include.dirs = F)[1]

msexp.sp <- read.table(mzml.sp, sep="\t", header = T)


####################################### READ FILE ###################################################
print("Isobaric Tagging Quantification")
setwd("~/Documents/GITHUB/P3/itraquant/Benchmark/itraq-benchmark/")

mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
                                full.names = F, recursive = F, ignore.case = T, include.dirs = F)[1]
mzml.files        <- list.files(path = ".", pattern ="mzML$|mzXML$", all.files = F, 
                                full.names = F, recursive = F, ignore.case = T, include.dirs = F)[1]

mzids.raw         <- mzID(mzid.files)
msexp.raw         <- readMSData(mzml.files)
mzid.flat         <- flatten(mzids.raw)


####################################### IDENTIFICATION & FILTERING ###################################################
print("Identifiying...")
msexp.id          <- addIdentificationData(msexp.raw, id = mzids.raw)
msexp.sp.id       <- merge(msexp.sp, mzid.flat, by.x = "ScanNumber", by.y="scan number(s)", all=T)

msexp.sp.filtered <- msexp.sp.id[!is.na(msexp.sp.id$accession),]
msexp.sp.filtered <- msexp.sp.filtered[duplicated(msexp.sp.filtered$accession) & msexp.sp.filtered$`ms-gf:evalue` <evalue_treshold,]
evalues           <- msexp.sp.filtered$`ms-gf:evalue`
####################################### IDENTIFICATION & FILTERING ###################################################

idSummary(msexp.id)
print("Filtering...")
k                 <- (fData(msexp.id)$'ms-gf:evalue'< evalue_treshold)
k[is.na(k)]       <- FALSE
msexp.filter1     <- removeNoId(msexp.id, keep=k)
msexp.filter2     <- removeMultipleAssignment(msexp.filter1)

qnt               <- quantify(msexp.filter2, method="max", reporters=iTRAQ4, strict=F, verbose=F)

qnt.id            <- merge(exprs(qnt),fData(qnt), by="row.names")

cols.sp <- c(cols, "Ion_114", "Ion_115", "Ion_116", "Ion_117")
cols.p3 <- c(cols, "iTRAQ4.114", "iTRAQ4.115", "iTRAQ4.116", "iTRAQ4.117")

msexp.sp.filtered$id_ <- paste(msexp.sp.filtered$ScanNumber, msexp.sp.filtered$pepseq, sep="_")
qnt.id$id_ <- paste(qnt.id$`scan number(s)`, qnt.id$pepseq, sep="_")
#qnt.merged        <- merge(qnt.id, msexp.sp.filtered[,c(cols.sp,"ScanNumber")], by.y = "ScanNumber", by.x = "scan number(s)", all.x =T, all.y = F)
qnt.merged <- merge(qnt.id, msexp.sp, by.x = "scan number(s)", by.y = "ScanNumber", all.x = T, all.y = F)
save.image("itraq-benchmark.RData")
