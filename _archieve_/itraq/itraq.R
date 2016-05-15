library(MSnbase)
library(mzID)
library(stringr)

setwd("/root/data")

####################################### READ FILE ###################################################
print("Isobaric Tagging Quantification")
mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)
mzml.files        <- list.files(path = ".", pattern ="mzML$|mzXML$", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

mzids.raw         <- mzID(mzid.files)
msexp.raw         <- readMSData(mzml.files)

####################################### IDENTIFICATION ###################################################
print("Identifiying...")
msexp             <- addIdentificationData(msexp.raw, id = mzids.raw)
idSummary(msexp)
msexp             <- removeNoId(msexp)
msexp             <- removeMultipleAssignment(msexp)

####################################### CLEAN UP ###################################################
rm(mzid.files)
rm(mzml.files)
rm(msexp.raw)
rm(mzids.raw)

####################################### QUANTIFICATION ###################################################
print("Quantifying...")
qnt               <- quantify(msexp, method="max", reporters=iTRAQ4, strict=F, verbose=F)
qnt               <- filterNA(qnt, pNA = 0)
agg               <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun="mean")

head(exprs(agg))
####################################### OUTPUT ###################################################
print("Writing the output...")
spectrum.count  <- as.data.frame(merge(fData(qnt)[,c("spectrum", "pepseq", "idFile", "ms-gf:evalue")], exprs(qnt), by="row.names"))
quantified      <- as.data.frame(cbind(Accession_ID=str_replace(row.names(agg),"ref\\|",""),exprs(agg)))
write.table(spectrum.count, quote=F, row.names=T, file="evalue.txt", sep ="\t")
write.table(quantified, row.names = F, quote=F, file="LabelledQuant.txt", sep = "\t")
#save.image(file="LabelledQuant.Rdata")


