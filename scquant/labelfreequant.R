library(MSnbase)
library(mzID)
library(stringr)

setwd("/root/data")

####################################### READ FILE ###################################################

mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
                                full.names = F, recursive = F, ignore.case = T, include.dirs = F)
mzml.files        <- list.files(path = ".", pattern ="mzML$|MzXML", all.files = F, 
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
qnt               <- quantify(msexp, method="count", verbose=T)
colnames(qnt)     <- "Count"
agg               <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun="sum")
head(exprs(agg))

####################################### OUTPUT ###################################################
print("Writing the output...")
spectrum.count  <- as.data.frame(merge(fData(qnt)[,c("spectrum", "pepseq", "idFile", "ms-gf:evalue")], exprs(qnt), by="row.names"))
quantified      <- as.data.frame(cbind(Accession_ID=str_replace(row.names(agg),"ref\\|",""),exprs(agg)))
write.table(spectrum.count, quote=F, row.names=T, file="evalue.txt", sep ="\t")
write.table(quantified, row.names = F, quote=F, file="LabelFreeQuant.txt", sep = "\t")
#save.image(file="LabelFreeQuant.Data")

