## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Docker Package: p3:itraquant
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)

library(MSnbase)
library(mzID)
library(stringr)


args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args)!=4) {
  stop("Areguments ", call.=FALSE)
} else {
  evalue_treshold  	= as.double(args[1])
  pNA 				= as.numeric(args[2])
  quant_method		= args[3]
  combine_by		= args[4]
}

print(paste("evalue_treshold:", evalue_treshold))
print(paste("pNA:", pNA))
print(paste("quant_method:", quant_method))
print(paste("combine_by:"), combine_by)

setwd("/root/data")

####################################### READ FILE ###################################################
print("Isobaric Tagging Quantification")
mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)
mzml.files        <- list.files(path = ".", pattern ="mzML$|mzXML$", all.files = F, 
                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

mzids.raw         <- mzID(mzid.files)
msexp.raw         <- readMSData(mzml.files)

####################################### IDENTIFICATION & FILTERING ###################################################
print("Identifiying...")
msexp.id          <- addIdentificationData(msexp.raw, id = mzids.raw)
idSummary(msexp.id)
k                 <- (fData(msexp.id)$'ms-gf:evalue'< evalue_treshold)
k[is.na(k)]       <- FALSE
msexp.filter1     <- removeNoId(msexp.id, keep=k)
msexp.filter2     <- removeMultipleAssignment(msexp.filter1)


####################################### QUANTIFICATION ###################################################
print("Quantifying...")
qnt               <- quantify(msexp.filter2, method=quant_method, reporters=iTRAQ4, strict=F, verbose=F)
qnt.filtered      <- filterNA(qnt, pNA = pNA)
result            <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun=combine_by)

head(exprs(result))
####################################### OUTPUT ###################################################
print("Writing the output...")
spectrum.count  <- as.data.frame(merge(fData(qnt)[,c("spectrum", "pepseq", "idFile", "ms-gf:evalue")], exprs(result), by="row.names"))
quantified      <- as.data.frame(cbind(Accession_ID=str_replace(row.names(result),"ref\\|",""),exprs(result)))
write.table(spectrum.count, quote=F, row.names=F, file="evalue.txt", sep ="\t")
rm(spectrum.count)
gc()
write.table(quantified, row.names = F, quote=F, file="LabelledQuant.txt", sep = "\t")
rm(quantified)
gc()
save(result, file = "msnset.rda")
save.image(file="itraq_results.RData")


