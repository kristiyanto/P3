## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)

library(MSnbase)
library(mzID)
library(stringr)

start.time = Sys.time()
args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args)!=6) {
  stop("Arguments ", call.=FALSE)
} else {
  mzml.files 		= args[1]
  mzid.files 		= args[2]
  evalue_treshold  	= as.double(args[3])
  pNA 				= as.numeric(args[4])
  quant_method		= args[5]
  out_file 			= args[6]
}

print("=======================================")
print(" iTRAQ4 Quantification")
print("=======================================")


print(paste("evalue_treshold:", evalue_treshold))
print(paste("pNA:", pNA))
print(paste("quant_method:", quant_method))


setwd("/root/data")

####################################### READ FILE ###################################################
print("Isobaric Tagging Quantification")
#mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
#                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)
#mzml.files        <- list.files(path = ".", pattern ="mzML$|mzXML$", all.files = F, 
#                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

mzids.raw         <- mzID(mzid.files)
msexp.raw         <- readMSData(mzml.files)

####################################### IDENTIFICATION & FILTERING ###################################################
print("Identifiying...")
msexp.id          <- addIdentificationData(msexp.raw, id = mzids.raw)
idSummary(msexp.id)
rm(mzid.files)
rm(mzml.files)
rm(mzids.raw)
rm(msexp.raw)
gc(verbose = FALSE)

print("Filtering...")
k                 <- (fData(msexp.id)$'ms-gf:evalue'< evalue_treshold)
k[is.na(k)]       <- FALSE
msexp.filter1     <- removeNoId(msexp.id, keep=k)
msexp.filter2     <- removeMultipleAssignment(msexp.filter1)

rm(msexp.id)
rm(msexp.filter1)
gc(verbose = FALSE)
####################################### QUANTIFICATION ###################################################
print("Quantifying...")
qnt               <- quantify(msexp.filter2, method=quant_method, reporters=iTRAQ4, strict=F, verbose=F)
qnt.filtered      <- filterNA(qnt, pNA = pNA)
rm(msexp.filter2)
rm(qnt)
gc(verbose = FALSE)
save(qnt.filtered, file=out_file)
#result            <- combineFeatures(qnt.filtered, groupBy = fData(qnt.filtered)$accession, fun=combine_by)

#head(exprs(result))
####################################### OUTPUT ###################################################
print("Writing the output...")
#evalue.table	  <- as.data.frame(merge(fData(qnt.filtered)[,c("accession","spectrum", "pepseq", "idFile", "ms-gf:evalue")], exprs(result), by.x ="accession", by.y ="row.names" ))
#quantified		  <- as.data.frame(cbind(Accession_ID=row.names(result),exprs(result)))
#write.table(evalue.table, quote=F, row.names=F, file="evalue.txt", sep ="\t")
#rm(evalue.table)
#gc()
#write.table(quantified, row.names = F, quote=F, file="LabelledQuant.txt", sep = "\t")
#save(result, file = "msnset.rda")
stop.time = Sys.time()
#rm(result)
print(paste("Start Time:", start.time," Stop Time:", stop.time, "Elapsed=",stop.time-start.time))

