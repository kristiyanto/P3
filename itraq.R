## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)
suppressMessages(library(mzID))
suppressMessages(library(MSnbase))
suppressMessages(library(stringr))
suppressMessages(library(BiocParallel))

BiocParallel::register(BiocParallel::SnowParam(2))
BiocParallel::register(BiocParallel::MulticoreParam(2))

start.time = Sys.time()
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Arguments.", call.=FALSE)
} else {
  mzml.files 		= args[1]
  mzid.files 		= args[2]
  specvalue_threshold  	= as.double(args[3])
  pNA 				  = as.numeric(args[4])
  quant_method	= args[5]
  combine_by    = args[6]
  out_file 			= args[7]
}


print("=======================================")
print(" Quantification")
if(quant_method=="count")
{
  print("Spectrum Count Quantification")
  if(specvalue_threshold > 0)
  {
    print(paste("Spec E-Value Threshold:", specvalue_threshold))
  }
  if(toupper(combine_by) != "SKIP")
  {
    print(paste("features combined using:", combine_by))
  }
}else{
  print("iTRAQ4 Quantification")
  if(specvalue_threshold > 0)
  {
    print(paste("Spec E-Value Threshold:", specvalue_threshold))
  }
  print(paste("pNA:", pNA))
  print(paste("quant_method:", quant_method))
  if(toupper(combine_by) != "SKIP")
  {
    print(paste("features combined using:", combine_by))
  }
}
print("=======================================")

setwd("/root/data")

####################################### READ FILE ###################################################
# mzid.files        <- list.files(path = ".", pattern ="mzid", all.files = F, 
#                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)
# mzml.files        <- list.files(path = ".", pattern ="mzML$|mzXML$", all.files = F, 
#                       full.names = F, recursive = F, ignore.case = T, include.dirs = F)

#mzids.raw         <- mzID(mzid.files)
msexp.raw         <- readMSData(mzml.files)

####################################### IDENTIFICATION & FILTERING ###################################################
print("Identifying...")
msexp.id          <- addIdentificationData(msexp.raw, id = mzid.files,verbose=F)
idSummary(msexp.id)
rm(mzid.files)
rm(mzml.files)
rm(msexp.raw)
gc(verbose = FALSE)

if(specvalue_threshold==0)
{
  msexp.filter1     <- msexp.id
} else {
  print("Filtering...")
  k                 <- (fData(msexp.id)$'ms-gf:specevalue'< specvalue_threshold)
  k[is.na(k)]       <- FALSE
  msexp.filter1     <- removeNoId(msexp.id, keep=k)
} 

rm(msexp.id)
gc(verbose = FALSE)
####################################### QUANTIFICATION ###################################################
print("Quantifying...")
if (quant_method=="count")
{
  qnt               <- quantify(msexp.filter1, method=quant_method)
  qnt.filtered      <- qnt
} else
{
  qnt               <- quantify(msexp.filter1, method=quant_method, reporters=iTRAQ4, strict=F, verbose=F)
  qnt.filtered      <- filterNA(qnt, pNA = pNA)
}

rm(msexp.filter1)
rm(qnt)
gc(verbose = FALSE)

if(toupper(combine_by) == "SKIP")
{
  result  <- qnt.filtered
} else  
{
  result  <- combineFeatures(qnt.filtered, groupBy = fData(qnt.filtered)$accession, fun=combine_by)
}
save(result, file=out_file)
####################################### OUTPUT ###################################################
print("Writing the output...")
if (quant_method=="count")
{
  Cnt=exprs(result)
  colnames(Cnt) = "Count"
  quantified		  <- as.data.frame(cbind(Scan=fData(result)$'scan number(s)', Spectrum=fData(result)$spectrum, Spec_Evalue=fData(result)$'ms-gf:specevalue', AccessionID=row.names(result), PepSeq=fData(result)$pepseq, Cnt))
  
} else
{
  quantified		  <- as.data.frame(cbind(Scan=fData(result)$'scan number(s)', Spectrum=fData(result)$spectrum, Spec_Evalue=fData(result)$'ms-gf:specevalue', AccessionID=row.names(result), PepSeq=fData(result)$pepseq, exprs(result)))
}

write.table(quantified, row.names = F, quote=F, file=str_replace(out_file,".rda",".txt"), sep = "\t")
stop.time = Sys.time()

print(paste("File:", out_file,"Start Time:", start.time," Stop Time:", stop.time, "Elapsed=",stop.time-start.time))
rm(result)
rm(qnt.filtered)
