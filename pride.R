## Download from ProteomeXchange Repository
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)
#######################


library("rpx")
print(args)

args = commandArgs(trailingOnly=TRUE)

px    <- PXDataset(args)
files  <- grep("mzXML|fasta|mzml|mzML", pxfiles(px), value = TRUE)

for(f in files)
{
  if(!file.exists(f)){
    file.create(paste0(f,".tmp"), showWarnings = F)
    pxget(px, f)
    file.remove(paste0(f,".tmp"))
  }
}