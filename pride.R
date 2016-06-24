## Download from ProteomeXchange Repository
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)
#######################


library("rpx")
args = commandArgs(trailingOnly=TRUE)

px    <- PXDataset(args)
print(paste("Species:",pxtax(px)))
strwrap(pxref(px))
pxurl(px)
cat(pxurl(px),file="pride_url.txt")

files  <- grep("fasta", pxfiles(px), value = TRUE)

for(f in files)
{
  if(!file.exists(f)){
    file.create(paste0(f,".tmp"), showWarnings = F)
    pxget(px, f)
    file.remove(paste0(f,".tmp"))
  }
}