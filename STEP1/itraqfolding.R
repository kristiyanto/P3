## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)
#######################

library(MSnbase)


setwd("/root/data")
msnset.files        <- list.files(path = ".", pattern ="rda", all.files = F, 
                      full.names = F, recursive = F, ignore.case = T, include.dirs = F)

for (i in 1:length(msnset.files))
{
    load(msnset.files[i])
    curr <- paste0(substr(msnset.files[i],1,4), i)
    #sampleNames(qnt.filtered) <- curr
    #qnt.filtered < updateFvarLabels(qnt.filtered, curr)
    dim(qnt.filtered)
    if (i == 1 )
    {
      msn.set <- qnt.filtered
      qnt.filtered <- NULL
    } else {
      msn.set <- combine(msn.set, qnt.filtered)
      qnt.filtered <- NULL
    }
}  

