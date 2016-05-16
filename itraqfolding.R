## Labelled (iTRAQ4) Mass Spectometry Protein Quantification
## Docker Package: p3
## Daniel Kristiyanto (daniel.kristiyanto@pnnl.gov)
#######################

library(MSnbase)


setwd("/root/data")
msnset.files        <- list.files(path = ".", pattern =".rda", all.files = F, 
                      full.names = F, recursive = F, ignore.case = T, include.dirs = F)

for (i in 1:length(msnset.files))
{
    load(msnset.files[i])
    if (i == 1 )
    {
      msn.set <- qnt.filtered
      qnt.filtered <- NULL
    } else {
      qnt.filtered <- updateFeatureNames(qnt.filtered)
      msn.set <- combine(msn.set, qnt.filtered)
      qnt.filtered <- NULL
    }
}  

