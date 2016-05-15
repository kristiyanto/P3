# Derrived from PNNL-Comp-Mass-Spec/celegans.aging


setwd("~/Documents/GITHUB/P3/scquant/Benchmark/")


################################# ORIGINAL SCRIPT ######################################################
library("rpx")
id <- "PXD002161"
px <- PXDataset(id)
try(setInternet2(FALSE),silent=TRUE)
library("jsonlite")
addr <- "http://www.ebi.ac.uk:80/pride/ws/archive/%s/list/project/%s"
files <- fromJSON(sprintf(addr, "file", id))$list
assays <- fromJSON(sprintf(addr, "assay", id))$list
files <- subset(files, grepl('mzid.gz', fileName),
                select=c("assayAccession","fileName"))
rownames(files) <- NULL
assays <- assays[,c("assayAccession",
                    "experimentalFactor",
                    "proteinCount",
                    "peptideCount",
                    "uniquePeptideCount",
                    "identifiedSpectrumCount",
                    "totalSpectrumCount")]
pttrn <- "Diet; fully fed, Mut: glp-4(bn2ts) daf-16(mgDf50); daf-2(e1370)"
assays <- subset(assays, grepl(pttrn, experimentalFactor, fixed=T))
assays <- with(assays, {data.frame(assayAccession,
                                   phenotype=sub("Age:\\s([^,]*).*", "\\1",
                                                 experimentalFactor),
                                   sampleName=sub(".*Name:\\s(\\w)-ctrl-FF.(\\d)", "\\1.\\2",
                                                  experimentalFactor),
                                   stringsAsFactors=F)})
files <- subset(files, assayAccession %in% assays$assayAccession)
pxget(px, files$fileName) # fetch them from PRIDE
files$datasetName <- sub('_msgfplus.mzid.gz','',files$fileName, fixed=T)
meta <- merge(files[,c("assayAccession","datasetName")], assays)
rownames(meta) <- meta$datasetName
meta <- meta[order(meta$sampleName),]
rownames(meta) <- NULL


library("MSnID")
prj <- MSnID()
prj <- read_mzIDs(object = prj, mzids = files$fileName)
prj <- assess_termini(prj)
prj <- assess_missed_cleavages(prj)
prj <- apply_filter(prj, "numIrregCleavages == 0")
prj <- apply_filter(prj, "numMissCleavages < 2")
prj$msmsScore <- -log10(prj$`MS-GF:SpecEValue`)
prj <- correct_peak_selection(prj)
prj$massError <- abs(mass_measurement_error(prj)) # ppm
fObj <- MSnIDFilter(prj)
fObj$msmsScore <- list(comparison=">", threshold=7.0)
fObj$massError <- list(comparison="<", threshold=20)
fObj.grid <- optimize_filter(fObj, prj, fdr.max=0.01,
                             method="Grid", level="peptide", n.iter=5000)
set.seed(0)
fObj.sann <- optimize_filter(fObj.grid, prj, fdr.max=0.01,
                             method="SANN", level="peptide", n.iter=5000)
prj <- apply_filter(prj, fObj.sann)


prj <- apply_filter(prj, "!grepl('Contaminant', accession)")
prj <- apply_filter(prj, "!grepl('XXX_', accession)")

###################### SKIP ENRICHMENT ###################
# 
# library("biomaRt")
# ens <- useMart("ensembl")
# ens <- useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
# ens.cel <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="celegans_gene_ensembl", host = "jul2015.archive.ensembl.org")
# 
# 
# ens.cel <- useDataset("celegans_gene_ensembl")
# wormpep2entrez    <- getBM(attributes=c("wormpep_id","entrezgene"),
#                            filters=c("wormpep_id"),
#                            values=proteins(prj), mart=ens.cel)
# # remove non-annotated proteins
# wormpep2entrez <- subset(wormpep2entrez, !is.na(entrezgene))
# # reorder and retain just the first EntrezGeneID if protein match multiple genes
# wormpep2entrez <- wormpep2entrez[with(wormpep2entrez,
#                                       order(wormpep_id, entrezgene)),]
# wormpep2entrez <- wormpep2entrez[with(wormpep2entrez,
#                                       !duplicated(wormpep_id)),]
# new.psms <- merge(psms(prj), wormpep2entrez, by.x="accession", by.y="wormpep_id")
# new.psms$wormpep_id <- new.psms$accession
# new.psms$accession <- new.psms$entrezgene
# psms(prj) <- new.psms




suppressPackageStartupMessages(library("MSnbase"))
mset <- as(prj, "MSnSet")
# take care of pheno data
sampleNames(mset) <- sub('.mzML', '', sampleNames(mset), fixed=T)
rownames(meta) <- meta$datasetName
pData(mset) <- meta[sampleNames(mset),]
sampleNames(mset) <- pData(mset)$sampleName
# roll-up to gene level
mset <- combineFeatures(mset, groupBy = fData(mset)$accession,
                        fun = 'sum', cv=F, redundancy.handler = 'unique')

Ground.Truth <- exprs(mset)

write.table(Ground.Truth, file="Celegans-Aging.txt", sep="\t")

################################# COMPARE IT WITH P3 RESULTS ######################################################

To.Compare   <- as.data.frame(read.table(file="SpectrumCount.txt",  stringsAsFactors=FALSE))
colnames(To.Compare) <- To.Compare[1,]
row.names(To.Compare) <- To.Compare[,1]
To.Compare <- To.Compare[-1,-1]

Is.Match <- Ground.Truth == To.Compare

print(paste("Total Unmatch:", toString(sum(Is.Match==F))))
