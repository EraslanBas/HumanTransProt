source('./PaperFigures/Conf.R')
source("./PaperFigures/Main.R")

numberOfTissues=29

###Order of tissues for proteinGroups_v4.txt
tissuesInOrder <- c("Tonsil", "Liver", "Spleen", "Stomach", 
                    "Brain", "Lung", "Testis", "Duodenum", 
                    "Smallintestine", "Urinarybladder", "Gallbladder", "Esophagus", 
                    "Heart", "Thyroid", "Endometrium", "Colon", 
                    "Fallopiantube", "Kidney", "Smoothmuscle", "BoneMarrow", 
                    "Prostate", "Appendices", "Pancreas", "Ovary", 
                    "Placenta", "Rectum", "Fat", "Lymphnode", 
                    "Salivarygland", "Adrenal")

protDF= as.data.frame(read.delim(paste0(tempFileDir, '/UnsharedData/proteinGroups_v4.txt'), 
                                 stringsAsFactors = FALSE))
                                
iBAQCols <- grep("iBAQ.", colnames(protDF))
ibaqDF <- protDF[,iBAQCols]
colnames(ibaqDF) <- paste0(tissuesInOrder,"_iBAQ")
ibaqDF <- ibaqDF[,sort(colnames(ibaqDF))]

### Bone marrow does not have a matching RNA sample, so it is left out
ibaqDF$BoneMarrow_iBAQ <- NULL
tissuesInOrder <- tissuesInOrder[-20]
ibaqDFRowSums <- apply(ibaqDF, 1, sum)

ibaqDF[,c("EnsemblGeneID","EnsemblProteinID","EnsemblTranscriptID")] <- protDF[,c("geneID","proteinID","transID")]


ibaqDF <- ibaqDF[ibaqDFRowSums != 0,]
ibaqDF_splitGeneID <- split(ibaqDF, ibaqDF$EnsemblGeneID)

ibaqDF_isoformsMerged <- lapply(ibaqDF_splitGeneID, function(gene){ 
  maxExpressedProtein = which.max(rowSums(gene[,1:numberOfTissues]))
  return(data.frame( gene[maxExpressedProtein,1:numberOfTissues], 
                     EnsemblGeneID=gene[maxExpressedProtein,(numberOfTissues+1)], 
                     EnsemblProteinID=gene[maxExpressedProtein,(numberOfTissues+2)], 
                     EnsemblTranscriptID=gene[maxExpressedProtein,(numberOfTissues+3)]))
})

ibaqDF_isoformsMerged <- do.call(rbind, ibaqDF_isoformsMerged)
ibaqDF_isoformsMerged[ibaqDF_isoformsMerged == 0] <- NA

temp <- log10(as.matrix(ibaqDF_isoformsMerged[,setdiff(colnames(ibaqDF_isoformsMerged),
                                                       c("EnsemblGeneID","EnsemblProteinID","EnsemblTranscriptID"))]))

#Convert IBAQ values to median centered IBAQ values(in log10) 
medianAll <- median(temp, na.rm = T)
gm <- apply(temp, 2, median, na.rm=T)
x <- sweep(temp, 2, gm, `-`)+medianAll
ibaqDF_isoformsMerged[,setdiff(colnames(ibaqDF_isoformsMerged),
                               c("EnsemblGeneID","EnsemblProteinID","EnsemblTranscriptID"))] <- x


saveRDS(ibaqDF_isoformsMerged, paste0(tempFileDir, "/UnsharedData/IBAQ_normalized.Rds"))
