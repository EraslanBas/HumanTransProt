source('./PaperFigures/Conf.R')
source("./PaperFigures/Main.R")

# gtf <- import(paste0(tempFileDir, "/UnsharedData/Homo_sapiens.GRCh38.83.gtf"), format="gtf")
# gtf <- gtf[mcols(gtf)$type=="exon",]
# 
# gr.list.gene <- split(gtf, as.factor(gtf$transcript_id))
# gr.list.gene.unified <- GRangesList(lapply(gr.list.gene, function(geneExons){
#   reducedExons <- reduce(geneExons, ignore.strand=FALSE)
#   return(reducedExons)
# }))
# 
# 
# unifiedEWidth <- as.list(width(gr.list.gene.unified))
# unifiedEWidths <- unlist(lapply(unifiedEWidth, sum))
# unifiedEWidths <- data.frame(totalExonLength = unifiedEWidths, transID= names(unifiedEWidths))
# saveRDS(unifiedEWidths, paste0(tempFileDir,"/UnsharedData/unifiedEWidthsTranscript_Human.Rds"))



###################################################
### Get exonic read counts only
allTissuesExon = list()
for (tissue in TissueNames){
  allTissuesExon[[tissue]] = updateObject(readRDS(paste0(tempFileDir,"/UnsharedData/TranscriptExonCounts/", 
                                                         tissue, "UniqueTranscriptReadCounts.Rds")))
}

allExonCountsCombined <- data.frame(assay(allTissuesExon[[1]]), 
                                    transID=rownames(assay(allTissuesExon[[1]])), 
                                    stringsAsFactors = FALSE)

for( elem in names(allTissuesExon[-1])){
  temp= data.frame(assay(allTissuesExon[[elem]]),
                   transID=rownames(assay(allTissuesExon[[elem]])),
                   stringsAsFactors = FALSE)
  allExonCountsCombined <- merge(allExonCountsCombined, temp, by="transID")
} 

rownames(allExonCountsCombined) <- allExonCountsCombined$transID
colnames(allExonCountsCombined)[-1] <- paste0(unlist(strsplit(colnames(allExonCountsCombined)[-1], "_.bam")), 
                                              "_exonCount")


#### Estimate size factors from the exon counts
sf <- estimateSizeFactorsForMatrix( allExonCountsCombined[,-1])

###################################################
#### Get whole transcript read counts

allTissuesTranscript = list()
for (tissue in TissueNames){
  allTissuesTranscript[[tissue]] = updateObject(readRDS(paste0(tempFileDir, "/UnsharedData/TranscriptCounts/",
                                                               tissue, "AllTranscriptReadCounts.Rds")))
}

allTranscriptCountsCombined <- data.frame(assay(allTissuesTranscript[[1]]), 
                                          transID=rownames(assay(allTissuesTranscript[[1]])),
                                          stringsAsFactors = FALSE)

for( elem in names(allTissuesTranscript[-1])){
  temp= data.frame(assay(allTissuesTranscript[[elem]]),
                   transID=rownames(assay(allTissuesTranscript[[elem]])),
                   stringsAsFactors = FALSE)
  allTranscriptCountsCombined <- merge(allTranscriptCountsCombined, temp, by="transID")
}

rownames(allTranscriptCountsCombined) <- allTranscriptCountsCombined$transID
colnames(allTranscriptCountsCombined)[-1] <- paste0(unlist(strsplit(colnames(allTranscriptCountsCombined)[-1], "_.bam")),
                                                    "_transcriptCount")


####################################################
#### Get intron counts
allTranscriptCountsCombined <- merge(allTranscriptCountsCombined, allExonCountsCombined, by='transID')

replicateTissueColNames = c()

for(elem in colnames(allExonCountsCombined)[-1]){
  tissueN <- strsplit(elem, "_")[[1]][1]
  replicateN <- strsplit(elem, "_")[[1]][2]
  tmp <- (allTranscriptCountsCombined[[paste0(tissueN,"_",replicateN,"_transcriptCount")]] - 
            allTranscriptCountsCombined[[paste0(tissueN,"_",replicateN,"_exonCount")]])
  tmp[tmp < 1] <- 0
  tmp[is.na(tmp)] <- 0
  allTranscriptCountsCombined[[paste0(tissueN,"_",replicateN,"_intronCount")]]<- tmp
  replicateTissueColNames <- c(replicateTissueColNames, paste0(tissueN,"_",replicateN))
}

###########################################################
##################################################
## Get intron total lengths per transcript

transcriptGtf <- import(paste0(tempFileDir, "/UnsharedData/Homo_sapiens.GRCh38.83.gtf"), format="gtf")
transcriptGtf <- transcriptGtf[mcols(transcriptGtf)$type=="transcript",]

#tp <- data.frame(GeneName=mcols(transcriptGtf)$gene_name,EnsemblGeneID=mcols(transcriptGtf)$gene_id,
#           EnsemblTranscriptID=mcols(transcriptGtf)$transcript_id)
#saveRDS(tp, paste0(tempFileDir,"Ensembl_geneNames_ids_GRCh38.83.rds"))

transcriptWidths <- data.frame(transcriptsWidth = width(transcriptGtf), 
                               transID= mcols(transcriptGtf)$transcript_id)

unifiedEWidths <- readRDS(paste0(tempFileDir,"/UnsharedData/unifiedEWidthsTranscript_Human.Rds"))
unifiedEWidths <- merge(unifiedEWidths, transcriptWidths, by="transID")
unifiedEWidths$intronLength <- (unifiedEWidths$transcriptsWidth - unifiedEWidths$totalExonLength)

allTranscriptCountsCombined <- merge(allTranscriptCountsCombined, unifiedEWidths, by="transID")

############################################################
##Get exonic length normalized counts
exonCountCols <- grep("_exonCount", colnames(allTranscriptCountsCombined))
allTranscriptCountsCombined[,exonCountCols] <- apply(allTranscriptCountsCombined[,exonCountCols], 2,  
                                                     function(x){ return((x * 10^3) / allTranscriptCountsCombined$totalExonLength)})


##Get intronic length normalized counts
intronCountCols <- grep("_intronCount", colnames(allTranscriptCountsCombined))
allTranscriptCountsCombined[,intronCountCols] <- apply(allTranscriptCountsCombined[,intronCountCols], 2,  
                                                       function(x){ return((x * 10^3)/ allTranscriptCountsCombined$intronLength)})
allTranscriptCountsCombined[is.na(allTranscriptCountsCombined)] <- NA

saveRDS(allTranscriptCountsCombined, paste0(tempFileDir,"/UnsharedData/allTranscriptCountsCombined.rds"))
# ################################################################
# Normalize exonic counts wrt intronic counts to get the mature mRNA counts
for(elem in replicateTissueColNames){

  intronNARows <- which(is.na(allTranscriptCountsCombined[, paste0(elem, "_intronCount")]))
  allTranscriptCountsCombined[intronNARows , paste0(elem, "_intronCount")] <- 0
  exon0Rows <- which(allTranscriptCountsCombined[, paste0(elem, "_exonCount")]  == 0)
  allTranscriptCountsCombined[exon0Rows, paste0(elem, "_exonCount")] <- NA
  
  allTranscriptCountsCombined[, paste0(elem, "_exonCount")] <- allTranscriptCountsCombined[, paste0(elem, "_exonCount")]- 
                                                                       allTranscriptCountsCombined[, paste0(elem, "_intronCount")]
}

############################################################

exonCols <- grep("_exonCount", colnames(allTranscriptCountsCombined))

## Normalize read counts wrt library size
allTranscriptCountsCombined[,exonCols] <- log10(as.data.frame((t(t(as.matrix(allTranscriptCountsCombined[,exonCols]))/sf ))))

onlyExon <- allTranscriptCountsCombined[,exonCols]


##Get median exonic value per tissue
for(elem in TissueNames){
  replicateCols <- grep(elem, colnames(onlyExon))
  allTranscriptCountsCombined[[paste0(elem,"_medianExon")]] <- rowMedians(as.matrix(onlyExon[,replicateCols]), na.rm = TRUE)
}


medianCols <- grep("_medianExon", colnames(allTranscriptCountsCombined))
allCountsCombined <- allTranscriptCountsCombined[,medianCols] 
allCountsCombined$EnsemblTranscriptID <-  allTranscriptCountsCombined$transID
allCountsCombined$totalExonLength <-  allTranscriptCountsCombined$totalExonLength

saveRDS(allCountsCombined, paste0(tempFileDir,"/UnsharedData/mRNA_normalized.Rds"))
