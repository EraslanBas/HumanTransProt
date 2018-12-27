  source('./PaperFigures/Conf.R')
  source("./PaperFigures/Main.R")
  
  
  #Read log10 transformed median centered IBAQ values
  numberofTissues <- 29
  IBAQ_normValues <- readRDS(paste0(tempFileDir, "/UnsharedData/IBAQ_normalized.Rds"))
  IBAQ_normValues$numNAProtein <- apply(IBAQ_normValues[,1:numberofTissues],1, 
                                        function(x)sum(is.na(x)))
  
  ################################################################################
  # Reaf log10 normalized mRNA count values (median of technical replicates)
  replicateCombinedValuesDF <- readRDS(paste0(tempFileDir,"/UnsharedData/mRNA_normalized.Rds"))
  RNA_Protein_Combined <- merge(replicateCombinedValuesDF, IBAQ_normValues, by="EnsemblTranscriptID")
  
  for (tissue in TissueNames){
    noisyOnes <- which(RNA_Protein_Combined[[paste0(tissue, "_medianExon")]] <= 1)
    RNA_Protein_Combined[noisyOnes, paste0(tissue, "_medianExon")] <- NA
    RNA_Protein_Combined[noisyOnes,paste0(tissue,"_iBAQ")] <- NA
  }
  
  for (tissue in TissueNames){
    RNA_Protein_Combined[[paste0(tissue,"_PTR")]] = (RNA_Protein_Combined[[paste0(tissue,"_iBAQ")]]-
                                                       (RNA_Protein_Combined[[paste0(tissue, "_medianExon")]] ))
  }
  
  ptrCols <- grep("_PTR", colnames(RNA_Protein_Combined))
  isAllPTRNA <- apply(RNA_Protein_Combined[,ptrCols], 1, function(x){ all(is.na(x))}  )
  RNA_Protein_Combined <- RNA_Protein_Combined[which(!isAllPTRNA),]
  
  RNA_Protein_Combined$rowMedianPTR <- rowMedians(as.matrix(RNA_Protein_Combined[,ptrCols]), na.rm = TRUE)
  RNA_Protein_Combined$numNAPTR <- apply(RNA_Protein_Combined[,ptrCols],1, function(x)sum(!is.na(x)))
  
  FPKMCols <- grep("_medianExon", colnames(RNA_Protein_Combined))
  IBAQCols <- grep("_iBAQ", colnames(RNA_Protein_Combined))
  
  for(tissue in TissueNames){
    PTRNA <- which(is.na(RNA_Protein_Combined[,paste0(tissue,"_PTR")]))
    RNA_Protein_Combined[PTRNA,paste0(tissue, "_medianExon")] <- NA
    RNA_Protein_Combined[PTRNA,paste0(tissue,"_iBAQ")] <- NA
  }
  
  RNA_Protein_Combined$medianExon <- rowMedians(as.matrix(RNA_Protein_Combined[,FPKMCols]), na.rm = T)
  RNA_Protein_Combined$medianIBAQ <- rowMedians(as.matrix(RNA_Protein_Combined[,IBAQCols]), na.rm = T)
  
  ################ Get CDS sequences
  CDSSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_CDSSeq.txt"),
                             seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  CDS_transIDs <- unlist(lapply(names(CDSSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  CDS_seqs <- as.vector(unlist(lapply(CDSSequences, function(x) {x[1]} )))
  CDS_seqs <- as.vector(unlist(lapply(CDS_seqs, function(x){gsub("T","U",x)})))
  CDS_length <- unlist(lapply(CDS_seqs, nchar))
  
  CDSSeqsDF <- data.frame(CDSSeq=as.character(CDS_seqs), 
                          CDS_length=CDS_length,
                          EnsemblTranscriptID= CDS_transIDs, stringsAsFactors = FALSE)
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, CDSSeqsDF, by="EnsemblTranscriptID")
  RNA_Protein_Combined <- RNA_Protein_Combined[which(RNA_Protein_Combined$CDSSeq != "Sequence unavailable"),]
  
  
  #############################################
  ## Get UTR3 sequences
  UTR3Seqs <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_3UTRSeq.txt"),
                         seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  UTR3_transID <- unlist(lapply(names(UTR3Seqs), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  UTR3_seqs <- as.vector(unlist(lapply(UTR3Seqs, function(x) {x[1]} )))
  UTR3_seqs <- as.vector(unlist(lapply(UTR3_seqs, function(x){gsub("T","U",x)})))
  
  
  UTR3_length <- unlist(lapply(UTR3_seqs, nchar))
  UTR3SeqsDF <- data.frame(UTR3Seq=as.character(UTR3_seqs),
                           UTR3_length=UTR3_length,
                           EnsemblTranscriptID= UTR3_transID, stringsAsFactors = FALSE)
  UTR3SeqsDF <- UTR3SeqsDF[UTR3SeqsDF$UTR3Seq != "Sequence unavailable",]
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, UTR3SeqsDF, by="EnsemblTranscriptID")
  
  #############################################
  ## Get UTR5 sequences
  UTR5Seqs <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_5UTRSeq.txt"), 
                         seqtype = "DNA",
                         as.string=TRUE ,
                         forceDNAtolower=FALSE)
  UTR5_transID <- unlist(lapply(names(UTR5Seqs), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  UTR5_seqs <- as.vector(unlist(lapply(UTR5Seqs, function(x) {x[1]} )))
  UTR5_seqs <- as.vector(unlist(lapply(UTR5_seqs, function(x){gsub("T","U",x)})))
  
  UTR5_length <- unlist(lapply(UTR5_seqs, nchar))
  UTR5SeqsDF <- data.frame(UTR5Seq=as.character(UTR5_seqs),
                           UTR5_length=UTR5_length,
                           EnsemblTranscriptID= UTR5_transID,stringsAsFactors = FALSE)
  UTR5SeqsDF <- UTR5SeqsDF[UTR5SeqsDF$UTR5Seq != "Sequence unavailable",]
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, UTR5SeqsDF, by="EnsemblTranscriptID")
  
  ################################################
  ## Get promoter sequence
  # PromoterSeqs <- read.fasta(paste0(rdataFolder, "Human_Promoter_1000.txt"), seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  # PromoterSeqs_seqs <- as.vector(unlist(lapply(PromoterSeqs, function(x) {x[1]} )))
  # PromoterSeqs_transID <- unlist(lapply(names(PromoterSeqs), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  # PromoterDF <- data.frame(PromoterSeqsSeq=PromoterSeqs_seqs, EnsemblTranscriptID= PromoterSeqs_transID, stringsAsFactors = FALSE)
  # RNA_Protein_Combined <- merge(RNA_Protein_Combined, PromoterDF, by="EnsemblTranscriptID")
  # saveRDS(RNA_Protein_Combined, paste0(tempFileDir, "PromoterMotifSearchDF.Rds"))
  
  
  ###################################################
  ## Get +500 PolyA sequence
  # UTR3Down500 <- read.fasta(paste0(dataDir, "/TxtFiles/", "Human_3UTRDownstream500.txt"), seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  # UTR3Down500_seqs <- as.vector(unlist(lapply(UTR3Down500, function(x) {x[1]} )))
  # UTR3Down500_seqs <- as.vector(unlist(lapply(UTR3Down500_seqs, function(x){gsub("T","U",x)})))
  # 
  # UTR3Down500_transID <- unlist(lapply(names(UTR3Down500), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  # UTR3DownDF <- data.frame(UTR3DownSeq=UTR3Down500_seqs, EnsemblTranscriptID= UTR3Down500_transID, stringsAsFactors = FALSE)
  # RNA_Protein_Combined <- merge(RNA_Protein_Combined, UTR3DownDF, by="EnsemblTranscriptID")
  # 
  # RNA_Protein_Combined$polyADown <- apply(RNA_Protein_Combined, 1, function(x){ paste0(substr(x["UTR3Seq"], nchar(x["UTR3Seq"])-100,nchar(x["UTR3Seq"])),x["UTR3DownSeq"])})
  # 
  # RNA_Protein_Combined$polyADownLen <- unlist(lapply(RNA_Protein_Combined$polyADown, nchar))
  # RNA_Protein_Combined$UTR3DownSeq <- NULL
  
  ####Add protein sequence
  proSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_proteinSeq.fasta"),
                             seqtype = "AA", as.string=TRUE , forceDNAtolower=FALSE)
  pro_transIDs <- unlist(lapply(names(proSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  pro_seqs <- as.vector(unlist(lapply(proSequences, function(x) {substr(x[1], 1, nchar(x[1])-1)} )))
  proSeq_length <- unlist(lapply(pro_seqs, nchar))
  ProSeqsDF <- data.frame(ProSeq=as.character(pro_seqs), 
                          Pro_length=proSeq_length,
                          EnsemblTranscriptID= pro_transIDs, stringsAsFactors = FALSE)
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, ProSeqsDF, by="EnsemblTranscriptID" )
  
  
  RNA_Protein_Combined$CDSNtermSeq <- unlist(lapply(RNA_Protein_Combined$CDSSeq, function(x){ substr(x, 1,300 ) }))
  RNA_Protein_Combined$CDSNterm_length <- unlist(lapply(RNA_Protein_Combined$CDSNtermSeq, nchar))
  
  RNA_Protein_Combined$CDSCtermSeq <- unlist(lapply(RNA_Protein_Combined$CDSSeq, function(x){ substr(x, nchar(x)-300,nchar(x)) }))
  RNA_Protein_Combined$CDSCterm_length <- unlist(lapply(RNA_Protein_Combined$CDSCtermSeq, nchar))
  
  RNA_Protein_Combined$EnsemblGeneID <- NULL
  
  transIDGeneID <- readRDS(paste0(tempFileDir,"Ensembl_geneNames_ids_GRCh38.83.rds"))
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, transIDGeneID, by="EnsemblTranscriptID")
  
  saveRDS(RNA_Protein_Combined, paste0(tempFileDir, "PTR_DF_UTR5_UTR3UniqueTranscript.Rds"))
