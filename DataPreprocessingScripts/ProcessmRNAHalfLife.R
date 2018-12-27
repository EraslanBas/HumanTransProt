  source('./PaperFigures/Conf.R')
  source("./PaperFigures/Main.R")
  
  ##### Hela Tet-off cells mRNA half-life
  ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3337439/
  
  refSeqEnsembl <- as.data.frame(read.table(paste0(tempFileDir,"/UnsharedData/RefSeq_Ensembl.txt"),
                                            sep = "\t", header = T, stringsAsFactors = F))
  colnames(refSeqEnsembl) <- c("EnsemblGeneID", "EnsemblTranscriptID", "RefSeqID")
  
  helaMRNAHL <- read.csv(paste0(tempFileDir,"Human_HeLamRNA_HL.csv"), stringsAsFactors = F)
  helaMRNAHL$RepName <- sapply(helaMRNAHL$RepName, function(x){strsplit(x,",")[[1]][1]})
  helaMRNAHL<- helaMRNAHL[,c(1,3,4)]
  colnames(helaMRNAHL) <- c("RefSeqID", "RPKM","HeLamRNAHL")
  
  helaMRNAHL <- merge(helaMRNAHL, refSeqEnsembl, by="RefSeqID")
  helaMRNAHL <- helaMRNAHL[helaMRNAHL$HeLamRNAHL != "N.D.",]
  helaMRNAHL[helaMRNAHL$HeLamRNAHL == ">24", "HeLamRNAHL"] <- 24
  helaMRNAHL$HeLamRNAHL <- as.numeric(helaMRNAHL$HeLamRNAHL)
  
  helaMRNAHL$HeLamRNAHL <- log10(helaMRNAHL$HeLamRNAHL*60)
  helaMRNAHL$RPKM <- as.numeric(helaMRNAHL$RPKM)
  helaMRNAHL$RPKM <- log10(helaMRNAHL$RPKM+1)
  
  CDSSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_CDSSeq.txt"),
                             seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  CDS_transIDs <- unlist(lapply(names(CDSSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  CDS_seqs <- as.vector(unlist(lapply(CDSSequences, function(x) {x[1]} )))
  CDS_seqs <- as.vector(unlist(lapply(CDS_seqs, function(x){gsub("T","U",x)})))
  CDS_length <- unlist(lapply(CDS_seqs, nchar))
  
  CDSSeqsDF <- data.frame(CDSSeq=as.character(CDS_seqs), 
                          CDS_length=CDS_length,
                          EnsemblTranscriptID= CDS_transIDs, stringsAsFactors = FALSE)
  
  helaMRNAHL <- merge(helaMRNAHL, CDSSeqsDF, by="EnsemblTranscriptID")
  helaMRNAHL <- helaMRNAHL[which(helaMRNAHL$CDSSeq != "Sequence unavailable"),]
  
  proSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_proteinSeq.fasta"),
                             seqtype = "AA", as.string=TRUE , forceDNAtolower=FALSE)
  pro_transIDs <- unlist(lapply(names(proSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  pro_seqs <- as.vector(unlist(lapply(proSequences, function(x) {substr(x[1], 1, nchar(x[1])-1)} )))
  proSeq_length <- unlist(lapply(pro_seqs, nchar))
  ProSeqsDF <- data.frame(ProSeq=as.character(pro_seqs), 
                          Pro_length=proSeq_length,
                          EnsemblTranscriptID= pro_transIDs, stringsAsFactors = FALSE)
  
  helaMRNAHL <- merge(helaMRNAHL, ProSeqsDF, by="EnsemblTranscriptID")
  
  
  CodonFreq <- generateCodonFrequency(seqs=helaMRNAHL$CDSSeq)
  colnames(CodonFreq) <- paste0("Codon_", colnames(CodonFreq))
  helaMRNAHL <- cbind(helaMRNAHL, CodonFreq)
  
  
  AAFreq <- generateAminoAcidFrequency(seqs=helaMRNAHL$ProSeq)
  colnames(AAFreq) <- paste0("AminoAcid_", colnames(AAFreq))
  helaMRNAHL <- cbind(helaMRNAHL, AAFreq)
  
  
  rbpPeaksFile <- as.data.frame(readRDS(paste0(tempFileDir,"/peak_center-gene_mapping.rds")))
  rbpPeaksFile <- rbpPeaksFile[rbpPeaksFile$binding_site,]
  rbpPeaksFile <- separate(rbpPeaksFile,
                           col="gene_id",
                           into=c("geneID"),
                           sep=c("\\."),
                           extra = "drop")
  
  rbpPeaksFileSplit <- split(rbpPeaksFile, f = rbpPeaksFile$rbp )
  rbpNames <- names(rbpPeaksFileSplit)
  
  
  for (i in 1:length(rbpPeaksFileSplit)){
    helaMRNAHL[,paste0("RBP_",rbpNames[i])] <- 0
    helaMRNAHL[helaMRNAHL$EnsemblGeneID %in% rbpPeaksFileSplit[[i]]$geneID,paste0("RBP_",rbpNames[i])] <- 1
  }
  rbpPeaksFile<- NULL
  rbpPeaksFileSplit <- NULL
  gc()
  
  geneNameID <- readRDS(paste0(tempFileDir,
                               "/UnsharedData/geneID_vs_Name.rds"))
  colnames(geneNameID) <- c("EnsemblGeneID", "GeneName")
  helaMRNAHL <- merge(helaMRNAHL, geneNameID, by="EnsemblGeneID")
  
  
  miRNAtargetGeneCutOff=200
  
  miRNATargetSites <- read.csv(file=paste0(tempFileDir, 'hsa_MTI.csv'), header=TRUE, sep=",")
  miRNATargetSites[,c("Experiments",
                      "Support.Type",
                      "References..PMID.",
                      "Species..miRNA.",
                      "Target.Gene..Entrez.Gene.ID.",
                      "Species..Target.Gene.")] <- NULL
  miRNATargetSites <-unique(miRNATargetSites)
  
  expVarmiRNAs <- split(miRNATargetSites, miRNATargetSites$miRNA)
  names(expVarmiRNAs) <- tolower(names(expVarmiRNAs))
  
  numberOfTargets <- sapply(expVarmiRNAs, nrow)
  expVarmiRNAs <- expVarmiRNAs[which(numberOfTargets>miRNAtargetGeneCutOff)]
  
  for(elem in expVarmiRNAs){
    miRNA <- elem[1,"miRNA"]
    miRNA <- gsub("/", "_", miRNA)
    miRNA <- gsub("-", "_", miRNA)
    #miRNA <- gsub("[.]", "", miRNA)
    helaMRNAHL[,paste0("sgn_MiRNA_",miRNA)] <- 0
    helaMRNAHL[helaMRNAHL$GeneName %in% elem$Target.Gene, 
                paste0("sgn_MiRNA_",miRNA)] <- 1
  }
  
  miRNATargetSites <- NULL
  expVarmiRNAs <- NULL
  gc()
  
  
  saveRDS(helaMRNAHL, paste0(tempFileDir, "/UnsharedData/humanHELAMRNAHL.rds"))
  
  #####################
  #### K562 cells mRNA half-life
  #### https://www.ncbi.nlm.nih.gov/pubmed/27257258
  
  k562_halfLife <- as.data.frame(readRDS(paste0(tempFileDir,
                                                "/UnsharedData/k562_features.rds")))
  
  x <- grep("Codon|AminoAcid|K562_FPKM", colnames(k562_halfLife))
  k562_halfLife[,x] <-NULL
  
  # geneExp <- read.table(paste0(tempFileDir,"K562_geneExpression.tsv"), header = T)
  # geneExp <- geneExp[,c("gene_id","FPKM")]
  # colnames(geneExp) <- c("EnsemblGeneID", "K562_FPKM")
  # geneExp$EnsemblGeneID <- sapply(geneExp$EnsemblGeneID, function(x){strsplit(x,"\\.")[[1]][1]})
  # geneExp$K562_FPKM <- log10(geneExp$K562_FPKM+1)
  # 
  # k562_halfLife <- merge(k562_halfLife, geneExp, by="EnsemblGeneID")
  # 
  
  CodonFreq <- generateCodonFrequency(seqs=k562_halfLife$CDSSeq)
  colnames(CodonFreq) <- paste0("Codon_", colnames(CodonFreq))
  k562_halfLife <- cbind(k562_halfLife, CodonFreq)
  
  AAFreq <- generateAminoAcidFrequency(seqs=k562_halfLife$ProSeq)
  colnames(AAFreq) <- paste0("AminoAcid_", colnames(AAFreq))
  k562_halfLife <- cbind(k562_halfLife, AAFreq)
  
  saveRDS(k562_halfLife, paste0(tempFileDir,
                                "/UnsharedData/k562_features.rds"))
  
  
  ##################
  #### Hek293 cells mRNA half-life
  ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49831
  
  Hek293_mRNAHL <- read.table(paste0(tempFileDir, 
                                     "GSE49831_HEK293_halflives.txt"), header = T)
  Hek293_mRNAHL$GeneName <- rownames(Hek293_mRNAHL)
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "PTR_DF_UTR5_UTR3UniqueTranscript.Rds"))
  
  Hek293_mRNAHL <- merge(Hek293_mRNAHL, RNA_Protein_Combined, by="GeneName")
  Hek293_mRNAHL$Hek293Cells_mRNAHL <- log10((Hek293_mRNAHL$HEK293_RNA.half.life..min..rep1 + Hek293_mRNAHL$HEK293_RNA.half.life..min..rep2)/2)
  
  CodonFreq <- generateCodonFrequency(seqs=Hek293_mRNAHL$CDSSeq)
  colnames(CodonFreq) <- paste0("Codon_", colnames(CodonFreq))
  Hek293_mRNAHL <- cbind(Hek293_mRNAHL, CodonFreq)
  
  
  AAFreq <- generateAminoAcidFrequency(seqs=Hek293_mRNAHL$ProSeq)
  colnames(AAFreq) <- paste0("AminoAcid_", colnames(AAFreq))
  Hek293_mRNAHL <- cbind(Hek293_mRNAHL, AAFreq)
  
  
  saveRDS(Hek293_mRNAHL, paste0(tempFileDir,
                                "/UnsharedData/Hek293_mRNAHL.rds"))
  
  ##############################
  lHL <- read.csv(paste0(tempFileDir,"LeoMRNAHL.tsv"), sep = "\t")
  lHL$EnsemblTranscriptID <- sapply(lHL$transcript_id , function(x){strsplit(x, "\\.")[[1]][1]})
  
  CDSSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_CDSSeq.txt"),
                             seqtype = "DNA", as.string=TRUE , forceDNAtolower=FALSE)
  CDS_transIDs <- unlist(lapply(names(CDSSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  CDS_seqs <- as.vector(unlist(lapply(CDSSequences, function(x) {x[1]} )))
  CDS_seqs <- as.vector(unlist(lapply(CDS_seqs, function(x){gsub("T","U",x)})))
  CDS_length <- unlist(lapply(CDS_seqs, nchar))
  
  CDSSeqsDF <- data.frame(CDSSeq=as.character(CDS_seqs), 
                          CDS_length=CDS_length,
                          EnsemblTranscriptID= CDS_transIDs, stringsAsFactors = FALSE)
  
  lHL <- merge(lHL, CDSSeqsDF, by="EnsemblTranscriptID")
  lHL <- lHL[which(lHL$CDSSeq != "Sequence unavailable"),]
  
  proSequences <- read.fasta(paste0(tempFileDir, "/UnsharedData/Human_proteinSeq.fasta"),
                             seqtype = "AA", as.string=TRUE , forceDNAtolower=FALSE)
  pro_transIDs <- unlist(lapply(names(proSequences), function(x) {regmatches(x, regexec('(ENST[0-9]+)',x))[[1]][1]} ))
  pro_geneIDs <- unlist(lapply(names(proSequences), function(x) {regmatches(x, regexec('(ENSG[0-9]+)',x))[[1]][1]} ))
  
  pro_seqs <- as.vector(unlist(lapply(proSequences, function(x) {substr(x[1], 1, nchar(x[1])-1)} )))
  proSeq_length <- unlist(lapply(pro_seqs, nchar))
  ProSeqsDF <- data.frame(ProSeq=as.character(pro_seqs), 
                          Pro_length=proSeq_length,
                          EnsemblTranscriptID= pro_transIDs, stringsAsFactors = FALSE)
  
  lHL <- merge(lHL, ProSeqsDF, by="EnsemblTranscriptID")
  
  CodonFreq <- generateCodonFrequency(seqs=lHL$CDSSeq)
  colnames(CodonFreq) <- paste0("Codon_", colnames(CodonFreq))
  lHL <- cbind(lHL, CodonFreq)
  
  
  AAFreq <- generateAminoAcidFrequency(seqs=lHL$ProSeq)
  colnames(AAFreq) <- paste0("AminoAcid_", colnames(AAFreq))
  lHL <- cbind(lHL, AAFreq)
  
  lHL$mRNAHalfLife <- log10(lHL$half.life..min.)
  
  saveRDS(lHL, paste0(tempFileDir, "LeoMRNAHalfLife.rds"))
  
  
  