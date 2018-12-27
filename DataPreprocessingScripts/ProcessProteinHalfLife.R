  source('./PaperFigures/Conf.R')
  source('./PaperFigures/Main.R')
  
  
  pHL <- read.csv(paste0(tempFileDir,"ProteinHalfLifes.csv"), stringsAsFactors = F)
  
  hlCols <- grep("half_life", colnames(pHL))
  pHL <- pHL[,c(1,hlCols)]
  pHL[pHL=="#N/A"] <- NA
  
  cCType=grep("Bcells|NK\\.cells|Hepatocytes|Monocytes", colnames(pHL))
  pHL <- pHL[,c(1,cCType)]
  
  
  pHL[,-1] <- lapply(pHL[,-1], as.numeric)
  pHL[pHL==Inf] <- 20000
  
  for(elem in c("Bcells", "NK", "Hepatocytes", "Monocytes")){
    eCols <- grep(elem, colnames(pHL))
    print(eCols)
    pHL[, paste0(elem, "_HL")] <- log10(rowMedians(as.matrix(pHL[,eCols]), na.rm = T))
  }
  
  repCols <- grep("replicate", colnames(pHL))
  pHL[,repCols] <- NULL
  colnames(pHL) <- c("GeneName", "BCells_ProHL", "NKCells_ProHL", "Hepatocytes_ProHL", "Monocytes_ProHL" )
  
  geneIDName <- readRDS(paste0(tempFileDir,"/UnsharedData/geneID_vs_Name.rds"))
  colnames(geneIDName) <- c("EnsemblGeneID", "GeneName")
  
  pHL <- merge(pHL, geneIDName, by="GeneName")
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, 
                                         "PTR_DF_UTR5_UTR3UniqueTranscript.Rds"))
  
  pHL <- merge(pHL, RNA_Protein_Combined, by="EnsemblGeneID")
  
  CodonFreq <- generateCodonFrequency(seqs=pHL$CDSSeq)
  colnames(CodonFreq) <- paste0("Codon_", colnames(CodonFreq))
  pHL <- cbind(pHL, CodonFreq)
  
  AAFreq <- generateAminoAcidFrequency(seqs=pHL$ProSeq)
  colnames(AAFreq) <- paste0("AminoAcid_", colnames(AAFreq))
  pHL <- cbind(pHL, AAFreq)
  
  
  saveRDS(pHL, paste0(tempFileDir, "MathiesonProteinHLData.rds"))

  
  #####################################
  #### HeLa cells protein turnover.
  
  kusterProteinHalfLife <- read.csv(paste0(tempFileDir,
                                           "KusterProteinTurnover.csv"),
                                    header = T,stringsAsFactors = F)
  
  kusterProteinHalfLife <- kusterProteinHalfLife[, c("UniProt.identifier.s.",
                                                     "Gene.name.s.",
                                                     "K..h.1.",
                                                     "k..h.1.",
                                                     "T50...h.",
                                                     "T1.2..h.")]
  colnames(kusterProteinHalfLife) <- c("UniProtID_T",
                                       "GeneName",
                                       "K",
                                       "k",
                                       "T50",
                                       "ProteinHalfLife")
  
  kusterProteinHalfLife$ProteinHalfLife <- log10(as.numeric(kusterProteinHalfLife$ProteinHalfLife))
  
  kusterProteinHalfLife$UniProtID <- sapply(kusterProteinHalfLife$UniProtID_T, function(x){
    strsplit(strsplit(x,"-")[[1]][1], ";")[[1]][1]
  })
  
  
  uniprotProtSeq <- read.csv(paste0(tempFileDir,"/UnsharedData/UniprotProteinSeq.txt"),sep = "\t",
                             stringsAsFactors = F)
  uniprotProtSeq <-uniprotProtSeq[,c("Entry", "Sequence")]
  colnames(uniprotProtSeq) <- c("UniProtID", "ProteinSeq")
  
  
  kusterProteinHalfLife <- merge(kusterProteinHalfLife,uniprotProtSeq,by="UniProtID")
  kusterProteinHalfLife$proteinLength <- sapply(kusterProteinHalfLife$ProteinSeq, 
                                                function(x){nchar(x)})
  
  aaFreqDFFraction <- generateAminoAcidFrequency(seqs=kusterProteinHalfLife$ProteinSeq)
  colnames(aaFreqDFFraction) <- paste0("AminoAcid_", colnames(aaFreqDFFraction))
  kusterProteinHalfLife <- cbind(kusterProteinHalfLife, aaFreqDFFraction)
  
  geneIDName <- readRDS(paste0(tempFileDir,"/UnsharedData/geneID_vs_Name.rds"))
  colnames(geneIDName) <- c("EnsemblGeneID", "GeneName")
  kusterProteinHalfLife <- merge(kusterProteinHalfLife, geneIDName, by="GeneName")
  
  ### Downloaded from https://egg2.wustl.edu/roadmap/data/byDataType/rna/expression/
  ### E117 is HeLa
  #rpkmFile <- read.table(paste0(tempFileDir,"/57epigenomes.RPKM.tsv"), header = T)
  #rpkmFile <- rpkmFile[,c("gene_id","E117")]
  #colnames(rpkmFile) <- c("EnsemblGeneID", "HeLa_RPKM")
  #kusterProteinHalfLife <- merge(kusterProteinHalfLife, rpkmFile, by="EnsemblGeneID")
  #kusterProteinHalfLife$HeLa_RPKM <- log10(kusterProteinHalfLife$HeLa_RPKM+1)
  
  kusterProteinHalfLife <- kusterProteinHalfLife[unique(kusterProteinHalfLife$EnsemblGeneID),]
  saveRDS(kusterProteinHalfLife,
          paste0(tempFileDir,"KusterProteinHalfLife.rds"))
  