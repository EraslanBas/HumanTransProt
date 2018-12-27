  source('./PaperFigures/Conf.R')
  source("./PaperFigures/Main.R")
  
  ############################################
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "PTR_DF_UTR5_UTR3UniqueTranscript.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[RNA_Protein_Combined$UTR3_length>6 &
                                               RNA_Protein_Combined$UTR5_length>6,]
  #### Region GC contents
  for(elem in mRNARegions){
    RNA_Protein_Combined[,paste0(elem,"GCPercent")] <- unlist(sapply(RNA_Protein_Combined[,paste0(elem,"Seq")],
                                                                     function(x){return(letterFrequency(RNAString(x),
                                                                                                        letters = "GC", as.prob = T)) }))
  }
  #################################################################################
  ##Get Kozak seq. content
  lroffset=6
  RNA_Protein_Combined$kozakSeq  <- apply(RNA_Protein_Combined, 1, 
                                          function(x){ paste0(substr(x["UTR5Seq"], nchar(x["UTR5Seq"])-lroffset+1,nchar(x["UTR5Seq"])), 
                                                          substr(x["CDSSeq"], 1, lroffset+3 )  )})
  
  isatg <- unlist(lapply(RNA_Protein_Combined$kozakSeq, 
                         function(x){substr(x,lroffset+1,lroffset+3)}))
  RNA_Protein_Combined <- RNA_Protein_Combined[which(isatg == "AUG"),]
  
  dmI = createDesignMatrixWithGivenSeqs(RNA_Protein_Combined$kozakSeq,
                                        lroffset,
                                        msize=3)
  kozakSeqMatrix <- as.data.frame(dmI$dmatrix.c)
  kozakSeqMatrix <- kozakSeqMatrix[,-1]
  
  motif = RNAStringSet(unlist(dmI$seq.c))
  background = apply(consensusMatrix(motif,as.prob=TRUE)[1:4,],2,which.max)
  ind_rm = sapply(1:length(background),function(x){(x-1)*4+background[x]})
  lm_pos = dmI$all_pos[-ind_rm]
  
  ###### Remove the consensus columns
  kozakSeqMatrix <- kozakSeqMatrix[,lm_pos]
  
  ###### Remove ATG cols
  #atgCols <- getColumnsOfSubstr(kozakSeqMatrix, inputStr=c("P1\\.","P2\\.","P3\\.") , exactMatch=T)
  #kozakSeqMatrix<-kozakSeqMatrix[,-atgCols]
  
  colnames(kozakSeqMatrix) <- paste0("Kozak_", colnames(kozakSeqMatrix))
  kozakSeqMatrix$EnsemblTranscriptID <- RNA_Protein_Combined$EnsemblTranscriptID
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined,
                                kozakSeqMatrix,
                                by="EnsemblTranscriptID",
                                all.x =T)
  
  
  #########################################################################
  ########
  
  restCodons <- sapply(RNA_Protein_Combined$CDSSeq,function(x){substr(x,10,nchar(x))})
  restCodonFreq <- generateCodonFrequency(seqs=restCodons)
  
  colnames(restCodonFreq) <- paste0("Codon_", colnames(restCodonFreq))
  RNA_Protein_Combined <- cbind(RNA_Protein_Combined, restCodonFreq)
  
  
  ###########################################################################
  #########
  
  restAminoAcids <- sapply(RNA_Protein_Combined$ProSeq,function(x){substr(x,4,nchar(x))})
  restAAFreq <- generateAminoAcidFrequency(seqs=restAminoAcids)
  
  colnames(restAAFreq) <- paste0("AminoAcid_", colnames(restAAFreq))
  RNA_Protein_Combined <- cbind(RNA_Protein_Combined, restAAFreq)
  
  
  ################################################################
  #### get codon pair bias
  codonPairList <- lapply(RNA_Protein_Combined$CDSSeq, function(x)
  {oligonucleotideFrequency(RNAString(x), 6, step=3)})
  
  codonPairBias <- as.data.frame(do.call(rbind, codonPairList), stringsAsFactors=F)
  codonPairBias[,which(colSums(codonPairBias)==0)] <- NULL
  codonPairBias <- as.data.frame(as.matrix(codonPairBias)+1)
  codonPairBias <- sweep(codonPairBias, 1, rowSums(codonPairBias), `/`)
  
  codonPairBias <- as.data.frame(myPrcomp(as.matrix(codonPairBias),n=2))
  colnames(codonPairBias) <- paste0("CPair_", 1:2)
  
  RNA_Protein_Combined <- cbind(RNA_Protein_Combined, codonPairBias)
  
  #########################################################################
  ################## Get 5' UTR folding energy values
  UTR5_tmp <- readRDS(paste0(tempFileDir, "UTR5_FE.rds"))
  UTR5_tmp[,-1] <- log2(abs(UTR5_tmp[-1])+1)
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined,UTR5_tmp, by="EnsemblTranscriptID", all.x =T)
  
  #####################################################################
  lroffset=6
  ##Get stop codon content
  RNA_Protein_Combined$sCodContext <- apply(RNA_Protein_Combined, 
                                                1,
                                                function(x){ paste0(substr(x["CDSSeq"], nchar(x["CDSSeq"])-lroffset-2,nchar(x["CDSSeq"])), 
                                                                    substr(x["UTR3Seq"], 1, lroffset )  )})
  
  isT<- unlist(lapply(RNA_Protein_Combined$sCodContext, 
                      function(x){substr(x,lroffset+1,lroffset+1)}))
  RNA_Protein_Combined <- RNA_Protein_Combined[which(isT == "U"),]
  
  dmI = createDesignMatrixWithGivenSeqs(RNA_Protein_Combined$sCodContext, 
                                        lroffset,
                                        msize=3)
  stopCodMatrix <- as.data.frame(dmI$dmatrix.c)
  stopCodMatrix <- stopCodMatrix[,-1]
  
  motif = RNAStringSet(unlist(dmI$seq.c))
  background = apply(consensusMatrix(motif,as.prob=TRUE)[1:4,],2,which.max)
  ind_rm = sapply(1:length(background),function(x){(x-1)*4+background[x]})
  lm_pos = dmI$all_pos[-ind_rm]
  
  ###### Remove the consensus columns
  stopCodMatrix <- stopCodMatrix[,lm_pos]
  colnames(stopCodMatrix) <- paste0("stop_", colnames(stopCodMatrix))
  stopCodMatrix$EnsemblTranscriptID <- RNA_Protein_Combined$EnsemblTranscriptID
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, stopCodMatrix,
                                by="EnsemblTranscriptID",
                                all.x =T)
  
  ############################################
  ##### Get N-end AA
  #RNA_Protein_Combined$nEndA_ <- substr(RNA_Protein_Combined$ProSeq,2,2)
  #RNA_Protein_Combined <- cbind(RNA_Protein_Combined, as.data.frame(model.matrix(~nEndA_, RNA_Protein_Combined)[,-1]))
  #RNA_Protein_Combined$nEndA_ <- NULL
  
  ############################################
  ##############################################
  ####### Add protein complex formation info
  allProteinComplexes <- readRDS(paste0(tempFileDir,"allProteinComplexes.Rds"))
  allProteinComplexes$complexID <- 1:dim(allProteinComplexes)[[1]]
  allTranscriptsInComplexes <- unlist(allProteinComplexes$transID)
  
  RNA_Protein_Combined$ProteinComplex <- 0
  RNA_Protein_Combined[which(RNA_Protein_Combined$EnsemblTranscriptID %in% allTranscriptsInComplexes),"ProteinComplex"] <- 1
  
  ###################################################
  ########## Add pest sequence info
  pestOuts <- as.data.frame(readRDS(paste0(tempFileDir,"pestOuts.rds")), stringsAsFactors=F)
  pestOuts$poorPESTNo <- NULL
  colnames(pestOuts) <- c("EnsemblTranscriptID","PESTNoPotential")
  
  pestOuts$PESTNoPotential <- unlist(pestOuts$PESTNoPotential)
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, pestOuts, by="EnsemblTranscriptID")
  
  ### Explained variance is higher when PEST regions are encoded as binary
  RNA_Protein_Combined[RNA_Protein_Combined$PESTNoPotential>0,"PESTNoPotential"] <- 1
  
  #############################################
  ######## Add mRNA m6A methylation peaks for HepG2 cell line
  ### https://www.nature.com/articles/nature11112
  m6AData <- data.frame(read.xls(paste0(tempFileDir, "/nature11112-s3.xls"),  header = TRUE), stringsAsFactors = F)
  m6AData$geneSymbol <- toupper(as.character(m6AData$geneSymbol))
  m6AData <- m6AData[!duplicated(m6AData$geneSymbol),]
  rnaMthData <- m6AData[,c("HumanHEPG2_peaksnum", "geneSymbol")]
  colnames(rnaMthData) <- c("RNAm6A", "GeneName")
  
  #### Keep the genes which do not have RNA modification peaks
  RNA_Protein_Combined <- merge(x=RNA_Protein_Combined, y=rnaMthData, by="GeneName", all.x =T)
  RNA_Protein_Combined[is.na(RNA_Protein_Combined$RNAm6A),"RNAm6A"] <- 0
  
  ############################################
  ######## Add miRNA data
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
    RNA_Protein_Combined[,paste0("sgn_MiRNA_",miRNA)] <- 0
    RNA_Protein_Combined[RNA_Protein_Combined$GeneName %in% elem$Target.Gene, paste0("sgn_MiRNA_",miRNA)] <- 1
  }
  
  miRNATargetSites <- NULL
  expVarmiRNAs <- NULL
  gc()
  
  k <- grep("sgn_MiRNA_", colnames(RNA_Protein_Combined))
  X <- RNA_Protein_Combined[,k]
  RNA_Protein_Combined[,k] <- NULL
  
  X[,which(colSums(X) < miRNAtargetGeneCutOff)] <- NULL
  rownames(X) <- RNA_Protein_Combined$EnsemblGeneID
  saveRDS(X, paste0(tempFileDir,"miRNATargetMatrix.rds"))
  
  miRNAs <- as.data.frame(myPrcomp(as.matrix(X),n=150))
  colnames(miRNAs) <- paste0("sgn_MiRNA_PC",1:ncol(miRNAs))
  RNA_Protein_Combined <- cbind(RNA_Protein_Combined, miRNAs)
  
  
  ########## RBPs
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
    RNA_Protein_Combined[,paste0("RBP_",rbpNames[i])] <- 0
    RNA_Protein_Combined[RNA_Protein_Combined$EnsemblGeneID %in% rbpPeaksFileSplit[[i]]$geneID,paste0("RBP_",rbpNames[i])] <- 1
  }
  rbpPeaksFile<- NULL
  rbpPeaksFileSplit <- NULL
  gc()
  
  
  ################################################################
  #### Add protein isoelectric point data
  
  tmp <- readRDS(paste0(tempFileDir,"proteinIsoelectricPoint.rds"))
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, tmp, by="EnsemblTranscriptID")
  
  ################################################################
  #### Add average hydrophobicity of first 15 amino acids
  NtermProteinHalfLife <- data.frame(read.table(paste0(tempFileDir,"N_term_halfLife.tsv"),
                                                header = T), stringsAsFactors = F)
  NtermProteinHalfLife$Hydrophobicity_values <- as.numeric(NtermProteinHalfLife$Hydrophobicity_values)
  NtermProteinHalfLife$aa <- as.character(NtermProteinHalfLife$aa)
  
  aaSeq  <- apply(RNA_Protein_Combined, 1, 
                                       function(x){ substr(x["ProSeq"], 2,16 )})
  
  aaIndv <- lapply(aaSeq, s2c)
  aaIndv <- data.frame(do.call(rbind, aaIndv), stringsAsFactors = F)
  
  for ( i in 1:20){
    aaIndv[aaIndv==NtermProteinHalfLife[i,"aa"]]<- as.numeric(NtermProteinHalfLife[i,"Hydrophobicity_values"])
  }
  
  aaIndv[] <- lapply(aaIndv, function(x) as.numeric(as.character(x)))
  
  RNA_Protein_Combined$ntermHydro <- apply(aaIndv,1, function(x){mean(x, na.rm=T)})
  
  #################################################################
  RNA_Protein_Combined <- addMotifsToDF(DF=RNA_Protein_Combined, 
                                        motifDB=motifDB[motifDB$consensusMotif != "AUG",],
                                        consOrUpdated="consensusMotif",
                                        allCount=F)
  
  ######## Only consider out-frame AUGs as upstream AUGs
  atgCount <- unlist(vcountPattern("AUG", RNA_Protein_Combined$UTR5Seq))
  atgPos <- start(vmatchPattern("AUG", RNA_Protein_Combined$UTR5Seq))
  
  frameInfo <- c()
  for(i in 1:length(atgPos)){
    if(length(atgPos[[i]])>0){
      atgDist <- RNA_Protein_Combined[i, "UTR5_length"] - atgPos[[i]]+1
      
      if(any((atgDist %% 3) != 0)){
        frameInfo <- c(frameInfo, 1)
      }else{
        frameInfo <- c(frameInfo, 0)
      }
      
    }else{
      frameInfo <- c(frameInfo, 0)
    }
  }
  
  RNA_Protein_Combined$Motif_UTR5_AUG <- frameInfo
  
  
  ##################################################################################
  
  RNA_Protein_Combined <- addMotifCounts(RNA_Protein_Combined, 
                                         AA_motifs, 
                                         RNA_Protein_Combined$ProSeq,
                                         "Motif_AA_", 
                                         binary=T, 
                                         fixedP=T,
                                         isAA=T,
                                         maxMismatch=0)
  
  CDS_Cterm_AASeq = apply(RNA_Protein_Combined,
                          1, 
                          function(x){ substr(x["ProSeq"], nchar(x["ProSeq"])-50,nchar(x["ProSeq"]))  })
  RNA_Protein_Combined <- addMotifCounts(RNA_Protein_Combined, 
                                         CDS_CtermAAMotifs,
                                         CDS_Cterm_AASeq,
                                         "Motif_AA_",
                                         binary=T,
                                         fixedP=T,
                                         isAA=T,
                                         maxMismatch=0)
  
  #####################################################################################
  ########### Include polyP-stretches
  #RNA_Protein_Combined$polyPstretches <- unlist(lapply(RNA_Protein_Combined$ProSeq, function(x){k = lengths(regmatches(x, gregexpr("P{3,}", x)) )}))
  #RNA_Protein_Combined[RNA_Protein_Combined$polyPstretches>1,"polyPstretches"]<-1
  
  ######################################################################################
  elmClas <- fread(paste0(tempFileDir,"elm_classes_all.tsv"), header = T)
  elmIdent <- gsub(" ", "_", elmClas$ELMIdentifier)
  elmIdent <- gsub("-", "_", elmIdent)
  elmClas$ELMIdentifier <- elmIdent
  signMotifs <- readRDS(paste0(tempFileDir,"significantLinearProteinMotifs.Rds"))
  
  elmClas <- elmClas[elmClas$ELMIdentifier %in% signMotifs,]
  
  
  for(i in 1:dim(elmClas)[[1]]){
    elem <- elmClas[i, "Regex"]
    Mname <- elmClas[i, "ELMIdentifier"]
    searchRes <- regexpr(elem, RNA_Protein_Combined$ProSeq , perl=TRUE)
    searchRes[searchRes==-1] <- 0
    searchRes[searchRes > 1] <- 1
    RNA_Protein_Combined[[paste0("Motif_LinAA_",Mname)]] = searchRes
  }
  
  ################################################################
  ####### Add protein post-translational modifications
  td <- readRDS(paste0(tempFileDir, "ProteinPostTranslationalMod.rds"))
  
  
  RNA_Protein_Combined <- merge(RNA_Protein_Combined, td , by="EnsemblGeneID",all.x =T)
  
  for(i in colnames(td)[-1])
  {
    RNA_Protein_Combined[is.na(RNA_Protein_Combined[,i]),i] <- 0
    RNA_Protein_Combined[,i] <- log2(RNA_Protein_Combined[,i]+1)
  }
  
  RNA_Protein_Combined[, c("CDS_length","UTR3_length","UTR5_length")] <- log2(RNA_Protein_Combined[, c("CDS_length","UTR3_length","UTR5_length")])
  
  
  ####################################
  ## Not sequence features but metadata
  gb <- read.csv(paste0(tempFileDir, "Gerstberger_RBP_data.csv"), header = T, stringsAsFactors = F)
  rownames(gb) <- gb$gene.name
  availableRBPS <- RNA_Protein_Combined[RNA_Protein_Combined$GeneName %in% gb$gene.name, "GeneName"]
  
  RNA_Protein_Combined$isRBP <- F
  RNA_Protein_Combined[RNA_Protein_Combined$GeneName %in% gb$gene.name, "isRBP"] <- T
  RNA_Protein_Combined[RNA_Protein_Combined$GeneName %in% gb$gene.name, "RBPBindingCategory"] <- gb[availableRBPS, "consensus.RNA.target"]
  
  
  
  RNA_Protein_Combined$mRNATisSScore <- calculateGerstbergerScore(as.matrix(RNA_Protein_Combined[,paste0(TissueNames,"_medianExon")]))
  RNA_Protein_Combined$proteinTisSScore <- calculateGerstbergerScore(RNA_Protein_Combined[,paste0(TissueNames,"_iBAQ")])
  RNA_Protein_Combined$PTRTisSScore <- calculateGerstbergerScore(RNA_Protein_Combined[,paste0(TissueNames,"_PTR")])
  
  
  essentialGeneList <- read.table(paste0(tempFileDir,"/HumanEssentialGenes.txt"), header = F)
  essentialGeneList <- essentialGeneList[essentialGeneList$V5=="E",]
  RNA_Protein_Combined$essentialGene <- F
  RNA_Protein_Combined[RNA_Protein_Combined$EnsemblGeneID %in% essentialGeneList$V4,"essentialGene"] <- T
  
  for(elem in TissueNames){
    fit <- lm(RNA_Protein_Combined[,paste0(elem,"_PTR")]~RNA_Protein_Combined[,paste0(elem,"_medianExon")], 
              na.action = na.exclude)
    RNA_Protein_Combined[,paste0(elem,"_mRNAINDP")] <- residuals(fit)
    RNA_Protein_Combined[,paste0(elem,"_mRNADP")] <- RNA_Protein_Combined[,paste0(elem,"_PTR")]-RNA_Protein_Combined[,paste0(elem,"_mRNAINDP")]
  }

  houseKeepingGenes <- read.csv(paste0(tempFileDir,"tissue_specificity_rna_any_expressed.tsv"), sep = '\t')
  RNA_Protein_Combined$housekeepingGene <- F
  RNA_Protein_Combined[RNA_Protein_Combined$EnsemblGeneID %in% houseKeepingGenes$Ensembl,"housekeepingGene"] <- T
  
  saveRDS(RNA_Protein_Combined, paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
