  source('./Conf.R')
  source("./Main.R")
  library("ggrepel")
  
  ########### EV3 A
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, 
                                         "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  
  codonUsageAll <- c(allfeatures$codonUsage, allfeatures$codonPairBias)
  
  codonExpVars <- c()
  
  X <- RNA_Protein_Combined[,codonUsageAll]
  for(tissue in TissueNames){
    res_lm <- compute10FCV_mixedEffectModel(X = X, 
                                            y = RNA_Protein_Combined[,paste0(tissue, "_PTR")], 
                                            rnaCol=NULL, 
                                            fixedEffectCols=colnames(X), 
                                            randomEffectCol=NULL, 
                                            isMixedEffect=F, 
                                            flds = NULL)
    
    codonExpVars <- c(codonExpVars, res_lm$cvR2)
  }
  
  k <- grep("AminoAcid_", colnames(RNA_Protein_Combined))
  X <- RNA_Protein_Combined[,k]
  
  
  aaExpVars <- c()
  
  for(tissue in TissueNames){
    res_lm <- compute10FCV_mixedEffectModel(X = X, 
                                            y = RNA_Protein_Combined[,paste0(tissue, "_PTR")], 
                                            rnaCol=NULL, 
                                            fixedEffectCols=colnames(X), 
                                            randomEffectCol=NULL, 
                                            isMixedEffect=F, 
                                            flds = NULL)
    
    aaExpVars <- c(aaExpVars, res_lm$cvR2)
  }
  
  tmp <- data.frame(codonExpVars=codonExpVars, 
                    aaExpVars=aaExpVars,
                    TissueNames=factor(TissueNamesPrint))
  
  tmpMelted <- melt(tmp)
  Fig_EV3B <- ggplot(tmpMelted, aes(TissueNames, y=100*value)) +   
    geom_bar(aes(fill = variable), position = "dodge", stat="identity")+
    scale_x_discrete(limits = rev(levels(tmpMelted$TissueNames)))+
    coord_flip()+
    ylab("Explained variance in PTR ratio (%)")+
    theme(legend.position = "top",
          axis.title.y =element_blank())+
    scale_fill_discrete(name="",
                        breaks=c("codonExpVars", "aaExpVars"),
                        labels=c("Codon frequency", "Amino acid frequency"))
  
  
  ########### EV3 C-D
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  RNA_Protein_Combined <- RNA_Protein_Combined[2^RNA_Protein_Combined$CDS_length>460,]
  getCodonFrequency <-  function(cdsSeq){
    seqLength <- sapply(cdsSeq , nchar)
    codonOccList <- lapply(cdsSeq, function(x)
    {trinucleotideFrequency(RNAString(x), step=3)
    })
    
    codonFreqDFNoBias <- as.data.frame(do.call(rbind, codonOccList), stringsAsFactors=F)
    codonFreqDFNoBias$UAA <- NULL
    codonFreqDFNoBias$UGA <- NULL
    codonFreqDFNoBias$UAG <- NULL
    
    codonFreqDFNoBias <- as.data.frame(as.matrix(codonFreqDFNoBias)+1)
    codonFreqDFNoBias <- sweep(codonFreqDFNoBias, 1, RNA_Protein_Combined$CDS_length, `/`)
    codonFreqDFNoBias <- log2(codonFreqDFNoBias)
    
    colnames(codonFreqDFNoBias) <- paste0("Codon_", colnames(codonFreqDFNoBias))
    
    return(codonFreqDFNoBias)
  }
  
  getWindowCodonEffectsOnPTR <- function(windowsStart, windowsEnd){
    
    RNA_Protein_CombinedTemp <- RNA_Protein_Combined
    allfeatures <- getAllFeatureColumns(RNA_Protein_CombinedTemp)
    
    RNA_Protein_CombinedTemp[,allfeatures$codonUsage] <- NULL
    
    
    if(!is.null(windowsEnd)){
      RNA_Protein_CombinedTemp$windowCodons  <- apply(RNA_Protein_CombinedTemp, 1, 
                                                      function(x){ substr(x["CDSSeq"], windowsStart, windowsEnd )})
    }else{
      RNA_Protein_CombinedTemp$windowCodons  <- apply(RNA_Protein_CombinedTemp, 1, 
                                                      function(x){ substr(x["CDSSeq"], windowsStart, nchar(x["CDSSeq"]) )})
    }
    
    
    codonFreqOfWindow <- getCodonFrequency(cdsSeq=RNA_Protein_CombinedTemp$windowCodons)
    codonMedianFreqOfWindow <- apply(codonFreqOfWindow,2,median)
    
    RNA_Protein_CombinedTemp <- cbind(RNA_Protein_CombinedTemp, codonFreqOfWindow)
    allfeatures <- getAllFeatureColumns(RNA_Protein_CombinedTemp)
    
    allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_CombinedTemp,
                                                allFeatureNames = colnames(RNA_Protein_CombinedTemp[,unlist(allfeatures)]),
                                                responseStr=c("_PTR"))
    codonFit <- getBetasOfFeature(allFeaturesFit, featureSubStr="Codon_")
    codonFit_median <- apply(as.data.frame(10^codonFit$betas),2, median)
    
    return(list(codonMedianFreqOfWindow=codonMedianFreqOfWindow, codonFit_median=codonFit_median))
  }
  
  
  tp <- data.frame(windowStarts=seq(4,454,45),
                   windowEnds=seq(48, 499, 45))
  
  allWindowCodonMeds <- apply(tp, 1, function(x){return(getWindowCodonEffectsOnPTR(x["windowStarts"], 
                                                                                   x["windowEnds"]))})
  
  medianFreqs <- list()
  for(elem in allWindowCodonMeds){
    medianFreqs <- lappend(medianFreqs, elem$codonMedianFreqOfWindow)
  }
  
  medianFreqs <- do.call(rbind, medianFreqs)
  
  
  medianFreqs <- t(as.matrix(medianFreqs))
  colnames(medianFreqs) <- paste0("codons ", paste0(round(((tp$windowStarts / 3)+ 1)), " to ", (tp$windowEnds / 3)))
  rownames(medianFreqs) <- unique(unlist(strsplit(rownames(medianFreqs), "Codon_")))[-1]
  
  
  medianFreqs <- t(scale(t(medianFreqs), center = T))
  
  ### Fig_EV3 C
  pheatmap(medianFreqs,
           cluster_cols = F,
           fontsize_row = 5,
           fontsize_col = 7,
           treeheight_row = 0)

  
  ### EV3 D
  firstWindow <- getWindowCodonEffectsOnPTR(4, 48)
  restWindow <- getWindowCodonEffectsOnPTR(48, NULL )
  
  k <- unlist(strsplit(names(firstWindow$codonFit_median), "Codon_"))
  k <- k[k!=""]
  tp <- data.frame(firstRegionEffect = firstWindow$codonFit_median, secondRegionEffect=restWindow$codonFit_median, codons=k)
  
  t <- cor.test(tp$firstRegionEffect, tp$secondRegionEffect, method="spearman")
  
  Fig_EV3D <- ggplot(tp, aes(x=firstRegionEffect, y=secondRegionEffect) ) + 
    geom_point(size=0.5, color="red")+
    geom_text(aes(label=as.character(codons)),
              hjust=0,
              vjust=0,
              size=2,
              check_overlap = F)+
    xlim(c(0.8, 1.3))+
    ylim(c(0.8, 1.3))+
    scale_x_continuous(trans="log2",breaks = seq(0.8,1.3,0.1), labels=seq(0.8,1.3,0.1))+
    scale_y_continuous(trans="log2",breaks = seq(0.8,1.3,0.1), labels=seq(0.8,1.3,0.1))+
    xlab(paste0("Effect of 2 fold codon frequency increase in \n proximal 15 codons on PTR ratio"))+
    ylab("Effect of 2 fold codon frequency increase in \n rest of the coding region on PTR ratio")+
    geom_abline(slope = 1, intercept = 0)+
    ggtitle(paste0("rho = ", round(t$estimate, digits = 2), " , P < 0.0024 "))
    