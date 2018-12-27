  source('./Conf.R')
  source("./Main.R")

  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)]),
                                               responseStr=c("_PTR"))
  
  
  codonBetaValues <- getBetasOfFeature(allFeaturesFit, featureSubStr="Codon_")
  medianBetas <- 10^apply(codonBetaValues$betas,2,median)
  names(medianBetas) <- sapply(colnames(codonBetaValues$betas), 
                               function(x){strsplit(x,"_")[[1]][2]})
  medianBetasDF <- data.frame(medianBetas=medianBetas, Codon=names(medianBetas))
  
  #### T. Tuller's decoding time estimates
  tTuller <- as.data.frame(read.table(paste0(tempFileDir, "/Tamir_Tuller_decoding_time.txt"), 
                                      stringsAsFactors=F,
                                      header = T))
  tTuller <- tTuller[,c("Codon.Organism.stageortissue.", "H.sapiens5HEK293")]
  colnames(tTuller) <- c("Codon", "Tuller_5HEK293")
  tTuller$Codon <- gsub("T","U", tTuller$Codon)
  
  ### RUST decoding time estimates
  fileNames <- list.files(path = paste0(tempFileDir, "/Riboseq_codon/"), 
                          all.files = T,
                          no.. = T,
                          pattern = "*.txt")
  
  tAll <- data.table(Codon=names(GENETIC_CODE))
  
  for(elem in fileNames){
    t <- read.table(paste0(tempFileDir, "/Riboseq_codon/", elem), stringsAsFactors=F)[,c(1,2,43)]
    t[,2] <- round(as.numeric(t[,2])/as.numeric(t[,3]), digits = 3)
    t[,3] <- NULL
    colnames(t) <- c("Codon", paste0("RUST_",str_split(elem, ".txt")[[1]][1]))
    tAll <- merge(t, tAll, by="Codon")
  }
  tAll$Codon <- gsub("T","U", tAll$Codon)
  t <- merge(tTuller, tAll, by="Codon")
  
  tSup <- t
  tSup <- merge(tSup, medianBetasDF, by="Codon")
  ###########  Kremer et al. Fibroblast cell line predictions
  PPR_fibroblast <- as.data.frame(readRDS(paste0(tempFileDir, "/ppr_ratio_fibrobased_length_weight.RDS")))
  
  selectedColNames <-  c("65126","73804", "78661","80248","80254","81273") 
  PPR_fibroblast <- PPR_fibroblast[selectedColNames]
  PPR_fibroblast$GeneName <- rownames(PPR_fibroblast)
  PPR_fibroblast$fibroblastMedianPTR <- rowMedians(as.matrix(PPR_fibroblast[,selectedColNames]), na.rm = TRUE)
  PPR_fibroblast <- PPR_fibroblast[!is.na(PPR_fibroblast$fibroblastMedianPTR),]
  PPR_fibroblast <- merge(PPR_fibroblast, RNA_Protein_Combined, by="GeneName")
  allfeatures <- getAllFeatureColumns(PPR_fibroblast)
  allFeatureNames <- colnames(PPR_fibroblast[,unlist(allfeatures)])
  
  fit_fibro <- lm(paste0(" fibroblastMedianPTR ~ ", paste(allFeatureNames, collapse = "+")), data = PPR_fibroblast)
  fibroFit <- as.data.frame(coef(summary(fit_fibro)))
  fibroCodonCols <- fibroFit[grepl("Codon_",rownames(fibroFit)) ,]
  fibroCodonCols$Codon <- substr( rownames(fibroCodonCols), 7, 10)
  fibroCodonCols <- fibroCodonCols[, c("Estimate", "Codon")]
  colnames(fibroCodonCols) <- c("Kremer_Fibroblast_PTR", "Codon")
  fibroCodonCols$Kremer_Fibroblast_PTR <- 10^fibroCodonCols$Kremer_Fibroblast_PTR
  
  tSup <- merge(tSup, fibroCodonCols, by="Codon")
  
  extarnalDMelted <- data.table(melt(tSup[,colnames(tSup) %ni% c("Kremer_Fibroblast_PTR")], id.vars=c("Codon", "medianBetas")))
  colnames(extarnalDMelted) <- c("Codon", "MedianCodonEffect", "variable", "value")
  extarnalDMelted[,spCor := round(cor.test(MedianCodonEffect, value, method = "spearman")$estimate, digits = 2),
                  by=variable]
  extarnalDMelted[,pVal := round(cor.test(MedianCodonEffect, value, method = "spearman")$p.value, digits = 3),
                  by=variable]
  
  extarnalDMelted[,variable := paste0(variable, ", R= ", spCor, ", P= ", pVal),]
  
  
  Fig_EV4A <- ggplot(extarnalDMelted, aes(x=value, y=MedianCodonEffect) ) + 
    geom_text(aes(label=Codon),hjust=0, vjust=0,  size=1.5, check_overlap = F)+
    facet_wrap(~variable, ncol=4, scales = "free_x")+
    ylab(paste0("Median 2 fold codon frequency increase \n effect on PTR ratio"))+
    xlab("Expected dwelling time for codon")+
    theme_bw()+
    theme( axis.title.x=element_text(size=14), 
           axis.title.y=element_text(size=14),
           axis.text.x=element_text(size=14),
           axis.text.y=element_text(size=14),
           plot.title = element_text(size=15, hjust = 0.5),
           strip.text.x = element_text(size = 8, colour = "black"))
  