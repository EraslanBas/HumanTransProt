  source('./Conf.R')
  source('./Main.R')
  
  
  ############ Get median effect estimates (across tissues) of our data
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  allFeatureNames = colnames(RNA_Protein_Combined[,unlist(allfeatures)])
  allFeaturesFit <- selectBetasAccrossTissues( inputDF=RNA_Protein_Combined, 
                                               allFeatureNames,
                                               responseStr=c("_PTR"))
  
  Motif_TV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Motif_UTR")
  pvals <- Motif_TV$pvalues
  
  Motifsign <- apply(pvals,2, function(x){length(which(x < 0.1))})
  signMotifs <- names(Motifsign[Motifsign>5])
  
  medianMotifEffect <- apply(Motif_TV$betas,2,median)
  medianMotifStdEr <- apply(Motif_TV$stdError,2,min)
  
  codon_TV <- getBetasOfFeature(allFeaturesFit, featureSubStr="Codon")
  medianCodonEffect <- apply(codon_TV$betas,2,median)
  
  ourData <- c(medianCodonEffect,medianMotifEffect)
  ourEstimates <- data.frame(feature= names(ourData),
                             EstimatesOur = 10^ourData,
                             data="tissueProteome")
  ############# Get fibroblast data estimates
  PPR_fibroblast <- as.data.frame(readRDS(paste0(tempFileDir, "/ppr_ratio_fibrobased_length_weight.RDS")))
  
  selectedColNames <-  c("65126","73804", "78661","80248","80254","81273") 
  PPR_fibroblast <- PPR_fibroblast[selectedColNames]
  PPR_fibroblast$GeneName <- rownames(PPR_fibroblast)
  PPR_fibroblast$fibroblastMedianPTR <- rowMedians(as.matrix(PPR_fibroblast[,selectedColNames]), na.rm = TRUE)
  PPR_fibroblast <- PPR_fibroblast[!is.na(PPR_fibroblast$fibroblastMedianPTR),]
  PPR_fibroblast <- merge(PPR_fibroblast, RNA_Protein_Combined, by="GeneName")
  
  fit_fibro <- lm(paste0(" fibroblastMedianPTR ~ ", paste(allFeatureNames, collapse = "+")), data = PPR_fibroblast)
  fibroFit <- as.data.frame(coef(summary(fit_fibro)))
  codonOrMotifCols <- grepl("Codon_",rownames(fibroFit)) | grepl("Motif_UTR",rownames(fibroFit))
  
  fibroEstimates <- data.frame(feature= rownames(fibroFit[codonOrMotifCols,]),
                               EstimatesFibro = 10^fibroFit[codonOrMotifCols,  "Estimate"],
                               data="Fibroblast")
  
  
  ################################
  ###### Plot estimates of our data vs. estimates of fibroblast data
  
  df <- merge(ourEstimates,fibroEstimates, by="feature")
  df <- df[df$feature %in% c(colnames(codon_TV$betas), signMotifs),]
  
  df$featureType <- "Codon"
  df[grepl("Motif_UTR3", df$feature), "featureType"]<-"Motif_UTR3"
  df[grepl("Motif_UTR5", df$feature), "featureType"]<-"Motif_UTR5"
  
  df$printLabel <- apply(df,1, 
                         function(x){strsplit(x["feature"], paste0(x["featureType"],"_"))[[1]][2]})
  df[df$featureType=="Codon","printLabel"] <- ""
  corRes <- cor.test(df$EstimatesOur,df$EstimatesFibro, method = "spearman")
  
  Figure6A <- ggplot(df, aes(x=EstimatesOur, y=EstimatesFibro, color=featureType,label=featureType))+
    geom_point(size=0.5) +
    scale_y_continuous(trans='log2',
                       breaks = seq(0.1,2.1,0.2) ,
                       labels = seq(0.1,2.1,0.2), 
                       limits = c(0.3,2.1))+
    scale_x_continuous(trans='log2', 
                       breaks = seq(0.1,2.1,0.1) ,
                       labels = seq(0.1,2.1,0.1),
                       limits = c(0.7,1.7))+
    ylab("Effect in fibroblast cells")+
    xlab("Median effect across 29 tissues")+
    geom_abline(color="darkgrey", slope = 1,intercept = 0)+
    geom_hline(yintercept = 1, linetype="dashed",color="darkgrey")+
    geom_vline(xintercept = 1, linetype="dashed",color="darkgrey")+
    geom_text_repel(
      aes(label = printLabel),
      size = 2,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.3, "lines")
    )+ourTheme+
    theme( legend.position = "bottom",
           legend.title = element_blank())+
    ggtitle(paste0(expression(rho)," = ",round(corRes$estimate, digits = 2),", P < 2.2e-16 "))
  

  ###################################
  ##### correct for experiment level and then do the nested ANOVA test
  
  utr3Data <- readRDS(paste0(tempFileDir, "utr3MotifValData.rds"))
  utr5Data <- readRDS(paste0(tempFileDir, "utr5MotifValData.rds"))
  
  fit <- lm(utr3Data$intensity ~ utr3Data$Experiment)
  utr3Data$intensity <- residuals(fit) + coef(summary(fit))[1,1]
  
  fit <- lm(utr5Data$intensity ~ utr5Data$Experiment)
  utr5Data$intensity <- residuals(fit) + coef(summary(fit))[1,1]
  
  
  
  valMotifs <- motifDB[motifDB$experimentsDone,c("consensusMotif",
                                                 "motifType",
                                                 "sequenceUsedInExperiment")]
  
  valMotifs <- valMotifs[valMotifs$consensusMotif != "AUG",]
  
  tp <- data.frame(consensusMotif= sapply(names(medianMotifEffect), function(x){ strsplit(x,"_")[[1]][3]}),
                   linearModelEstimate=10^medianMotifEffect, side="greater", stringsAsFactors = F)
  tp[tp$linearModelEstimate < 1, "side"] <- "less"
  
  valMotifs <- merge(valMotifs, tp, by="consensusMotif")
  
  getMotifSignificanceNestedAnova<- function(dfTemp, motif){
    dfTempK <- dfTemp[dfTemp$motif %in% c(paste0("SCRAMBLED_",motif), motif),]
    dfTempK$motif <- drop.levels(dfTempK$motif)
    dfTempK$motif <- factor(dfTempK$motif, levels = c(paste0("SCRAMBLED_",motif), motif))
    dfTempK <- dfTempK[,c("motif", 
                          "Experiment",
                          "TB", 
                          "intensity",
                          "time.point")]
    dfTempK$TB <- paste0(dfTempK$TB, "_", dfTempK$motif, "_", dfTempK$Experiment)
    
    model = lme(intensity ~ motif, random=~1|TB, 
                data=dfTempK, 
                method="REML")
    anovaRes <- anova.lme(model, 
                          type="sequential", 
                          adjustSigma = FALSE)
    
    return(list(pVal = anovaRes$`p-value`[2], fc=10^model$coefficients$fixed[2]))
  }
  
  
  
  testPValues <- list()
  foldChanges <- list()
  for(i in 1:nrow(valMotifs)){
    if(valMotifs[i,"motifType"]=="UTR3"){
      k <- getMotifSignificanceNestedAnova(dfTemp=utr3Data, 
                                           motif=as.character(valMotifs[i,"sequenceUsedInExperiment"]))
    }else{
      k <- getMotifSignificanceNestedAnova(dfTemp=utr5Data, 
                                           motif=as.character(valMotifs[i,"sequenceUsedInExperiment"]))
    }
    
    
    testPValues <- lappend(testPValues, k[[1]])
    foldChanges <- lappend(foldChanges, k[[2]])
  }
  
  
  
  testPValues <- as.data.frame(do.call(rbind,testPValues))
  colnames(testPValues) <- "PVAL"
  foldChanges <- as.data.frame(do.call(rbind, foldChanges))
  colnames(foldChanges) <- "Median_FC"
  
  
  valMotifs <- cbind(valMotifs, foldChanges)
  valMotifs <- cbind(valMotifs, testPValues)
  
  
  valMotifs$FDR <- p.adjust(valMotifs$PVAL, method = "BH")
  valMotifs$fcSame <- "greater"
  
  valMotifs[valMotifs$Median_FC < 1, "fcSame"] <- "less"
  
  valMotifs$significant <- "No"
  valMotifs[valMotifs$FDR < 0.1 & (valMotifs$side == valMotifs$fcSame), "significant"] <- "Yes"
  
  valMotifs$significant <- factor(valMotifs$significant, levels=c("Yes", "No"))
  
  ##### Box plots of validated motifs
  utr5Data <- readRDS(paste0(tempFileDir, "utr5MotifValData.rds"))
  augData <- readRDS(paste0(tempFileDir, "augValData.rds"))
  
  fit <- lm(utr5Data$intensity ~ utr5Data$Experiment)
  utr5Data$intensity <- residuals(fit) + coef(summary(fit))[1,1]
  
  fit <- lm(augData$intensity ~ augData$Experiment)
  augData$intensity <- residuals(fit) + coef(summary(fit))[1,1]
  
  utr5Data$intensity <- 10^utr5Data$intensity
  augData$intensity <- 10^augData$intensity
  
  plotMotifValidationBoxplots <- function(df, selectedLevels, legendLabels, breakWidth=2, plotXLabel=F, setBreaks=F, breaksSeq, heightL){
    dfTempK <- df[df$motif %in% selectedLevels,]
    dfTempK$motif <- drop.levels(dfTempK$motif)
    dfTempK$motif <- factor(dfTempK$motif, levels = selectedLevels)
    
    p <- ggplot(dfTempK, aes(y=intensity, fill=motif, x=time.point))+
      geom_boxplot(alpha = 0.80, lwd=0.1, outlier.shape = NA) +
      geom_point(aes(fill = motif), size = 0.5, shape = 21, position = position_jitterdodge()) +
      ylab("GLuc / SEAP \n intensity")+
      scale_fill_discrete(name = "Inserted motif", labels = legendLabels)+
      theme(text = element_text(family = fontFamily), 
            axis.title.x=element_blank(), 
            axis.title.y=element_text(family = fontFamily, size=7),
            axis.text.x=element_text(family = fontFamily, size=7),
            axis.text.y=element_text(family = fontFamily, size=7),
            plot.title = element_blank(),
            legend.text = element_text(family = fontFamily, size=7),
            legend.title = element_text(family = fontFamily, size=7),
            legend.position = c(0.53,0.98),
            #legend.position = "right",
            legend.direction = "horizontal",
            legend.box.margin=margin(-2,-2,-2,-2),
            legend.box.spacing=unit(0.1, 'cm'))
    
    if(plotXLabel){
      p <- p+ xlab("Time (minutes)")+theme(axis.title.x=element_text(family = fontFamily, size=7))
    }
    
    if(!setBreaks){
      p <- p+ scale_y_continuous(trans="log10",
                         breaks = seq(round((min(dfTempK$intensity)-1)),round((max(dfTempK$intensity)+1)),breakWidth))
    }else{
      p <- p+ scale_y_continuous(trans="log10",
                                 breaks = breaksSeq)
    }
    
    print(p)
    
  }
  
  
  plotMotifValidationBoxplots(df=utr5Data, 
                              selectedLevels=c("CUGUCCU", "SCRAMBLED_CUGUCCU"), 
                              legendLabels=c("CUGUCCU", "Scrambled CUGUCCU"),
                              breakWidth=4,
                              plotXLabel=T,
                              heightL=1.5)
  plotMotifValidationBoxplots(df=utr5Data, 
                              selectedLevels=c("UUCCGUUCCG", "SCRAMBLED_UUCCGUUCCG"), 
                              legendLabels=c("UUCCG", "Scrambled UUCCG"),
                              plotXLabel=T,
                              breakWidth=4,
                              heightL=1.5)
  
  augData$motif <- factor(augData$motif, levels = c("no motif", "AUG", "AUG+SP"))
  plotMotifValidationBoxplots(df=augData, 
                              selectedLevels=c("no motif", "AUG", "AUG+SP"), 
                              legendLabels=c("No \n AUG", "Out-frame \n uAUG", "Out-frame \n uORF"),
                              breakWidth=0.8,
                              plotXLabel=T,
                              setBreaks=T,
                              breaksSeq=c(c(0, 0.1, 0.2, 0.4, 0.6, 1, 1.8), seq(2.4, 7, 1.8)),
                              heightL=1.8)
  