  source('./Conf.R')
  source("./Main.R")
  
  #### Fig EV5A
  aaProp <-  read.table(file = paste0(tempFileDir,"N_term_halfLife.tsv"), header = T)
  rownames(aaProp) <- aaProp$AA
  
  proHL <- readRDS(paste0(tempFileDir, "MathiesonProteinHLData.rds"))
  allfeaturesproHL <- getAllFeatureColumns(proHL)
  
  proHLEf <- getCodonAAVals(proHL,
                            responseStr="_ProHL",
                            rLab=" protein half-life",
                            allfeatures=allfeaturesproHL)
  
  pH2 <- proHLEf$aaBetas
  pH2 <- data.frame(meanProtHLAAEffects=sapply(pH2, median),
                    aa=colnames(pH2))
  
  rownames(pH2) <- toupper(aaProp[pH2$aa,"N_end_aa"])
  
  
  popDB <- data.table(read.table(file = file.path(tempFileDir,"PoP_DB.txt")))
  popDB[,c("WildAA","MutatedAA"):=tstrsplit(V1,"-"),]
  popDB[,WildAA := substring(WildAA, nchar(WildAA)-2,nchar(WildAA)),]
  popDB[,MedianChange := median(V2),by=c("WildAA","MutatedAA")]
  popDBMedian <- unique(popDB[,c("WildAA","MutatedAA","MedianChange")])
  
  phalfLifeEf <- data.frame(gtools::permutations(n=nrow(pH2),
                                                 r=2,
                                                 v=rownames(pH2),
                                                 repeats.allowed=F), stringsAsFactors = F)
  colnames(phalfLifeEf) <- c("WildAA","MutatedAA")
  phalfLifeEf$halfLifeEfChange <- pH2[phalfLifeEf$WildAA,"meanProtHLAAEffects"]- pH2[phalfLifeEf$MutatedAA,"meanProtHLAAEffects"]
  
  popDBMedian <- merge(popDBMedian,
                       phalfLifeEf,
                       by=c("WildAA","MutatedAA"))
  rownames(aaProp) <- toupper(aaProp$N_end_aa)
  popDBMedian$WildAA <- aaProp[popDBMedian$WildAA,"aa"]
  popDBMedian$MutatedAA <- aaProp[popDBMedian$MutatedAA,"aa"]
  popDBMedian$label <- paste0(popDBMedian$WildAA,"->",popDBMedian$MutatedAA)
  
  Fig_EV5A <- plotCors(inDF=as.data.frame(popDBMedian),
                x="halfLifeEfChange",
                y="MedianChange",
                xLab="Amino acid effect difference on protein half-life",
                yLab="Median amino acid mutation effects \n on protein stability",
                labelStr="label")
  
  
  #### Fig EV5B
  
  t <- as.data.frame(merge( pH2, aaProp, by="aa"))
  
  Fig_EV5B <- plotCors(inDF=t,
                x="Hydrophobicity_values",
                y="meanProtHLAAEffects",
                xLab="Amino acid hydrophobicity",
                yLab="2 fold amino acid frequency increase \n effect on protein half-life",
                labelStr="aa")
  
  
  
  #### Fig EV5C
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allFeatures <- getAllFeatureColumnsWithProteinPTM(RNA_Protein_Combined)
  
  PTRexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined,
                                            isAAIncluded=T,
                                            responseSubstring="_PTR",
                                            ylabStr=" PTR ratio",
                                            allFeatures)
    
  
  ptr<- plotResultFigures(expVarsOfFeatures=PTRexpVar$expVarsOfFeatures,
                          mainStr="Explained variance in protein to mRNA ratio (%)")
  Fig_EV5C <- ptr$p+ scale_y_continuous(breaks = seq(0,30,5), labels = seq(0,30,5))
  
  
  #### Fig EV5 D
  proHL <- readRDS(paste0(tempFileDir, "MathiesonProteinHLData.rds"))
  rownames(proHL) <- proHL$EnsemblGeneID
  c <- grep("Motif_LinAA_", colnames(proHL))
  proHL[,c] <- NULL
  
  kusterProteinHalfLife <- readRDS(paste0(tempFileDir, "/KusterProteinHalfLife.rds"))
  proHL <- merge(kusterProteinHalfLife, proHL, by="EnsemblGeneID", all.x=T, all.y=T)
  
  a <- grep("_ProHL|ProteinHalfLife|EnsemblGeneID", colnames(proHL))
  proHL <- proHL[,a]
  
  colnames(proHL) <- c("EnsemblGeneID", "BCells_ProHL", "Hepatocytes_ProHL", "Monocytes_ProHL")
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  
  proHL <- merge(RNA_Protein_Combined, proHL, by="EnsemblGeneID")
  allFeatures <- getAllFeatureColumnsWithProteinPTM(proHL)
  
  ProHLexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined = proHL,
                                           isAAIncluded=T,
                                           responseSubstring=c("BCells_ProHL",
                                                               "Hepatocytes_ProHL",
                                                               "Monocytes_ProHL"),
                                           ylabStr=" protein half-life",
                                           allFeatures,
                                           isTissues=F)
  
  
  pro <- plotResultFigures(expVarsOfFeatures=ProHLexpVar$expVarsOfFeatures,
                           mainStr="Explained variance in protein half-life (%)")
  
  Fig_EV5D <- pro$p+ scale_y_continuous(breaks = seq(0,25,5), labels = seq(0,25,5))
  
  #### Fig EV5E
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  allFeatures <- getAllFeatureColumnsWithProteinPTM(RNA_Protein_Combined)
  
  mRNAexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined,
                                           isAAIncluded=T,
                                           responseSubstring="_medianExon",
                                           ylabStr=" mRNA",
                                           allFeatures)
  
  
  mRNA <- plotResultFigures(expVarsOfFeatures=mRNAexpVar$expVarsOfFeatures,
                            mainStr="Explained variance in mRNA levels (%)")
  Fig_EV5E <- mRNA$p+ scale_y_continuous(breaks = seq(0,40,5), labels = seq(0,40,5))
  
  
  #### Fig EV5F
  
  #Hek293Cells_mRNAHL
  Hek293MHL <- readRDS(paste0(tempFileDir,"/Hek293_mRNAHL.rds"))
  Hek293MHL <- Hek293MHL[,c("Hek293Cells_mRNAHL", "EnsemblTranscriptID", "EnsemblGeneID")]
  
  #HeLa cells
  helaMHL <- readRDS(paste0(tempFileDir,"/humanHELAMRNAHL.rds"))
  helaMHL <- helaMHL[,c("HeLamRNAHL","EnsemblTranscriptID", "EnsemblGeneID")]
  colnames(helaMHL) <- c("HeLa_mRNAHL","EnsemblTranscriptID", "EnsemblGeneID")
  
  
  #K562 cells
  k562 <- readRDS(paste0(tempFileDir,"/k562_features.rds"))
  k562 <- k562[, c("EnsemblGeneID", "half_life_mean")]
  colnames(k562) <- c("EnsemblGeneID", "K562_mRNAHL")
  
  mHL <- merge(helaMHL, Hek293MHL, by=c("EnsemblGeneID","EnsemblTranscriptID"), all.x=T)
  mHL <- merge(mHL, k562, by="EnsemblGeneID", all.x=T)
  mHL <- mHL[!duplicated(mHL$EnsemblGeneID),]
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, "RNA_Protein_CombinedWithFeatures.Rds"))
  mHL <- merge(mHL, RNA_Protein_Combined, by=c("EnsemblGeneID","EnsemblTranscriptID"))
  
  allFeatures <- getAllFeatureColumnsWithProteinPTM(mHL)
  
  mRNAHLexpVar <- plotFeaturesExplainedVariance(RNA_Protein_Combined = mHL,
                                           isAAIncluded=T,
                                           responseSubstring=c("Hek293Cells_mRNAHL",
                                                               "HeLa_mRNAHL"),
                                           ylabStr=" mRNA half-life",
                                           allFeatures,
                                           isTissues=F)
  
  
  mrna <- plotResultFigures(expVarsOfFeatures=mRNAHLexpVar$expVarsOfFeatures,
                            mainStr="Explained variance in mRNA half-life (%)")
  
  Fig_EV5F <- mrna$p+ scale_y_continuous(breaks = seq(0,35,5), labels = seq(0,35,5))
