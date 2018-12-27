  source('./Conf.R')
  source("./Main.R")
  
  ######################## 3A, 3B
  
  RNA_Protein_Combined <- readRDS(paste0(tempFileDir, 
                                         "RNA_Protein_CombinedWithFeatures.Rds"))
  
  allfeatures <- getAllFeatureColumns(RNA_Protein_Combined)
  
  ptrEf <- getCodonAAVals(RNA_Protein_Combined,
                          responseStr="_PTR",
                          rLab=" protein-to-mRNA ratio",
                          includeAA=T,
                          allfeatures)
  Figure3A <- ptrEf$p1
  
  Figure3B <- ptrEf$p2 + theme(axis.text.x = element_text(size=6, hjust = 1, vjust=0.5))
  
  
  ######################## 3C
  
  humanCodonFreq <- read.csv(paste0(tempFileDir, "/Human_codon_frequency_tAI.csv"))
  t <- readRDS(paste0(tempFileDir, "/AllCodonDwellingTimes.rds"))
  cCols <- grep("RUST_|Tuller_", colnames(t))
  t$medianCodonDwel <- apply(t[,cCols],1,median)
  
  humanCodonFreq <- merge(humanCodonFreq, t, by="Codon")
  humanCodonFreq <- as.data.table(humanCodonFreq)
  
  codonDwellTimes <- subset(humanCodonFreq, select = grep("RUST", names(humanCodonFreq))) 
  codonDwellTimes[,AA := humanCodonFreq$AA,][,genomicCodonFreq := humanCodonFreq$genomicCodonFreq,]
  
  codDwell <- colnames(codonDwellTimes[,grep("RUST", colnames(codonDwellTimes)), with=F])
  cdDwellByAA <- codonDwellTimes[,lapply(.SD,function(x){sum(x*genomicCodonFreq)/sum(genomicCodonFreq)}),by=AA, .SDcols=codDwell]
  humanCodonFreq[,medianCodonDwellByAA := median(medianCodonDwel),by=AA]
  humanCodonFreq[,meanCodonDwellByAA := sum(medianCodonDwel*genomicCodonFreq)/sum(genomicCodonFreq),by=AA]
  
  humanCodonFreq <- humanCodonFreq[order(humanCodonFreq$meanCodonDwellByAA),]
  humanCodonFreq$AA <- factor(humanCodonFreq$AA, levels = unique(humanCodonFreq$AA))
  
  Figure3C <- ggplot(humanCodonFreq, aes(x=AA, y=medianCodonDwel)) + 
    geom_point( cex=0.5)+
    geom_text(aes(label=Codon),
              hjust=+0.6, 
              vjust=0,
              size=2,
              check_overlap = F)+
    geom_point(aes(x=AA, y=meanCodonDwellByAA,col="red"))+
    theme(legend.position = "None")+
    xlab("")+
    ylab("Median codon dwelling time based on \n 17 ribosome profiling data sets (z-score)")
  
  
  
  ######################## 3D
  
  s <- data.frame(Codon=names(ptrEf$codonMedianBetas),
                  PTRCodonEffects=ptrEf$codonMedianBetas)
  
  humanCodonFreq <- read.csv(paste0(tempFileDir, "/Human_codon_frequency_tAI.csv"))
  t <- readRDS(paste0(tempFileDir, "/AllCodonDwellingTimes.rds"))
  humanCodonFreq <- merge(humanCodonFreq, t, by="Codon")
  
  allCodon <- merge(humanCodonFreq, s, by="Codon")
  
  Figure3D <- plotCors(allCodon,
                x="medianCodonDwel",
                y="PTRCodonEffects",
                xLab="Median codon dwelling time \n across 17 Ribo-Seq data estimates (z-score)",
                yLab="2 fold codon-frequency increase \n effect on PTR ratio",
                labelStr = "Codon")
  
  
  ######################## 3E
  
  kusterProteinHalfLife <- readRDS(paste0(tempFileDir, "/KusterProteinHalfLife.rds"))
  aaCols <- grep("AminoAcid_", colnames(kusterProteinHalfLife))
  
  halfLifeEf <- plotAAAndCodonEffect(inputDF=kusterProteinHalfLife,
                                     aaCols,
                                     aaCols,
                                     y=kusterProteinHalfLife$ProteinHalfLife,
                                     responseLabel= " HeLa cells protein half-life")
  
  proHL <- readRDS(paste0(tempFileDir, "MathiesonProteinHLData.rds"))
  allfeaturesproHL <- getAllFeatureColumns(proHL)
  
  proHLEf <- getCodonAAVals(proHL,
                            responseStr="_ProHL",
                            rLab=" protein half-life",
                            allfeatures=allfeaturesproHL)
  
  pH2 <- proHLEf$aaBetas
  a <- halfLifeEf$aaRes$AAeffects
  names(a) <- halfLifeEf$aaRes$aa
  a <- a[colnames(pH2)]
  pH2 <- rbind(pH2, Hela_ProHL=a)
  
  pH2 <- pH2[3:5,]
  
  pH2Med <- apply(pH2,2,median)
  pH2Med <-pH2Med[order(pH2Med)]
  aaOrder <- names(pH2Med)
  
  pH2Melted <- melt(pH2)
  pH2Melted$variable <- factor(pH2Melted$variable, levels = aaOrder)
  
  Figure3E <- ggplot(pH2Melted, aes(x=variable,y=value, fill=variable))+
    geom_boxplot(outlier.shape = NA, width=0.5, lwd=0.3) + 
    scale_fill_manual(values=rep(c("#999999","#E69F00"),10))+
    ylab(paste0("2 fold amino acid frequency increase \n effect on protein half-life"))+
    xlab("")+
    theme(panel.spacing = unit(0, "lines"), 
          axis.text.x = element_text(size=7, hjust = 1, vjust=0.5), 
          plot.title = element_text(size=7, hjust = 0.5), 
          legend.position = "None", 
          strip.text = element_text(size=7))+
    geom_hline(yintercept=1, color="darkgrey")+
    scale_y_continuous(trans='log2', breaks = seq(0.2,2.0,0.1) , labels = seq(0.2,2,0.1))
  
  
  
  ######################## 3F
  pH2_2 <- data.frame(meanProtHLAAEffects=sapply(pH2, median),
                    aa=colnames(pH2))
  
  tp = data.frame(aa=names(ptrEf$aaMedianBetas), PTRAAEffects=ptrEf$aaMedianBetas)
  
  allAA <- merge(tp,pH2_2, by="aa")
  
  Figure3F <- plotCors(as.data.frame(allAA),
                x="meanProtHLAAEffects",
                y="PTRAAEffects",
                xLab="2 fold amino acid frequency increase \n effect on protein half-life",
                yLab="2 fold amino acid frequency increase \n effect on PTR ratio",
                labelStr = "aa")
  
  
  
  ######################## 3H
  
  Figure3H <- plotCors(allCodon,
                x="tAI",
                y="PTRCodonEffects",
                xLab="Codon tRNA adaptiveness",
                yLab="2 fold codon frequency increase \n effect on PTR ratio",
                labelStr = "Codon")
  
 