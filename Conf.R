peptideFileName="peptides.txt"
tempFileDir="./DATA/"

TissueNames = c("Adrenal", "Appendices", "Brain", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopiantube", "Fat", "Gallbladder", "Heart", "Kidney", "Liver", "Lung", "Lymphnode", "Ovary", "Pancreas", "Placenta","Prostate", "Rectum", "Salivarygland", "Smallintestine", "Smoothmuscle", "Spleen", "Stomach", "Testis", "Thyroid", "Tonsil", "Urinarybladder")
TissueNamesPrint = c("Adrenal gland", "Appendix", "Brain", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian tube", "Fat", "Gallbladder", "Heart", "Kidney", "Liver", "Lung", "Lymph node", "Ovary", "Pancreas", "Placenta","Prostate", "Rectum", "Salivary gland", "Small intestine", "Smooth muscle", "Spleen", "Stomach", "Testis", "Thyroid", "Tonsil", "Urinary bladder")

annotationFileName=paste0(tempFileDir,"/Homo_sapiens.GRCh38.83.gtf")

motifDB <- data.frame(read.csv(paste0(tempFileDir,"motifDB.csv")), stringsAsFactors = F)
mRNARegions <- c("UTR3","UTR5","CDS")
AA_motifs=c("KRR","NS")
CDS_CtermAAMotifs=c("CG")
stopCodons=c("TAG", "TAA", "TGA")

## Dynamic range for RNA and proteins 
Dynamic_Range = c(0.1, 0.9)
ColSequential="PuRd"
ColDiverging="PiYG"
ColPink="#fde0ef"
ColLight="#e6f5d0"

panelSpacingUnit = 0.5
LABFCTDIR="/Users/basak/WORKSPACE/gagneurlab_shared/r/"
fontFamily="Helvetica"
GoDir="/Users/basak/WORKSPACE/gagneurlab_shared/GO/"

tmpnull= sapply(list.files(file.path(LABFCTDIR,'/plotters/'), full.names = T), source)
tmpnull= sapply(list.files(file.path(LABFCTDIR,'/knitr_helper/'), full.names = T), source)
tmpnull= sapply(list.files(file.path(LABFCTDIR,'/go_enrichment/'), full.names = T), source)
tmpnull= sapply(list.files(file.path(LABFCTDIR,'/disease/'), full.names = T), source)

MGSA_GO_FULL <- load_mgsaset_for_organism(organism='human', folderPath=GoDir)

featureCols <- data.frame(read.csv(paste0(tempFileDir, 
                                          "FeatureColors.csv"), sep = ","), 
                          stringsAsFactors = F)
rownames(featureCols) <- featureCols$FeatureID
