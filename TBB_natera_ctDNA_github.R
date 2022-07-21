# set directory
setwd("~/Documents/Clinical_topics/TBB/TBB_data/natera/ctDNA_analysis")

#load data
signatera <- read.csv("SignateraRUO_Report_Stanford_BeyondBRCA_29SEP2020.csv", header = T)
pwes <- read.csv("20201229_pWES_all_data.csv", header=T)
Sumpwes <- read.csv("20201229_pWES_summary.csv", header=T)


# clean up notation
# note: TP=timepoint (TP1 = baseline, TP2 = progression)
library(stringr)
TP1 <- str_replace(pwes$X..Plasma.TP1.VAF, pattern="%", "")
pwes$X..Plasma.TP1.VAF <- as.numeric(TP1)/100

TP2 <- str_replace(pwes$X..Plasma.TP2.VAF, patter="%", "")
pwes$X..Plasma.TP2.VAF <- as.numeric(TP2)/100

# change NAs to small value for visualization
library(tidyr)
pwes2 <- pwes
TP1_NAreplace <- replace_na(data=pwes2$X..Plasma.TP1.VAF, replace = 0.009)
pwes2$X..Plasma.TP1.VAF <- TP1_NAreplace

TP2_NAreplace <- replace_na(data=pwes2$X..Plasma.TP2.VAF, replace = 0.009)
pwes2$X..Plasma.TP2.VAF <- TP2_NAreplace

#remove TBB-B-010
head(pwes2)
pwes2 <- pwes2[which(pwes2$Primary.Institute.Patient.ID != "TBB-B-010"), ]
Sumpwes <- Sumpwes[which(Sumpwes$Primary.Institute.Patient.ID != "TBB-B-010"), ]
write.csv(pwes2, file="pwes2.csv")

# set colors
library(RColorBrewer)
library(scales)
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
palette(mycolors)
dim(Sumpwes)
palette()

# plot pWES SNVs
# set up pchs
pwesPlot <- pwes
levels(pwesPlot$Primary.Institute.Patient.ID) <- rep(15:20,5)
pchs<-as.numeric(as.character(pwesPlot$Primary.Institute.Patient.ID))
# construct plot
dev.off()
par(mar=c(5.1, 4.1, 2.1, 13.1), xpd=NA)
plot(pwes2$X..Plasma.TP2.VAF, pwes2$X..Plasma.TP1.VAF, 
     col=pwes$Primary.Institute.Patient.ID, log="xy", pch=pchs,
     xlab="plasma VAF at progression", 
     ylab="plasma VAF at baseline", 
     main="pWES SNVs")
abline(0, 1, xpd=F, lty=2, lwd=3, col="grey")
abline(0.5, 1, xpd=F, lty=2)
abline(-0.5, 1, xpd=F, lty=2)
anno <- paste(levels(pwes2$Primary.Institute.Patient.ID), Sumpwes$genetics, sep = ", ")
legend(0.85, 0.70, legend=Sumpwes$genetics, col=mycolors, border=NA,
       fill=NA, title="Patient Mutation(s)", pch=c(15:20))

# select the variants associated with progression (5-fold or 10-fold increased) & NAs
pwes2$logdiff <- log10(pwes2$X..Plasma.TP2.VAF) - log10(pwes2$X..Plasma.TP1.VAF)
pwes2$logdiff
tenFoldwNAs <- pwes2[which(pwes2$logdiff > 1), ]
tenFoldwNAs
dim(tenFoldwNAs)
fiveFoldwNAs <- pwes2[which(pwes2$logdiff > 0.5), ]
dim(fiveFoldwNAs)
negFiveFoldwNAs <- pwes2[which(pwes2$logdiff < -0.5), ]
dim(negFiveFoldwNAs)
write.csv(tenFoldwNAs, file = "tenfoldwNAs.csv")
write.csv(fiveFoldwNAs, file = "fivefoldwNAs.csv")
write.csv(negFiveFoldwNAs, file = "negfivefoldwNAs.csv")


#plot responders SNVs (TBB 3, 6, 7, 16, 19, 22), (gPALB2, gPALB2, gCHEK2/FANCA/sPTEN, gPALB2, gPALB2, gPALB2 -- respectively)
responders <- read.csv("fivefoldwNAs_responders_deleterious.csv", header=T)
tissue <- str_replace(responders$X..Tissue.VAF, patter="%", "")
responders$X..Tissue.VAF <- as.numeric(tissue)/100
set1cols <- c("#4178AE", "#4AA956", "#678C69", "#AF5B39", "#F182BC")
palette(set1cols)
par(mar=c(5.1, 4.1, 2.1, 13.1), xpd=NA)
plot(responders$X..Plasma.TP2.VAF, pch=1, col=as.factor(responders$Natera.Patient.ID), 
     xlab="mutation", ylab="VAF", log="y", ylim=c(0.01, 0.8), 
     main="Deleterious SNVs enriched at progression in responders")
points(responders$X..Plasma.TP1.VAF, pch=2, col=as.factor(responders$Natera.Patient.ID), 
     add=T)
points(responders$X..Tissue.VAF, pch=15, col=as.factor(responders$Natera.Patient.ID), 
     add=T)
text(x=c(50, 51), y=c())
legend(145, 0.95, legend=c("gPALB2", "gPALB2 *", "gCHEK2/gFANCA/sPTEN", 
"gPALB2", "gPALB2"), fill=set1cols, title="Patient Mutation(s)")
legend(145, 0.2, legend=c("tumor biopsy", 
                           "ctDNA-progression", "ctDNA-baseline"), 
       pch=c(15, 1, 2))

# same analysis for non-responders 
nonres <- read.csv("fivefoldwNAs_nonresponders_deleterious.csv", header=T)
nrtissue <- str_replace(nonres$X..Tissue.VAF, patter="%", "")
nonres$X..Tissue.VAF <- as.numeric(nrtissue)/100
setNRcols <- brewer.pal(9, "Set1")
palette(setNRcols)
par(mar=c(5.1, 4.1, 2.1, 13.1), xpd=NA)
plot(nonres$X..Plasma.TP2.VAF, pch=1, col=as.factor(nonres$Natera.Patient.ID), 
     xlab="mutation", ylab="VAF", log="y", ylim=c(0.01, 0.8), 
     main="Deleterious SNVs enriched at progression in non-responders")
points(nonres$X..Plasma.TP1.VAF, pch=2, col=as.factor(nonres$Natera.Patient.ID), 
     add=T)
points(nonres$X..Tissue.VAF, pch=15, col=as.factor(nonres$Natera.Patient.ID), 
     add=T)
text(x=c(50, 51), y=c())
legend(105, 0.95, legend=c("gBRIP1", "sPTEN", "gPALB2/gBRIP1", 
"gCHEK2", "sPTEN", "sPTEN", "sRAD50", "gATM", "gATM"), fill=setNRcols, title="Patient Mutation(s)")
legend(105, 0.075, legend=c("tumor biopsy", 
                           "ctDNA-progression", "ctDNA-baseline"), 
       pch=c(15, 1, 2))

# plot signatera results
sigVars <- signatera$Variant_AF_In_Plasma....
signatera$Variant_AF_In_Plasma.... <- as.numeric(sigVars)/100
signatTP1 <- signatera[which(signatera$Timepoint_ID == 1), ]
signatTP1 <- signatTP1[which(signatTP1$Primary_Institute_Patient_ID != "TBB-B-010"), ]
signatTP2 <- signatera[which(signatera$Timepoint_ID == 2), ]
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
palette(mycolors[c(2:4, 6:7, 10:18)])
par(mar=c(5.1, 4.1, 2.1, 14.1), xpd=NA)
options(scipen=999)
levels(signatTP1$Primary.Institute.Patient.ID)
levels(signatTP1$Primary_Institute_Patient_ID) <- rep(15:20,5)
SigPchs<- as.numeric(as.character(signatTP1$Primary_Institute_Patient_ID))
class(SigPchs)
plot(y = signatTP1$Variant_AF_In_Plasma....+0.0009, 
     x = signatTP2$Variant_AF_In_Plasma....+0.0009, log="xy", 
     xlab = "plasma VAF at progression", 
     ylab = "plasma VAF at baseline", 
     col=signatTP1$Primary_Institute_Sample_ID, 
     pch=SigPchs, main="Signatera variants", 
     xlim=c(0.001, 1), ylim=c(0.001,1))
legend(1.4, 1.3, legend=Sumpwes$genetics[c(2:4, 6:7, 10:18)], fill=NA, border=NA,
       col=mycolors[c(2:4, 6:7, 10:18)], title="Patient Mutation(s)", pch=c(15:20))
abline(0, 1, xpd=F, lty=2, lwd=3, col="grey")
abline(0.5, 1, xpd=F, lty=2)
abline(-0.5, 1, xpd=F, lty=2)



# read in indels
tissueindels <- read.csv("indels_tissue.csv", header=T)
pindels1 <- read.csv("indels_plasma_TP1.csv", header=T)
pindels2 <- read.csv("indels_plasma_TP2.csv", header=T)

# merge indel files based on Gene
tissueindels$newID <- paste(tissueindels$Stanford.ID, tissueindels$Gene.ID, sep=".")
pindels1$newID <- paste(pindels1$Stanford.ID, pindels1$Gene.ID, sep=".")
pindels2$newID <- paste(pindels2$Stanford.ID, pindels2$Gene.ID, sep=".")
pindelsmerge <- merge(x=pindels1, y=pindels2, by="newID", all=T, sort=F, 
                      suffixes=c(".1", ".2"))
#save merged file
write.csv(pindelsmerge, "pindelsmergeByGene.csv")

# convert VAF % to fraction
TP1 <- str_replace(pwes$X..Plasma.TP1.VAF, pattern="%", "")
pindelsmerge$pWES.TP1.VAF <- as.numeric(str_replace(pindelsmerge$pWES.TP1.VAF, pattern="%", ""))/100
pindelsmerge$pWES.TP2.VAF <- as.numeric(str_replace(pindelsmerge$pWES.TP2.VAF, pattern="%", ""))/100

#NA replace
pindels1_NAreplace <- replace_na(data=pindelsmerge$pWES.TP1.VAF, replace = 0)
pindelsmerge$pWES.TP1.VAF <- pindels1_NAreplace

pindels2_NAreplace <- replace_na(data=pindelsmerge$pWES.TP2.VAF, replace = 0)
pindelsmerge$pWES.TP2.VAF <- pindels2_NAreplace

# set up color scheme
pindelsmerge$Stanford.ID.1 <- as.character(pindelsmerge$Stanford.ID.1)
pindelsmerge$Stanford.ID.2 <- as.character(pindelsmerge$Stanford.ID.2)
pindelsmerge$Stanford.ID.1[201:381] <- pindelsmerge$Stanford.ID.2[201:381]
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
palette(mycolors)

#remove TBB-B-010
pindelsmerge <- pindelsmerge[which(pindelsmerge$Stanford.ID.1 != "TBB-B-010"), ]
pindelsmerge <- pindelsmerge[which(pindelsmerge$Stanford.ID.1 != "TBB-B-004"), ]

# set Pch variable
IndelPchs <- as.factor(pindelsmerge$Stanford.ID.1)
levels(IndelPchs) <- rep(15:20,5)
IndelPchs<- as.numeric(as.character(IndelPchs))

# scatterplot indels
plot(y = pindelsmerge$pWES.TP1.VAF+0.01, 
     x = pindelsmerge$pWES.TP2.VAF+0.01,
     xlab = "plasma VAF at progression", log="xy",
     ylab = "plasma VAF at baseline", 
     col=as.factor(pindelsmerge$Stanford.ID.1), 
     pch=IndelPchs, main="pWES indels",
     xlim=c(0.01, 1), ylim=c(0.01,1))
legend(1.3, 1.2, legend=Sumpwes$genetics, fill=NA, border=NA,
       col=mycolors, title="Patient Mutation(s)", pch=c(15:20))
abline(0, 1, xpd=F, lty=2, lwd=3, col="grey")
abline(0.5, 1, xpd=F, lty=2)
abline(-0.5, 1, xpd=F, lty=2)

# graph SNVs and indels differing between baseline and progression
# focus on breast gPALB2
negFFsnvgP <- negFiveFoldwNAs[which(negFiveFoldwNAs$Primary.Institute.Patient.ID %in%
                                           c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
FFsnvgP    <- fiveFoldwNAs[which(fiveFoldwNAs$Primary.Institute.Patient.ID %in% 
                                           c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
negFFindelsgP <- negfiveFoldindels[which(negfiveFoldindels$Stanford.ID.1 %in% 
                                                  c("TBB-B-003", "TBB-B-006", 
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
FFindelsgP    <- fiveFoldindels[which(fiveFoldindels$Stanford.ID.1 %in% 
                                                c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]


# non-Breast non-gPALB2
negFFsnvNgP <- negFiveFoldwNAs[which(!negFiveFoldwNAs$Primary.Institute.Patient.ID %in%
                                           c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
FFsnvNgP    <- fiveFoldwNAs[which(!fiveFoldwNAs$Primary.Institute.Patient.ID %in% 
                                           c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
negFFindelsNgP <- negfiveFoldindels[which(!negfiveFoldindels$Stanford.ID.1 %in% 
                                                  c("TBB-B-003", "TBB-B-006", 
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]
FFindelsNgP    <- fiveFoldindels[which(!fiveFoldindels$Stanford.ID.1 %in% 
                                                c("TBB-B-003", "TBB-B-006",  
                                             "TBB-B-016", "TBB-B-019", "TBB-B-022")), ]


#extract Breast/gPALB2 data for prism

       i= as.data.frame(table(droplevels(negFFsnvNgP$Primary.Institute.Patient.ID))) 
       j= as.data.frame(table(droplevels(FFsnvNgP$Primary.Institute.Patient.ID)))
       k= as.data.frame(table(negFFindelsNgP$Stanford.ID.1))
       l= as.data.frame(table(FFindelsNgP$Stanford.ID.1))

NgPsnvTable1 <- merge(i, j, by="Var1", all=T)
NgPsnvTable2 <- merge(k, l, by="Var1", all=T)
NgPsnvTable <- merge(NgPsnvTable1, NgPsnvTable2, by="Var1", all=T)
colnames(NgPsnvTable) <- c("patients", "negFFsnvNgP", "FFsnvNgP", "negFFindelsNgP", "FFindelsNgP")

        m= as.data.frame(table(droplevels(negFFsnvgP$Primary.Institute.Patient.ID)))
        n= as.data.frame(table(droplevels(FFsnvgP$Primary.Institute.Patient.ID)))
        o= as.data.frame(table(negFFindelsgP$Stanford.ID.1))
        p= as.data.frame(table(FFindelsgP$Stanford.ID.1))
gPSnvTable1 <- merge(m, n, by="Var1", all=T)
gPSnvTable2 <- merge(o, p, by="Var1", all=T)
gPSnvTable  <- merge(gPSnvTable1, gPSnvTable2, by="Var1", all=T)
colnames(gPSnvTable) <- c("patients", "negFFsnvgP", "FFsnvgP", "negFFindelsgP", "FFindelsgP")
Breast_gPALB2_SnvTable <- merge(NgPsnvTable, gPSnvTable, by="patients", all=T)
write.csv(Breast_gPALB2_SnvTable, "Breast_gPALB2_variantsTable.csv")
# further statistical testing performed in prism

# ctDNA HRD with scar
scar <- read.csv("scarHRD scores.csv", header=T)
scarTumor <- scar[which(scar$Timepoint == 0), ]
scarTP1 <- scar[which(scar$Timepoint == 1), ]
scarTP2 <- scar[which(scar$Timepoint == 2), ]

# look at how scarHRD compares to myriad HRD on tumors
plot(x=scarTumor$HRD.sum, y=scarTumor$hrd_score_prim)
cor.test(x=scarTumor$HRD.sum, y=scarTumor$hrd_score_prim, method="pearson",
         use="complete.obs") # r=0.85, p=0.004
# conclusion: scarHRD on tumor is closely related to myriadHRD on tumor

# match the natera sample to the sample used for myriad HRD
plot(x=scarTumor$HRD.sum, y=scarTumor$myriad_matchToNatera, xlim=c(0, 70), ylim=c(0,70))
cor.test(x=scarTumor$HRD.sum, y=scarTumor$myriad_matchToNatera, 
         method="pearson", use="complete.obs") # r=0.58, p=0.028

# plot pWES-HRD versus clinical outcomes
# first ensure complete cases (remove TBB004, 010, 011)
scarTP1 <- scarTP1[-which(scarTP1$Patient %in% c("CI004", "CI010", "CI011")), ]
plot(scarTP1$HRD.sum, scarTP2$HRD.sum, xlim=c(0, 70), ylim=c(0,70))
cor.test(scarTP1$HRD.sum, scarTP2$HRD.sum, method="pearson", 
         use="complete.obs") # r= 0.808, p = 0.00008661


# certain ctDNA samples have low complexity (low VAF) = TBB-B-008, 020, 013, 010, 002, 005
# aka CI008, CI019, CI013, CI010, CI002, CI005
# remove these samples
scarTP1best <- scarTP1[-which(scarTP1$Patient %in% c("CI008", "CI019", "CI013", 
                                                     "CI002", "CI005")), ]

# plot ctDNA HRD versus TTP
par(xpd=T)
plot(scarTP1best$HRD.sum, scarTP1best$TTP, col=as.character(scarTP1best$colorhex), pch=16, cex=2.2, 
     xlab = "plasma WES HRD score (ctHRD)", ylab= "Time to Progression (days)")
legend(x=62, y=399, legend=c("breast", "colon", "pancreas", "parotid", "testicular", "uterine"), 
       fill=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F" ), 
       cex=1)
cor.test(scarTP1best$HRD.sum, scarTP1best$TTP, method="pearson", 
         use="complete.obs")  # r = 0.61, p= 0.037
text(x=15, y=375, labels=" r=0.61")
text(x=16, y=360, labels="p=0.037")
text(x=scarTP1best$HRD.sum, y=scarTP1best$TTP, labels=as.character(scarTP1best$genetics), 
     pos=   c(1, 3, 1, 2, 1, 4, 3, 4, 2, 1, 4, 1), 
     offset=rep(0.75, 12))
lm(scarTP1best$TTP ~ scarTP1best$HRD.sum)
par(xpd=F)
abline(79.397, 4.192, xlim=c(10, 60), lty=2, col="gray")





# create Venn Diagrams of pWES versus tWES mutations
install.packages("BioVenn")
library(BioVenn)
colnames
twesMuts <- pwes[which(is.na(pwes$X..Tissue.VAF) == F), ]
pwesTP1  <- pwes[which(is.na(pwes$X..Plasma.TP1.VAF) == F), ]
pwesTP2  <- pwes[which(is.na(pwes$X..Plasma.TP2.VAF) == F), ]

pts <- as.vector(unique(pwes$Primary.Institute.Patient.ID))
Sumpwes2 <- read.csv("20201229_pWES_summary.csv", header=T) #include sample 10
ptGen <- Sumpwes2$genetics

#this will generate 19 plots as pdf files
for (i in 1:19) {
    draw.venn(twesMuts[which(twesMuts$Primary.Institute.Patient.ID == pts[i]), ]$Mutation, 
          pwesTP1[which(pwesTP1$Primary.Institute.Patient.ID == pts[i]), ]$Mutation, 
          pwesTP2[which(pwesTP2$Primary.Institute.Patient.ID == pts[i]), ]$Mutation, 
          title = pts[i],t_fb=1, st_s=1.5, st_fb=1, nr_s=1.25, t_f="", nr_f="", nr_fb=1,
          xtitle = "", ytitle = "", ztitle = "", output="pdf", width=450, height=450, 
          filename=pts[i], subtitle=ptGen[i], st_f="")
}

# analyze mutational signatures
install.packages("mutSignatures")
install.packages("kableExtra")
library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
install.packages("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") 
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
hg19 <- BSgenome.Hsapiens.UCSC.hg19

# import plasma WES (pwes) data
ss <- separate(data=pwes,col="Mutation", into=c("CHROM", "POS", "REF", "ALT"), sep="_" )
ss <- filterSNV(dataSet = ss, seq_colNames = c("REF", "ALT"))
ss <- attachContext(mutData = ss,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)

ss <- removeMismatchMut(mutData = ss,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")  
ss <- attachMutType(mutData = ss,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 

ss.counts <- countMutTypes(mutTable = ss,
                             mutType_colName = "mutType",
                             sample_colName = "Primary.Institute.Patient.ID")

num.sign <- 2
ss.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 20,               # bootstrapping: usually 500-1000
    num_parallelCores = 4)                # total num of cores to use (parallelization)
ss.analysis <- 
  decipherMutationalProcesses(input = ss.counts,
                              params = ss.params)
head(ss.analysis)

# Retrieve signatures (results)
ss.sig <- ss.analysis$Results$signatures
as.data.frame(ss.sig)

# Retrieve exposures (results)
ss.exp <- ss.analysis$Results$exposures

# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(ss.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(ss.sig, signature = 2, ylim = c(0, 0.05), title= "")  # this is similar to HRD signature

# Export Signatures as data.frame
ss.xprt <- coerceObj(x = ss.sig, to = "data.frame") 
ss.exp.xprt <- as.matrix(coerceObj(x=ss.exp, to = "data.frame"))

#change colnames to genetics
colnames(ss.exp.xprt) <- ptGen

# regenerate export with TBB names
ss.exp.tbb <- as.matrix(coerceObj(x=ss.exp, to = "data.frame"))
   # then transform to percentage
   tbb_percentage <- apply(ss.exp.tbb, 2, function(x){x*100/sum(x,na.rm=T)})

# Transform this data to %
data_percentage <- apply(ss.exp.xprt, 2, function(x){x*100/sum(x,na.rm=T)})
d <- rbind(ss.exp.xprt, data_percentage)
rownames(d) <- c("Sign.01", "Sign.02", "s1", "s2")
d
e <- as.data.frame(t(d))

f <- e[order(-e$s2), ]
g <- as.data.frame(t(f))
g <- rename(g, "gCHEK2..gFANCA..sPTEN"="gCHEK2/gFANCA/sPTEN")
g <- rename(g, "gPALB2..gBRIP1" = "gPALB2/gBRIP1")
g
newcols <- gsub("\\..","",colnames(g))
colnames(g) <- newcols
g
g <- as.matrix(g)

dev.off()
par(mar=c(12,5,1,6), oma=c(1,1,0,0), xpd=NA)
barplot(g[3:4, ], 
        col=c("#ff7f00", "#a6cee3") , 
        border="white", 
        space=0.04, 
        font.axis=2, 
        ylab="signature percentage", las=2)
legend(x =20.5, y=80, title="signatures",
       c("gPALB2", "other"), fill=c("#a6cee3", "#ff7f00"), horiz=F, cex=1)

# use the HRD percentage to correlate with response variables
gmerge <- merge(x=t(tbb_percentage), y=Sumpwes2, by.x=0, by.y="Primary.Institute.Patient.ID")

# plot correlation
dev.off()
plot(x=gmerge$Sign.02, y=gmerge$best_change_SLD, xlab = "Percentage gPALB2 signature", 
     ylab = "Best Change in SLD by RECIST (%)", pch=19, ylim = rev(range(gmerge$best_change_SLD)))
cor.test(gmerge$Sign.02, -gmerge$best_change_SLD, method="pearson", 
         use="complete.obs")  # r = 0.78, p= 0.00007681 (for Spearman r=0.857, p=0.000002741)
lm(gmerge$best_change_SLD ~ gmerge$Sign.02)
par(xpd=F)
abline(a=40.8444, b=-0.8874, lty=2, col="gray")
text(x=5, y=-95, labels="      r=0.78")
text(x=12, y=-85, labels="p=0.000087")
# there is a highly significant correlation between mutational signatures (tWES+pWES) and response


# repeat signature analysis on tumor WES data 
tt <- separate(data=twesMuts,col="Mutation", into=c("CHROM", "POS", "REF", "ALT"), sep="_" )
head(tt)
tt <- filterSNV(dataSet = tt, seq_colNames = c("REF", "ALT"))
tt <- attachContext(mutData = tt,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = hg19)

tt <- removeMismatchMut(mutData = tt,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")  
tt <- attachMutType(mutData = tt,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context") 

tt.counts <- countMutTypes(mutTable = tt,
                             mutType_colName = "mutType",
                             sample_colName = "Primary.Institute.Patient.ID")

num.sign <- 2
tt.params <- 
  mutSignatures::setMutClusterParams( 
    num_processesToExtract = num.sign,    # num signatures to extract
    num_totIterations = 20,               # bootstrapping: usually 500-1000
    num_parallelCores = 4)                # total num of cores to use (parallelization)
tt.analysis <- 
  decipherMutationalProcesses(input = tt.counts,
                              params = tt.params)

# Retrieve signatures (results)
tt.sig <- tt.analysis$Results$signatures

# Retrieve exposures (results)
tt.exp <- tt.analysis$Results$exposures

# Plot signature 1, 2
msigPlot(tt.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(tt.sig, signature = 2, ylim = c(0, 0.10))

# Export Signatures as data.frame
tt.xprt <- coerceObj(x = tt.sig, to = "data.frame") 
tt.exp.xprt <- as.matrix(coerceObj(x=tt.exp, to = "data.frame"))

# Transform this data in %
data_percentage <- apply(tt.exp.xprt, 2, function(x){x*100/sum(x,na.rm=T)})
tt.d <- rbind(tt.exp.xprt, data_percentage)
rownames(tt.d) <- c("Sign.01", "Sign.02", "s1", "s2")
tt.d
tt.e <- as.data.frame(t(tt.d))
tt.f <- tt.e[order(-tt.e$s2), ]
tt.g <- as.data.frame(t(tt.f))
tt.g <- as.matrix(tt.g)


# use the HRD percentage to correlate with response variables
tt.gmerge <- merge(x=t(tt.g[3:4, ]), y=Sumpwes2, by.x=0, by.y="Primary.Institute.Patient.ID")

dev.off()
plot(tt.gmerge$s2, -tt.gmerge$best_change_SLD, xlab = "Percentage HRD signature", 
     ylab = "Best Change in SLD by RECIST (%)", pch=19)
cor.test(tt.gmerge$s2, -tt.gmerge$best_change_SLD, method="pearson", 
         use="complete.obs")  # r = 0.7856098 , p= 0.000518 

# Repeat correlation with gPALB2 samples removed
tt.gmergeRemoved <- tt.gmerge[which(tt.gmerge$genetics != "gPALB2"), ]
cor.test(tt.gmergeRemoved$s2, -tt.gmergeRemoved$best_change_SLD, method="pearson", 
         use="complete.obs")

