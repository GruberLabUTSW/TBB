# set directory
setwd("~/Documents/Clinical_topics/TBB/TBB_data")
install.packages("zoo")

# load data
tbb <- read.csv("TBB_waterfall_data.csv", header=T)
tbb
tbb <- tbb[order(tbb$best_change_SLD, decreasing = T), ]
tbb
tbb$type

# set color scheme
library(RColorBrewer)
levels(tbb$type)
brewer.pal(n = 6, name = "Set2")
palette(brewer.pal(n = 6, name = "Set2")) 

# waterfall plot
dev.off()
bp <- barplot(tbb$best_change_SLD, col=as.factor(tbb$type), ylab="Best Change in SLD by RECIST (%)", 
              ylim=c(-100,100))
legend(x=1, y=-40, legend=levels(as.factor(tbb$type)), 
       fill=brewer.pal(n = 6, name = "Set2"))
abline(h=-30, lty=2, col="grey", lwd=2)
text(x=bp+0.25, y=tbb$label2+12, labels=tbb$genetics, cex=1, pos=3, srt=90, adj=0)


# Swimmer's plots
par(mar=c(5, 12, 4, 2))
tbb
tbb2 <- tbb[order(tbb$last.first, decreasing = F), ]
tbb2
sp <- barplot(tbb2$last.first, col=as.factor(tbb2$type), horiz=T, xlab="Duration of Therapy (days)", 
              names.arg=tbb2$genetics, las=1)
points(x=tbb2$SD1, y=sp, pch=5, cex=1.5, lwd=2)
points(x=tbb2$SD2, y=sp, pch=5, cex=1.5, lwd=2)
points(x=tbb2$SD3, y=sp, pch=5, cex=1.5, lwd=2)
points(x=tbb2$SD4, y=sp, pch=5, cex=1.5, lwd=2)
points(x=tbb2$PR1, y=sp, pch=8, cex=1.5, lwd=2)
points(x=tbb2$PR2, y=sp, pch=8, cex=1.5, lwd=2)
points(x=tbb2$PR3, y=sp, pch=8, cex=1.5, lwd=2)
points(x=tbb2$PD1, y=sp, pch=15, cex=2, lwd=2)
points(x=tbb2$PD2, y=sp, pch=15, cex=2, lwd=2)
legend(x=200, y=sp[[6]], legend=c("SD", "PR", "PD"),
       pch=c(5,8, 15), 
       #col=c(4,3,2,"black", "black"), 
       pt.cex=c(1.5,1.5,2))
legend(x=250, y=sp[[7]], legend=levels(as.factor(tbb$type)), 
       fill=brewer.pal(n = 6, name = "Set2"))

# plot HRD scores
dev.off()
tbb2
dim(tbb2)
tbb3 <- tbb2[, c(1, 2, 4, 7, 25, 29, 30, 32, 33, 34)]
tbb3
tbb3 <- tbb3[complete.cases(tbb3), ]
plot(y=tbb3$best_change_SLD, x=tbb3$HRDbest, col=as.factor(tbb3$type),  
            ylab="Best Change in SLD by RECIST (%)", 
            pch=19, cex=2.5,
            las=1, xlab="HRD score", 
            ylim = c(50, -100),
            #ylim = rev(range(tbb3$best_change_SLD)), 
            xlim = c(0, 80))
cor.test(x=tbb3$best_change_SLD, y=tbb3$HRDbest, method="pearson")
lm(tbb3$best_change_SLD ~ tbb3$HRDbest)
abline(35.28, -1.23)
abline(v=42, lty=2, lwd=2, col="grey")
abline(v=33, lty=2, lwd=2, col="grey")
text(y=tbb3$best_change_SLD, x=tbb3$HRDbest, labels=tbb3$genetics, pos=tbb3$pos3)
legend(x=65, y=0, legend=levels(as.factor(tbb$type)), 
       fill=brewer.pal(n = 6, name = "Set2"))
text(x=8, y=-100, labels=c("Pearson's r=0.64"))
text(x=4, y=-94, labels=c("p=0.008"))


# correlation if remove gPALB2s
othersHRD <- c(NA,40,NA,5,24,NA,36,26,24,15,9,12,26,25)
othersSLD <- c(0,18.5,6.9,0,-35.3,-1.6,0,0,29.1,46.2,-11.7,20.5,15,25)
cor.test(othersHRD, othersSLD, method="pearson")

# calculate confidence intervals for ORRs
BiocManager::install("Hmisc")
library(Hmisc)
#breast ORR
binconf(x = 4, n = 13, method = "exact")

# breast CBR
binconf(x = 7, n=13, method = "exact")

# Non-breast ORR
binconf(x = 0, n = 7, method = "exact")

# Non-breast CBR 
binconf(x = 2, n=7, method = "exact")

# combined ORR
binconf(x = 4, n=20, method = "exact")

# combined CBR
binconf(x = 9, n = 20, method = "exact")


