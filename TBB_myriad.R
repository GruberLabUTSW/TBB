# set directory
setwd("~/Documents/Clinical_topics/TBB/TBB_data")

# load data
myr <- read.csv("20200527_TTBB_MyriadMatrix.csv", header=T)
myr

library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

# load data matrices
myra <- read.csv("20200527_TTBB_MyriadMatrix_anno.csv", header=T, 
                 stringsAsFactors=F, check.names=F,
                 row.names=1)
gline <- read.csv("20200606_TBB_entryGermlineMutations.csv", header=T,
                  stringsAsFactors=F, check.names=F, 
                  row.names=1)
somatic <- read.csv("20200606_TBB_entrySomaticMutations.csv", header=T,
                    stringsAsFactors=F, check.names=F,
                    row.names=1)
gline2 <- read.csv("20200606_TBB_entryGermlineMutations2.csv", header=T,
                  stringsAsFactors=F, check.names=F, 
                  row.names=1, colClasses="character")


myra_mat <- as.matrix(myra)
myra_mat
col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
?grid.text
?sprintf

# set color scheme
library(RColorBrewer)
palette(brewer.pal(n = 6, name = "Set2")) 
display.brewer.pal(n = 6, name = "Set2")
display.brewer.pal(n=8, name="Set2")
library(circlize)
col_fun = circlize::colorRamp2(c(0,1), c("beige", "cornflowerblue"))
column_ha = HeatmapAnnotation(tumor = colnames(myra_mat), 
                              annotation_name_side = "left",
                              col=list(tumor=c("breast"=palette()[1], 
                                               "colon"=palette()[2],
                                               "pancreas"=palette()[3],
                                               "parotid"=palette()[4],
                                               "testicular"=palette()[5],
                                               "uterine"=palette()[6])))
row_ha = rowAnnotation(n = anno_barplot(myra_mat), gp = gpar(col = "grey"))
row_ha2 = rowAnnotation(n=anno_barplot(rowSums(myra_mat)))

h=Heatmap(myra_mat, column_title="TBB cohort B tumor mutations", show_column_names = F,
        cluster_rows=T, show_column_dend=T, col=col_fun,
        border=T, 
        rect_gp=gpar(col = "white", lwd = 2),
        row_dend_side = "left", column_dend_side = "top", 
        column_dend_height = unit(1, "cm"), 
        row_dend_reorder = T, 
        top_annotation = column_ha, 
        right_annotation = row_ha2, 
        heatmap_legend_param = list(
                title = "mutation", at = 0:1, 
                labels = c("absent", "present"), border="black", 
                legend_gp = gpar(fill = 0:2),
                direction="horizontal", legend_width = unit(1.5, "cm"), 
                labels_rot = 60),
        cell_fun = function(j, i, x, y, width, height, fill) {
                if(gline2[i,j] == "g")
                        grid.text("g", x, y, 
                           gp = gpar(fontsize = 10, fontface="bold"))
                if(somatic[i,j] == 1)
                        grid.text("s", x, y,
                           gp = gpar(fontsize = 10, fontface="bold"))
                }
        )

draw(h, heatmap_legend_side = "right", annotation_legend_side = "right")
