setwd("~/Documents/Clinical_topics/TBB/TBB_data/myriad")
df.fig2 <- read.table('Fig2.tsv', header = T, sep = '\t')

require(ggplot2)
pd <- position_dodge(0.15)
pal.loc <- c(breast="#66C2A5",pancreas="#8DA0CB", uterine="#FFD92F", colon="#FC8D62",  parotid="#E78AC3")  # testicular="#A6D854", 

level_order <- c('Primary', 'Metastasis')
ggplot(df.fig2, aes(x=factor(location, level=level_order),y=HRD.score,group=Patient,col=type)) + 
  xlab('Location') + ylab('HRD score') + 
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(), panel.background=element_blank())+
  theme(axis.text = element_text(size=16), axis.title = element_text(size=16))+
  theme(legend.title = element_blank(), legend.text = element_text(size=14))+
  geom_hline(yintercept=c(33,42), linetype='dashed', color='grey')+
  geom_line(position=pd, size=2, alpha=0.6)+
  scale_color_manual(values = pal.loc)+
  geom_point(position=pd, size=4, alpha=0.6, shape=19)
