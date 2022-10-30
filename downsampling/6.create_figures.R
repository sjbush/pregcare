library(tidyverse)
library(ggpubr)
library(dplyr)
library(scales)
theme_set(theme_bw())

## FAM27 ##

df<-read.csv('C:/Users/sbush/Desktop/pregcare/FAM27-downsampled.csv',header=T)
df.sub<-subset(df,(df$group == '27-F-3' | df$group == '27-F-4' | df$group == '27-F-6' | df$group == '27-Co'))
df.sub$group <- recode_factor(df.sub$group, "27-Co" = "control", "27-F-3" = "F3 (blood)", "27-F-4" = "F4 (saliva)", "27-F-6" = "F6 (semen)")
df.sub$group <- factor(df.sub$group, levels=c('control','F3 (blood)','F4 (saliva)','F6 (semen)'))
a<-ggplot(df.sub,aes(x=as.factor(avg_depth_to_which_reads_downsampled),y=DEL.ALL,fill=group)) + geom_boxplot(outlier.alpha=0.05, na.rm = TRUE) + labs(x='Downsampled depth (x-fold)', y='% of reads calling variant', title='A', fill='Sample') + scale_y_continuous(limits = c(0,0.6),labels = number_format(accuracy = 0.01)) + scale_fill_manual(values=c('grey','blue','darkblue','purple')) + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1))

df<-read.csv('C:/Users/sbush/Desktop/pregcare/FAM27-downsampled-corrected.csv',header=T)
df.sub<-subset(df,(df$group == '27-F-3' | df$group == '27-F-4' | df$group == '27-F-6' | df$group == '27-Co'))
df.sub$group <- recode_factor(df.sub$group, "27-Co" = "control", "27-F-3" = "F3 (blood)", "27-F-4" = "F4 (saliva)", "27-F-6" = "F6 (semen)")
df.sub$group <- factor(df.sub$group, levels=c('control','F3 (blood)','F4 (saliva)','F6 (semen)'))
b<-ggplot(df.sub,aes(x=as.factor(avg_depth_to_which_reads_downsampled),y=(vaf_corrected*100),fill=group)) + geom_boxplot(outlier.alpha=0.05, na.rm = TRUE) + labs(x='Downsampled depth (x-fold)', y='Corrected VAF (%)', title='B', fill='Sample') + scale_y_continuous(limits = c(0,0.6),labels = number_format(accuracy = 0.01)) + scale_fill_manual(values=c('grey','blue','darkblue','purple')) + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1))

fam27<-ggarrange(a,b,ncol=2,nrow=1,common.legend=TRUE,legend='bottom')

## FAM34 ##

df<-read.csv('C:/Users/sbush/Desktop/pregcare/FAM34-downsampled.csv',header=T)
df.sub<-subset(df,(df$group == '34-F-1' | df$group == '34-F-2' | df$group == '34-Co'))
df.sub$group <- recode_factor(df.sub$group, "34-Co" = "control", "34-F-1" = "F1 (left buccal swab)", "34-F-2" = "F2 (right buccal swab)")
df.sub$group <- factor(df.sub$group, levels=c('control','F1 (left buccal swab)','F2 (right buccal swab)'))
a<-ggplot(df.sub,aes(x=as.factor(avg_depth_to_which_reads_downsampled),y=DEL.ALL,fill=group)) + geom_boxplot(outlier.alpha=0.05, na.rm = TRUE) + labs(x='Downsampled depth (x-fold)', y='% of reads calling variant', title='C', fill='Sample') + scale_y_continuous(limits = c(0,1),labels = number_format(accuracy = 0.01)) + scale_fill_manual(values=c('grey','lightblue','cadetblue')) + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1))

df<-read.csv('C:/Users/sbush/Desktop/pregcare/FAM34-downsampled-corrected.csv',header=T)
df.sub<-subset(df,(df$group == '34-F-1' | df$group == '34-F-2' | df$group == '34-Co'))
df.sub$group <- recode_factor(df.sub$group, "34-Co" = "control", "34-F-1" = "F1 (left buccal swab)", "34-F-2" = "F2 (right buccal swab)")
df.sub$group <- factor(df.sub$group, levels=c('control','F1 (left buccal swab)','F2 (right buccal swab)'))
b<-ggplot(df.sub,aes(x=as.factor(avg_depth_to_which_reads_downsampled),y=(vaf_corrected*100),fill=group)) + geom_boxplot(outlier.alpha=0.05, na.rm = TRUE) + labs(x='Downsampled depth (x-fold)', y='Corrected VAF (%)', title='D', fill='Sample') + scale_y_continuous(limits = c(0,1),labels = number_format(accuracy = 0.01)) + scale_fill_manual(values=c('grey','lightblue','cadetblue')) + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1))

fam34<-ggarrange(a,b,ncol=2,nrow=1,common.legend=TRUE,legend='bottom')

## COMBINED FIGURE ##

pdf('C:/Users/sbush/Desktop/pregcare/SupplementaryFigureS3.pdf',width=10,height=12) # PDFs are 7x7 inches by default
ggarrange(fam27,fam34,ncol=1,nrow=2,common.legend=FALSE)
dev.off()
