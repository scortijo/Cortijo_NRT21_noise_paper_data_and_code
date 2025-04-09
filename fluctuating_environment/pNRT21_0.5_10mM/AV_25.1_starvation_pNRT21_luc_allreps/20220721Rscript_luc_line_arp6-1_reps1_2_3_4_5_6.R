
setwd("/Users/vettoral/Documents/luciferase_imaging_horizontal/2022/2022_arp6-1_pNRT21_LUC_compilation_rep/arp6_vs_pNRT21_compilation/")

library(tidyverse)
library(scales)
library(ggpubr)
library(RColorBrewer)

#

luc_reps <- read.table("NRT21_apr6-1_luc_analysis_ALLreps.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("1mM","10mM")),
         genotype=factor(genotype, levels = c("pNRT2.1:LUC","arp6-1/pNRT2.1:LUC")))


luc_reps_stats <- filter(luc_reps, seedling!="background") %>% 
  group_by(genotype, KNO3, replicate)%>%
  summarise(average.NormRawIntDen_normLen = mean(RawIntDen_normBackg_normLen, 
                                                 na.rm = TRUE),
            n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/average.NormRawIntDen_normLen) %>% 
  unite(Sample, genotype, KNO3, replicate, remove=FALSE) %>% 
  mutate(n=paste("n =",n))



luc_reps12345678 <- read.table("NRT21_apr6-1_luc_analysis_ALLreps.txt",
                       header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("1mM","10mM")),
         genotype=factor(genotype, levels = c("pNRT2.1:LUC","arp6-1/pNRT2.1:LUC"))) %>% 
  filter(replicate%in%c("rep1","rep2","rep3","rep4","rep6", "rep7", "rep8"))


luc_reps12345678_stats <- filter(luc_reps12345678, seedling!="background") %>% 
  group_by(genotype, KNO3, replicate)%>%
  summarise(average.NormRawIntDen_normLen = mean(RawIntDen_normBackg_normLen, 
                                                 na.rm = TRUE),
            n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/average.NormRawIntDen_normLen) %>% 
  unite(Sample, genotype, KNO3, replicate, remove=FALSE) %>% 
  mutate(n=paste("n =",n))


write_tsv(luc_reps12345678_stats, file="STATS_NRT21_arp6-1_luc_analysis_reps12345678.txt")

luc_stats_reps <- read.table("STATS_NRT21_arp6-1_luc_analysis_reps12345678.txt",
                       header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("1mM","10mM")),
         genotype=factor(genotype, levels = c("pNRT2.1:LUC","arp6-1/pNRT2.1:LUC")))

ggplot(luc_stats_reps, aes(x=genotype, y=CV)) +
  geom_point(aes( color=replicate), size=8) +
  ggtitle("Comparison CV arp6-1/pNRT2.1:LUC reps12345678 ")+
  geom_line(aes(group=replicate), col="grey", linetype="dotted") +
  facet_grid(~KNO3)+
  scale_color_brewer(palette="Dark2", labels=c("rep1"="rep1 AV_01","rep2"="rep2 AV_08",
                                               "rep3"="rep3 SC_13","rep4"="rep4 SC_15",
                                               "rep6"="rep6 AV_16",
                                               "rep7"="rep7 AV_18", "rep8"="rep8 SC_21")) +
  ylim(0,0.6) +
  ylab(expression(paste("CV for ", italic("pNRT2.1:LUC"), " signal"))) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20), 
        axis.title=element_text(size=22),
        plot.title = element_text(size=20,face="bold", hjust=0,5),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5,
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text.x = element_text(size = 20)) +
  stat_compare_means(label="p.format", method = "wilcox.test",  paired = TRUE,
                     label.x = 1.4, 
                     label.y = 0.55) +
  scale_x_discrete(labels=c("pNRT2.1:LUC" = "WT", "arp6-1/pNRT2.1:LUC" = expression(italic("arp6-1"))))



ggsave("Comparison_CV_reps12345678_arp6-1_nrt21.jpg")

###########

filter(luc_stats_reps, KNO3=="1mM") %>%
  ggplot( aes(x=genotype, y=CV))+
  geom_point(aes( color=replicate), size=8) +
  ggtitle("Comparison CV arp6-1/pNRT2.1:luc")+
  geom_line(aes(group=replicate), col="grey", linetype="dotted") +
  facet_grid(~KNO3)+
  scale_color_brewer(palette="Dark2", labels=c("rep1"="rep1","rep2"="rep2",
                                               "rep3"="rep3","rep4"="rep4",
                                               "rep6"="rep6",
                                               "rep7"="rep7", "rep8"="rep8")) +
  ylim(0,0.6) +
  ylab(expression(paste("CV for ", italic("pNRT2.1:LUC"), " signal"))) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=22, face="bold"),
        axis.text=element_text(size=20), 
        axis.title=element_text(size=22),
        plot.title = element_text(size=20,face="bold", hjust=0,5),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5,
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text.x = element_text(size = 20)) +
  stat_compare_means( method = "wilcox.test",  paired = TRUE,
                     label.x = 1.3, 
                     label.y = 0.55) +
  scale_x_discrete(labels=c("pNRT2.1:LUC" = "WT", "arp6-1/pNRT2.1:LUC" = expression(italic("arp6-1"))))

ggsave("Comparison_CV_allreps_arp6-1_1mM.jpg")

########

ggplot(luc_reps1234567_stats, aes(x=genotype, y=CV, group=replicate, color=replicate)) +
  geom_point(size=6) +
  facet_grid(~KNO3)+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5,
        panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.text.x = element_text(size = 16)) +
  stat_compare_means(label="p.format", method = "wilcox.test", paired = TRUE,
                     label.x = 0.2,
                     label.y = 0.8, comparisons = list(c("arp6-1/pNRT2.1:LUC", "pNRT2.1:LUC")))






###################### all length root size




library(ggpubr)

luc <- read.table("NRT21_apr6-1_luc_analysis_ALLreps.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("1mM","10mM")))


luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC", "arp6-1/pNRT2.1:LUC" ))

levels(luc$genotype)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=KNO3.x , y=length_mm, col=genotype.x))+
  geom_violin(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 0.5, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB", "#CB2027"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC")), 
                             expression(italic("arp6/pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB", "#CB2027"))+
  geom_text(aes(label=n , x=KNO3.x, y=20 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  stat_compare_means(label="p.format", method = "wilcox.test",
                     label.x = 20,
                     label.y = 60) +
  xlab(expression(KNO[3]~concentration)) +
  ylab("length roots") +
  guides(col="none")


ggsave("length_arp6_replicate.jpg")



####### all reps arp6 signal #####


luc <- read.table("NRT21_apr6-1_luc_analysis_ALLreps.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("1mM","10mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3, replicate)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            QCD=(as.numeric(quantile(RawIntDen_normBackg_normLen, 0.75, na.rm = TRUE))-
                   as.numeric(quantile(RawIntDen_normBackg_normLen, 0.25, na.rm = TRUE)))/
              (as.numeric(quantile(RawIntDen_normBackg_normLen, 0.75, na.rm = TRUE))+
                 as.numeric(quantile(RawIntDen_normBackg_normLen, 0.25, na.rm = TRUE)))) %>% 
  unite(Sample, genotype, KNO3, replicate, remove=FALSE) %>% 
  mutate(n=paste("n",n),
         CV=paste("", round(CV,2)),
         QCD=paste("", round(QCD,2)))


luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC", "arp6/pNRT2.1:LUC" ))

levels(luc$genotype)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, replicate, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=replicate.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB", "#CB2027"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC")), 
                             expression(italic("arp6/pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB", "#CB2027"))+
  geom_text(aes(label=n , x=replicate.x, y=-1000 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=replicate.x, y=18000 ),position=position_dodge(0.8)) +
  geom_text(aes(label=QCD , x=replicate.x, y=17000 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Luciferase signal \n (Raw integrated density normalisez by leng root)") +
  guides(col="none")


ggsave("RawIntDen_normBackg_normLen_arp_all_replicate.jpg")






