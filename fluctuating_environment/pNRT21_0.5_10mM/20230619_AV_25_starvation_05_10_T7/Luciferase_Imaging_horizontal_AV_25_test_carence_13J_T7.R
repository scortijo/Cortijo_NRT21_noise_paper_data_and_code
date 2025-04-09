####script propre#####

setwd("/Users/alex6/Desktop/AV_25.1_pNRT21_carence/20230619_AV_25_starvation_05_10_T7/")
getwd()

library(tidyverse)
library(scales)

#### test effet box ####

testbox <- read.table("AV_25_test_carence_13J_T7.txt",
                     header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

bind_rows(testbox) %>% filter(seedling!="background") %>% 
  ggplot( aes(x=genotype, y=RawIntDen_normBackg_normLen, col=PLATE))+
  geom_boxplot( width = 0.5, position = position_dodge(0.8),
                outlier.shape = NA) + 
  geom_dotplot(aes(fill = PLATE), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3, scales= "free")+
  theme_bw()+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"))+
  xlab("genotype") +
  ylab("RawIntDen_normBackg_normLen")

ggsave("compilation_DIFFERENCE_PLATE_AV_25.1_13j_T7.jpg")

######### Luciferase signal Raw integrated density #########


luc <- read.table("AV_25_test_carence_13J_T7.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=KNO3.x , y=RawIntDen.norm_bckg, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 5000, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=KNO3.x, y=-10000 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=KNO3.x, y=250000 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Luciferase signal \n (Raw integrated density)") +
  guides(col="none")

ggsave("Luciferase_signal_AV_25.1_13J_T7.jpg")


################ Mean norm data #################


#Do the same as above but with data also normalised by the root length#

All_meanForNorm <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(average.NormRawIntDen = mean(RawIntDen.norm_bckg, na.rm = TRUE),
            n=n(),
            CV=sd(RawIntDen.norm_bckg, 
                  na.rm = TRUE)/average.NormRawIntDen) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)))

All_MeanNormalised <- filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(., All_meanForNorm, by="Sample") %>% 
  mutate(MeanNormalised_RawIntDenSignal=log2(RawIntDen.norm_bckg/average.NormRawIntDen))

ggplot(All_MeanNormalised, aes(x=KNO3.x , y=MeanNormalised_RawIntDenSignal, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 0.1, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=KNO3.x, y=-3 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=KNO3.x, y=1.5 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Mean normalised raw integrated density") +
  guides(col="none")

ggsave("Mean_normalized.jpg")


############# Mean norm data for length norm signal ##############


All_meanForNorm <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(average.NormRawIntDen_normLen = mean(RawIntDen_normBackg_normLen, 
                                                 na.rm = TRUE),
            n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/average.NormRawIntDen_normLen) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)))

All_MeanNormalised <- filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(., All_meanForNorm, by="Sample") %>% 
  mutate(MeanNormalised_RawIntDenSignal_normLen=
           log2(RawIntDen_normBackg_normLen/average.NormRawIntDen_normLen))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

ggplot(All_MeanNormalised, aes(x=KNO3.x , y=MeanNormalised_RawIntDenSignal_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 0.12, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=KNO3.x, y=-4 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=KNO3.x, y=2 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Mean normalised raw integrated density \n (length normalised)") +
  guides(col="none")

ggsave("Meannorm_data_length_norm_AV_25.1_13j_T7.jpg")


############## RawIntDen_normBackg_normLen ############


luc <- read.table("AV_25_test_carence_13J_T7.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            QCD=(as.numeric(quantile(RawIntDen_normBackg_normLen, 0.75, na.rm = TRUE))-
                   as.numeric(quantile(RawIntDen_normBackg_normLen, 0.25, na.rm = TRUE)))/
              (as.numeric(quantile(RawIntDen_normBackg_normLen, 0.75, na.rm = TRUE))+
                 as.numeric(quantile(RawIntDen_normBackg_normLen, 0.25, na.rm = TRUE)))) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)),
         QCD=paste("QCD =", round(QCD,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=KNO3.x , y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=KNO3.x, y=-900 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=KNO3.x, y=6000 ),position=position_dodge(0.8)) +
  geom_text(aes(label=QCD , x=KNO3.x, y=5000 ),position=position_dodge(0.8)) +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Luciferase signal \n (Raw integrated density)") +
  guides(col="none")

  guides(col="none")+stat_compare_means(label="p.format", method = "wilcox.test",
                                        label.x = 1.4,
                                        label.y = 8000)

ggsave("RawIntDen_normBackg_normLen_AV_25.1_13J_T7.jpg")


###########lenght_roots############

library(ggpubr)

luc <- read.table("AV_25_test_carence_13J_T7.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  mutate(n=paste("n =",n),
         CV=paste("CV =", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

my_comparisons <- list( c("0.5mM->10mM", "0.5mM->0.5mM"),
                        c("0.5mM->10mM", "10mM->0.5mM"),
                        c("0.5mM->10mM", "10mM->10mM"),
                        c("0.5mM->0.5mM", "10mM->0.5mM"),
                        c("0.5mM->0.5mM", "10mM->10mM"),
                        c("10mM->0.5mM", "10mM->10mM"))

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=KNO3.x , y=length_mm, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 0.7, position = position_dodge(0.8))+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=KNO3.x, y=35),position=position_dodge(0.8)) +
    theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  stat_compare_means(comparisons = my_comparisons, label="p.format", method = "wilcox.test")+
  xlab(expression(KNO[3]~concentration)) +
  ylab("length roots") +
  guides(col="none")


ggsave("length_AV_25.1_13J_T7.jpg")

