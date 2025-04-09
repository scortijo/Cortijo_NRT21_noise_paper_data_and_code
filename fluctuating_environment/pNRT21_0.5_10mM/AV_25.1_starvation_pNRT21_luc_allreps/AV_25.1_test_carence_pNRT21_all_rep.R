
setwd("/Users/vettoral/Documents/luciferase_imaging_horizontal/AV_25.1_pNRT21_nitrate_starvation/AV_25.1_starvation_pNRT21_luc_allreps/")
getwd()

library(tidyverse)
library(scales)


######## all rep #########"

luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=DAY.x, y=-100 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=DAY.x, y=8000 ),position=position_dodge(0.8)) +
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


ggsave("Luciferase_signal_starvation_AV_25.1_all_rep.jpg")




filter(luc_stats) %>%
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1","T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=CV.x, col=KNO3.x))+
  geom_line(aes(group=KNO3.x), linetype="dashed", size=1)+
  geom_point(size=3)+
  scale_color_manual(values = c("#00AFBB", "#CB2027", "#66CC00", "#9900CC"), 
                     name="genotype", 
                     labels=c(expression(italic("0.5mM->0.5mM")),
                              expression(italic("0.5mM->10mM")),
                              expression(italic("10mM->10mM")),
                              expression(italic("10mM->0.5mM"))))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("Time after transfert (T0)") +
  ylab("CV") + 
  ggtitle(expression(""))


ggsave("AV_25.1_nitrate_starvation_CV_all_rep.jpg")


############### only 0.5mM->10mM pNRT21 ########

luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="10mM->0.5mM", KNO3!="10mM->10mM") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="10mM->0.5mM", KNO3!="10mM->10mM") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=DAY.x, y=-100 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=DAY.x, y=8000 ),position=position_dodge(0.8)) +
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


ggsave("Luc_signal_starvation_25.1_pNRT21_luc_0.5mM_10mM.jpg")


filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", 
       KNO3!="10mM->0.5mM", KNO3!="10mM->10mM",
       PLATE!="plate3", PLATE!="plate4",) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 0.5mM->10mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_1and2.jpg")

filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", 
       KNO3!="10mM->0.5mM", KNO3!="10mM->10mM",
       PLATE!="plate1", PLATE!="plate2",) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 0.5mM->10mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_3and4.jpg")




filter(luc_stats) %>%
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1","T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=CV.x, col=KNO3.x))+
  geom_line(aes(group=KNO3.x), linetype="dashed", size=1)+
  geom_point(aes(fill = KNO3.x, color = KNO3.x), size=3)+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("0.5mM->10mM"))))+
  scale_color_manual(values = c("#00AFBB"))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("Time after transfert (T0)") +
  ylab("CV") + 
  ggtitle(expression(""))



ggsave("starvation_CV_0and1mM_pNRT21.jpg")





############### only 0.5mM->0.5mM pNRT21 ########

luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", KNO3!="10mM->10mM") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", KNO3!="10mM->10mM") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=DAY.x, y=-100 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=DAY.x, y=8000 ),position=position_dodge(0.8)) +
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

ggsave("Luc_signal_starvation_25.1_pNRT21_luc_0.5mM_0.5mM.jpg")





filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", 
       KNO3!="10mM->10mM", PLATE!="plate7", PLATE!="plate8") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 0.5mM->0.5mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_5and6.jpg")


filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", 
       KNO3!="10mM->10mM", PLATE!="plate5", PLATE!="plate6",) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 0.5mM->0.5mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_7and8.jpg")



### selection seedling for transgeneration ###


filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", 
       KNO3!="10mM->10mM", seedling!="seedling35", seedling!="seedling45",
       seedling!="seedling46", seedling!="seedling49", seedling!="seedling54",
       seedling!="seedling60", DAY!="T0.1") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 50, position = position_dodge(0))+
   theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling on 0.5mM->0.5mM"))+
  facet_grid(.~KNO3, scales = "free") +
  guides(col="none")


ggsave("seedling_signal_selection_for_transgeneration.jpg")



filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="10mM->0.5mM", 
       KNO3!="10mM->10mM", seedling!="seedling35", seedling!="seedling45",
       seedling!="seedling46", seedling!="seedling49", seedling!="seedling54",
       seedling!="seedling60", DAY!="T0.1", DAY!="T0", DAY!="T1", DAY!="T2",
       DAY!="T3", DAY!="T6", DAY!="T5", DAY!="T7") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 50, position = position_dodge(0))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling on 0.5mM->0.5mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_selection_for_transgeneration_T4_only.jpg")



####

luc <- read.table("AV_25_T4_selection_trans.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->0.5mM")))

luc_stats <- filter(luc, seedling!="background") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$signal <- fct_relevel(luc$signal, c("low", "medium","high"))

levels(luc$signal)

filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=signal))+
   geom_dotplot(aes(fill = signal, color = signal), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 60, position = position_dodge(0))+
  scale_fill_manual(values = c("#00AFBB","#FFCC00", "#CB2027"), name="signal", 
                    labels=c(expression(italic("low")),
                            expression(italic("medium")),
                              expression(italic("high"))))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling on 0.5mM->0.5mM"))+
  facet_grid(.~DAY, scales = "free") +
  scale_color_manual(values = c("#00AFBB","#FFCC00", "#CB2027"))+
  guides(col="none")

ggsave("signal_selection_T4_3_level.jpg")


###point ###



filter(luc, seedling!="background") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T4")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 50, position = position_dodge(0))+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling on 0.5mM->0.5mM"))+
  facet_grid(.~KNO3, scales = "free") +
  guides(col="none")



############### only 10mM->0.5mM pNRT21 ########

luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="0.5mM->10mM", KNO3!="10mM->10mM") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="0.5mM->10mM", KNO3!="10mM->10mM") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=DAY.x, y=-100 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=DAY.x, y=8000 ),position=position_dodge(0.8)) +
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

ggsave("Luc_signal_starvation_25.1_pNRT21_luc_10mM_0.5mM.jpg")





filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="0.5mM->10mM",
       KNO3!="10mM->10mM", PLATE!="plate11", PLATE!="plate12", PLATE!="plate9",
       DAY!="T0.1") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 10mM->0.5mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_10_10vers0.5.jpg")



filter(luc, seedling!="background", KNO3!="0.5mM->0.5mM", KNO3!="0.5mM->10mM",
       KNO3!="10mM->10mM", PLATE!="plate9", PLATE!="plate10",) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 150, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 10mM->0.5mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_11and12.jpg")



############### only 10mM->10mM pNRT21 ########

luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM", "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="0.5mM->0.5mM", KNO3!="10mM->0.5mM") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="0.5mM->0.5mM", KNO3!="10mM->0.5mM") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=RawIntDen_normBackg_normLen, col=genotype.x))+
  geom_boxplot(aes(color = genotype.x), width = 0.5, position = position_dodge(0.8),
               outlier.shape = NA) +
  geom_dotplot(aes(fill = genotype.x, color = genotype.x), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 100, position = position_dodge(0.8))+
  facet_grid(~KNO3.x, scales= "free")+
  scale_fill_manual(values = c("#00AFBB"), name="genotype", 
                    labels=c(expression(italic("pNRT2.1:LUC"))))+
  scale_color_manual(values = c("#00AFBB"))+
  geom_text(aes(label=n , x=DAY.x, y=-100 ),position=position_dodge(0.8)) +
  geom_text(aes(label=CV , x=DAY.x, y=8000 ),position=position_dodge(0.8)) +
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

ggsave("Luc_signal_starvation_25.1_pNRT21_luc_10mM_10mM.jpg")





filter(luc, seedling!="background", KNO3!="0.5mM->10mM", KNO3!="0.5mM->0.5mM", 
       KNO3!="10mM->0.5mM", PLATE!="plate15", PLATE!="plate16",) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 120, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 10mM->10mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_13and14.jpg")


filter(luc, seedling!="background", DAY!="T0.1", KNO3!="0.5mM->10mM", KNO3!="0.5mM->0.5mM", 
       KNO3!="10mM->0.5mM", PLATE!="plate13", PLATE!="plate14", PLATE!="plate16") %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=seedling))+
  geom_line(aes(group=seedling), linetype="dashed", 
            size=1)+
  geom_dotplot(aes(fill = seedling, color = seedling), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 120, position = position_dodge(0))+
  facet_wrap(~KNO3, scales= "free")+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 10mM->10mM"))+
  facet_grid(.~PLATE, scales = "free") +
  guides(col="none")

ggsave("seedling_signal_PLATE_15and16.jpg")


#### CV ####


luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM",
                                      "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", DAY!='T0.1') %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)

filter(luc, seedling!="background", DAY!='T0.1') %>%
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", 
                                    "T5", "T6", "T7")))  %>%
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=CV, col=KNO3.x))+
  geom_line(aes(group=KNO3.x), linetype="dashed", size=1)+
  geom_point(aes(fill = KNO3.x, color = KNO3.x), size=3)+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("Time") +
  ylab("CV") + 
  ggtitle(expression(""))





######## mean signal ########


luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM",
                                      "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", DAY!='T0.1') %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(signalmean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)

filter(luc, seedling!="background", DAY!="T0.1") %>%
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", 
                                    "T5", "T6", "T7")))  %>%
  full_join(luc_stats, by="Sample") %>%
  ggplot( aes(x=DAY.x, y=signalmean, col=KNO3.x))+
  geom_line(aes(group=KNO3.x), linetype="dashed", size=1)+
  geom_point(aes(fill = KNO3.x, color = KNO3.x), size=3)+
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab("Time") +
  ylab("mean_signal_luc") + 
  ggtitle(expression(""))





###### seedling filter  ###



luc <- read.table("AV_25.1_test_carence_pNRT2_all_reps_V2.txt",
                  header=TRUE, sep='\t') %>% 
  mutate(KNO3=factor(KNO3, levels = c("0.5mM->10mM","0.5mM->0.5mM", "10mM->0.5mM",
                                      "10mM->10mM")))

luc_stats <- filter(luc, seedling!="background", KNO3!="0.5mM->10mM",
                    KNO3!="10mM->10mM", KNO3!="10mM->0.5mM") %>% 
  group_by(genotype, DAY, KNO3)%>%
  summarise(n=n(),
            CV=sd(RawIntDen_normBackg_normLen, 
                  na.rm = TRUE)/ mean(RawIntDen_normBackg_normLen, na.rm = TRUE),
            mean=mean(RawIntDen_normBackg_normLen, na.rm = TRUE)) %>% 
  unite(Sample, genotype, KNO3, DAY, remove=FALSE) %>% 
  mutate(n=paste("",n),
         CV=paste("", round(CV,2)))

luc$genotype <- fct_relevel(luc$genotype, c("pNRT2.1:LUC"))

levels(luc$genotype)

filter(luc, seedling!="background", DAY!="T0.1", KNO3!="0.5mM->10mM", 
       KNO3!="0.5mM->0.5mM", 
        PLATE!="plate13", PLATE!="plate14", PLATE!="plate16",
       PLATE!="plate10", PLATE!="plate11", PLATE!="plate12",
       DAY!="T7") %>%
  unite(Sample, genotype, KNO3, DAY, remove=FALSE)  %>%
  mutate(DAY=factor(DAY, levels = c("T0", "T0.1", "T1", "T2", "T3", "T4", 
                                    "T5", "T6", "T7")))  %>%
  ggplot( aes(x=DAY, y=RawIntDen_normBackg_normLen, col=KNO3))+
  geom_line(aes(group=seedling), linetype="dashed", size=1)+
  geom_dotplot(aes(fill = KNO3, color = KNO3), binaxis='y', stackdir='center', 
               dotsize = 1, binwidth = 120, position = position_dodge(0))+
    theme(legend.text=element_text(size=16),
        legend.title = element_text(size=18, face="bold"), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18),
        plot.title = element_text(size=20,face="bold"),
        legend.key = element_rect(fill = "white"),
        legend.title.align=0.5, 
        panel.background = element_rect(fill = "white", colour = "grey50"))+
  xlab(expression(KNO[3]~concentration)) +
  ylab("Normalised luciferase signal") + 
  ggtitle(expression("pNRT2.1 luciferase signal for seedling by plate on 10mM->10mM"))+
  facet_grid(.~genotype, scales = "free") +
  guides(col="none")

ggsave("Luc_signal_starvation_25.1_pNRT21_luc_.jpg")

