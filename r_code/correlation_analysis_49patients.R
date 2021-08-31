## 49 patients
library(readxl)
library(dplyr)
library(ggsci)
library(ggplot2)
library(psych)
library(cowplot)
library(survival)
library(survminer)
library(MASS)
library(ltm)

dat_all <- read_excel("data/210605_FINAL VERSION_patient file_49patients.xlsx",
                      sheet = "PATIENT FILE_DISS_FCZ ", skip = 1)

dat_all <- dat_all %>%
  mutate(
    Protocol = factor(Protocol),
    Sex = factor(Sex),
    Stage = factor(Stage),
    Progression_y_n = as.numeric(Progression_y_n),
    log_GTV_absolut = log(GTV_absolut),
    log_GTV_rel_weight = log(GTV_rel_weight),
    log_GTV_rel_nasal_cavity = log(GTV_rel_nasal_cavity),
    log_GTV_rel_BSA = log(GTV_rel_BSA),
    stage4_status = factor(if_else(Stage == 4, 1, 0)),
    stage34_status = factor(if_else(Stage == 4|Stage ==3, 1, 0))
  )



## Correlation analysis
covariate_list <- c("Weight",
                    "body_surface_area", 
                    "Nasal_Cavity", 
                    "Stage", 
                    "GTV_absolut",
                    "GTV_rel_weight", 
                    "GTV_rel_nasal_cavity", 
                    "GTV_rel_BSA")

# correlation between body sizes
cor.test(dat_all$Weight,dat_all$Nasal_Cavity, method="pearson", exact = F) 
cor.test(dat_all$body_surface_area,dat_all$Nasal_Cavity, method="pearson", exact = F) 


# GTVs vs body sizes
all_GTV <- list(dat_all$GTV_absolut, dat_all$GTV_rel_weight,  dat_all$GTV_rel_BSA, dat_all$GTV_rel_nasal_cavity)
lapply(all_GTV, FUN = cor.test, y = dat_all$Weight)
lapply(all_GTV, FUN = cor.test, y = dat_all$body_surface_area)
lapply(all_GTV, FUN = cor.test, y = dat_all$Nasal_Cavity)

#plots
#=======
give.n <- function(x){
  return(data.frame(y =0.4*median(x)+max(x) , label = paste("n =",length(x))))
}

plot_stages <- function(data, stage_type, GTVtype){
  arg <- match.call()
  ggplot(aes(x= eval(arg$stage_type), 
             y = eval(arg$GTVtype), 
             fill = eval(arg$stage_type)),
         data = data) +
    geom_violin(trim = TRUE) +
    scale_fill_brewer(palette="Greens") +
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() +  
    stat_summary(fun.data = give.n, geom = "text")+
    xlab("")+ theme(legend.position = "none")
}

# with stage1 vs stage 2 vs stage 3 vs stage 4
label_xaxis <- scale_x_discrete(labels=c("Stage 1",
                                           "Stage 2",
                                           "Stage 3", 
                                           "Stage 4"))  
g1 <- plot_stages(dat_all, Stage, GTV_absolut) +
  label_xaxis +
  ylab("GTV")
g2 <- plot_stages(dat_all, Stage, GTV_rel_weight)  +
  label_xaxis +
  ylab("GTV rel weight")
g3 <- plot_stages(dat_all, Stage, GTV_rel_BSA)  +
  label_xaxis +
  ylab("GTV rel BSA")
g4 <- plot_stages(dat_all, Stage, GTV_rel_nasal_cavity)  +
  label_xaxis +
  ylab("GTV rel nasal cavity")
plot_grid(g1,g2,g3,g4)


# with stage1-2 vs stage 3-4
label_xaxis2 <- scale_x_discrete(labels=c("Stage 1-2",
                                         "Stage 3-4"))  
h1 <- plot_stages(dat_all, stage34_status, GTV_absolut) +
  label_xaxis2 +
  ylab("GTV")
h2 <- plot_stages(dat_all, stage34_status, GTV_rel_weight)  +
  label_xaxis2 +
  ylab("GTV rel weight")
h3 <- plot_stages(dat_all, stage34_status, GTV_rel_BSA)  +
  label_xaxis2 +
  ylab("GTV rel BSA")
h4 <- plot_stages(dat_all, stage34_status, GTV_rel_nasal_cavity)  +
  label_xaxis2 +
  ylab("GTV rel nasal cavity")
plot_grid(h1,h2,h3,h4)


# with stage1-2 vs stage 3-4
label_xaxis3 <- scale_x_discrete(labels=c("Stage 1-3",
                                          "Stage 4"))  
j1 <- plot_stages(dat_all, stage4_status, GTV_absolut) +
  label_xaxis3 +
  ylab("GTV")
j2 <- plot_stages(dat_all, stage4_status, GTV_rel_weight)  +
  label_xaxis3 +
  ylab("GTV rel weight")
j3 <- plot_stages(dat_all, stage4_status, GTV_rel_BSA)  +
  label_xaxis3 +
  ylab("GTV rel BSA")
j4 <- plot_stages(dat_all, stage4_status, GTV_rel_nasal_cavity)  +
  label_xaxis3 +
  ylab("GTV rel nasal cavity")
plot_grid(j1,j2,j3,j4)



# Statistical Tests
#==================
# with stage1 vs stage 2 vs stage 3 vs stage 4------
kruskal.test(GTV_absolut ~ Stage, data = dat_all)
kruskal.test(GTV_rel_weight ~ Stage, data = dat_all)
kruskal.test(GTV_rel_BSA ~ Stage, data = dat_all)
kruskal.test(GTV_rel_nasal_cavity ~ Stage, data = dat_all)

pairwise.wilcox.test(dat_all$GTV_rel_nasal_cavity, dat_all$Stage,
                     p.adjust.method = "BH")
#between stage 1 and stage 3, 
#stage 1 and stage 4, 
#stage 2 and stage 3, 
#stage 2 and stage 4

# compare stage1-2 vs stage3-4------
wilcox.test(GTV_absolut ~ stage34_status, data = dat_all)
wilcox.test(GTV_rel_weight ~ stage34_status, data = dat_all) #0.02777
wilcox.test(GTV_rel_BSA ~ stage34_status, data = dat_all) #0.02476
wilcox.test(GTV_rel_nasal_cavity ~ stage34_status, data = dat_all) #0.000236



# compare stage1-3 vs stage4------
wilcox.test(GTV_absolut ~ stage4_status, data = dat_all)
wilcox.test(GTV_rel_weight ~ stage4_status, data = dat_all)
wilcox.test(GTV_rel_BSA ~ stage4_status, data = dat_all) 
wilcox.test(GTV_rel_nasal_cavity ~ stage4_status, data = dat_all) #0.04554



