## 34 patients
library(readxl)
library(dplyr)
library(ggsci)
library(ggplot2)
library(psych)
library(cowplot)
library(survival)
library(survminer)
library(MASS)
library(pROC)
library(tableone)
library(asbio)

dat <- read_excel("data/210605_FINAL VERSION_outcome file_34patients.xlsx", skip = 1)
tumor_codes <- read_excel("data/210605_FINAL VERSION_patient file_49patients.xlsx",sheet = "PATIENT FILE_DISS_FCZ ", skip = 1) %>%
  data.frame() %>%
  dplyr::select(KG_Nr, TumorCode_simplyfied)


dat <- dat %>%
  mutate(
    Protocol = factor(Protocol),
    Sex = factor(Sex),
    is_stage4 = if_else(Stage =="4", 1,0),
    is_stage34 = if_else(Stage =="3"|Stage =="4", 1,0),
    Stage = factor(Stage),
    Progression_y_n = as.numeric(Progression_y_n),
    log_GTV_absolut = log(GTV_absolut),
    log_GTV_rel_weight = log(GTV_rel_weight),
    log_GTV_rel_nasal_cavity = log(GTV_rel_nasal_cavity),
    log_GTV_rel_BSA = log(GTV_rel_BSA)
  ) %>%
  left_join(tumor_codes)

median_TTP <- dat %>%
  group_by(Protocol) %>%
  summarize(mean = mean(PFS),
            sd = sd(PFS),
            median = median(PFS)) %>%
  data.frame()

dat_regular <- dat %>% filter(Protocol==0)
dat_boost <- dat %>% filter(Protocol==1)

ci.median(dat_regular$PFS, conf = 0.95)
ci.median(dat_boost$PFS, conf = 0.95)

kruskal.test(PFS ~ Protocol, data = dat)



# Kapler-mier curves
KM_protcol <- survfit(Surv(PFS,Progression_y_n) ~ Protocol , data = dat)

KMplots <- list()
KMplots[[1]] <- ggsurvplot(KM_protcol,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = pal_jco()(4),
           legend.labs = c("regular protocol", "boost protocol"))


KM_stage <- survfit(Surv(PFS,Progression_y_n) ~ is_stage4 , data = dat)

KMplots[[2]] <- ggsurvplot(KM_stage,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = pal_jco()(4),
           legend.labs = c("Stage 1-3", "Stage 4"))

arrange_ggsurvplots(KMplots, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.3)


KM_stage34 <- survfit(Surv(PFS,Progression_y_n) ~ is_stage34 , data = dat)

KMplot_stage34 <- ggsurvplot(KM_stage34,
                           pval = TRUE, conf.int = TRUE,
                           risk.table = TRUE, # Add risk table
                           risk.table.col = "strata", # Change risk table color by groups
                           linetype = "strata", # Change line type by groups
                           surv.median.line = "hv", # Specify median survival
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           palette = pal_jco()(4),
                           legend.labs = c("Stage 1-2", "Stage 3-4"))


KMplot_stage34


# Univariate analysis stratified by protocol

#GTV_absolut
#------------------
m_abs<- coxph(Surv(PFS,Progression_y_n) ~ GTV_absolut + strata(Protocol) , data = dat)
summary(m_abs)
cox.zph(m_abs)
#GTV_rel weight
#------------------
m_weight<- coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_weight + strata(Protocol) , data = dat)
summary(m_weight)
cox.zph(m_weight)
#GTV_rel_BSA
#------------------
m_bsa<- coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_BSA + strata(Protocol) , data = dat)
summary(m_bsa)
cox.zph(m_bsa)

#GTV_rel_nasal_cavity
#------------------
m_nasal<- coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_nasal_cavity + strata(Protocol) , data = dat)
summary(m_nasal)
cox.zph(m_nasal)



# Multivariate analysis stratified by protocol adjusted for stage 34 status
m_nasal_stage34<- coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_nasal_cavity + is_stage34+  strata(Protocol) , data = dat)
summary(m_nasal_stage34)

summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_absolut + is_stage34+  strata(Protocol) , data = dat))
summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_weight + is_stage34+  strata(Protocol) , data = dat))
summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_BSA + is_stage34+  strata(Protocol) , data = dat))

# Multivariate analysis stratified by protocol adjusted for stage 34 status
m_nasal_stage4<- coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_nasal_cavity + is_stage4+  strata(Protocol) , data = dat)
summary(m_nasal_stage4)



summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_absolut + is_stage4+  strata(Protocol) , data = dat))
summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_weight + is_stage4+  strata(Protocol) , data = dat))
summary(coxph(Surv(PFS,Progression_y_n) ~ GTV_rel_BSA + is_stage4+  strata(Protocol) , data = dat))

