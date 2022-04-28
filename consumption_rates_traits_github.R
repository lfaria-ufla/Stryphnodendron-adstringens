###3 sheets are necessary to run the code
## consumption_defense.csv, consumption.csv; consumption_quantity.csv

rm(list=ls())
setwd("~/estatistica/R/juliana")
dado=read.delim("consumption_defense.csv", sep= ";", header = TRUE)
dados=read.delim(file= "consumption.csv", sep= ";", header = T)
names(dados)

library(tidyverse)
library(MuMIn)
library(glmmTMB)
library(data.table)

#combining data
rate = dados%>%
  select(taxa_pred_par, taxa_pred_herb, area, cut, fert)
quim = dado%>%
  select(total_totalfenols, total_hydrolizable, total_condensed)

data = as.data.frame(cbind(rate, quim))
str(data)
data$area = as.factor(data$area)
data$cut = as.factor(data$cut)
data$fert = as.factor(data$fert)

##exploiting data
plot(taxa_pred_herb ~ total_totalfenols, data = data)
plot(taxa_pred_herb ~ total_hydrolizable, data = data)
plot(taxa_pred_herb ~ total_condensed, data = data)

plot(taxa_pred_par ~ total_totalfenols, data = data)
plot(taxa_pred_par ~ total_hydrolizable, data = data)
plot(taxa_pred_par ~ total_condensed, data = data)


###multimodels GLMM - model inference
names(data)
m1 = glmer.nb(taxa_pred_herb ~ total_totalfenols + 
                total_hydrolizable + total_condensed+(1| area), 
              family = "binomial"(link = "logit"), data = data, 
              na.action = "na.fail")
summary(m1)

####selecting variables and collecting data####
my_sel=dredge(m1)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) 
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

### plotting herbivory rate on defense chemicals (not shown in the manuscript)
mhc=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="blue", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour = "blue")+theme_classic(base_size = 20)+ ylab("")+
  xlab("Parameters")
mhc = mhc+ scale_x_discrete(labels=c("total_totalfenols" = "Total Fenols", 
                                     "total_hydrolizable" = "Hydrolyzable tannins",
                                     "total_condensed" = "Condensed tannins"))
mhc


####now model for parasitoid rates
### multimodel inference (GLMM)
m2 = glmer.nb(taxa_pred_par ~ total_totalfenols + total_hydrolizable + total_condensed + (1| area), 
              family = "binomial"(link = "logit"), 
              data = data, na.action = "na.fail")
summary(m2)

####Selecting variables and collecting data ####
my_sel=dredge(m2)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) #2 modelos
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

##plotting parasitoid rate on defense chemicals (not shown in the manuscript)
mpc=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="firebrick", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour="firebrick")+theme_classic(base_size = 20)+ ylab("")+
  xlab("")
mpc = mpc+ scale_x_discrete(labels=c("total_totalfenols" = "Total Fenols", 
                                     "total_hydrolizable" = "Hydrolyzable tannins",
                                     "total_condensed" = "Condensed tannins"))
mpc



######################################
#####now lets fits quality traits ###
names(dados)
quality = dados%>%
  select(taxa_pred_par, taxa_pred_herb, area, cut, fert,media_carbs_totais, 
         media_prots_totais, gords_totais)
str(quality)

##exploiting data
plot(taxa_pred_herb ~ media_carbs_totais, data = quality)
plot(taxa_pred_herb ~ media_prots_totais, data = quality)
plot(taxa_pred_herb ~ gords_totais, data = quality)

plot(taxa_pred_par ~ media_carbs_totais, data = quality)
plot(taxa_pred_par ~ media_prots_totais, data = quality)
plot(taxa_pred_par ~ gords_totais, data = quality)


##fiting model herbivory
names(data)
m1 = glmer.nb(taxa_pred_herb ~ media_carbs_totais + 
                media_prots_totais + gords_totais+(1| area), 
              family = "binomial"(link = "logit"), data = quality, 
              na.action = "na.fail")


####Selecao variavel####
my_sel=dredge(m1)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) 
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

##plotting figure 6A
mhq=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="blue", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour = "blue")+theme_classic(base_size = 20)+ 
  ylab("Coefficient values [95%CI]")+
  xlab("Parameters")
mhq = mhq+ scale_x_discrete(labels=c("media_carbs_totais" = "Carbohydrate", 
                                     "media_prots_totais" = "Protein",
                                     "gords_totais" = "Lipid"))
mhq



##fiting model parasitoid
m2 = glmer.nb(taxa_pred_par ~ media_carbs_totais + 
                media_prots_totais + gords_totais + (1| area), 
              family = "binomial"(link = "logit"), 
              data = quality, na.action = "na.fail")


####selecting variables and collecting data
my_sel=dredge(m2)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) 
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

##plotting figure 6B
mpq=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="firebrick", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour="firebrick")+theme_classic(base_size = 20)+ 
  ylab("Coefficient values [95%CI]")+
  xlab("")
mpq = mpq+ scale_x_discrete(labels=c("media_carbs_totais" = "Carbohydrate", 
                                     "media_prots_totais" = "Protein",
                                     "gords_totais" = "Lipid"))
mpq

######plotting figures in one page (figure 6A-B)
library(ggpubr)
tiff("figure_rate_traits.tiff", units="in", width=12, height=5, res=300)
ggarrange(mhq, mpq, labels = c("A", "B"),
          ncol = 2, nrow = 1, common.legend = T, align = "h")
dev.off()





######################################################
########now lets fits quantity traits ###
dado=read.delim("consumption_quantity.csv", sep= ";", header = TRUE)
names(dado)
quatity = dado%>%
  select(taxa_pred_par, taxa_pred_herb, area, area_semente, 
         biofruto, number_seed_avalability)
str(quatity)
quatity$area = as.factor(quatity$area)

###exploiting data
plot(taxa_pred_herb ~ area_semente, data = quatity)
plot(taxa_pred_herb ~ biofruto, data = quatity)
plot(taxa_pred_herb ~ number_seed_avalability, data = quatity)

plot(taxa_pred_par ~ area_semente, data = quatity)
plot(taxa_pred_par ~ biofruto, data = quatity)
plot(taxa_pred_par ~ number_seed_avalability, data = quatity)


###fitting model on herbivory
names(data)
m1 = glmer.nb(taxa_pred_herb ~ area_semente + 
                biofruto + number_seed_avalability+(1| area), 
              family = "binomial"(link = "logit"), data = quatity, 
              na.action = "na.fail")
summary(m1)

####Selecting variable and collecting data####
my_sel=dredge(m1)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) 
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

##plotting herbivory rate on resources quantity (not shown in the manuscript)
mht=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="blue", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour = "blue")+theme_classic(base_size = 20)+
  ylab("Coefficient values [95%CI]")+  xlab("Parameters")
mht = mht+ scale_x_discrete(labels=c("area_semente" = "Seeds size (mm²)", 
                                     "biofruto" = "Fruit biomass (g)",
                                     "number_seed_avalability" = "Seeds number"))
mht


### fitting model on parasitoid rate
m2 = glmer.nb(taxa_pred_par ~ number_seed_avalability + 
                area_semente + biofruto + (1| area), 
              family = "binomial"(link = "logit"), 
              data = quatity, na.action = "na.fail")
summary(m2)

####Selecting variable####
my_sel=dredge(m2)
my_sel_d2=get.models(my_sel,cumsum(weight) <= .95) #2 modelos
my_sel_d2
avg_model=model.avg(my_sel_d2)
avg_model
m_final=summary(avg_model) 
m_final
confint(avg_model,full = F)
importance(subset(model.sel(my_sel, rank = AICc),
                  cumsum(weight) <= .95))

df1 = as.data.frame(m_final$coefmat.subset) 
CI = as.data.frame(confint(avg_model, full=F)) 
df1$CI.min = CI$`2.5 %` 
df1$CI.max = CI$`97.5 %`
setDT(df1, keep.rownames = "coefficient") 
names(df1) = gsub(" ", "", names(df1))

##plotting parasitoid rate on defense chemicals (not shown in the manuscript)
mpt=ggplot(data=df1[2:4,], aes(x=coefficient, y=Estimate))+ 
  geom_hline(yintercept=0, color = "gray",linetype="dashed", lwd=1.5)+ 
  geom_errorbar(aes(ymin=CI.min, ymax=CI.max), colour="firebrick", 
                width=0, lwd=1.5) +
  coord_flip()+ 
  geom_point(size=5, colour="firebrick")+theme_classic(base_size = 20)+ ylab("Coefficient values [95%CI]")+
  xlab("")
mpt = mpt+ scale_x_discrete(labels=c("area_semente" = "Seeds size (mm²)", 
                                     "biofruto" = "Fruit biomass (g)",
                                     "number_seed_avalability" = "Seeds number"))
mpt


###end code






