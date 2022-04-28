rm(list=ls())
setwd("~/estatistica/R/juliana")
dados=read.delim("general_data.csv", sep= ";", header = TRUE)
names(dados)
str(dados)
library(MuMIn)
library(lme4)
library(lmerTest)
library(Rcpp)

####################################################
####start fitting GLMM for resources quantity traits
#organizing quantity traits data
names(dados)
quantity=dados[,c("area", "trat_cut", "trat_fert",
                 "biofruto", "number_seed_avalability", "area_semente")]

###biomass
biomassa=subset(quantity, quantity$biofruto>0)
biomassas=na.exclude(biomassa)
hist(biomassa$biofruto)


####fitting GLMM
m1=lmer(biofruto~trat_cut*trat_fert+(1|area), data=biomassas)
shapiro.test(resid(m1))
summary(m1)


###plotting figure 1
library(sjlabelled)
library(sjmisc)
library(snakecase)
library(ggeffects)
library(ggplot2)
library(sjPlot)

##figure 1D
b = plot_model(m1, type = "pred", terms =  c("trat_cut", "trat_fert"), 
                sort.est = T, title = "", legend.title = "Fertilization", 
                colors = c("firebrick", "blue"))
b1 = b + labs(x = "Clipping", y = "Fruits Biomass (g)")+
  theme_classic(base_size = 12) + 
  scale_x_discrete(limits =c("No", "Yes"))
b1


##seed numbers
number=subset(quantity, quantity$number_seed_avalability>0)
number=na.exclude(number)
hist(number$number_seed_avalability)

### fitting model GLMM
m2=glmer.nb(number_seed_avalability~trat_cut*trat_fert+(1|area), 
         family = "poisson", data=number)
shapiro.test(resid(m2))
summary(m2)

###plotting figure 1E
n = plot_model(m2, type = "pred", terms =  c("trat_cut", "trat_fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"))
n1 = n + labs(x = "Clipping", y = "No. of Seeds")+
  theme_classic(base_size = 12) + 
  scale_x_discrete(limits =c("No", "Yes"))
n1

##seed size
area=subset(quantity, quantity$area_semente>0)
area=na.exclude(area)
hist(area$area_semente)

###fitting GLMM
m3=lmer(number_seed_avalability~trat_cut*trat_fert+(1|area), data=area)
shapiro.test(resid(m3))
summary(m3)

###plotting figure 1F
s = plot_model(m3, type = "pred", terms =  c("trat_cut", "trat_fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"))
s1 = s + labs(x = "Clipping", y = "Seeds Size (mm²)")+
  theme_classic(base_size = 12) + 
  scale_x_discrete(limits =c("No", "Yes"))
s1


##################################################################
###next we fit GLMM for defense chemicals in fruits
dados=read.delim("consumption_defense.csv", sep= ";", header = TRUE)
names(dados)

str(dados)
dados$trat_cut = as.factor(dados$trat_cut)
dados$trat_fert = as.factor(dados$trat_fert)
dados$area = as.factor(dados$area)


##fitting fenols fruits GLMM
hist(dados$total_totalfenols)
m1=glmer.nb(total_totalfenols~trat_cut*trat_fert+(1|area), data=dados)
shapiro.test(resid(m1))
summary(m1)


###plotting figure 1C
p = plot_model(m1, type = "pred", terms = c("trat_cut", "trat_fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"), grid = F)
p1 = p+ labs(x = "Clipping", y = "Total Phenols (Fruits)")+
  theme_classic(base_size = 12) +
  scale_x_discrete(limits =c("No", "Yes"))
p1


##total_hydrolizable fruits GLMM
hist(dados$total_hydrolizable)
m2=glmer.nb(total_hydrolizable~trat_cut*trat_fert+(1|area), data=dados)
shapiro.test(resid(m2))
summary(m2)

###plotting figure 1A
h = plot_model(m2, type = "pred", terms = c("trat_cut", "trat_fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"), grid = F)
h1 = h+ labs(x = "Clipping", y = "Hydrolyzable tannins (Fruits)")+
  theme_classic(base_size = 12) +
  scale_x_discrete(limits =c("No", "Yes"))
h1


##condense fruit GLMM
hist(dados$total_condensed)
m3=glmer.nb(total_condensed~trat_cut*trat_fert+(1|area), data=dados)
shapiro.test(resid(m3))
summary(m3)

###plotting figure 1B
c = plot_model(m3, type = "pred", terms = c("trat_cut", "trat_fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"), grid = F)
c1 = c+ labs(x = "Clipping", y = "Condensed tannins (Fruits)")+
  theme_classic(base_size = 12) +
  scale_x_discrete(limits =c("No", "Yes"))
c1

############################
###chemical defense on seeds
dados=read.delim("seed_defense.csv", sep= ";", header = TRUE)
names(dados)

##fenols seed (no data from cut treatment) GLMM
hist(dados$total_totalfenols)
m4=glmer.nb(total_totalfenols~trat_fert+(1|area), data=dados)
shapiro.test(resid(m4))
summary(m4)

####plotting figure 1G
ts = plot_model(m4, type = "pred", terms =  "trat_fert", 
                sort.est = T, title = "", legend.title = "Fertilization", 
                colors = "blue")
ts1 = ts + labs(x = "Fertilization", y = "Total phenols (Seeds)")+
  theme_classic(base_size = 12) + 
  scale_x_discrete(limits =c("No", "Yes"))
ts1

##condensed tannins GLMM
hist(dados$total_condensed)
m5=glmer.nb(total_condensed~trat_fert+(1|area),data=dados)
shapiro.test(resid(m5))
summary(m5)


##hydrolizable tannins GLMM
hist(dados$total_hydrolizable)
m6=glmer.nb(total_hydrolizable~trat_fert+(1|area), data=dados)
shapiro.test(resid(m6))
summary(m6)



######plotting figures in one page
library(ggpubr)
tiff("figure1.tiff", units="in", width=18, height=15, res=300)
leg = get_legend(h1)
h1 = h1+theme(legend.position = "none")
c1 = c1+theme(legend.position = "none")
p1 = p1+theme(legend.position = "none")
n1 = n1+theme(legend.position = "none")
s1 = s1+theme(legend.position = "none")
b1 = b1+theme(legend.position = "none")
ggarrange(h1,c1,p1,b1,n1,s1,ts1,leg,  labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 3, nrow = 3)
dev.off() 



################################################################
###finally we fit GLMM for resource quality traits for fruits and seeds
#Organizing data by quality variables for fruits data
quality=dados[,c("area", "trat_cut", "trat_fert",
                 "media_prots_totais", "media_carbs_totais", "gords_totais")]

##lipids
lipideo=subset(quality, quality$gords_totais>0)
lipideo=na.exclude(lipideo)
hist(dados$gords_totais, breaks = 3)

###fitting GLMM for lipids (we had hard time fitting lipids)
m1=glmer.nb(gords_totais~trat_cut*trat_fert+(1|area), data=dados)
shapiro.test(resid(m1))
summary(m1)

##fitting carbs 
hist(dados$media_carbs_totais, breaks = 5)
m2=lmer(media_carbs_totais~trat_cut*trat_fert+(1|area), data=dados)
summary(m2)

##fitting proteins
hist(dados$media_prots_totais)
m3=lmer(media_prots_totais~trat_cut*trat_fert+(1|area), data=dados)
shapiro.test(resid(m3))
summary(m3)

###################
###fitting GLMM for resource quality regarding seeds
###seeds#######
dados<-read.table("seeds_resource_quality.csv",sep=";", h=T)
names(dados)

###GLMM carbs
hist(dados$media_carbs_totais)
m4=lmer(media_carbs_totais~trat_fert+(1|area),data=dados)
summary(m4)

##protein
hist(dados$media_prots_totais)
m5=lmer(media_prots_totais~trat_fert+(1|area),data=dados)
summary(m5)

##lipid (we had hard time fitting lipids)
hist(dados$gords_totais)
m6=glmer.nb(gords_totais~trat_fert+(1|area), data=dados)
shapiro.test(resid(m6))
summary(m6)


##end code