#cabeçalho
rm(list=ls())
setwd("~/estatistica/R/juliana")
dados=read.delim(file= "consumption.csv", sep= ";", header = T)
names(dados)
str(dados)
dados$area = as.factor(dados$area)
dados$fert = as.factor(dados$fert)
dados$cut = as.factor(dados$cut)
str(dados)
summary(dados)

mean(dados$taxa_pred_par[dados$fert == "yes" & dados$cut == "yes"])
sd(dados$taxa_pred_par[dados$fert == "yes" & dados$cut == "yes"])

mean(dados$taxa_pred_herb[dados$fert == "yes" & dados$cut == "no"])
sd(dados$taxa_pred_herb[dados$fert == "yes" & dados$cut == "no"])


###models
library(lme4)
library(parameters)
library(performance)

##model for herbivore consumption
m1 = glmer.nb(taxa_pred_herb ~ cut * fert + (1| area), 
           family = "binomial"(link = "logit"), data = dados)
summary(m1)
model_parameters(m1)
check_model(m1)

####overdispersion test####
overdisp_fun = function(model) {
   vpars = function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df = sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf = nrow(model.frame(model))-model.df
  rp = residuals(model,type="pearson")
  Pearson.chisq = sum(rp^2)
  prat = Pearson.chisq/rdf
  pval = pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(m1)


##model for parasitoid consumption
m2 = glmer.nb(taxa_pred_par ~ cut * fert + (1| area), 
              family = "binomial"(link = "logit"), data = dados)
summary(m2)
model_parameters(m2)
check_model(m2)

m2n = glmer.nb(taxa_pred_par ~ 1+ (1| area), 
              family = "binomial"(link = "logit"), data = dados)
anova(m2,m2n, test = "LTR")

###plot figure 5
library(sjlabelled)
library(sjmisc)
library(snakecase)
library(ggeffects)
library(ggplot2)
library(sjPlot)

##figure 5A
h = plot_model(m1, type = "pred", terms =  c("cut", "fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"))
h1 = h + labs(x = "", y = "Proportion of seed attacked by herbivores")+
  theme_classic(base_size = 15) + 
  scale_x_discrete(limits =c("No", "Yes"))
h1

##figure 5B
p = plot_model(m2, type = "pred", terms =  c("cut", "fert"), 
               sort.est = T, title = "", legend.title = "Fertilization", 
               colors = c("firebrick", "blue"))
p1 = p + labs(x = "Clipping", y = "Proportion of seed attacked by parasitoids")+
  theme_classic(base_size = 15) + 
  scale_x_discrete(limits =c("No", "Yes"))
p1


######plotting figures in one page
library(ggpubr)
tiff("figure_rate.tiff", units="in", width=10, height=10, res=300)
ggarrange(h1,p1,  labels = c("A", "B"),
          ncol = 1, nrow = 2, common.legend = T)
dev.off() 



###end code
