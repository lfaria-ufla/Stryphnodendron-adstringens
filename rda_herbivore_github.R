rm(list=ls())
setwd("~/estatistica/R/juliana")
library(vegan)
library(tidyverse)

dado=read.delim("nmds_rda_herbivore.csv", sep= ";", header = TRUE)
names(dado)
dados=na.exclude(dado)
summary(dados)


######organizing data in two sheets (community species and traits)
tracos=dados%>%
 select(area,tratamento, condensed_g_1g_plant, totalphenols_g_1g_plant, hydrolyzable_g_1g_plant,
        media_carbs_totais, media_prots_totais, gords_totais, 
        biofruto, area_semente, number_seed_avalability)
write.csv(tracos, "tracos.csv")

composicao = dados %>%
  select(bruchinae:blattaria_sp1)
write.csv(composicao, "composic.csv")

composicao=read.delim("composic.csv", sep= ",", header = TRUE)
tracos=read.delim("tracos.csv", sep= ",", header = TRUE)

tracos=tracos%>%
  select(area,tratamento, condensed_g_1g_plant, totalphenols_g_1g_plant, hydrolyzable_g_1g_plant,
         media_carbs_totais, media_prots_totais, gords_totais, 
         biofruto, area_semente, number_seed_avalability)

###getting names to figure
tracos= rename(tracos, "fruit biomass" = biofruto, "seed size" = area_semente, 
              "seed numbers" = number_seed_avalability,
             "total carbs" = media_carbs_totais)
composicao = composicao %>%
  select(bruchinae:blattaria_sp1)

summary(tracos)
summary(composicao)
str(composicao)

###testing for rda
### first null model
mod0 <- dbrda(composicao ~ 1, tracos)
mod0
plot(mod0)
## Now all traits (Full model)
modgeral <- dbrda(composicao ~ ., tracos)
modgeral
vif.cca(modgeral)
plot(modgeral)

## Automatic selection of variables by permutation P-values
mod <- ordistep(mod0, scope=formula(modgeral))
mod$anova
plot(mod)
## Permutation test for all variables
anova(mod)
## Permutation test or significance when a term
## is added to the model after all other terms
anova.cca(mod, by = "margin")

##final model
m1 = dbrda(composicao~`seed numbers`, data=tracos)
rda = summary(m1)
rda

###pplotting figure 3B
with(tracos, levels(as.factor(tratamento)))
scl <- 2 ## scaling = 2
colvec <- c("red2", "mediumblue", "green")
bp=scores(mod, display = "bp")

tiff("traits_herb.rda.tiff", units="in", width=10, height=8, res=300)
ordiplot(mod,scaling = scl, type = "n", xlim=c(-10,4.0), ylim = c(-4,4.0), 
         xlab = "dbRDA (100%)", ylab = "MDS1")
with(mod,points(mod,display = "sites",
                pch = c(3 ,4, 8) [as.numeric(as.factor(tracos$tratamento))], 
                col = c("blue","red") [as.numeric(as.factor(tracos$area))], 
                bg = c("black", "black", "green", "red"), 
                cex = 1.5))
text(mod, display = "bp", col = "black", cex = 1.3)
with(mod, legend("topleft", legend = c(levels(as.factor(tracos$tratamento)),
                                          levels(as.factor(tracos$area))), 
                 pch = c(3 ,4, 8, 16, 16), title = c("Treatment/Area(colored)"),
                 col = c("black","black", "black", "blue","red"),
                 bty = "n")) # displays symbol and colour legend
with(mod, legend("topright", legend = "B", bty = "n"))
with(tracos, ordiellipse(m1, tratamento, label = F, 
      col=c("purple","orange", "darkgreen"), lty = c(1, 1, 1)),
      pch = c(1,1, 1) (as.factor(tracos$tratamento)), 
      col = c("blue","red", "green") (as.factor(tracos$tratamento)), 
      cex = 2)
dev.off()


###end code



