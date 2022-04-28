rm(list=ls())
setwd("~/estatistica/R/juliana")
library(vegan)
library(tidyverse)

dado=read.delim("nmds_rda_paras.csv", sep= ";", header = TRUE)
names(dado)
dados=na.exclude(dado)
dados=dados[-1,]
summary(dados)
###arrumando os dados em dois conjuntos (traços e composição)
tracos=dados%>%
 select(area,tratamento, condensed_g_1g_plant, totalphenols_g_1g_plant, hydrolyzable_g_1g_plant,
        media_carbs_totais, media_prots_totais, gords_totais, 
        biofruto, area_semente, number_seed_avalability)

write.csv(tracos, "tracos.csv")

composicao = dados %>%
  select(hy_sp2:hy_sp17)

write.csv(composicao, "composic.csv")


composicao=read.delim("composic.csv", sep= ",", header = TRUE)
tracos=read.delim("tracos.csv", sep= ",", header = TRUE)

tracos=tracos%>%
  select(area,tratamento, condensed_g_1g_plant, totalphenols_g_1g_plant, hydrolyzable_g_1g_plant,
         media_carbs_totais, media_prots_totais, gords_totais, 
         biofruto, area_semente, number_seed_avalability)
tracos= rename(tracos, "seed size" = area_semente,"total carbs" = media_carbs_totais, 
              "seed numbers" = number_seed_avalability,
             "total proteins" = media_prots_totais)

composicao = composicao %>%
  select(hy_sp2:hy_sp17)

composicao= rename(composicao, "Prodecatoma sp. 1" = hy_sp3, 
                    "Torymus sp. 1" = hy_sp4, 
                    "P. alvarengai" = hy_sp5)
summary(tracos)
summary(composicao)
str(composicao)

###testing for rda
### first null model
mod0 <- dbrda(composicao ~ 1, tracos)
mod0
plot(mod0)
## Full model
modgeral <- dbrda(composicao ~ ., tracos)
modgeral
vif.cca(modgeral)
plot(modgeral)
## Automatic selection of variables by permutation P-values
mod <- ordistep(mod0, scope=formula(modgeral))
mod$anova
plot(mod)
## Permutation test for all variables
set.seed(123)
anova.cca(mod, by = "term", permutations = 999)

#final model
m1 = dbrda(composicao~`seed size`+`total carbs`+`total proteins`, data=tracos)
rda = summary(m1)

### plotting figure 4B
with(tracos, levels(as.factor(tratamento)))
scl <- 2 
colvec <- c("red2", "mediumblue", "green")
bp=scores(mod, display = "bp")

tiff("traits_paras.rda.tiff", units="in", width=10, height=8, res=300)
ordiplot(mod,scaling = scl, type = "n", xlim=c(-5,5.0), ylim = c(-5,5.0), 
         xlab = "dbRDA (84%)", ylab = "dbRDA (11%)")
with(mod,points(mod,display = "sites",
                pch = c(3 ,4, 8) [as.numeric(as.factor(tracos$tratamento))], 
                col = c("blue","red") [as.numeric(as.factor(tracos$area))], 
                bg = c("black", "black", "blue", "red"),
                cex = 1.5))
text(mod, display = "bp", col = "black", cex = 2)
with(mod, legend("topleft", legend = c(levels(as.factor(tracos$tratamento)),
                                          levels(as.factor(tracos$area))), 
                 pch = c(3 ,4, 8, 16, 16), title = c("Treatment/Area(colored)"),
                 col = c("black","black", "black", "blue","red"),
                 bty = "n", cex = 1)) 
with(mod, legend("topright", legend = "B", bty = "n"))
with(tracos, ordiellipse(m1, tratamento, label = FALSE, 
      col=c("purple","orange", "darkgreen")),
      pch = c(1, 1, 1) (as.factor(tracos$tratamento)), 
      col = c("blue","red", "green") (as.factor(tracos$tratamento)), 
      cex = 2)
dev.off()



##end code



