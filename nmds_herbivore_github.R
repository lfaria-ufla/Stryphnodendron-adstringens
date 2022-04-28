rm(list=ls())
setwd("~/estatistica/R/juliana")
dados=read.delim("nmds_rda_herbivore.csv", sep= ";", header = TRUE)
names(dados)
summary(dados)

library(vegan)
library(tidyverse)

###organizing data in two sheets (community species and treatments)
dados = dados%>%
  select(tratamento, area, bruchinae:blattaria_sp1)%>%
  drop_na()
summary(dados)

##treatments
tracos=dados%>%
  select(area,tratamento)
write.csv(tracos, "tracos.csv")
tracos=read.delim("tracos.csv", sep= ",", header = TRUE)
tracos=tracos%>%
  select(area,tratamento)

#composition
composicao = dados %>%
  select(bruchinae:blattaria_sp1)
write.csv(composicao, "composic.csv")
composicao=read.delim("composic.csv", sep= ",", header = TRUE)
composicao = composicao %>%
  select(bruchinae:blattaria_sp1)

composicao= rename(composicao, "A. gregorioi" = bruchinae, 
                   "Lepidoptera sp.1" = lep_sp1, 
                   "Lepidoptera sp.2" = lep_sp2, 
                   "Lepidoptera sp.3" = lep_sp3, 
                   "Diptera sp.1" = diptera_sp1)

summary(tracos)
summary(composicao)
str(composicao)

###decomposing
compos.rel = decostand(composicao, method="hellinger")

##eval matriz of distance
compos.mat = vegdist(composicao, method = "bray")
compos.mat = as.matrix(compos.mat, labels=T) 
write.csv(compos.mat, "compos.mat.csv")

##nmds
nmd=metaMDS(compos.mat, k = 2,maxit = 999, trymax=500)
nmd
barba.nmds=metaMDS(compos.mat, k = 2,maxit = 999, previous.best = nmd)
barba.nmds

##species score
barba.nmds.scr =  `sppscores<-`(barba.nmds, composicao)

##species correlation
barba.nmds.cor = cor(composicao, barba.nmds$points,
                     use="complete.obs",
                     method = "pearson")
barba.nmds.cor
write.csv(barba.nmds.cor, file = "barba.nmds.PearsonCor.csv")

#Shepards test/goodness of fit
gof = goodness(barba.nmds)
gof
stressplot(barba.nmds)

##testing species importance
set.seed(123) 
barba.spp.fit = envfit(barba.nmds, composicao, permutations = 999, 
                       strata=tracos$tratamento)
barba.spp.fit

###plotting figure 3A
tiff("nmds_herb.tiff", units="in", width=10, height=8, res=300)
plot(barba.nmds, type = "n")
with(tracos,points(barba.nmds,display = "sites", cex = 1.5,
                   pch = c(21,22,23,24) [as.numeric(as.factor(tracos$tratamento))], 
                   col = c("blue","red") [as.numeric(as.factor(tracos$area))]))
with(tracos, legend("topright", legend = "A", bty = "n"))

ordihull(barba.nmds,
         tracos$tratamento,
         display = "sites",
         draw = "polygon",
         col = NULL,
         border = c("purple", "firebrick", "orange", "darkgreen"),
         lty = c(1, 1, 1, 1),
         lwd = 2, 
         cex = 2)

legend("bottomleft", legend = c(levels(as.factor(tracos$tratamento)),
                                 levels(as.factor(tracos$area))), 
       pch = c(21,22,23,24, 19,19), title = c("Treatment/Area(colored)"),
       col = c("black","black", "black", "black", "blue", "red"),
       bty = "n", cex = 1) # displays symbol and colour legend
legend("topleft", "stress = 0.1402", bty = "n", cex = 1) 

scrs = scores(barba.nmds, display = "sites", "species")
cent = aggregate(scrs ~ tratamento, data = tracos, FUN = "mean")
names(cent) [-1] <- colnames(scrs)
points(cent [,-1],
       pch = c(16,15,18,17),
       col = c("purple", "firebrick", "orange", "darkgreen"),
       bg = c("black"),
       lwd = 2.0,
       cex = 3.0 # Plots centroids as points on ordination
)
plot(barba.spp.fit, p.max = 0.05, col = "black", cex = 1)
dev.off()

###testing area and treatments
set.seed(123)
area = anosim(composicao, tracos$area, distance = "bray", 
              permutations = 999)
area
set.seed(123) 
tratamento = anosim(composicao, tracos$tratamento, distance = "bray", 
                    permutations = 999)
tratamento

###end code