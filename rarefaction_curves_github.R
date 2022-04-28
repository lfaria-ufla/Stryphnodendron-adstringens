rm(list=ls())
setwd("~/estatistica/R/juliana")
dado=read.delim("general_community.csv", sep= ";", header = TRUE)
names(dado)
dados=na.exclude(dado)
summary(dados)

library(vegan)
library(tidyverse)
library(iNEXT)

###selecting variables and organizing data by treatments and herbivores
datum=dados %>% 
  select(Tratamento, bruchinae:hy_sp1)

data_h=datum %>% 
  select(Tratamento, bruchinae:hy_sp1)%>%
  gather(key = "species", value = "abundance", bruchinae:hy_sp1)%>%
  group_by(Tratamento, species)%>%
  summarise_if(is.numeric, sum, na.rm = TRUE)%>%
  spread(key = Tratamento, value = abundance)%>%
  select(Co, C, F, CF)

data_h=as.data.frame(data_h)

###rarefaction curve herbivores by treatments
out = iNEXT(data_h, q=0, datatype = "abundance")
ar=ggiNEXT(out, type=1) +theme_classic(base_size=10)+
  theme(legend.title=element_blank())+
  labs(y="Expected no. of herbivores", x = "")+
  theme(axis.title = element_text(size = 15))
ar


###selecting variables and organizing data by treatments and parasitoids
datu= dados%>%
  select(Tratamento, hy_sp2:hy_sp17)

data_p=datu %>% 
  select(Tratamento,hy_sp2:hy_sp17)%>%
  gather(key = "species", value = "abundance", hy_sp2:hy_sp17)%>%
  group_by(Tratamento, species)%>%
  summarise_if(is.numeric, sum, na.rm = TRUE)%>%
  spread(key = Tratamento, value = abundance)%>%
  select(Co, C, F, CF)


data_p=as.data.frame(data_p)

###rarefaction curve parasitoid by treatments
out = iNEXT(data_p, q=0, datatype = "abundance")
ap=ggiNEXT(out, type=1) +theme_classic(base_size=10)+
  theme(legend.title=element_blank())+
  labs(y="Expected no. of parasitoids", x = "No. of individuals sampled")+
  theme(axis.title = element_text(size = 15))
ap

######plotting figures in one page
library(ggpubr)
tiff("figure_2.tiff", units="in", width=10, height=15, res=100)
ggarrange(ar,ap,  labels = c("A", "B"), legend = "right",
          ncol = 1, nrow = 2, common.legend = T, legend.grob = NULL)
dev.off()


###end code