dat <- read.table('/home/bertrand/Cours Bioinfo/Projet_3/Fst_1-2.txt')
dat1 <- read.table('/home/bertrand/Cours Bioinfo/Projet_3/Fst_1-3.txt')
dat2 <- read.table('/home/bertrand/Cours Bioinfo/Projet_3/Fst_1-4.txt')

colnames(dat)<- c('chro','pos','Fst')
colnames(dat1)<- c('chro','pos','Fst')
colnames(dat2)<- c('chro','pos','Fst')

library(dplyr)
visuFstChro <- function(Chro) {
datAnd <- dat%>%filter(chro==Chro) 
datAnd$File <- "G0-G12_cerise"
datAnd1 <- dat1%>%filter(chro==Chro) 
datAnd1$File <- "G0-G12_cran"
datAnd2 <- dat2%>%filter(chro==Chro) 
datAnd2$File <- "G0-G12_fraise"
datF <- rbind(datAnd,datAnd1,datAnd2)
library(ggplot2)
ggplot(datAnd, aes(x=pos, y=Fst)) +   geom_line()
ggplot(datF, aes(x=pos, y=Fst, color=File)) + geom_line() + 
  xlab(Chro) + facet_grid(File ~ .) + theme(legend.position = "none")
}
visuFstChro('Andromeda')
