################################################################################
# Article: Coupling microbial incubations and FTIR-ATR to assess organic matter 
#compounds persistence in salt marsh soils
# 3/3
# V 1.2
# Author: Nerea Piñeiro-Juncal (@NPJuncal)
# Year: 2025
################################################################################


setwd( "C:/Users/npjun/Dropbox/Proxectos/Rematados/Margarita Salas 2021/Ecoplates 2022/Ecoplates2022")# Set working directory

library(ggplot2)
library(dplyr)
library(psych)
library(corrplot)
library(matrixStats)
library(TTR)
library(reshape)
library(chemometrics)
library(ggrepel)



# PCA transpuesta + Ecoplates ---------------------------------------------


File<-"data/Spectra.Ecop.csv"

A<-read.csv(File, header=T, sep=";", dec=".")
A<-as.data.frame(A)

head(A)

pca<-A[,-1]

PCA<-principal(pca, nfactors=20, residual=T, rotate="varimax", covar=F)
print(PCA)

# loading #
loa <-loadings(PCA)
loa <-as.data.frame(loa[,1:ncol(loa)])
write.csv(loa,"results/loa_t_Ecop.csv",sep=";", dec=",")

# scores #
sco<-PCA$scores
sco<-as.data.frame(sco)
sco<-cbind(A[,1],sco)
sco[, 1] <- sapply(sco[, 1], as.numeric)



# plots
colnames(sco)<-c("WN","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9"
                 ,"PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

colnames(loa)<-c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9"
                 ,"PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

loa.s<-loa[c(1:24),]
loa.s$Depth<-c(4,16,26,42,4,16,26,42,54,66,86,4,16,26,42,54,66,76,4,16,26,42,54,64)
loa.s$Station<-c("BB","BB","BB","BB","BOL","BOL","BOL","BOL","BOL","BOL","BOL","ATR","ATR","ATR","ATR","ATR","ATR","ATR","J","J","J","J","J","J")



write.csv(sco,"results/Sco_t_Ecop.csv", sep=";", dec=".")


Spec<-ggplot(A, aes(WL,runMean(scale(F4.1))))+xlab("Wavenumber (cm-1)")+ylab("Scaled abs and scores")+
  geom_line(size=1)+
  #geom_line(aes(WL,runMean(scale(C4.1))),col="grey50",size=1)+
  geom_line(aes(WL,runMean(scale(sco$PC17))),col="blue",size=1)+
  annotate(geom="text", x=3000, y=2, label="PC17",color="blue")+
  annotate(geom="text", x=3000, y=3, label="Glycyl-L-Glutamic Acid",color="black")+
  #annotate(geom="text", x=3000, y=5, label="L-Phenilalanine",color="grey50")+
  scale_x_reverse()


Loadings<-ggplot(loa.s,aes(Depth,PC17))+ ylab("Loadings")+
  geom_point(aes(color=Station), size=3)+
  geom_line(aes(color=Station), size=1)+
  scale_color_manual(values=c("blue","red", "orange", "black"))

  
g<-gridExtra::grid.arrange(Spec, Loadings, ncol=2)

ggsave("results/Eco.PC17.jpg", g, width = 19, height =5, units = "cm")  
  



# PCA loadings + geochemistry  --------------------------------------------


File<-"data/FTIR_GE.csv"

FTIR_GE<-read.csv(File, header=T, sep=";", dec=".")
FTIR_GE<-as.data.frame(FTIR_GE)


#transform close data and scale
E<-FTIR_GE[,c(21:22)]# clr transformation of OM and carbonates data
E_clr<-clr(E)
FTIR_GE[,c(21:22)]<-E_clr


Gr<-FTIR_GE[,c(23:29)]# clr transformation of granulometric data
Gr<-Gr+1
Gr_clr<-clr(Gr)
FTIR_GE[,c(23:29)]<-Gr_clr


#delet those PCs that have negative loadings in the samples (compunds not present in the soil): 
#α-Ketobutyric Acid (PC4), L-Phenilalanine (PC5), Itaconic Acid (PC7), 2-Hydroxy Benzoic Acid (PC10), and Glycyl-L-Glutamic Acid (PC17)
pca<-scale(FTIR_GE[,c(4:6,9,11,12,14:19,21:33)])



# PCA with varimax rotation #
PCA<-principal(pca, nfactors=1, residual=T, rotate="varimax", covar=F) #Modify number of factor until you have the maximun number of components with more than 1 explained variance
AutVal<-sum(PCA$values>1)
PCA<-principal(pca, nfactors=AutVal, residual=T, rotate="varimax", covar=F)
print(PCA)

fa.parallel(pca)


# loading #
loa <-loadings(PCA)
loa <-as.data.frame(loa[,1:ncol(loa)])

# comunality #

loa2 <- loa^2 #Cuadrados de los loadings
Com <-t(loa2) # Transpuesta

barplot(Com, col=c("skyblue1","tan1","springgreen2","lightgoldenrod2","violetred3","aquamarine", "chartreuse4"),
        ylim=c(0,1), axisnames = T, las=2, cex.names=0.8,xaxt="s",main="Fraccionamiento comunalidad")# barplot

write.csv(loa,"results/FTIR_GE_loa.csv",sep=";", dec=",")

# scores #
sco<-PCA$scores
sco<-as.data.frame(sco)
sco<-cbind(FTIR_GE[,1],sco)
write.csv(sco,"results/FTIR_GE_sco.csv",sep=";", dec=",")

### biplots
colnames(sco)<-c("Sample","PC1","PC2","PC3","PC4","PC5","PC6")
colnames(loa)<-c("PC1","PC2","PC3","PC4","PC5","PC6")
fit<- as.data.frame(sco)
Station<-FTIR_GE[,2]
Depth<-FTIR_GE[,3]
sco$depth<-Depth


p1<-  ggplot(data=loa, aes(x=PC1,y=PC2))+
  #geom_point(data=fit,aes(x=PC1,y=PC2,colour=Station), size=2)+
  geom_text(data=fit,aes(x=PC1,y=PC2,label=sco$depth,colour=Station)
            ,fontface = "bold",position=position_jitter(width=1,height=1))+
  #coord_fixed(ratio=1)+ 
  theme_bw() +
  guides(color = "none")+
  xlab("PC1 (26%)")+ylab("PC2 (22%)")+ #theme(legend.position = "none")+
  #xlim(c(-4.8,4.8))+
  # ylim(c(-3,5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")

p1


p2<-ggplot(data=loa, aes(x=PC3,y=PC4)) +
  #geom_point(data=fit,aes(x=PC3,y=PC4,colour=Station), size=2)+ 
  geom_text(data=fit,aes(x=PC3,y=PC4,label=sco$depth,colour=Station)
            ,fontface = "bold",position=position_jitter(width=1,height=1))+
  #coord_fixed(ratio=1)+
  theme_bw() +
  guides(color = "none")+
  xlab("PC3 (11%)")+ylab("PC4 (11%)")+ #theme(legend.position = "none")+
  #xlim(c(-4.8,4.8))+
  # ylim(c(-3,5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")

p3<-ggplot(data=loa, aes(x=PC5,y=PC6)) +
  #geom_point(data=fit,aes(x=PC5,y=PC6,colour=Station), size=2)+ 
  geom_text(data=fit,aes(x=PC5,y=PC6,label=sco$depth,colour=Station)
            ,fontface = "bold",position=position_jitter(width=1,height=1))+
  #coord_fixed(ratio=1)+ 
  theme_bw() +
  guides(color = "none")+
  xlab("PC5 (7%)")+ylab("PC6 (7%)")+ #theme(legend.position = "none")+
  #xlim(c(-4.8,4.8))+
  # ylim(c(-3,5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")


loa.p1<-subset(loa, abs(PC1) > 0.5 |abs(PC2) > 0.5 )
loa.p2<-subset(loa, abs(PC3) > 0.5 |abs(PC4) > 0.5 )
loa.p3<-subset(loa, abs(PC5) > 0.5 |abs(PC6) > 0.5 )

p4<-ggplot(data=loa.p1, aes(x=PC1,y=PC2)) +theme_bw() +
  geom_text_repel(x=(loa.p1$PC1), y=(loa.p1$PC2), label=row.names(loa.p1), col="red3", size= 4)+
  #annotate("text", x=(loa.p1$PC1*8.8), y=(loa.p1$PC2*8.8), label=row.names(loa.p1), col="red3", size= 5)+
  xlab("PC1 (26%)")+ylab("PC2 (22%)")+ #theme(legend.position = "none")+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")+
  geom_segment(data=loa.p1, aes(x =0 , y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(1/2, 'picas')), color = "black", linetype='solid', size=0.5)

p4


p5<-ggplot(data=loa.p2, aes(x=PC3,y=PC4)) +theme_bw() +
  geom_text_repel(x=(loa.p2$PC3), y=(loa.p2$PC4), label=row.names(loa.p2), col="red3", size= 4)+
  #annotate("text", x=(loa.p2$PC1*8.8), y=(loa.p2$PC2*8.8), label=row.names(loa.p2), col="red3", size= 5)+
  xlab("PC3 (11%)")+ylab("PC4 (11%)")+ #theme(legend.position = "none")+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")+
  geom_segment(data=loa.p2, aes(x =0 , y = 0, xend = PC3, yend = PC4), arrow = arrow(length = unit(1/2, 'picas')), color = "black", linetype='solid', size=0.5)

p6<-ggplot(data=loa.p3, aes(x=PC5,y=PC6)) +theme_bw() +
  geom_text_repel(x=(loa.p3$PC5), y=(loa.p3$PC6), label=row.names(loa.p3), col="red3", size= 4)+
  #annotate("text", x=(loa.p3$PC1*8.8), y=(loa.p3$PC2*8.8), label=row.names(loa.p3), col="red3", size= 5)+
  xlab("PC5 (7%)")+ylab("PC6 (7%)")+ #theme(legend.position = "none")+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")+
  geom_segment(data=loa.p3, aes(x =0 , y = 0, xend = PC5, yend = PC6), arrow = arrow(length = unit(1/2, 'picas')), color = "black", linetype='solid', size=0.5)


library(cowplot)
Biplot<-plot_grid(p1, p4, p2,p5, p3, p6, ncol=2, labels = c("A", "B","C","D","E","F"), align = "v")
Biplot

ggsave("results/Biplot.jpg", Biplot, width = 19, height =20, units = "cm")  


bplabel<-ggplot(data=loa, aes(x=PC1,y=PC2)) +geom_point(data=fit,aes(x=PC1,y=PC2,colour=Station), size=2)+ 
  #coord_fixed(ratio=1)+ 
  theme_bw() +
  xlab("PC1 (26%)")+ylab("PC2 (22%)")+ #theme(legend.position = "none")+
  #xlim(c(-4.8,4.8))+
  # ylim(c(-3,5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  theme(legend.title = element_text(colour = "white"),
        legend.position = "bottom")

ggsave("results/label.jpg", bplabel, width = 9, height =5, units = "cm")  


## significant diferences among station for pca scores

sco$station<-c("BB", "BB", "BB",  "BOL", "BOL", "BOL", "BOL","BOL",
               "ATR", "ATR","ATR","ATR","ATR","ATR","ATR"
               ,"J","J","J","J","J","J")

pairwise.wilcox.test(sco$PC1, sco$station,
                           p.adjust.method = "BH")

vegan::adonis2(
  sco[ , c("PC1", "PC2")] ~ station,
  data = sco,
  method = "euc"
)



pairwise_permanova <- function(sp_matrix, group_var, dist = "bray", adj = "fdr", perm = 10000) {
  
  require(vegan)
  
  ## list contrasts
  group_var <- as.character(group_var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA
  )
  
  for (i in seq(nrow(contrasts))) {
    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i] 
    contrast_matrix <- sp_matrix[sp_subset,]
    
    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist, 
      perm = perm
    )
    
    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
  }
  
  ## adjust p-values for multiple comparisons
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)
  
  return(list(
    contrasts = contrasts, 
    "p-value adjustment" = adj, 
    permutations = perm
  ))
}

groups <- as.character(sco$station) # Apparently, my function doesn't like factors
pairwise_permanova(sco[,c(2:3)], groups)

