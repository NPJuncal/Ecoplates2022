################################################################################
# Article: Coupling microbial incubations and FTIR-ATR to assess organic matter 
#compounds persistence in salt marsh soils
# 2/3
# Author: Nerea Piñeiro-Juncal (@NPJuncal)
# V 1.0
# Year: 2025
################################################################################

setwd( "C:/Users/npjun/Dropbox/Proxectos/Rematados/Margarita Salas 2021/Ecoplates 2022/Ecoplates2022")# Set working directory

# Load data -------------------------------------------------------

Data<-read.csv("data/Cores.csv", header=T, sep=";", dec=".")
Data<-as.data.frame(Data)

Data2<-subset(Data, Data$Core_ID == "ATR_C" | Data$Core_ID == "BOL_A" | Data$Core_ID == "BB_A" | Data$Core_ID == "J_C")


# OM and CaCO3 -------------------------------------------------------

# signifficative diferences

shapiro.test(Data$OM)### normality (>0.05 normal, <0.05 no normal) 
shapiro.test(log(Data$OM))

anova<-(aov(log(OM) ~ Station, Data)) # ANOVA
summary(anova)
tukey<-TukeyHSD(aov(log(OM) ~ Station, Data))
print(tukey)


## collecting the groups information for the boxplot
library(multcompView)
let<-multcompLetters4(anova,tukey)
print(let)

Tk<-group_by(Data,Station)%>%
  summarise(mean=mean(OM,,na.rm=T),quant=quantile(OM, probs=0.75,na.rm=T)) %>%
  arrange(desc(mean))

let<-as.data.frame.list(let$Station)
Tk$let<-let$Letters
print(Tk)


#CaCO3
shapiro.test(Data$CaCO3) ### normality (>0.05 normal, <0.05 no normal) 
shapiro.test(log(Data$CaCO3))

#summary(aov(CaCO3 ~ Station, Data)) # ANOVA
#TukeyHSD(aov(CaCO3 ~ Station, Data))


pair<-pairwise.wilcox.test(Data$CaCO3, Data$Station,
                           p.adjust.method = "BH")
print(pair)

#manually adding the letter for the boxplots
Tk2<-group_by(Data,Station)%>%
  summarise(mean=mean(CaCO3,,na.rm=T),quant=quantile(OM, probs=0.75,na.rm=T)) %>%
  arrange(desc(mean))
Tk2$let<-c("b","b","b","a")

# plots

Data2$Station = factor(Data2$Station, levels = c("BB","BOL","ATR","J"))

P1<-ggplot(subset(Data2, !is.na(Data2$OM)),aes(Depth,OM))+
  ggtitle("A") +
  xlab("Depth (cm)")+ ylab("%")+
  geom_point(aes( shape="Organic matter"), color= "black")+
  geom_line( color= "black")+
  
  geom_point(aes(Depth,CaCO3, shape="Carbonates"), color= "grey56")+
  geom_line( aes(Depth,CaCO3), color= "grey56")+
  scale_x_reverse()+
  coord_flip()+
  facet_grid( ~ Station, scale ='fixed')+
  theme(plot.title = element_text(hjust = 0)
        ,legend.title=element_blank()
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed"))

P1

Data$Station = factor(Data$Station, levels = c("BB","BOL","ATR","J"))

P2<-ggplot(Data,aes(x=factor(Station), OM)) + ylab("Organic Matter % (LOI 450ºC)")+xlab("Station")+
  ggtitle("B") +
  geom_boxplot(aes(x=factor(Station), OM),fill="grey",alpha=0.5)+
  #geom_jitter()+
  geom_text(data=Tk, aes(label=let, x=Station, y=50), size=4)+
  theme(panel.background = element_rect(fill="white", color = "black"
                                        #, colour, size, linetype, color
  )
  ,panel.grid.major = element_line(colour="grey",linetype="dashed"))
#hjust=-1,

P2

P3<-ggplot(Data,aes(x=factor(Station), CaCO3)) + labs(y = expression("CaCO"["3"] ~ "% (LOI 950ºC)"))+xlab("Station")+
  ggtitle("C") +
  geom_boxplot(aes(x=factor(Station), CaCO3),fill="grey",alpha=0.5)+
  #geom_jitter()+
  geom_text(data=Tk2, aes(label=let, x=Station, y=17), size=4)+
  theme(panel.background = element_rect(fill="white", color = "black"
                                        #, colour, size, linetype, color
  )
  ,panel.grid.major = element_line(colour="grey",linetype="dashed"))

P3

OM_CaCO3<-grid.arrange(P1, P2, P3, nrow=2 ,ncol = 2, layout_matrix = cbind(c(1,2),c(1,3)))

ggsave("results/OM-CaCO3.jpg",OM_CaCO3, units="cm", width = 19, height = 15)


aggregate(Data[, c(2,13,14)], list(Data$Station), mean, na.rm = TRUE)
aggregate(Data[, c(2,13,14)], list(Data$Station), sd, na.rm = TRUE)


# grain size -------------------------------------------------------

DataGr <- melt(Data2[,c(2,5:12)], id.vars = c("Depth","Station"))
DataGr$variable <- factor(DataGr$variable,levels = c("Silt_Clay","VF_Sand","F_Sand","M_Sand","C_Sand","VC_Sand","Gravel"))
DataGr$value<-DataGr$value*100

P1<-
  
  ggplot(DataGr, aes(fill=variable, y=value, x=Depth)) + ylab("Proportion of grain size particles")+
  coord_flip()+
  ggtitle("A") +
  scale_x_reverse()+
  scale_y_continuous(labels = function(x) {
    x[seq_along(x) %% 2 == 1] <- ""
    x
  })+
  facet_grid( ~ Station, scale ='fixed')+
  geom_bar(position="fill", stat="identity")+
  #scale_fill_grey(start = 0, end = .9)+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0))


P2<-ggplot(Data,aes(x=factor(Station), Silt_Clay)) + ylab("Silt and Clay (<63 µm) % ")+
  geom_boxplot(aes(x=factor(Station), Silt_Clay),fill="grey",alpha=0.5)+
  geom_jitter()+
  ggtitle("B") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed")
  )

P3<-ggplot(Data,aes(x=factor(Station), VF_Sand)) + ylab("Very fine Sand (63-125 µm) %")+xlab("Station")+
  ggtitle("C") +
  geom_boxplot(aes(x=factor(Station), VF_Sand),fill="grey",alpha=0.5)+
  geom_jitter()+
  theme(panel.background = element_rect(fill="white", color = "black"
                                        #, colour, size, linetype, color
  )
  ,panel.grid.major = element_line(colour="grey",linetype="dashed"))



GranDis<-grid.arrange(P1, P2, P3, ncol = 3, 
                      layout_matrix = cbind(c(1,1), c(1,1),c(2,3)))

ggsave("results/GranDis.jpg",GranDis, units="cm", width = 19, height = 15)

