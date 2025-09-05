################################################################################
# Article: Coupling microbial incubations and FTIR-ATR to assess organic matter 
#compounds persistence in salt marsh soils
# 1/2
# Author: Nerea Piñeiro-Juncal (@NPJuncal)
# V 1.0
# Year: 2025
################################################################################


setwd( "C:/Users/npjun/Dropbox/Proxectos/Rematados/Margarita Salas 2021/Ecoplates 2022/Ecoplates2022")# Set working directory

dir.create("results")

#load libraries
library(plyr) #aaply function
library(vegan) #Shannon diversity
library(ggplot2) #figures
library(reshape) #function melt
library(gridExtra) #function grid.arrange
library(dplyr) #%>%
library(chemometrics) # clr function
library(psych) # function principal for PCA with varimax rotation

# load and clean data -----------------------------------------------------

File<-"data/T0.csv"
A<-read.csv(File, header=T, sep=";", dec=".")
T0<-as.data.frame(A)


load_and_clean<-function(name=NULL){
  
  File<-paste("data/",name,".csv", sep = "")
  df<-read.csv(File, header=T, sep=";", dec=".")
  df<-as.data.frame(df)
  
  Fdf<-cbind(df[,c(1:5)],(df[,c(6:101)])-(T0[,c(6:101)]))
  
  #some of the values for the control wells A0,A5,A9 and other wells were there was no consumption are negative
  #measure method error? we check how negative the values
  
  min(na.omit(Fdf[,c(6:101)]))
  
  #substitute negative values for 0
  Fdf[Fdf < 0] <- 0 
  
  write.csv(Fdf,paste("data/F",name, ".csv", sep=""),sep=";", dec=",")
  
  # estimation of the Coefficient of variation among plate replicates
  
  FR<-Fdf[,c(6:37)]#First replicate
  SR<-Fdf[,c(38:69)]#second replicate
  TR<-Fdf[,c(70:101)]#third replicate
  
  all.dat <- list(dat1=FR, dat2=SR, dat3=TR)
  FLR = aaply(laply(all.dat, as.matrix), c(2, 3), mean) # media de las tres replicas
  FLR.sd = aaply(laply(all.dat, as.matrix), c(2, 3), sd) # estandar deviation of three replicates
  FLR.cv<-FLR.sd/FLR # var coef of three replicates
  
  #There are some value in our plates that are very different from their replicates
  #For example HAL C 40-42 ANA T0 (A1: 0.5, A5: 0.1, A9: 0.1)
  #In this case we want to identify the deviated value and eliminate it
  # in some cases, the values are so low that a small deviation results in a cv > 0.3
  # we do not eliminate those values <0.5 in mean absorbance
  #From FLR.cv we are going to get a list of compounds and plates were cv is higher than 0.3
  
  
  # add condition mean > 0.5
  coords<-as.data.frame(which(FLR.cv > 0.3 & FLR > 0.5, arr.ind = TRUE))
  
  
  CV<- data.frame(R1=numeric(), 
                  R2=numeric(),
                  R3=numeric()) #dataframe values to check
  
  
  for(i in 1:nrow(coords)) {
    
    CV[i,1]<-FR[coords[i,1],coords[i,2]]
    CV[i,2]<-SR[coords[i,1],coords[i,2]]
    CV[i,3]<-TR[coords[i,1],coords[i,2]]
    
  }
  
  CV<-cbind(coords,CV) # visual identification of outlayer
  
  print(nrow(CV))
  
  CV<-select_outlier(CV)
  
  Fdf2<-Fdf[,c(6:101)]
  
  for (i in 1:nrow(Fdf2)) {
    
    Fdf2[CV[i,1],CV[i,9]]<-NA
    
  }
  
  Fdf3<-cbind(Fdf[,c(1:5)], Fdf2)
  
  return(Fdf3)
  
}

select_outlier<-function(CV){
  
  CV$d1<-rowMeans(CV[,c(3:5)])-CV$R1
  CV$d2<-rowMeans(CV[,c(3:5)])-CV$R2
  CV$d3<-rowMeans(CV[,c(3:5)])-CV$R3
  
  for (i in 1:nrow(CV)) {
    
    if (max(abs(CV[i,c(6:8)]))==abs(CV[i,6])) {CV[i,"X2F"]<-CV[i,"X2"]} 
    if (max(abs(CV[i,c(6:8)]))==abs(CV[i,7])) {CV[i,"X2F"]<-CV[i,"X2"]+32}
    if (max(abs(CV[i,c(6:8)]))==abs(CV[i,8])) {CV[i,"X2F"]<-CV[i,"X2"]+64}

  }
  return(CV)
  }

T1<-load_and_clean(name="T1")
T2<-load_and_clean(name="T2")
T3<-load_and_clean(name="T3")
T4<-load_and_clean(name="T4")

# Mean consumption and index ---------------------------------------------

# estimation of the mean value per plate and compound

estimate_mean_abs<-function (df) {
  
  FR<-df[,c(6:37)]#First replicate
  SR<-df[,c(38:69)]#second replicate
  TR<-df[,c(70:101)]#third replicate
  
  all.dat <- list(dat1=FR, dat2=SR, dat3=TR)
  plate_m = aaply(laply(all.dat, as.matrix), c(2, 3), mean, na.rm=TRUE)
  Fdf<-cbind(df[,c(1:5)],plate_m)
  
  return (Fdf)
}

estimate_AWCD_H<- function (df) {
  
  for (i in 1:nrow(df)) {
    
    vec<-as.numeric(df[i,c(7:36)])
    A1<-df[i,6]
    df[i,"AWCD"]<- (sum(vec-A1)/31)
    df[i,"H"]<-diversity(df[i,c(7:36)], "shannon")}
  
  return (df)
}

T1<-estimate_AWCD_H(estimate_mean_abs(T1))
T2<-estimate_AWCD_H(estimate_mean_abs(T2))
T3<-estimate_AWCD_H(estimate_mean_abs(T3))
T4<-estimate_AWCD_H(estimate_mean_abs(T4))


# Comparison among replicates --------------------------------------------

mean_se_rep<-function (df){

  dfA<-subset(df, df$Station == "ATR" & df$Met == "AER" & df$Core == "C")
  dfB<-subset(df, df$Station =="BOL" & df$Met == "AER" & df$Core == "A")
  dfJ<-subset(df, df$Station =="J" & df$Met == "AER" & df$Core == "C")
  dfBB<-subset(df, df$Station =="BB" & df$Met == "AER" & df$Core == "A")
  
  dfAm<-subset(df, df$Station == "ATR" & df$Met == "MOX" & df$Core == "C")
  dfBm<-subset(df, df$Station =="BOL" & df$Met == "MOX" & df$Core == "A")
  dfJm<-subset(df, df$Station =="J" & df$Met == "MOX" & df$Core == "C")
  dfBBm<-subset(df, df$Station =="BB" & df$Met == "MOX" & df$Core == "A")
  
  df2<-as.data.frame(rbind(dfA,dfB,dfJ,dfBB,dfAm,dfBm,dfJm,dfBBm))
  df_se<-df2
  df_se[,c(6:39)]<-NA
  
  ATR4<-subset(df,df$Depth == 4 & df$Met == "AER" & df$Station == "ATR")
  df2[which(df2$Station == "ATR" & df2$Met == "AER" & df2$Depth == 4),c(6:39)]<-sapply(ATR4[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "ATR" & df_se$Met == "AER" & df_se$Depth == 4),c(6:39)]<-sapply(ATR4[,c(6:39)],std.error)
  
  ATR4m<-subset(df,df$Depth == 4 & df$Met == "MOX" & df$Station == "ATR")
  df2[which(df2$Station == "ATR" & df2$Met == "MOX" & df2$Depth == 4),c(6:39)]<-sapply(ATR4m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "ATR" & df_se$Met == "MOX" & df_se$Depth == 4),c(6:39)]<-sapply(ATR4m[,c(6:39)],std.error)
  
  ATR42<-subset(df,df$Depth == 42 & df$Met == "AER" & df$Station == "ATR")
  df2[which(df2$Station == "ATR" & df2$Met == "AER" & df2$Depth == 42),c(6:39)]<-sapply(ATR42[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "ATR" & df_se$Met == "AER" & df_se$Depth == 42),c(6:39)]<-sapply(ATR42[,c(6:39)],std.error)
  
  ATR42m<-subset(df,df$Depth == 42 & df$Met == "MOX" & df$Station == "ATR")
  df2[which(df2$Station == "ATR" & df2$Met == "MOX" & df2$Depth == 42),c(6:39)]<-sapply(ATR42m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "ATR" & df_se$Met == "MOX" & df_se$Depth == 42),c(6:39)]<-sapply(ATR42m[,c(6:39)],std.error)
  
  BB4<-subset(df,df$Depth == 4 & df$Met == "AER" & df$Station == "BB")
  df2[which(df2$Station == "BB" & df2$Met == "AER" & df2$Depth == 4),c(6:39)]<-sapply(BB4[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BB" & df_se$Met == "AER" & df_se$Depth == 4),c(6:39)]<-sapply(BB4[,c(6:39)],std.error)
  
  BB4m<-subset(df,df$Depth == 4 & df$Met == "MOX" & df$Station == "BB")
  df2[which(df2$Station == "BB" & df2$Met == "MOX" & df2$Depth == 4),c(6:39)]<-sapply(BB4m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BB" & df_se$Met == "MOX" & df_se$Depth == 4),c(6:39)]<-sapply(BB4m[,c(6:39)],std.error)
  
  BB42<-subset(df,df$Depth == 42 & df$Met == "AER" & df$Station == "BB")
  df2[which(df2$Station == "BB" & df2$Met == "AER" & df2$Depth == 42),c(6:39)]<-sapply(BB42[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BB" & df_se$Met == "AER" & df_se$Depth == 42),c(6:39)]<-sapply(BB42[,c(6:39)],std.error)
  
  BB42m<-subset(df,df$Depth == 42 & df$Met == "MOX" & df$Station == "BB")
  df2[which(df2$Station == "BB" & df2$Met == "MOX" & df2$Depth == 42),c(6:39)]<-sapply(BB42m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BB" & df_se$Met == "MOX" & df_se$Depth == 42),c(6:39)]<-sapply(BB42m[,c(6:39)],std.error)
  
  BOL4<-subset(df,df$Depth == 4 & df$Met == "AER" & df$Station == "BOL")
  df2[which(df2$Station == "BOL" & df2$Met == "AER" & df2$Depth == 4),c(6:39)]<-sapply(BOL4[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BOL" & df_se$Met == "AER" & df_se$Depth == 4),c(6:39)]<-sapply(BOL4[,c(6:39)],std.error)
  
  BOL4m<-subset(df,df$Depth == 4 & df$Met == "MOX" & df$Station == "BOL")
  df2[which(df2$Station == "BOL" & df2$Met == "MOX" & df2$Depth == 4),c(6:39)]<-sapply(BOL4m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BOL" & df_se$Met == "MOX" & df_se$Depth == 4),c(6:39)]<-sapply(BOL4m[,c(6:39)],std.error)
  
  BOL42<-subset(df,df$Depth == 42 & df$Met == "AER" & df$Station == "BOL")
  df2[which(df2$Station == "BOL" & df2$Met == "AER" & df2$Depth == 42),c(6:39)]<-sapply(BOL42[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BOL" & df_se$Met == "AER" & df_se$Depth == 42),c(6:39)]<-sapply(BOL42[,c(6:39)],std.error)
  
  BOL42m<-subset(df,df$Depth == 42 & df$Met == "MOX" & df$Station == "BOL")
  df2[which(df2$Station == "BOL" & df2$Met == "MOX" & df2$Depth == 42),c(6:39)]<-sapply(BOL42m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "BOL" & df_se$Met == "MOX" & df_se$Depth == 42),c(6:39)]<-sapply(BOL42m[,c(6:39)],std.error)
  
  J4<-subset(df,df$Depth == 4 & df$Met == "AER" & df$Station == "J")
  df2[which(df2$Station == "J" & df2$Met == "AER" & df2$Depth == 4),c(6:39)]<-sapply(J4[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "J" & df_se$Met == "AER" & df_se$Depth == 4),c(6:39)]<-sapply(J4[,c(6:39)],std.error)
  
  J4m<-subset(df,df$Depth == 4 & df$Met == "MOX" & df$Station == "J")
  df2[which(df2$Station == "J" & df2$Met == "MOX" & df2$Depth == 4),c(6:39)]<-sapply(J4m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "J" & df_se$Met == "MOX" & df_se$Depth == 4),c(6:39)]<-sapply(J4m[,c(6:39)],std.error)
  
  J42<-subset(df,df$Depth == 42 & df$Met == "AER" & df$Station == "J")
  df2[which(df2$Station == "J" & df2$Met == "AER" & df2$Depth == 42),c(6:39)]<-sapply(J42[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "J" & df_se$Met == "AER" & df_se$Depth == 42),c(6:39)]<-sapply(J42[,c(6:39)],std.error)
  
  J42m<-subset(df,df$Depth == 42 & df$Met == "MOX" & df$Station == "J")
  df2[which(df2$Station == "J" & df2$Met == "MOX" & df2$Depth == 42),c(6:39)]<-sapply(J42m[,c(6:39)],mean, na.rm=TRUE)
  df_se[which(df_se$Station == "J" & df_se$Met == "MOX" & df_se$Depth == 42),c(6:39)]<-sapply(J42m[,c(6:39)],std.error)
  
  results<-list(df2, df_se)
  
  return(results)
    
}

std.error <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

T1F<-mean_se_rep(T1)[[1]]
T1F_se<-mean_se_rep(T1)[[2]]
T2F<-mean_se_rep(T2)[[1]]
T2F_se<-mean_se_rep(T2)[[2]]
T3F<-mean_se_rep(T3)[[1]]
T3F_se<-mean_se_rep(T3)[[2]]
T4F<-mean_se_rep(T4)[[1]]
T4F_se<-mean_se_rep(T4)[[2]]



# Temporal differences -----------------------------------------------------


plot(T1F$AWCD, T4F$AWCD)

#temp<-cbind(T1F[,c(1:5, 38,39)], T2F[,c(38,39)], T3F[,c(38,39)], T4F[,c(38,39)])
#colnames(temp)<-c("sample", "station", "core", "depth", "met", "T1_AWCH", "T1_H",
#                  "T2_AWCH", "T2_H", "T3_AWCH", "T3_H", "T4_AWCH", "T4_H")

T1F$time<-62
T2F$time<-110
T3F$time<-158
T4F$time<-234

temp<-rbind(T1F[,c(1:5, 38:40)], T2F[,c(1:5,38:40)], T3F[,c(1:5,38:40)], T4F[,c(1:5,38:40)])

pairwise.wilcox.test(temp$AWCD, temp$time, # sig diff <0.05
                     p.adjust.method = "BH")

pairwise.wilcox.test(temp$H, temp$time, # sig diff <0.05
                     p.adjust.method = "BH")


ggplot(T1F, aes(T1F$AWCD, T4F$AWCD))+
  geom_point(aes(color=Station, shape=Met))+
  xlim(0.5, 1.8) + ylim(0.5,1.8)


# Average compound consumption --------------------------------------------


ATR_AERmean<-as.data.frame(colMeans(subset(T4F, Station=="ATR" & Met == "AER")[,c(7:37)]))
ATR_AERse<-as.data.frame(sapply(subset(T4F, Station=="ATR" & Met == "AER")[,c(7:37)],std.error))
ATR_MOXmean<-as.data.frame(colMeans(subset(T4F, Station=="ATR" & Met == "MOX")[,c(7:37)]))
ATR_MOXse<-as.data.frame(sapply(subset(T4F, Station=="ATR" & Met == "MOX")[,c(7:37)],std.error))

BOL_AERmean<-as.data.frame(colMeans(subset(T4F, Station=="BOL" & Met == "AER")[,c(7:37)]))
BOL_AERse<-as.data.frame(sapply(subset(T4F, Station=="BOL" & Met == "AER")[,c(7:37)],std.error))
BOL_MOXmean<-as.data.frame(colMeans(subset(T4F, Station=="BOL" & Met == "MOX")[,c(7:37)]))
BOL_MOXse<-as.data.frame(sapply(subset(T4F, Station=="BOL" & Met == "MOX")[,c(7:37)],std.error))

J_AERmean<-as.data.frame(colMeans(subset(T4F, Station=="J" & Met == "AER")[,c(7:37)]))
J_AERse<-as.data.frame(sapply(subset(T4F, Station=="J" & Met == "AER")[,c(7:37)],std.error))
J_MOXmean<-as.data.frame(colMeans(subset(T4F, Station=="J" & Met == "MOX")[,c(7:37)]))
J_MOXse<-as.data.frame(sapply(subset(T4F, Station=="J" & Met == "MOX")[,c(7:37)],std.error))

BB_AERmean<-as.data.frame(colMeans(subset(T4F, Station=="BB" & Met == "AER")[,c(7:37)]))
BB_AERse<-as.data.frame(sapply(subset(T4F, Station=="BB" & Met == "AER")[,c(7:37)],std.error))
BB_MOXmean<-as.data.frame(colMeans(subset(T4F, Station=="BB" & Met == "MOX")[,c(7:37)]))
BB_MOXse<-as.data.frame(sapply(subset(T4F, Station=="BB" & Met == "MOX")[,c(7:37)],std.error))

station_AER_mean<-cbind(ATR_AERmean, BOL_AERmean, J_AERmean, BB_AERmean)
colnames(station_AER_mean)<-c("ATR", "BOL", "J", "BB")
station_AER_mean$Compound<-rownames(station_AER_mean)
station_AER_se<-cbind(ATR_AERse, BOL_AERse, J_AERse, BB_AERse)
colnames(station_AER_se)<-c("ATR", "BOL", "J", "BB")
station_AER_se$Compound<-rownames(station_AER_se)
station_AER_mean$Guilt<-NA
station_AER_mean$Met<-"AER"
station_AER_se$Met<-"AER"

station_AER_mean$Guilt[station_AER_mean$Compound == "A1"]<-"Control"
station_AER_mean$Guilt[station_AER_mean$Compound == "B1"|station_AER_mean$Compound == "G2"|station_AER_mean$Compound == "H2"] <- "Misc. "
station_AER_mean$Guilt[station_AER_mean$Compound == "C1"|station_AER_mean$Compound == "D1"|station_AER_mean$Compound == "F1"|station_AER_mean$Compound == "E1"] <- "Polymers"
station_AER_mean$Guilt[station_AER_mean$Compound == "G1"|station_AER_mean$Compound =="H1"|station_AER_mean$Compound == "A2"| station_AER_mean$Compound =="B2"| station_AER_mean$Compound =="C2"| station_AER_mean$Compound =="D2"| station_AER_mean$Compound =="E2"] <- "Carbohydrates"
station_AER_mean$Guilt[station_AER_mean$Compound == "F2"|station_AER_mean$Compound =="A3"|station_AER_mean$Compound == "B3"| station_AER_mean$Compound =="C3"| station_AER_mean$Compound =="D3"| station_AER_mean$Compound =="E3"| station_AER_mean$Compound =="F3"| station_AER_mean$Compound =="G3"| station_AER_mean$Compound =="H3"] <- "Carboxylic acids"
station_AER_mean$Guilt[station_AER_mean$Compound == "A4"|station_AER_mean$Compound =="B4"|station_AER_mean$Compound == "C4"| station_AER_mean$Compound =="D4"| station_AER_mean$Compound =="E4"| station_AER_mean$Compound =="F4"] <- "Amino acids"
station_AER_mean$Guilt[station_AER_mean$Compound == "G4"|station_AER_mean$Compound == "H4"] <- "Amines"

station_MOX_mean<-cbind(ATR_MOXmean, BOL_MOXmean, J_MOXmean, BB_MOXmean)
colnames(station_MOX_mean)<-c("ATR", "BOL", "J", "BB")
station_MOX_mean$Compound<-rownames(station_MOX_mean)
station_MOX_se<-cbind(ATR_MOXse, BOL_MOXse, J_MOXse, BB_MOXse)
colnames(station_MOX_se)<-c("ATR", "BOL", "J", "BB")
station_MOX_se$Compound<-rownames(station_MOX_se)
station_MOX_mean$Guilt<-NA
station_MOX_mean$Met<-"MOX"
station_MOX_se$Met<-"MOX"

station_MOX_mean$Guilt[station_MOX_mean$Compound == "A1"]<-"Control"
station_MOX_mean$Guilt[station_MOX_mean$Compound == "B1"|station_MOX_mean$Compound == "G2"|station_MOX_mean$Compound == "H2"] <- "Misc. "
station_MOX_mean$Guilt[station_MOX_mean$Compound == "C1"|station_MOX_mean$Compound == "D1"|station_MOX_mean$Compound == "F1"|station_MOX_mean$Compound == "E1"] <- "Polymers"
station_MOX_mean$Guilt[station_MOX_mean$Compound == "G1"|station_MOX_mean$Compound =="H1"|station_MOX_mean$Compound == "A2"| station_MOX_mean$Compound =="B2"| station_MOX_mean$Compound =="C2"| station_MOX_mean$Compound =="D2"| station_MOX_mean$Compound =="E2"] <- "Carbohydrates"
station_MOX_mean$Guilt[station_MOX_mean$Compound == "F2"|station_MOX_mean$Compound =="A3"|station_MOX_mean$Compound == "B3"| station_MOX_mean$Compound =="C3"| station_MOX_mean$Compound =="D3"| station_MOX_mean$Compound =="E3"| station_MOX_mean$Compound =="F3"| station_MOX_mean$Compound =="G3"| station_MOX_mean$Compound =="H3"] <- "Carboxylic acids"
station_MOX_mean$Guilt[station_MOX_mean$Compound == "A4"|station_MOX_mean$Compound =="B4"|station_MOX_mean$Compound == "C4"| station_MOX_mean$Compound =="D4"| station_MOX_mean$Compound =="E4"| station_MOX_mean$Compound =="F4"] <- "Amino acids"
station_MOX_mean$Guilt[station_MOX_mean$Compound == "G4"|station_MOX_mean$Compound == "H4"] <- "Amines"


station_mean<-rbind(station_AER_mean, station_MOX_mean)
station_se<-rbind(station_AER_se, station_MOX_se)

mstation_mean<-melt(station_mean, ID= c("Guilt", "Compound", "Met"))
colnames(mstation_mean)<-c("Compound","Guilt",  "Met", "Station", "Mean")
mstation_se<-melt(station_se, ID= c("Compound",  "Met"))
colnames(mstation_se)<-c("Compound","Met", "Station", "SE")
mstation<-cbind(mstation_mean, mstation_se[, "SE"])
names(mstation)[names(mstation) == 'mstation_se[, "SE"]'] <- 'SE'

plot_consumption<-function(guilt=NULL){
  
  p<-ggplot(subset(mstation, Guilt== guilt), aes( Mean, Station))+ ggtitle(guilt)+
    geom_point(aes(color=Compound, shape=Met))+
    geom_errorbarh(aes(xmax = Mean + SE, xmin = Mean - SE, color=Compound, height = 0))+
    labs(color=NULL)+
    xlim(-0.1, 2.5)+
    geom_vline(xintercept = mean(as.numeric(T4F[,"A1"])))+
    geom_vline(xintercept = mean(as.numeric(T4F[,"A1"]))-sd(as.numeric(T4F[,"A1"])),linetype = "dashed")+
    geom_vline(xintercept = mean(as.numeric(T4F[,"A1"]))+sd(as.numeric(T4F[,"A1"])),linetype = "dashed")+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}

Amines<-plot_consumption(guilt="Amines")
AA<-plot_consumption(guilt="Amino acids")
Ch<-plot_consumption(guilt="Carbohydrates")
CA<-plot_consumption(guilt="Carboxylic acids")
P<-plot_consumption(guilt="Polymers")
M<-plot_consumption(guilt="Misc. ")

CsC<-grid.arrange(Amines, AA, Ch, CA, M, P, nrow=3 ,ncol = 2, left = "Station", bottom= "Average Color Development at T4")

ggsave("results/Carbon sources comp.jpg",CsC, units="cm", width = 19, height = 19)


#### heat map ####

#AER 

temp<-melt(station_AER_mean[,c(1:6)], id=c('Compound',"Guilt"))

A1_lim<-mean(as.numeric(T4F[,"A1"]))+sd(as.numeric(T4F[,"A1"]))

#temp$value <- replace(temp$value, temp$value < A1_lim, 0) 

temp$value<-temp$value-A1_lim
temp$value <- replace(temp$value, temp$value < 0, 0) 


p1<-ggplot(temp, aes(Compound, variable, fill = value))+
  ylab("Aerobic Consumption")+
  #xlab("Average Color Development at T4 minus the control")+
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red")+
  facet_grid(~ Guilt, scales = "free", space = "free_x")

p1

# MOX

temp<-melt(station_MOX_mean[,c(1:6)], id=c('Compound',"Guilt"))

A1_lim<-mean(as.numeric(T4F[,"A1"]))+sd(as.numeric(T4F[,"A1"]))

temp$value<-temp$value-A1_lim
temp$value <- replace(temp$value, temp$value < 0, 0) 


p2<-ggplot(temp, aes(Compound, variable, fill = value))+
  ylab("Microoxic Consumption")+
  xlab("Average Color Development at T4 minus the control (A1)")+
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red")+
  facet_grid(~ Guilt, scales = "free", space = "free_x")+


p2

heat_map<-grid.arrange(p1, p2, ncol = 1)
ggsave("results/ecoplates_heat_map.jpg",heat_map, units="cm", width = 19, height = 15)

# Comparison among stations ----------------------------------------------


box_station<-function(variable, y_title) {

AER<-subset(T4F, Met == "AER")
AER_se<-subset(T4F_se, Met == "AER")
MOX<-subset(T4F, Met == "MOX")
MOX_se<-subset(T4F_se, Met == "MOX")

p1<-
  ggplot(AER, aes(x=factor(Station), AER[,variable]))+ ylab(y_title)+ ggtitle("Aerobic")+
  geom_boxplot(aes(x=factor(Station), AER[,variable]),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station))+
  #scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=AER[,variable]-AER_se[,variable], ymax=AER[,variable]+AER_se[,variable], color=Station),
    width = 0.1
    #,linetype = "dotted"
    #,position=position_dodge(width=0.5)
    )

p2<-
  ggplot(MOX, aes(x=factor(Station), MOX[,variable]))+ ylab(y_title)+  ggtitle("Microoxic")+
  geom_boxplot(aes(x=factor(Station), MOX[,variable]),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station))+
  #scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=MOX[,variable]-MOX_se[,variable], ymax=MOX[,variable]+MOX_se[,variable], color=Station),
    width = 0.1
    #,linetype = "dotted"
    #,position=position_dodge(width=0.5)
  )
  
p<-grid.arrange(p1, p2, ncol=2)


return(p)
}

box_station(variable="AWCD", y_title = "AWCD")
box_station(variable="H", y_title = "H")
box_station(variable="A1" , y_title = "A1")

pairwise.wilcox.test(T4F$G4, T4F$Station,#### are significantly different (p < 0.05)
                     p.adjust.method = "BH")


# create extra columnd with avergae consumtion per guilt
T4F_g<-T4F
T4F_g$Amines<-(T4F$H4+T4F$G4)/2
T4F_g$Miscellaneous<-(T4F$B1+T4F$G2+T4F$H2)/3
T4F_g$Polymers<-(T4F$C1+T4F$D1+T4F$F1+T4F$E1)/4
T4F_g$Carbohydrates<-(T4F$G1+T4F$H1+T4F$A2+T4F$B2+T4F$C2+T4F$D2+T4F$E2)/7
T4F_g$CarboxylicAcids<-(T4F$F2+T4F$A3+T4F$B3+T4F$C3+T4F$D3+T4F$E3+T4F$F3+T4F$G3+T4F$H3)/9
T4F_g$AminoAcids<-(T4F$A4+T4F$B4+T4F$C4+T4F$D4+T4F$E4+T4F$F4)/6


#subset per incubation condition
AER<-subset(T4F_g, Met == "AER")
pairwise.wilcox.test(AER$G4, AER$Station,#### are significantly different (p < 0.05)
                     p.adjust.method = "BH")

MOX<-subset(T4F_g, Met == "MOX")


diff_aer<-function (var) {
  
  test<-pairwise.wilcox.test(var, AER$Station,#### are significantly different (p < 0.05)
                       p.adjust.method = "BH")
   if(any(as.vector(test$p.value)<0.05, na.rm = TRUE)==TRUE) {print(test)}
}
diff_mox<-function (var) {
  
  test<-pairwise.wilcox.test(var, MOX$Station,#### are significantly different (p < 0.05)
                             p.adjust.method = "BH")
  if(any(as.vector(test$p.value)<0.05, na.rm = TRUE)==TRUE) {print(test)}
}




lapply(T4F_g[,c(6:45)], diff_aer)
lapply(T4F_g[,c(6:45)], diff_mox)

box_station(variable="D3", y_title= "Average Well Color Development")

p1<-box_station(variable="AWCD", y_title= "Average Well Color Development")
p2<-box_station(variable="H", y_title= "Shannon-Weber index")
p3<-box_station(variable="G4", y_title= "Phenil ethyl amine")
p4<-box_station(variable="E2", y_title= "N-Acetyl-Glucosamine")

sig_dif_m<-grid.arrange(p1, p2, ncol = 1)
sig_dif<-grid.arrange(p1, p2, p3, p4, ncol = 2)

ggsave("results/site_dif_m.jpg",sig_dif_m, units="cm", width = 25, height = 7)


# Comparison between metabolism ------------------------------------------

C<- T4F[-c(13,14,42),]
C_se<- T4F_se[-c(13,14,42),] #data.frame with only samples where we have data for both AER and ANA 

#AWCD

AWCD.AER<-C$AWCD[C$Met == "AER"]
AWCD.MOX<-C$AWCD[C$Met == "MOX"]

shapiro.test(AWCD.AER) ### normality (>0.05 normal, <0.05 no normal)
shapiro.test(AWCD.MOX) ### normality (>0.05 normal, <0.05 no normal) 

#F-test requires normality

var.test(AWCD ~ Met, data = C) #>0.05 no significant difference between the two variances

t.test(AWCD.AER, AWCD.MOX) # p <0.05 significant differences

Bxp.AWCD<-ggplot(C,aes(x=factor(Met), AWCD)) + ggtitle("Average Well Color Development (AWCD)")+
  geom_boxplot(aes(x=factor(Met), AWCD),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station),position=position_dodge(width=0.5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=AWCD-C_se$AWCD, ymax=AWCD+C_se$AWCD, color=Station),
    width = 0.1,
    #linetype = "dotted",
    position=position_dodge(width=0.5))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none", 
        plot.title = element_text(size = 11)
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed")
        )


station_AWCD<-ggplot(C,aes(x=factor(Met), AWCD)) + ylab("Average Well Color Development (AWCD)")+
  geom_boxplot(aes(x=factor(Met), AWCD),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station),position=position_dodge(width=0.5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=AWCD-C_se$AWCD, ymax=AWCD+C_se$AWCD, color=Station),
    width = 0.1,
    #linetype = "dotted",
    position=position_dodge(width=0.5))+
  facet_wrap(~Station, scale="fixed")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none"
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed"))

C$Depth<-as.numeric(C$Depth)

ggplot(C,aes(Depth,AWCD))+ ylab("Average Well Color Development (AWCD)")+
  geom_point(aes(color=Met))+
  coord_flip()+
  scale_x_reverse()+
  facet_wrap(~Station, scale="fixed")


####Shannon-Wiener Diversity

H.AER<-C$H[C$Met == "AER"]
H.MOX<-C$H[C$Met == "MOX"]

H.AER<-H.AER[-21]#delete outlier
H.MOX<-H.MOX[-21]#in the other met as well to compensate the groups

shapiro.test(H.AER) ### normality (>0.05 normal, <0.05 no normal) 

#F-test requires normality

var.test(H ~ Met, data = C[-c(39,42),]) #>0.05 no significant difference between the two variances
#C[-c(39,42),] because of the outlier
t.test(H.AER, H.MOX)


Bxp.H<-ggplot(C,aes(x=factor(Met), H)) + ggtitle("Shannon-Wiener diversity index (H')")+
  geom_boxplot(aes(x=factor(Met), H),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station),position=position_dodge(width=0.5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=H-C_se$H, ymax=H+C_se$H, color=Station),
    width = 0.1,
    #linetype = "dotted",
    position=position_dodge(width=0.5))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        plot.title = element_text(size = 11)
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed"))

station_H<-ggplot(C,aes(x=factor(Met), H)) + ylab("Shannon-Wiener diversity index (H')")+
  geom_boxplot(aes(x=factor(Met), H),fill="grey",alpha=0.5)+
  geom_point(aes(color=Station),position=position_dodge(width=0.5))+
  scale_colour_manual(values = c("blue","red","orange","black"))+
  geom_errorbar(
    aes(ymin=H-C_se$H, ymax=H+C_se$H, color=Station),
    width = 0.1,
    #linetype = "dotted",
    position=position_dodge(width=0.5))+
  facet_wrap(~Station, scale="fixed")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()
        ,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
        )
        ,panel.grid.major = element_line(colour="grey",linetype="dashed"))


Aer.vs.ana<-cowplot::plot_grid(Bxp.AWCD,Bxp.H, station_AWCD, station_H, align = "h", nrow = 2, rel_widths = c(5/6, 1),labels=c('A', 'B', "C", "D"))

ggsave("results/Aer.vs.ana.jpg",Aer.vs.ana, units="cm", width = 19, height = 20)


test_station_dif <- function (var, station) {
  
  data<-subset(C, Station==station)

  data.AER<- subset(data, Met == "AER")[,var]
  data.MOX<- subset(data, Met == "MOX")[,var]
  
  print(shapiro.test(data.AER)) 
  print(shapiro.test(data.MOX))  

  print(t.test(data.AER, data.MOX))
  
}

test_station_dif(var="AWCD", station="J")


# Consumption with depth --------------------------------------------------

ggplot(T4F_g, aes(Depth, A1))+
  geom_point(aes(color=Met))+
  coord_flip()+
  scale_x_reverse() +
  facet_wrap(~Station)
  

p1<-ggplot(T4F_g, aes(Depth, AWCD))+
  xlab("Depth (cm)")+
  ylab("Average Well Color Development")+
  geom_point(aes(color=Met))+
  geom_line(aes(color=Met))+
  coord_flip()+
  scale_x_reverse() +
  facet_wrap(~Station)+
  theme(,panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
  )
  ,panel.grid.major = element_line(colour="grey",linetype="dashed"))

p2<-ggplot(T4F_g, aes(Depth, H))+
  xlab("Depth (cm)")+
  ylab("Shannon-Wiener diversity index (H)")+
  geom_point(aes(color=Met))+
  geom_line(aes(color=Met))+
  coord_flip()+
  scale_x_reverse() +
  facet_wrap(~Station)+
  theme(axis.text.x = element_text(angle=45),
    panel.background = element_rect(fill="white", color = "black"
                                         #, colour, size, linetype, color
  )
  ,panel.grid.major = element_line(colour="grey",linetype="dashed"))


eco_depth<-grid.arrange(p1, p2,ncol = 2
                      #,layout_matrix = cbind(c(1,1), c(1,1),c(2,3))
                      )

ggsave("results/eco_depth.jpg",eco_depth, units="cm", width = 19, height = 15)


