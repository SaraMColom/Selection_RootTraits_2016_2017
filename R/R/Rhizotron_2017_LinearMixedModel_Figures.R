#######################################################
### R Code for publication in prep
### Title :  Belowground competition can influence the evolution of root traits 

        # Sara Colom and Regina Baucom #

#     2016 Rhizotron common garden experiment 
#   Linear mixed models for fixed and random effects
#         Least square means of root traits 
#######################################################
### *Note that this code does not show all preliminary#
### analysis                                        ###
#######################################################

#######################################################
### Getting started									###
#######################################################

#Load libraries
library(tidyr)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(PerformanceAnalytics)

source("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/R/SummarySE.R")

#Save theme text setting
Tx<-theme(axis.text.x = element_text(face="bold",  
                                     size=14),
          axis.text.y = element_text(face="bold", 
                                     size=14))

# Root architecture trait data from rhizotron experiment
totaldat<-read.csv("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/CleanData/Rhiz_root_traits.csv")
#Check structure and correct structure
str(totaldat)
totaldat$Id<-as.factor(totaldat$Id)
totaldat$Experiment<-as.factor(totaldat$Experiment)

# Count species
table(totaldat$Species)
# 196 Ihed; 125 Ip

# Count species table(totaldat$Species,totaldat$Experiment)
table(totaldat$Species,totaldat$Experiment)
# 102 Ihed in Exp 1 and 94  in Exp 2; 67 Ip in Exp 1 and 58 in Exp 2

#####################################################################
#################### Linear mixed model analysis ###############
#                                 &
#                 Least square means of root traits
#####################################################################

#Average root angle
mod1.angleAvAng<-lmer(AvAng~Species+Experiment+Population+
                          (1|ML),totaldat)
anova(mod1.angleAvAng) # F-statistics for fixed terms
ranova(mod1.angleAvAng) # Log-liklihood for random terms
emmeans(mod1.angleAvAng,"Species")  #Least square means by species

#Primary root length
mod1.primary<-lmer(PrimaryRootLength~Species+Population+Experiment+(1|ML),totaldat)
anova(mod1.primary) # F-statistics for fixed terms
ranova(mod1.primary) # Log-liklihood for random terms
emmeans(mod1.primary, "Species")#Least square means by species

#Root system width
mod1.RSW<-lmer(RootSystemWidth~Species+Population+Experiment+(1|ML),totaldat)
anova(mod1.RSW) # F-statistics for fixed terms
ranova(mod1.RSW) # Log-liklihood for random terms
emmeans(mod1.RSW, "Species") #Least square means by species 

#Total root surface area
mod1.TA<-lmer(Total.Area~Species+Population+Experiment+(1|ML),totaldat)
anova(mod1.TA) # F-statistics for fixed terms
ranova(mod1.TA) # Log-liklihood for random terms
emmeans(mod1.TA, "Species") #Least square means by species 


#Mixed models on on traits within Ihed and Ip
Ip1<-droplevels(subset(totaldat,Species=="Ip"))
Ihed1<-droplevels(subset(totaldat,Species=="Ihed"))

# First I.hederacea
mod1.angleAvH<-lmer(AvAng~Experiment+Population+
                      (1|ML),Ihed1)
anova(mod1.angleAvH) #F statistics for fixed effects
ranova(mod1.angleAvH) #LRT statistics for ranodm effects

mod1.primaryH<-lmer(PrimaryRootLength~Experiment+Population+(1|ML),Ihed1)
anova(mod1.primaryH) #F statistics for fixed effects
ranova(mod1.primaryH) #LRT statistics for random effects

mod1.RSW<-lmer(RootSystemWidth~Experiment+Population+(1|ML),Ihed1)
anova(mod1.RSW) #F statistics for fixed effects
ranova(mod1.RSW) #LRT statistics for random effects

mod1.TRSh<-lmer(Total.Area~Population+Experiment+(1|ML),Ihed1)
anova(mod1.TRSh)
ranova(mod1.TRSh)

### Now Ip
#Average root angle
mod1.angleAVP<-lmer(AvAng~Population+Experiment+
                      (1|ML),Ip1)
anova(mod1.angleAVP) #F statistics for fixed effects
ranova(mod1.angleAVP) #LRT statistics for random effects

#Primary root length
mod1.primaryP<-lmer(PrimaryRootLength~Population+Experiment+(1|ML),Ip1)
anova(mod1.primaryP) #F statistics for fixed effects
ranova(mod1.primaryP) #LRT statistics for random effects

mod1.RSWp<-lmer(RootSystemWidth~Population+Experiment+(1|ML),Ip1)
anova(mod1.RSWp) #F statistics for fixed effects
ranova(mod1.RSWp) #LRT statistics for random effects

mod1.TRSp<-lmer(Total.Area~Population+Experiment+(1|ML),Ip1)
anova(mod1.TRSp)
ranova(mod1.TRSp)

##############################################
#                  FIGURES 
###############################################

# Figure 2

allData<-na.omit(totaldat) #Merge

#PCA
var.sp<-allData[c("Species","Total.Area","AvAng","PrimaryRootLength","RootSystemWidth")]
res.pca<-PCA(var.sp,quali.sup=1,scale.unit=T,graph=T)#scale mean of zero and std of one
#Dim 1 versus Dim 2
plot.PCA(res.pca,axes=c(1,2),choix="ind",habillage=1,col.hab=c('#bdbdbd','black')) # 

p<-fviz_pca_ind(res.pca, label="none", habillage=var.sp$Species,
                addEllipses=TRUE, ellipse.level=0.95)+
  theme_bw()
p + scale_color_manual(values=c('#bdbdbd','black'))+
    scale_fill_manual(values=c('#bdbdbd','darkgrey'))+
    theme_classic()+Tx

# Bar graph of loading scores
Loads<-data.frame(res.pca$var$contrib)
Loads$Trait<-c("Root area","Root angle","primary root length","Root system width")
# Use position=position_dodge()
q<-ggplot(data=Loads, aes(x=Trait, y=Dim.1), fill="black") +
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()
q +Tx+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

r<-ggplot(data=Loads, aes(x=Trait, y=Dim.2), fill="black") +
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()

r +Tx+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

