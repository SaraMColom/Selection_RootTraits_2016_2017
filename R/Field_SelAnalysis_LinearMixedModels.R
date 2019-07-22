#######################################################
### R Code for publication in prep 
### Title :  Belowground competition can influence the evolution of root traits 

        
          # Sara Colom and Regina Baucom #

#     2017 Field experiment selection analysis &
#   Linear mixed models for fixed and random effects
#     Least square means of root and fitness traits 
#######################################################
### *Note that this code does not show all preliminary#
### analysis                                        ###
#######################################################

#######################################################
### Getting started									###
#######################################################

#Libraries loaded
library(dplyr)
library(tidyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggpubr)

source("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/R/SummarySE.R")

#Aesthetics
Tx<-theme(axis.text.x = element_text(face="bold",  
                                     size=15),#,color="black"),
          axis.text.y = element_text(face="bold", 
                                     size=15),#,color="black"),
          axis.title.x = element_text(size=20),
          axis.title.y = element_text(size=20))+
          theme(plot.title=element_text(face="italic",size=30,hjust=0.5),
          legend.text=element_text(size=15))

# Key:
# MaxWidh* = Root system width
# AreaConvexHull* = Root size
# AvAng* = Average root angle
# Of_exp* = Cohort; O*= cohort 1, R*= cohort 2
# ML* = Maternal line
# Combos* = species by maternal line by species maternal line competition pairing
# Important note: maternal lines are explicitly nested within population--they only occur within specific populations and species.


# Data prep 

### Load Data and manage data

      Fitness<-read.csv("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/CleanData/FitnessData.csv")
      BRT<-read.csv("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/CleanData/BRTdata.csv")
      
      ###### Examine data structure and adjust
      
      str(Fitness)
      str(BRT)
      
      #Adjust
      Fitness[c("Block","Position")]<-lapply(Fitness[c("Block","Position")],as.factor)
      BRT[c("Block","Position")]<-lapply(BRT[c("Block","Position")],as.factor)
      
# Count sample sizes
      
# Root phenotyping
table(BRT$Species) 
table(BRT$Species,BRT$Trt) 

# Fitness seed number
table(Fitness$Species)
table(Fitness$Species,Fitness$Of_exp)
table(Fitness$Species,Fitness$Of_exp,Fitness$Trt)

###################################################################
#################### Linear mixed model analysis ###############
#                                 &
#                 Least square means of root traits
###################################################################
      
## Both species and Block included
      
MW<-lmer(MaxWidth~Trt+Block+Species+(1|ML),BRT)
anova(MW) # Fixed effects
ranova(MW) # Random effects
emmeans(MW,~Species|Trt)# Calculate LS mean of Species by treatment
      
ACH<-lmer(AreaConvexHull~Trt+Block+Species+(1|ML),BRT)
anova(ACH) # Fixed effects
ranova(ACH) # Random effects
# Calculate LS mean of Species by treatment
emmeans(ACH,~Species|Trt)
      
ARA<-lmer(AvAng~Trt+Block+Species+(1|ML),BRT)
anova(ARA) # Fixed effects
ranova(ARA) # Random effects
# Calculate LS mean of Species by treatment
emmeans(ARA,~Species|Trt) 
      
SN<-lmer(Seed_Number~Trt+Block+Species+(1|ML),Fitness)
anova(SN) # Fixed effects
ranova(SN) # Random effects
# Calculate LS mean of Species by treatment
emmeans(SN,~Species|Trt)

## Within species and Block included

Ihed<-droplevels(subset(BRT,Species=="H"))
Ip<-droplevels(subset(BRT,Species=="P"))

anova(lmer(AvAng~Trt+Block+(1|ML),Ihed))
ranova(lmer(AvAng~Trt+Block+(1|ML),Ihed))

anova(lmer(MaxWidth~Trt+Block+(1|ML),Ihed))
ranova(lmer(MaxWidth~Trt+Block+(1|ML),Ihed))

anova(lmer(AreaConvexHull~Trt+Block+(1|ML),Ihed))
ranova(lmer(AreaConvexHull~Trt+Block+(1|ML),Ihed))

anova(lmer(AvAng~Trt+Block+(1|ML),Ip))
ranova(lmer(AvAng~Trt+Block+(1|ML),Ip))

anova(lmer(MaxWidth~Trt+Block+(1|ML),Ip))
ranova(lmer(MaxWidth~Trt+Block+(1|ML),Ip))

anova(lmer(AreaConvexHull~Trt+Block+(1|ML),Ip))
ranova(lmer(AreaConvexHull~Trt+Block+(1|ML),Ip))

## post-hoc 

# Tukey's Test

summary(glht(ACH, linfct = mcp(Trt = "Tukey")), test = adjusted("holm"))

# Least square means
difflsmeans(ACH) # Across species; Intra v Alone is marginally signif and Intra v Inter is signif
# Within species
difflsmeans(lmer(AreaConvexHull~Trt+Block+(1|ML),Ihed)) # Evidence for trt effect
difflsmeans(lmer(AreaConvexHull~Trt+Block+(1|ML),Ip)) # No Evidence for trt effect

# Evidence indicates that treatment effect is within I.hederacea when we compare inter to intra

#######################################################################################
#################### CALCULATING RELATIVE FITNESS & FAMILY MEAN BIOMASS ###############
#######################################################################################

# Count number of samples per group for fitness data 
SampleSize<-Fitness %>% group_by(Species,Trt,Combos,Of_exp) %>% count(ML)

## ## ## ## ## ##   ## ## ## ## ## ##  Calculate relative fitness ## ## ## ## ## ## ## ## ## ## ## ##

#Calculate mean seed number  by Species and Cohort these will be analyzed separatley later!
SdMnSpeciesCohort<- aggregate(Seed_Number ~ Trt+Species+Of_exp, Fitness, mean)

names(SdMnSpeciesCohort)[4]<- "MeanSdNm" #Rename coloumn for mean seed number

# Merge average fitness of species and cohort 
Fitmean<-merge(Fitness,SdMnSpeciesCohort,by=c("Of_exp","Species","Trt"))

# Calculate relative fitness as the mean family seed number by cohort and treatment and the total mean seed number of that species,treatment and cohort
Fitmean$Rel_Fit<-Fitmean$Seed_Number/Fitmean$MeanSdNm

######################################################################################
                    ### Calculate family means of root traits
#######################################################################################

# The following for loop calculates the means of root traits by ML x Trt x Species 
trait<-which(names(BRT)%in%c("AreaConvexHull", "MaxWidth", "AvAng"))

# Aggregate 
for(i in trait) {
  if(i==trait[1]){
    df5<-aggregate(BRT[,i] ~ML+Trt+Species, BRT, mean)
    colnames(df5)[4]<-paste(colnames(BRT)[i])
  }
  else{
    df10<-aggregate(BRT[,i] ~ML+Trt+Species, BRT, mean) 
    colnames(df10)[4]<-paste(colnames(BRT)[i])
    df5<-merge(df10,df5,by=c("ML","Trt","Species"))
  }
}

# Merge family root trait means to fitness data
df5a<-unique(merge(Fitmean,df5,by=c("ML","Trt","Species")))

# Standardize root trait values
df5a$StdAreaConvexHull<-(df5a$AreaConvexHull-mean(df5a$AreaConvexHull))/sd(df5a$AreaConvexHull)
df5a$StdMaxWidth<-(df5a$MaxWidth-mean(df5a$MaxWidth))/sd(df5a$MaxWidth)
df5a$StdAvAng<-(df5a$AvAng-mean(df5a$AvAng))/sd(df5a$AvAng)

# Square values
df5a$StdAreaConvexHull2<-(df5a$StdAreaConvexHull*df5a$StdAreaConvexHull)
df5a$StdMaxWidth2<-(df5a$StdMaxWidth*df5a$StdMaxWidth)
df5a$StdAvAng2<-(df5a$StdAvAng*df5a$StdAvAng)

######################################################################################
                                ### Selection analysis
#######################################################################################

### Instruction for selection analysis
# Next we took averages of standardized fitness and all traits according to maternal line, treatment, species, and cohort. Have to do this for genotypic selection analysis!
# Made this data with and without block, ie average within and across blocks

# Analysis decisions: 
# 1) ANOVA for selection, using average across blocks, so use Seln datasheet. 
# 2) ANCOVA using competition treatment and blocks in model, so Seln2 datasheet. 
# 3) Present selection analyses results with cohort combined. Although fitness very different between cohorts (ie lower in cohort 2), pattern of selection on root traits did not differ according to cohort. Thus our preliminary examination suggested pattern of selection similar regardless of cohort. Note: first relativized fitness within cohort before testing.

## Average across block
Index<-which(names(df5a) %in% c("Seed_Number","AG_Biomass","MeanSdNm","Rel_Fit","AvAng",
                                "MaxWidth","AreaConvexHull","StdAreaConvexHull","StdMaxWidth","StdAvAng",
                                "StdAreaConvexHull2","StdMaxWidth2","StdAvAng2"))
Seln<-aggregate(list(df5a[Index]), by = list(df5a$ML,df5a$Species,df5a$Trt,df5a$Of_exp), mean) #Average everything by maternal line, species, treatment and cohort
colnames(Seln)[1:4]<-c("ML","Species","Trt","Of_exp")

## Average within blocks
Seln2<-aggregate(list(df5a[Index]), by = list(df5a$Block,df5a$ML,df5a$Species,df5a$Trt,df5a$Of_exp), mean) #Average everything by block*, maternal line, species, treatment and cohort
colnames(Seln2)[1:5]<-c("Block","ML","Species","Trt","Of_exp")

# Subset data by species

purpurea <- subset(Seln, Species=="P")
hederacea <- subset(Seln, Species=="H")
purpurea2 <- subset(Seln2, Species=="P")
hederacea2 <- subset(Seln2, Species=="H")

### ** Table B3 ** Selection gradients only per species within each treatment*

purpureaAlone <- subset(purpurea, Trt=="Alone")
purpureaInter <- subset(purpurea, Trt=="Inter")

hederaceaAlone <- subset(hederacea, Trt=="Alone")
hederaceaInter <- subset(hederacea, Trt=="Inter")
hederaceaIntra <- subset(hederacea, Trt=="Intra")

#Ipurp
summary(lm(Rel_Fit ~  StdMaxWidth + StdAvAng +StdAreaConvexHull, data= purpureaAlone))
summary(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull , data= purpureaInter))

#Ihed
summary(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull  , data= hederaceaAlone))
summary(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull  , data= hederaceaInter))
summary(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull  , data= hederaceaIntra))


## **Table B4 ** F-statistics: Selection gradients using ANCOVA to test Trt x Trait * 
#                         **F statistics reported in Table 4**

# Selection gradients using ANCOVA to test trt:trait interactions per species, block:trait etc 
# While all effects are in model, we present a summary of the analysis in text
# For I.hederacea we did Alone and Inter, and Alone and Inter

hederacea2.inter<-subset(hederacea2, Trt!= "Intra")
hederacea2.intra<-subset(hederacea2, Trt!= "Inter")

anova(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull + Trt + Block + Block:Trt + Trt:StdMaxWidth +  Trt:StdAvAng + Trt:StdAreaConvexHull + Block:StdMaxWidth + Block:StdAvAng + Block:StdAreaConvexHull + Block:Trt:StdMaxWidth + Block:Trt:StdAvAng + Block:Trt:StdAreaConvexHull, data= purpurea2)) # Inter v Alone (Ip)
anova(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull + Trt + Block + Block:Trt + Trt:StdMaxWidth +  Trt:StdAvAng + Trt:StdAreaConvexHull + Block:StdMaxWidth + Block:StdAvAng + Block:StdAreaConvexHull + Block:Trt:StdMaxWidth + Block:Trt:StdAvAng + Block:Trt:StdAreaConvexHull, data= hederacea2.inter)) # Inter v Alone
anova(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull + Trt + Block + Block:Trt + Trt:StdMaxWidth +  Trt:StdAvAng + Trt:StdAreaConvexHull + Block:StdMaxWidth + Block:StdAvAng + Block:StdAreaConvexHull + Block:Trt:StdMaxWidth + Block:Trt:StdAvAng + Block:Trt:StdAreaConvexHull, data= hederacea2.intra)) # Intra v Alone
anova(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull + Trt + Block + Block:Trt + Trt:StdMaxWidth +  Trt:StdAvAng + Trt:StdAreaConvexHull + Block:StdMaxWidth + Block:StdAvAng + Block:StdAreaConvexHull + Block:Trt:StdMaxWidth + Block:Trt:StdAvAng + Block:Trt:StdAreaConvexHull, data= hederacea2%>%filter(Trt!="Alone"))) # Inter v Intra
anova(lm(Rel_Fit ~  StdMaxWidth + StdAvAng + StdAreaConvexHull + Trt + Block + Block:Trt + Trt:StdMaxWidth +  Trt:StdAvAng + Trt:StdAreaConvexHull + Block:StdMaxWidth + Block:StdAvAng + Block:StdAreaConvexHull + Block:Trt:StdMaxWidth + Block:Trt:StdAvAng + Block:Trt:StdAreaConvexHull, data= hederacea2)) # All

#######################################################
              ### Correlation analysis 
#######################################################

#Plotting Correlation Matrix
subvars_purp <- names(purpurea) %in% c("Rel_Fit", "AvAng","MaxWidth", "AreaConvexHull")
sub_purpurea <- purpurea[subvars_purp]

subvars_hed <- names(hederacea) %in% c("Rel_Fit", "AvAng","MaxWidth", "AreaConvexHull")
sub_hederacea <- hederacea[subvars_hed]

## Correlation table

# Each species across each treatment
pcor<-round(cor(sub_purpurea[2:4]),2)
lower_purp<-pcor
lower_purp[lower.tri(pcor)]<-""
lower_purp<-as.data.frame(lower_purp)
lower_purp

hcor<-round(cor(sub_hederacea[2:4]),2)
lower_hed<-hcor
lower_hed[lower.tri(pcor)]<-""
lower_hed<-as.data.frame(lower_hed)
lower_hed

# Correlation significance
cor.test(sub_purpurea$MaxWidth,sub_purpurea$AreaConvexHull) # Root size and root width in Ip
cor.test(sub_hederacea$MaxWidth,sub_hederacea$AreaConvexHull) # Root size and root width in Ihed
cor.test(sub_hederacea$MaxWidth,sub_hederacea$AvAng) # Root width and root angle in Ihed
cor.test(sub_hederacea$AreaConvexHull,sub_hederacea$AvAng) # Root size and root angle in Ihed

##############################################
#                  FIGURES 
###############################################

dev.off()

# Regression selectin plots

# (1) Assaign the beta values for ea regression slope according to selection gradient

Beta1Alone=paste("beta[Alone] == ",0.56)
Beta1Inter=paste("beta[Inter] == ",-0.15)

Beta2Alone=paste("beta[Alone] == ",0.01)
Beta2Inter=paste("beta[Inter] == ",0.23)


# (2) Plot regression for standard root size on rel fit for I.purp

plot1<-ggplot(data=Seln2%>%filter(Species=="P"),aes(x=StdAreaConvexHull,y=Rel_Fit))+  
  geom_point(aes(color=Trt))+
  scale_color_manual(values=c('grey','black'))+
  geom_abline(intercept = 1.07, slope = -0.15, color="black", 
              linetype="dashed", size=1.5)+
  geom_abline(intercept = 1.17, slope = 0.63, color="grey", 
              linetype="solid", size=1.5)+
  ylab("Relative fitness")+
  xlab("Standardized root size")+
  theme_classic()+
  Tx+
  ggtitle("I. purpurea")+
    ggplot2::annotate("text",label=as.character(Beta1Alone),parse=T,x=0.5,y=5,size=5)+
    ggplot2::annotate("text",label=paste("± 0.29^"),x=1.1,y=5,size=5)+
    ggplot2::annotate("text",label=as.character(Beta1Inter),parse=T,x=0.5,y=4.7,size=5)+
    ggplot2::annotate("text",label=paste("± 0.22"),x=1.1,y=4.7,size=5)

# (3) Plot regression for standard root angle on rel fit for I.hed

plot2<-ggplot(data=Seln2%>%filter(Species=="H",Trt!="Intra"),aes(x=StdAvAng,y=Rel_Fit))+  
  geom_point(aes(color=Trt))+
  scale_color_manual(values=c('grey','black'))+
  geom_abline(intercept = 0.97, slope = 0.22, color="black", 
              linetype="dashed", size=1.5)+
  geom_abline(intercept = 1.04, slope = 0.01, color="grey", 
              linetype="solid", size=1.5)+
  ylab("")+
  xlab("Standardized root angle")+
  theme_classic()+
  Tx+
  ggtitle("I. hederacea")+
    ggplot2::annotate("text",label=as.character(Beta2Alone),parse=T,x=2.5,y=3,size=5)+
    ggplot2::annotate("text",label=paste("± 0.08"),x=3.3,y=3,size=5)+
    ggplot2::annotate("text",label=expression(bold(paste(beta[Inter], " = 0.23"))),parse=T,x=2.5,y=2.8,size=5)+
    ggplot2::annotate("text",label=paste("± 0.09"),x=3.3,y=2.8,size=5,fontface =2)

# (3) Combine plots
ggarrange(plot1,plot2,common.legend = T)

# Supplimentary PCA figure 1
PCA.df<-na.omit(BRT[c("Species","AreaConvexHull","MaxWidth","AvAng")])
res.pca<-PCA(PCA.df,quali.sup=1,scale.unit=T,graph=T)#scale mean of zero and std of one

# Extract the results for variables
var <- get_pca_var(res.pca)
# Contribution of variables
ContribVar<-data.frame(var$contrib)

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2)

#Add color and ellipse
p<-fviz_pca_ind(res.pca, label="none", habillage=PCA.df$Species,
                addEllipses=TRUE, ellipse.level=0.95,palette = c("grey", "black"))+
  theme_classic()+Tx

p+scale_fill_manual(values = c("white","lightgrey"))

# Bar graph of loading scores
ContribVar$Trait<-row.names(ContribVar)
# Use position=position_dodge()
q<-ggplot(data=ContribVar, aes(x=Trait, y=Dim.1),fill="black") +
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()
q 

r<-ggplot(data=ContribVar, aes(x=Trait, y=Dim.2),fill="black") +
  geom_bar(stat="identity", position=position_dodge())+
  theme_classic()
r 

