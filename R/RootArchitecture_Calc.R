# Calculate root width, primary root length and root angle code

CoordDat<-read.csv("https://raw.githubusercontent.com/SaraMColom/Selection_RootTraits_2016_2017/master/RawData/RhizotronExpCoordinates.csv") # Read in data

head(CoordDat)

# Create empty data.frame
df<-data.frame()

###Root system measurements
for(x in 1:nrow(CoordDat)) {
  
  X.lat <-CoordDat[x,]$LatRightX  # the X value of coordinate of most horz.distal root
  Y.lat<- CoordDat[x,]$LatRightY  # the Y value of coordinate of most horz.distal root
  
  X.prim1<-CoordDat[x,]$PRX # X coordinate of point of primary root at the soil surface
  Y.prim1<-CoordDat[x,]$PRY # Y coordinate of point of primary root at the soil surface
  
  a1 =sqrt( (X.lat-X.prim1)^2 + (Y.lat-Y.prim1)^2) # Length of lateral Root  to top of primary root
    
#Repeat steps above with the min value of X to find the most distal lateralroot on the other side of indivdual
    
  X.lat.2 <-CoordDat[x,]$LatLeftX  # the X value of coordinate of most horz.distal root
  Y.lat.2<- CoordDat[x,]$LatLeftY  # the Y value of coordinate of most horz.distal root
   
  a2=sqrt( (X.lat.2-X.prim1)^2 + (Y.lat.2-Y.prim1)^2) #Length of lateral Root (LLR) to top of primary root
    
 ### Next value is the distance between the most distal root lateral roots and the root tip estimated as a line perpendicular to the soil surface
    
  Y.prim2<-CoordDat[x,]$PRYtip # Y coordinate of point of primary root at the soil surface
  
    b1=sqrt( (X.lat-X.prim1)^2 + (Y.lat-Y.prim2)^2)     #Length of lateral Root (LLR) to primary root tip for one side
    b2=sqrt( (X.lat.2-X.prim1)^2 + (Y.lat.2-Y.prim2)^2) #Length of lateral Root (LLR) to tip of primary root of the other side
    
### Calculate distance of primary root top to primary root tip (also length of primary root)

c=sqrt((Y.prim1-Y.prim2)^2) #Length of primary root top to primary root tip
    
##Angle between primary root top, most horizontally distal lateral root and primary root tip
    
    Angle_1=acos((a1^2+c^2-b1^2)/(2*a1*c))
    Angle_2=acos((a2^2+c^2-b2^2)/(2*a2*c))
    
    #### Call variables to a table
    df[x,"LLR_top1"]=a1
    df[x,"LLR_tip1"]=b1
    df[x,"LLR_top2"]=a2
    df[x,"LLR_tip2"]=b2
    df[x,"PrimaryRootLength"]=c
    df[x,"Angle_1"]=Angle_1
    df[x,"Angle_2"]=Angle_2
    df[x,"RootSystemWidth"]=sqrt( (X.lat.2-X.lat)^2 + (Y.lat.2-Y.lat)^2)#Next varialbe is the width of the root system, therfore the distance between the two most distal lateral root tips
    df[x,"Id"] = gsub('.txt','',CoordDat[x,"Id"])#add column with individual Id number
    df[x,"Experiment"]=CoordDat[x,"Experiment"]
    df[x,"Code"]=CoordDat[x,"Code"]
    df[x,"Population"]=CoordDat[x,"Population"]
    df[x,"Species"]=CoordDat[x,"Species"]
    df[x,"Comment"]=CoordDat[x,"Comment"]
  }


Root.measurements<-df

# R gives the arcos in radians to convert to degrees apply the following function:
rad2deg = function(rad) {
  return((180 * rad) / pi)
}

# Multiply by factor to properly scale.

# One final step is to calibrate the sizes of the lengths by a factor of 1.47
MultiplyFactor<-function(x){x/1.47}
Transform<-c("LLR_top1", "LLR_tip1", "LLR_top2", "LLR_tip2", "PrimaryRootLength","RootSystemWidth")
Root.measurements[which(colnames(Root.measurements)%in%Transform)]<-lapply(Root.measurements[which(colnames(Root.measurements)%in%Transform)],MultiplyFactor)


#Convert to degrees
Root.measurements[grep("Angle",colnames(Root.measurements))]<-rad2deg(Root.measurements[grep("Angle",colnames(Root.measurements))])

#Calculate mean root angle
Root.measurements$AvAng<-rowMeans(Root.measurements[,c("Angle_1", "Angle_2")], na.rm=TRUE) 

