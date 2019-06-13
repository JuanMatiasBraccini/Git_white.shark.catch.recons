###### SCRIPT FOR RECONSTRUCTING WHITE SHARK CATCHES FROM COMMERCIAL FISHING IN WA  ########

#NOTES: 1. Explore data
#       2. Fit a range of models to min.catch to define optimum model
#       3. Set up bootstrapping procedure for selecting a catch sample between min and max reported
#           Then fit optimum model to bootstrapped sample. Then predict catch of unserveyed fishers

library(ggplot2)
library(car)
library(lattice)
library(bbmle)
library(lme4)
library(MASS)
library(randomForest)
library(tree)
library(rpart.plot)  			# Enhanced tree plots


#DATA SECTION

setwd("M:/Fisheries Research/FinFish/Shark/Steve/Review of fisheries management/Survey/Analysis/R")
#Data=read.csv("Cpue.ws.csv")
Data=read.csv("Cpue2.ws.Min and Max.csv")
Data$Min.catch.period=round(Data$Min.catch*Data$Number.years)
Data$Max.catch.period=round(Data$Max.catch*Data$Number.years)

#change integers to factors
Data$Time=as.factor(Data$Time)


#Select gear
# Data.gn=subset(Data,Gear.type=="GN")
# Data.wl=subset(Data,Gear.type=="WL")
Data.gn=Data

#Select fisher type
Data.gn.pred=subset(Data.gn,Fisher.type==2)
Data.gn=subset(Data.gn,Fisher.type==1)

#Drop levels for Fisher ID factor
Fishers=sort(unique(Data.gn$Fisher.ID))
Data.gn$Fisher.ID=factor(Data.gn$Fisher.ID,levels=Fishers)

#Data.gn=Data.gn[order(Data.gn$Fisher.ID),]    
#Data.gn$Fisher.ID=as.factor(Data.gn$Fisher.ID)
#Data.gn.pred$Fisher.ID=as.factor(Data.gn.pred$Fisher.ID)



#1. EXPLORE DATA
str(Data.gn) #look at data structure


#do tables of predictors
cfac=function(x,breaks=NULL){
  if(is.null(breaks)) breaks=unique(quantile(x,na.rm=T))
  x=cut(x,breaks,include.lowest=T,right=F)
  levels(x)=paste(breaks[-length(breaks)],ifelse(diff(breaks)>1,
                                                 c(paste('-',breaks[-c(1,length(breaks))]-1,sep=''),'+'),''),sep='')
  return(x)
}

fn.explore=function(DAT,Response)
{
  id=match(Response,names(DAT))
  
  Tab1=with(DAT,table(Time,Fishing.region))
  
  par(mfcol=c(3,2),mai=c(.7,.8,.3,.4))
  
  #look at 0 inflation
  hist(DAT[,id],main="zero inflation?",xlab=Response)
  
  
  #look at relation between predictor and response
  plot(DAT[,id]~Fishing.region,data=DAT,ylab=Response)
  plot(DAT[,id]~Gear.type,data=DAT,ylab=Response)
  plot(DAT[,id]~Fisher.ID,data=DAT,ylab=Response)
  plot(DAT[,id]~Time,data=DAT,ylab=Response)
  plot(DAT[,id]~cfac(Effort),data=DAT,ylab=Response)  
  
  
  #see effects of all factors
  par(mfcol=c(1,1))
  plot.design(DAT[,id]~Fishing.region+Fisher.ID+Time+cfac(Effort),data=DAT,cex.lab=1.5)
  
  
  #see interaction of factors
  interaction.plot(cfac(DAT$Effort),DAT$Fishing.region,DAT[,id],ylab=Response,
                   xlab="Region")
  
  interaction.plot(DAT$Time,DAT$Fishing.region,DAT[,id],ylab=Response,
                   xlab="Region")
  
  interaction.plot(cfac(DAT$Effort),DAT$Time,DAT[,id],ylab=Response,
                   xlab="Region")
  
  return(list(Tab1=Tab1))
}


boxplot(Min.catch.period ~ cfac(Effort)*Fishing.region,col=1:4,Data.gn)


  #Min catch
GN.only.min=fn.explore(Data.gn,"Min.catch.period")


#Max catch
GN.only.max=fn.explore(Data.gn,"Max.catch.period")



#see effects of three factors
xyplot(Data$Min.catch.period~Effort|Fishing.region*Time,data=Data,ylab="Min.catch.period",
       xlab="Effort", panel=function(x,y){panel.grid(h=-1,v=2)
                                          panel.points(x,y,col=2)
                                          panel.loess(x,y,span=.3,col=3,lwd=2)})


#2. DEFINE OPTIMAL MODEL

#note: using min.catch

  #2.1. Test model structure

    #2.1.1 Building fixed effects model
  #Poisson
simple <- glm( Min.catch.period ~ 1, data = Data.gn,family =poisson)
simple.one <- glm( Min.catch.period ~ Fisher.ID, data = Data.gn,family =poisson)
simple.two <- glm( Min.catch.period ~ Fisher.ID+ Effort, data = Data.gn,family =poisson)
simple.three <- glm( Min.catch.period ~ Fisher.ID+ Effort+Fishing.region, data = Data.gn,family =poisson)
simple.four <- glm( Min.catch.period ~ Fisher.ID+ Effort+Fishing.region+Time , data = Data.gn,family =poisson)
simple.five <- glm( Min.catch.period ~ Fisher.ID+ Effort*Fishing.region+Time, data = Data.gn,family =poisson)
simple.six <- glm( Min.catch.period ~ Effort*Fishing.region+Time+Fishing.region/Fisher.ID, data = Data.gn,family =poisson)
simple.seven<- glm( Min.catch.period ~ Fisher.ID+ Effort*Fishing.region+Time*Fishing.region, data = Data.gn,family =poisson)

  #NB
simple.nb <- glm.nb( Min.catch.period ~ 1, data = Data.gn)
simple.one.nb <- glm.nb( Min.catch.period ~ Fisher.ID, data = Data.gn)
simple.two.nb <- glm.nb( Min.catch.period ~ Fisher.ID+ Effort, data = Data.gn)
simple.three.nb <- glm.nb( Min.catch.period ~ Fisher.ID+ Effort+Fishing.region, data = Data.gn)
simple.four.nb <- glm.nb( Min.catch.period ~ Fisher.ID+ Effort+Fishing.region+Time , data = Data.gn)
simple.five.nb<- glm.nb( Min.catch.period ~ Fisher.ID+ Effort*Fishing.region+Time, data = Data.gn)
simple.six.nb <- glm.nb( Min.catch.period ~ Effort*Fishing.region+Time+Fishing.region/Fisher.ID, data = Data.gn)
simple.seven.nb<- glm.nb( Min.catch.period ~ Fisher.ID+ Effort*Fishing.region+Time*Fishing.region, data = Data.gn)


AICtab(simple,simple.one,simple.two,simple.three,simple.four,simple.five,simple.six,simple.seven,simple.nb,
       simple.one.nb,simple.two.nb,simple.three.nb,simple.four.nb,simple.five.nb,simple.six.nb,simple.seven.nb)
plot(predict(simple.five,type='response'))
points(predict(simple.six,type='response'),col=2,pch=19)
points(predict(simple,type='response'),col=4,pch=19)
anova(update(simple.six,.~1),simple.six,test="Chisq")

#Percent of deviance explained
100*(simple.seven$null.deviance - simple.seven$deviance) / simple.seven$null.deviance

#Terms significance
#drop1(simple.seven,test="Chi")
anova(simple.seven,test="Chi")

    #2.1.2 Building mixed model
simple <- glmer( Min.catch.period ~ 1+ (1|Fisher.ID), data = Data.gn,family =poisson)
simple.one <- glmer( Min.catch.period ~ Effort+ (1|Fisher.ID), data = Data.gn,family =poisson)
simple.two <- glmer( Min.catch.period ~ Effort+Fishing.region +(1|Fisher.ID), data = Data.gn,family =poisson)
simple.three <- glmer( Min.catch.period ~ Effort+Fishing.region+Time +(1|Fisher.ID), data = Data.gn,family =poisson)
simple.four <- glmer( Min.catch.period ~ Effort+Fishing.region+Time +(1|Fisher.ID), data = Data.gn,family =poisson)
simple.five <- glmer( Min.catch.period ~ Effort*Fishing.region+Time+ (1|Fisher.ID), data = Data.gn,family =poisson)
simple.six <- glmer( Min.catch.period ~ Effort*Fishing.region+Time+ (1|Fishing.region:Fisher.ID), data = Data.gn,family =poisson)
AICtab(simple,simple.one,simple.two,simple.three,simple.four,simple.five,simple.six)
plot(predict(simple.five,type='response'))
points(predict(simple.six,type='response'),col=2,pch=19)
points(predict(simple.two,type='response'),col=4,pch=19)


# get model deviance
dev.simple.five <- -2*logLik(simple.five)[1]
dev.simple.null <- -2*logLik(simple)[1]



  #2.2. Test different error structures

  #FIXED TERMS
    #a. main effects only
  lm1 <- glm( Min.catch.period ~ Fisher.ID+Fishing.region+Time+Effort, data = Data.gn,family =poisson)
  slm1.min <- step(lm1)
  summary(slm1.min)
  
    #b. main effects and effort as offset
  lm2 <- glm( Min.catch.period ~ Fisher.ID+Fishing.region+Time+offset(Effort), data = Data.gn,family =poisson)
  slm2.min <- step(lm2)
  summary(slm2.min)
  
    #c. effort-region interaction
  #Poisson
  lm3 <- glm( Min.catch.period ~ Fisher.ID+Time+Effort*Fishing.region, data = Data.gn,family =poisson)
  slm3.min <- step(lm3)
  summary(slm3.min)
  
  #Negative Binomial
  lm4 <- glm.nb( Min.catch.period ~ Fisher.ID+Time+Effort*Fishing.region, data = Data.gn,link=log)
  slm4.min <- step(lm4)
  summary(slm4.min)

  #d.effort-region and time-fishing region interactions
  lm5 <- glm(Min.catch.period ~ Fisher.ID+Time*Fishing.region+Effort*Fishing.region, data = Data.gn, family=poisson)
  slm5.min <- step(lm5)
  summary(slm5.min)

  
  #MIXED MODEL (FIXED AND RANDOM TERMS) 

    #a. fisher as random and main effects
Boris.main<-glmer( Min.catch.period ~ Effort+Fishing.region+Time+(1|Fisher.ID), data = Data.gn,family =poisson)

    #b. fisher as random and effort-region interaction
Boris<-glmer( Min.catch.period ~ Effort*Fishing.region+Time+(1|Fisher.ID), data = Data.gn,family =poisson)

    #c. fisher as random and effort-region and Time-Fishing.region interaction
#NEW
Boris2<-glmer( Min.catch.period ~ Fishing.region+Time+(1|Fisher.ID),offset=Effort, data = Data.gn,family =poisson)
#Boris2<-glmer( Min.catch.period ~ Effort*Fishing.region+Time*Fishing.region+(1|Fisher.ID), data = Data.gn,family =poisson)


    #d.fisher as random nested in region  and effort-region interaction
Boris.nested <- glmer( Min.catch.period ~ Time+Effort*Fishing.region+ (1|Fishing.region:Fisher.ID), data = Data.gn,family =poisson)


#compare all selected models
AICtab(slm1.min,slm2.min,slm3.min,slm4.min,slm5.min,Boris.main,Boris,Boris2,Boris.nested)
BIC(slm5.min)

#updated by Steve on Sep 8
AIC(Boris.main,Boris,Boris2,Boris.nested)

BIC(Boris.main,Boris,Boris2,Boris.nested)


  #2.3. Diagnostics for selected model

#Summary of selected model
summary(Boris2)
anova(Boris2)

#get significance for final model
#NOTE: this doesn't work for mixed models!! do it manually...

#1.this tells if full model is Significant compared to null
# anova(update(slm3.min,.~1),slm3.min,test="Chisq") 
# # (or same thing by doing)
# sg1.min <- summary(slm3.min)
# devdiff.min <- with(sg1.min,null.deviance-deviance)
# dfdiff.min <- with(sg1.min,df.null-df.residual)
# pchisq(abs(devdiff.min),df=dfdiff.min,lower.tail=FALSE) 
# 
# #2. this tells significance of each term
# anova(slm3.min,test = "Chisq")  #you get the Chisq test by setting the test argument....

#plot fit
par(mfcol=c(2,2))
plot(slm3.min)

hist(residuals(slm3.min))

# Cook distance to see outliers or overdisperse data (if distance >1)
plot(slm3.min,which=4)


# Calculate overdispersion (if larger than 1 then overdispersed!)
E1.min=residuals(slm3.min)
Npars.min=length(coef(slm3.min))
Overdisp1.min=sum(E1.min^2)/(nrow(Data)-(Npars.min))


#proportion of variance explained (res deviance/null deviance)
#Prop.var.expl=62.946/476.781
Prop.var.expl=1-(slm3.min$dev/slm3.min$null)

#Summary of selected model updated by Steve to look at slm4
summary(slm4.min)

#get significance for final model

#1.this tells if full model is Significant compared to null
anova(update(slm4.min,.~1),slm4.min,test="Chisq") 
# (or same thing by doing)
sg1.min <- summary(slm4.min)
devdiff.min <- with(sg1.min,null.deviance-deviance)
dfdiff.min <- with(sg1.min,df.null-df.residual)
pchisq(abs(devdiff.min),df=dfdiff.min,lower.tail=FALSE) 

#2. this tells significance of each term
anova(slm4.min,test = "Chisq")  #you get the Chisq test by setting the test argument....

#Summary of selected model updated by Steve to look at slm5
summary(slm5.min)

#get significance for final model

#1.this tells if full model is Significant compared to null
anova(update(slm5.min,.~1),slm5.min,test="Chisq") 
# (or same thing by doing)
sg1.min <- summary(slm5.min)
devdiff.min <- with(sg1.min,null.deviance-deviance)
dfdiff.min <- with(sg1.min,df.null-df.residual)
pchisq(abs(devdiff.min),df=dfdiff.min,lower.tail=FALSE) 

#2. this tells significance of each term
anova(slm5.min,test = "Chisq")  #you get the Chisq test by setting the test argument....



#note: don't need to use max.catch because in the bootstrapping, the response var is a sample
#     between min and max catch....

##########################################
#1. Max Catch

#Define optimal model

#main effects only
lm1 <- glm( Max.catch.period ~ Fisher.ID+Fishing.region+Time+Effort, data = Data.gn,family =poisson)
slm1 <- step(lm1)
summary(slm1)

#main effects and effort as offset
lm2 <- glm( Max.catch.period ~ Fisher.ID+Fishing.region+Time+offset(Effort), data = Data.gn,family =poisson)
slm2 <- step(lm2)
summary(slm2)

#effort and region interaction
#Poisson
lm3 <- glm( Max.catch.period ~ Fisher.ID+Time+Effort*Fishing.region, data = Data.gn,family =poisson)
slm3 <- step(lm3)
summary(slm3)

#Negative Binomial
lm4 <- glm.nb( Max.catch.period ~ Fisher.ID+Time+Effort*Fishing.region, data = Data.gn,link=log)
slm4 <- step(lm4)
summary(slm4)

#Added by Steve Poisson with additional time*fishing region interaction
lm5 <- glm( Min.catch.period ~ Fisher.ID+Time*Fishing.region+Effort*Fishing.region, data = Data.gn, family=poisson)
slm5 <- step(lm5)
summary(slm5)
#compare all selected models
AICtab(slm1,slm2,slm3,slm4,slm5)

#MIXED MODEL (FIXED AND RANDOM TERMS) 

#a. fisher as random and main effects
Boris.mainmax<-glmer( Max.catch.period ~ Effort+Fishing.region+Time+(1|Fisher.ID), data = Data.gn,family =poisson)

#b. fisher as random and effort-region interaction
Borismax<-glmer( Max.catch.period ~ Effort*Fishing.region+Time+(1|Fisher.ID), data = Data.gn,family =poisson)

#c. fisher as random and effort-region and Time-Fishing.region interaction
Boris2max<-glmer( Max.catch.period ~ Effort*Fishing.region+Time*Fishing.region+(1|Fisher.ID), data = Data.gn,family =poisson)


#d.fisher as random nested in region  and effort-region interaction
Boris.nestedmax <- glmer( Max.catch.period ~ Time+Effort*Fishing.region+ (1|Fishing.region:Fisher.ID), data = Data.gn,family =poisson)


#compare all selected models
AICtab(slm1.min,slm2.min,slm3.min,slm4.min,slm5.min,Boris.mainmax,Borismax,Boris2max,Boris.nestedmax)
BIC(slm5.min)

#updated by Steve on Sep 8
AIC(Boris.mainmax,Borismax,Boris2max,Boris.nestedmax)

BIC(Boris.mainmax,Borismax,Boris2max,Boris.nestedmax)

#get significance for final model

#1.this tells if full model is Significant compared to null
anova(update(slm4,.~1),slm4,test="Chisq") 
stepAIC(slm4)


sg1 <- summary(slm4)
devdiff <- with(sg1,null.deviance-deviance)
dfdiff <- with(sg1,df.null-df.residual)
pchisq(abs(devdiff),df=dfdiff,lower.tail=FALSE) 

#2. this tells significance of each term
anova(slm4,test = "Chisq")  #you get the Chisq test by setting the test argument....


#plot fit
par(mfcol=c(2,2))
plot(slm4)

hist(residuals(slm4))

# Cook distance to see outliers or overdisperse data (if distance >1)
plot(slm4,which=4)


# Calculate overdispersion (if larger than 1 then overdispersed!)
E1=residuals(slm4)
Npars=length(coef(slm4))
Overdisp1=sum(E1^2)/(nrow(Data)-(Npars))

E4=resid(slm4,type="pearson")
Dispersion=sum(E4^2)/slm4$df.resid


#proportion of variance explained (res deviance/null deviance)



#Check observed-predicted

PRED.Min=predict(slm3.min,type="response",se.fit =T)
PRED.Max=predict(slm4,type="response",se.fit =T)

Data.gn$GLM.pred.min=round(PRED.Min$fit)
Data.gn$GLM.pred.min.SE=round(PRED.Min$se.fit)
Data.gn$GLM.pred.max=round(PRED.Max$fit)
Data.gn$GLM.pred.max.SE=round(PRED.Max$se.fit)

tiff(file="GLM.fit.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")

par(mfcol=c(2,1),mgp=c(2.5,1,0))
plot(Data.gn$Min.catch.period,Data.gn$GLM.pred.min,ylim=c(0,max(Data.gn$Min.catch.period)),
     ylab="Predicted catch", xlab="Observed catch",las=1,cex.axis=1.5,cex.lab=2,pch=19,cex=1.25,
     main="Min catch")
lines(Data.gn$Min.catch.period,Data.gn$Min.catch.period,col=2,lwd=2.5)

plot(Data.gn$Max.catch.period,Data.gn$GLM.pred.max,ylim=c(0,max(Data.gn$Min.catch.period)),
     ylab="Predicted catch", xlab="Observed catch",las=1,cex.axis=1.5,cex.lab=2,pch=19,cex=1.25,
     main="Max catch")
lines(Data.gn$Max.catch.period,Data.gn$Max.catch.period,col=2,lwd=2.5)

dev.off()

####################################



#3. BOOTSTRAP DATA WITHIN FISHER
Data.boot=read.csv("Cpue2.ws.Min and Max.csv")
Data.boot.non.int=subset(Data.boot,Fisher.type==2)

# extract total effort by region and period
Tot.Eff=aggregate(Effort~Fishing.region+Time,Data.boot,sum)
Interviewd.Eff=aggregate(Effort~Fishing.region+Time,subset(Data.boot,Fisher.type==1),sum)

tiff(file="Effort.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mai=c(1,1,.1,.35),oma=c(.1,.1,.1,.1),las=1,mgp=c(3,.5,0))
Tot.Eff.reshaped=reshape(Tot.Eff,v.names = "Effort", idvar = "Time",timevar = "Fishing.region", direction = "wide")
with(Tot.Eff.reshaped,plot(Time,Effort.SS1,type="b",ylab="Effort (1000 km gn d)",pch=19,cex=2,xaxt="n",
                           cex.lab=2,cex.axis=1.5,ylim=c(0,max(Effort.SS2))))
with(Tot.Eff.reshaped,points(Time,Effort.SS2,type="b",pch=19,cex=2,col="grey30"))
with(Tot.Eff.reshaped,points(Time,Effort.WC,type="b",pch=19,cex=2,col="grey60"))
legend("topright",c("Zn1","Zn2","WC"),bty="n",pch=19,lty=1,col=c("black","grey40","grey60"),cex=2)
axis(1,at=c(1,2,3),labels=c("1988-1996","1997-2004","2005-2013"),tck=-0.02,cex.axis=1.50)
dev.off()

names(Tot.Eff)[3]="Total.Effort"

#aggregate fisher ID with effort under threshold
Get.threshold=subset(Data.boot,Min.catch>0)
Get.Min=function(DAT) min(DAT$Effort)
Reg=unique(Data.boot$Fishing.region)
TIM=unique(Data.boot$Time)
Store.reg=matrix(nrow=length(Reg),ncol=length(TIM))
colnames(Store.reg)=TIM
rownames(Store.reg)=Reg
for (p in 1:length(Reg))
{
  Store.dum=rep(NA,length(TIM))
  for (x in 1:length(TIM)) Store.dum[x]=Get.Min(subset(Get.threshold,Fishing.region==Reg[p] & Time ==TIM[x]))
  Store.reg[p,]=Store.dum
}


Data.boot.Pred=subset(Data.boot,Fisher.type==2)
Data.boot.Pred$Aggregate=with(Data.boot.Pred,ifelse(Fishing.region=="WC" & Time==1 & Effort<Store.reg[1,1],"YES",
                      ifelse(Fishing.region=="WC" & Time==2 & Effort<Store.reg[1,2],"YES",
                      ifelse(Fishing.region=="WC" & Time==3 & Effort<Store.reg[1,3],"YES",
                      ifelse(Fishing.region=="SS2" & Time==1 & Effort<Store.reg[2,1],"YES",
                      ifelse(Fishing.region=="SS2" & Time==2 & Effort<Store.reg[2,2],"YES",
                      ifelse(Fishing.region=="SS2" & Time==3 & Effort<Store.reg[2,3],"YES",
                      ifelse(Fishing.region=="SS1" & Time==1 & Effort<Store.reg[3,1],"YES",
                      ifelse(Fishing.region=="SS1" & Time==2 & Effort<Store.reg[3,2],"YES",
                      ifelse(Fishing.region=="SS1" & Time==3 & Effort<Store.reg[3,3],"YES","NO"))))))))))

Data.boot.Pred=subset(Data.boot.Pred,Aggregate=="YES")
Get.Fisher.ID=unique(Data.boot.Pred$Fisher.ID)
Data.boot.Pred.eff.sum=aggregate(Effort~Fishing.region+Time,Data.boot.Pred,sum)

Data.boot.Pred.eff.sum$Fisher.ID=paste("Agg",paste(Data.boot.Pred.eff.sum$Fishing.region,Data.boot.Pred.eff.sum$Time))
Data.boot.Pred.eff.sum$Fisher.type=2
Data.boot.Pred.eff.sum$Min.catch=NA
Data.boot.Pred.eff.sum$Max.catch=NA
Data.boot.Pred.eff.sum$Number.years=NA
Data.boot.Pred.eff.sum$Gear.type="GN"

Data.boot=subset(Data.boot,!Fisher.ID%in%Get.Fisher.ID)

Data.boot.Pred.eff.sum=subset(Data.boot.Pred.eff.sum,select=names(Data.boot))

Data.boot=rbind(Data.boot,Data.boot.Pred.eff.sum)



  #3.1 change integers to factors
Data.boot$Time=as.factor(Data.boot$Time)

  #3.2 Drop levels for Fisher ID factor
Data.boot=Data.boot[order(Data.boot$Fisher.ID),]
Data.boot.non.sampled=subset(Data.boot,Fisher.type==2)
Data.boot=subset(Data.boot,Fisher.type==1)
Fishers=sort(unique(Data.boot$Fisher.ID))
Data.boot$Fisher.ID=factor(Data.boot$Fisher.ID,levels=Fishers)
Data.boot$Min.per=round(Data.boot$Min.catch*Data.boot$Number.years)
Data.boot$Max.per=round(Data.boot$Max.catch*Data.boot$Number.years)

N.boots=1000

Dat.list=Dat.pred=GLM.list=vector('list',length=N.boots)

#Either keeping aggregated fisher or not
#Keep.agg="NO"
Keep.agg="YES"
if(Keep.agg=="NO")
{
  Data.boot.non.sampled=subset(Data.boot.non.sampled,!Fisher.ID%in%c("Agg SS1 1","Agg SS1 2","Agg SS1 3",
            "Agg SS2 1","Agg SS2 2","Agg SS2 3","Agg WC 1","Agg WC 2","Agg WC 3"))
}

#number of fishers to infer by region
Table.Inf=with(Data.boot.non.sampled,table(Fishing.region))

#interviewed fishers by zone
Data.boot.WC=subset(Data.boot,Fishing.region=="WC")
Data.boot.Zn1=subset(Data.boot,Fishing.region=="SS1")
Data.boot.Zn2=subset(Data.boot,Fishing.region=="SS2")

#Number of interviewed and total fishers by region and time
Table.interviewed=table(Data.boot$Fishing.region,Data.boot$Time)
Table.non.interviewed=table(Data.boot.non.int$Fishing.region,Data.boot.non.int$Time)
Table.total=Table.interviewed+Table.non.interviewed


  #bootstrapping function for fisher
Resampl=function(DATs) DATs[sample(1:nrow(DATs),nrow(DATs),replace=T),]

  #bootstrapping function for catch
fn.boot=function(a,b) round(runif(1,a,b))

  #3.3 Bootstrap observed effort
system.time(for(j in 1:N.boots)
{
  #1.a. bootstrap interviewed fisher within each zone with replacement     #NEW
  Boot.WC=Resampl(Data.boot.WC)
  Boot.Zn1=Resampl(Data.boot.Zn1)
  Boot.Zn2=Resampl(Data.boot.Zn2)
  dummy=rbind(Boot.WC,Boot.Zn1,Boot.Zn2)
  dummy$Fisher.ID=dummy$Fisher.ID[, drop=TRUE]
  
  #1.b bootstap catch within reported range
  dummy$Catch.boot=mapply(fn.boot,dummy$Min.per,dummy$Max.per)                                 
  
  DAT=dummy[,match(c("Fisher.ID","Fishing.region","Time","Effort","Catch.boot"),names(dummy))]
  DAT=DAT[order(DAT$Fisher.ID),]
  DAT$Pred=DAT$Catch.boot
    
  #2. store stuff
    Dat.list[[j]]=DAT
})
  
  #3.4 Predict Total catch from total effort 
All.catch.folly=All.catch.average=vector('list',length=N.boots)
for(j in 1:N.boots)
{
  DaT=Dat.list[[j]]
  if(!is.null(DaT))
  {
    #Folly approach
    CATCH.obs=aggregate(Catch.boot~Fishing.region+Time,DaT,sum)
    EFF.obs=aggregate(Effort~Fishing.region+Time,DaT,sum)
    if(length(CATCH.obs[,3])==length(Tot.Eff$Total.Effort))
    {
      CATCH.obs=merge(CATCH.obs,EFF.obs,by=c("Fishing.region","Time"))
      CATCH.obs$CPUE=CATCH.obs$Catch.boot/CATCH.obs$Effort
      CATCH.obs=merge(CATCH.obs,Tot.Eff,by=c("Fishing.region","Time"))
      CATCH.obs$Total.Catch=(CATCH.obs$CPUE*CATCH.obs$Total.Effort)
      All.catch.folly[[j]]=CATCH.obs
    }
    
    #Average approach
    if(length(CATCH.obs[,3])==length(Tot.Eff$Total.Effort))
    {
      DaT$cpue=DaT$Catch.boot/DaT$Effort
      Mean.obs.cpue=aggregate(cpue~Fishing.region+Time,DaT,mean)
      Mean.obs.cpue=merge(Mean.obs.cpue,Tot.Eff,by=c("Fishing.region","Time"))
      Mean.obs.cpue$Total.Catch=(Mean.obs.cpue$cpue*Mean.obs.cpue$Total.Effort)
      All.catch.average[[j]]=Mean.obs.cpue
    }
  }

}


  #3.5 Get mean catch and CI
RS.All.catch.folly=RS.All.catch.average=matrix(nrow=nrow(Tot.Eff),ncol=N.boots)
RS.All.cpue=RS.All.catch.folly
for(j in 1:N.boots)
{
  if(!is.null(All.catch.folly[[j]]))RS.All.catch.folly[,j]=All.catch.folly[[j]]$Total.Catch
  if(!is.null(All.catch.average[[j]]))RS.All.catch.average[,j]=All.catch.average[[j]]$Total.Catch
  if(!is.null(All.catch.folly[[j]]))RS.All.cpue[,j]=All.catch.folly[[j]]$CPUE
}


Mean.Catch.folly=round(apply(RS.All.catch.folly, 1, mean,na.rm=T))
SD.Catch.folly=round(apply(RS.All.catch.folly, 1, sd,na.rm=T))
Low95.Catch.folly=round(apply(RS.All.catch.folly, 1, function(x)quantile(x,probs=0.025,na.rm=T)))
Up95.Catch.folly=round(apply(RS.All.catch.folly, 1, function(x)quantile(x,probs=0.975,na.rm=T)))
Catch.Folly=cbind(Mean.Catch.folly,Low95.Catch.folly,Up95.Catch.folly)
Tot.Eff=Tot.Eff[order(Tot.Eff$Fishing.region),]
Catch.Folly=cbind(Tot.Eff,Catch.Folly)

# Mean.Catch.avrg=round(apply(RS.All.catch.average, 1, mean,na.rm=T))
# Low95.Catch.avrg=round(apply(RS.All.catch.average, 1, function(x)quantile(x,probs=0.025,na.rm=T)))
# Up95.Catch.avrg=round(apply(RS.All.catch.average, 1, function(x)quantile(x,probs=0.975,na.rm=T)))
# Catch.avrg=cbind(Mean.Catch.avrg,Low95.Catch.avrg,Up95.Catch.avrg)
# Catch.avrg=cbind(Tot.Eff,Catch.avrg)


#correct confidence bounds

    #Finite population correction
# CI.bi.co=function(n,N) (1-(n/N))^0.5
# Prop.Eff=merge(Tot.Eff,Interviewd.Eff,by=c("Fishing.region","Time"))
# 
# Catch.Folly$Low95.Catch.folly.c=Catch.Folly$Low95.Catch.folly
# 
# for(i in 1:nrow(Catch.Folly))
#   {
#     Catch.Folly$Low95.Catch.folly.c[i]=Catch.Folly$Mean.Catch.folly[i]-1.96*
#             ((SD.Catch.folly[i])^0.5/(c(t(Table.interviewed))[i])^0.5)*
#             CI.bi.co(c(t(Table.interviewed))[i],c(t(Table.total))[i])
#     
#     Catch.Folly$Up95.Catch.folly.c[i]=Catch.Folly$Mean.Catch.folly[i]+1.96*
#       ((SD.Catch.folly[i])^0.5/(c(t(Table.interviewed))[i])^0.5)*
#       CI.bi.co(c(t(Table.interviewed))[i],c(t(Table.total))[i])
#   }



#3.6 Plot mean and CI
fn.CI=function(x,y,z,jit,COL)
{
  segments(x,y,x,z,col=COL)
  segments(x,y,x,y,col=COL)
  segments(x,z,x,z,col=COL)
}

fn.plot=function(DATA,YMAX)
{
  names(DATA)[4:6]=c("Mean","Low95","Up95")
  RS.DATA=reshape(DATA[,-match("Total.Effort",names(DATA))], v.names = c("Mean","Low95","Up95"),
                  idvar = "Time", timevar = "Fishing.region", direction = "wide")
  with(RS.DATA,
       {
         plot(Time,Mean.SS1,ylim=c(0,YMAX),xlab="",
              ylab="",cex.lab=2,xaxt="n",cex=1.75,pch=15,
              xlim=c(0.9,3.2),cex.axis=1.5)
         points(Time+.1,Mean.SS2,col=1,cex=1.75,pch=17)
         points(Time-.1,Mean.WC,col=1,cex=1.75,pch=16)
         
         fn.CI(Time,Low95.SS1,Up95.SS1,1.02,1)
         fn.CI(Time+.1,Low95.SS2,Up95.SS2,1.07,1)
         fn.CI(Time-.1,Low95.WC,Up95.WC,0.96,1)
         
         axis(1,at=c(1,2,3),labels=c("1988-1996","1997-2004","2005-2013"),tck=-0.02,cex.axis=1.50)
        })
}

  #Export
write.table(Catch.Folly,"Total.catch.csv",sep=",",row.names = F)



#Explore boostrapped catch distribution
LABS=paste(All.catch.folly[[1]]$Fishing.region,All.catch.folly[[1]]$Time)
par(mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1),las=1,mgp=c(3,.5,0))

tiff(file="Explore.Boot.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
par(mfcol=c(3,3))
for(i in 1:nrow(RS.All.catch.folly))hist(RS.All.catch.folly[i,],main=LABS[i],xlab="Catch (numbers)")
dev.off()


Mean.cpue.folly=apply(RS.All.cpue, 1, mean,na.rm=T)
SD.cpue.folly=apply(RS.All.cpue, 1, sd,na.rm=T)
Low95.cpue.folly=apply(RS.All.cpue, 1, function(x)quantile(x,probs=0.025,na.rm=T))
Up95.cpue.folly=apply(RS.All.cpue, 1, function(x)quantile(x,probs=0.975,na.rm=T))
CPUE.Folly=cbind(Mean.cpue.folly,SD.cpue.folly,Low95.cpue.folly,Up95.cpue.folly)
CPUE.Folly=cbind(Tot.Eff,CPUE.Folly)


#Export
write.table(CPUE.Folly,"CPUE.Folly.csv",sep=",",row.names = F)


  #Method 2: Mean cpue and CI
#note: this method uses zone 1 cpue to construct zone 2 cpue for time 1 and 2
All.catch.folly.Method2=vector('list',length=N.boots)
for(j in 1:N.boots)
{
  dat=All.catch.folly[[j]]
  if(!is.null(dat))
  {
    #first do time 2 zone 2
    A=dat$CPUE[dat$Fishing.region=="SS1" & dat$Time==2]/dat$CPUE[dat$Fishing.region=="SS1" & dat$Time==3]
    id=which(dat$Fishing.region=="SS2" & dat$Time==2)
    dat$CPUE[id]=A*dat$CPUE[dat$Fishing.region=="SS2" & dat$Time==3]
    
    #then do time 1 zone 2
    B=dat$CPUE[dat$Fishing.region=="SS1" & dat$Time==1]/dat$CPUE[dat$Fishing.region=="SS1" & dat$Time==2]
    id=which(dat$Fishing.region=="SS2" & dat$Time==1)
    dat$CPUE[id]=B*dat$CPUE[dat$Fishing.region=="SS2" & dat$Time==2]
    
    dat$Total.Catch=dat$CPUE*dat$Total.Effort
    
    All.catch.folly.Method2[[j]]=dat
    
  }
}


####
#Get mean catch and CI
RS.All.catch.folly.Method2=matrix(nrow=nrow(Tot.Eff),ncol=N.boots)
RS.All.cpue.Method2=RS.All.catch.folly
for(j in 1:N.boots)
{
  if(!is.null(All.catch.folly.Method2[[j]]))RS.All.catch.folly.Method2[,j]=All.catch.folly.Method2[[j]]$Total.Catch
  if(!is.null(All.catch.folly.Method2[[j]]))RS.All.cpue.Method2[,j]=All.catch.folly.Method2[[j]]$CPUE
}


Mean.Catch.folly.Method2=round(apply(RS.All.catch.folly.Method2, 1, mean,na.rm=T))
SD.Catch.folly.Method2=round(apply(RS.All.catch.folly.Method2, 1, sd,na.rm=T))
Low95.Catch.folly.Method2=round(apply(RS.All.catch.folly.Method2, 1, function(x)quantile(x,probs=0.025,na.rm=T)))
Up95.Catch.folly.Method2=round(apply(RS.All.catch.folly.Method2, 1, function(x)quantile(x,probs=0.975,na.rm=T)))
Catch.Folly.Method2=cbind(Mean.Catch.folly.Method2,Low95.Catch.folly.Method2,Up95.Catch.folly.Method2)
Catch.Folly.Method2=cbind(Tot.Eff,Catch.Folly.Method2)


#Export
write.table(Catch.Folly.Method2,"Total.catch.Method2.csv",sep=",",row.names = F)

tiff(file="Figure2AllFishers.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff
par(mfcol=c(2,1),mai=c(.6,1,.1,.1),oma=c(.75,.1,.5,.1),las=1,mgp=c(3,.5,0))
#par(mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1),las=1,mgp=c(3,.5,0))
fn.plot(Catch.Folly,YMAX=870)
legend("topright",c("Zone 1","Zone 2","West Coast"),bty='n',pch=c(15,17,16),
       col=1,bg=c(1,1,1),cex=1.5,pt.cex=2)
fn.plot(Catch.Folly.Method2,YMAX=870)
mtext("Mean number caught (+/- 95% CI)",2,line=-2,outer=T,cex=2,las=3)
mtext("Time period",1,line=-.5,outer=T,cex=2)
dev.off()



Mean.cpue.folly.Method2=apply(RS.All.cpue.Method2, 1, mean,na.rm=T)
SD.cpue.folly.Method2=apply(RS.All.cpue.Method2, 1, sd,na.rm=T)
Low95.cpue.folly.Method2=apply(RS.All.cpue.Method2, 1, function(x)quantile(x,probs=0.025,na.rm=T))
Up95.cpue.folly.Method2=apply(RS.All.cpue.Method2, 1, function(x)quantile(x,probs=0.975,na.rm=T))
CPUE.Folly.Method2=cbind(Mean.cpue.folly.Method2,SD.cpue.folly.Method2,
                         Low95.cpue.folly.Method2,Up95.cpue.folly.Method2)
CPUE.Folly.Method2=cbind(Tot.Eff,CPUE.Folly.Method2)


#Export
write.table(CPUE.Folly.Method2,"CPUE.Folly.Method2.csv",sep=",",row.names = F)


######################################
#NOTE USED

#3.3 Bootstrap observed effort
system.time(for(j in 1:N.boots)
{
  #1. fit selected model to bootstrapped data
  
  #1.a. bootstrap interviewed fisher within each zone with replacement     #NEW
  Boot.WC=Resampl(Data.boot.WC)
  Boot.Zn1=Resampl(Data.boot.Zn1)
  Boot.Zn2=Resampl(Data.boot.Zn2)
  dummy=rbind(Boot.WC,Boot.Zn1,Boot.Zn2)
  dummy$Fisher.ID=dummy$Fisher.ID[, drop=TRUE]
  
  #1.b bootstap catch within reported range
  dummy$Catch.boot=mapply(fn.boot,dummy$Min.per,dummy$Max.per)                                 
  
  tryCatch({
    fit=glm( Catch.boot ~ Fisher.ID+ Effort*Fishing.region+Time*Fishing.region, data = dummy,family =poisson)
    
    
    #  fit=glmer( Catch.boot ~ Fishing.region+Time+(1|Fisher.ID),offset=Effort, 
    #                       data = dummy,family =poisson)
    #   fit=glmer( Catch.boot ~ Effort*Fishing.region+Time*Fishing.region+(1|Fisher.ID), 
    #              data = dummy,family =poisson)
    #   
    
    #2. predict catches
    
    #2.1. observations (i.e. interviewed fishers)
    #RAND=unlist(ranef(fit)["Fisher.ID"])
    #FIX=fixef(fit)
    #Rand.Coef=data.frame(Fisher.ID=levels(dummy$Fisher.ID),Rand.coef=RAND)
    
    
    DAT=dummy[,match(c("Fisher.ID","Fishing.region","Time","Effort","Catch.boot"),names(dummy))]
    DAT=DAT[order(DAT$Fisher.ID),]
    
    ##: previous previuos
    #: previous
    
    ##mm = model.matrix(terms(fit),DAT)
    
    ##DAT=merge(DAT,Rand.Coef,by="Fisher.ID")
    ##DAT=DAT[order(DAT$Fisher.ID),]
    ##   DAT$Pred1 = mm %*% fixef(fit)
    #   
    ##   for (i in 1:length(DAT$Pred)) DAT$Pred1[i]=DAT$Pred1[i]+RAND[DAT$Fisher.ID[i]]
    # 
    ##   DAT$Pred1=exp(DAT$Pred1)
    ##   DAT$fitted=fitted(fit)
    
    #DAT$Pred1=round(predict(fit,type='response'))
    DAT$Pred=DAT$Catch.boot
    
    #2.2. non-interviewed fishers
    ##missing.fisher=length(unique(Data.boot.non.sampled$Fisher.ID))
    #     missing.zn1=Table.Inf[1]
    #     missing.zn2=Table.Inf[2]
    #     missing.WC=Table.Inf[3]
    #     
    ##Rand.Cof=subset(DAT,select=c(Fisher.ID,Fishing.region,Rand.coef))
    #     Rand.Cof=subset(DAT,select=c(Fisher.ID,Fishing.region))
    #     Rand.Cof=Rand.Cof[!duplicated(Rand.Cof$Fisher.ID),]
    #     Rand.Cof.WC=subset(Rand.Cof,Fishing.region=="WC")
    #     Rand.Cof.Zn1=subset(Rand.Cof,Fishing.region=="SS1")
    #     Rand.Cof.Zn2=subset(Rand.Cof,Fishing.region=="SS2")
    
    ##Fish.samp=sample(RAND,missing.fisher,replace=T) #take random sample of fisher coefficients   #NEW
    ##Fish.samp=rnorm(missing.fisher,mean(RAND),sd(RAND)) #take random sample of fisher coefficients
    
    ##Fish.samp.WC=rnorm(missing.WC,mean(Rand.Cof.WC$Rand.coef),sd(Rand.Cof.WC$Rand.coef)) 
    ##Fish.samp.Zn1=rnorm(missing.zn1,mean(Rand.Cof.Zn1$Rand.coef),sd(Rand.Cof.Zn1$Rand.coef))
    ##Fish.samp.Zn2=rnorm(missing.zn2,mean(Rand.Cof.Zn2$Rand.coef),sd(Rand.Cof.Zn2$Rand.coef))
    
    
    #Get mean vessel
    #ID=which(substr(names(coef(fit)),1,9)=="Fisher.ID")
    #Ves.coef=sort(coef(fit)[ID])
    #Mean.Ves=Ves.coef[round(length(Ves.coef)/2)]
    #Mean.Ves=substr(names(Mean.Ves),10,12)
    #Mean.Ves=factor(Mean.Ves,levels=levels(DAT$Fisher.ID))
    
    ##Fish.samp.WC=rep(Mean.Ves,missing.WC,replace=T)      #NEW
    ##Fish.samp.Zn1=rep(Mean.Ves,missing.zn1,replace=T)
    ##Fish.samp.Zn2=rep(Mean.Ves,missing.zn2,replace=T)
    
    
    #Fish.samp.WC=sample(Rand.Cof.WC$Fisher.ID,missing.WC,replace=T)      #NEW
    #Fish.samp.Zn1=sample(Rand.Cof.Zn1$Fisher.ID,missing.zn1,replace=T)
    #Fish.samp.Zn2=sample(Rand.Cof.Zn2$Fisher.ID,missing.zn2,replace=T)
    
    
    #Dat.WC=subset(Data.boot.non.sampled,Fishing.region=="WC",select=c(Fisher.ID))
    #Dat.Zn1=subset(Data.boot.non.sampled,Fishing.region=="SS1",select=c(Fisher.ID))
    #Dat.Zn2=subset(Data.boot.non.sampled,Fishing.region=="SS2",select=c(Fisher.ID))
    
    #Dat.WC$Fisher.ID.new=Fish.samp.WC
    #Dat.Zn1$Fisher.ID.new=Fish.samp.Zn1
    #Dat.Zn2$Fisher.ID.new=Fish.samp.Zn2
    
    ##  Rand.Coef.miss=data.frame(Fisher.ID=unique(Data.boot.non.sampled$Fisher.ID),Rand.coef=Fish.samp)
    #Rand.Coef.miss=rbind(Dat.WC,Dat.Zn1,Dat.Zn2)
    
    
    
    #Data.boot.non.sampled$Catch.boot=0
    #DAT.pred=Data.boot.non.sampled[,match(c("Fisher.ID","Fishing.region","Time","Effort","Catch.boot"),
    #                                      names(Data.boot.non.sampled))]
    #DAT.pred=DAT.pred[order(DAT.pred$Fisher.ID),]
    
    ##datnew=DAT.pred[,match(c("Fishing.region","Time","Effort","Catch.boot"),names(DAT.pred))]
    ##mm = model.matrix(terms(fit),DAT.pred)
    
    #DAT.pred=merge(DAT.pred,Rand.Coef.miss,by="Fisher.ID")
    ##DAT.pred=DAT.pred[order(DAT.pred$Fisher.ID),]
    ##DAT.pred$Pred = mm %*% fixef(fit)
    ##DAT.pred$Pred=exp(DAT.pred$Rand.coef+DAT.pred$Pred)
    
    #DAT.pred$Fisher.ID=DAT.pred$Fisher.ID.new
    #NEWDATA=DAT.pred[,match(c("Fisher.ID","Fishing.region","Time","Effort"),names(DAT.pred))]
    #DAT.pred$Pred=round(predict(fit,newdata=NEWDATA,type="response"))
    
    #3. store stuff
    Dat.list[[j]]=DAT
    GLM.list[[j]]=fit
    #Dat.pred[[j]]=DAT.pred
  }, error = function(e) {
    
  })
})
#NOT USED

# #Check fitted and predicted
# par(mfcol=c(4,4))
# for (i in 1:16)
# {
#   DAT=Dat.list[[i]]
#   plot(DAT$Catch.boot,ylim=c(0,91))
#   points(DAT$fitted,col=2)
#   points(DAT$Pred,col=3,pch=19,cex=.9)
# }


  #3.4 Put bootstrapped results in nice shape

  #3.4.1 Predictions
This.stuff=c("Fishing.region","Time","Effort")
ID=match(This.stuff,names(Data.boot.non.sampled))
source.Dat=Data.boot.non.sampled[,ID]
source.Dat=source.Dat[order(source.Dat$Fishing.region,source.Dat$Time),]  #sort new data to match FisherID order
Store.App1=expand.grid(Fishing.region=unique(source.Dat$Fishing.region),Time=unique(source.Dat$Time))
Store.App1=Store.App1[order(Store.App1$Time,Store.App1$Fishing.region),]
#Store.App1=subset(Store.App1,!(Fishing.region=="SS1" &Time==2))
Store.App1.final=as.data.frame(matrix(,nrow=nrow(Store.App1),ncol=N.boots))
Store.App1.final=cbind(Store.App1,Store.App1.final)

for(i in 1:N.boots)
{
  dat=Dat.pred[[i]]
  SUM=aggregate(Pred~Fishing.region+Time,dat,sum)
  Store.App1.final[,2+i]=SUM[,3]
}

    #Get mean and 95% CI
Store.App1$Mean.Catch.App1=round(apply(as.matrix(Store.App1.final[3:(N.boots+2)]), 1, mean))
Store.App1$Low95.Catch.App1=round(apply(as.matrix(Store.App1.final[3:(N.boots+2)]), 1, function(x)quantile(x,probs=0.025)))
Store.App1$Up95.Catch.App1=round(apply(as.matrix(Store.App1.final[3:(N.boots+2)]), 1, function(x)quantile(x,probs=0.975)))
Pred.Catch.Missing.Fish.App1=Store.App1


  #3.4.2 Bootstrapped observations
    #Get mean and 95% CI for observed (boostrapped) data
Store.obs.App1=expand.grid(Fishing.region=unique(source.Dat$Fishing.region),Time=unique(source.Dat$Time))
Store.obs.App1=Store.obs.App1[order(Store.obs.App1$Time,Store.obs.App1$Fishing.region),]

#first sum by region-period
Ag.sum=matrix(nrow=nrow(Store.obs.App1),ncol=N.boots)
for(i in 1:N.boots)
{
  dat=Dat.list[[i]]
  SUM=aggregate(Catch.boot~Fishing.region+Time,dat,sum)
  if(length(SUM[,3])==length(Ag.sum[,i]))Ag.sum[,i]=SUM[,3]
}

#then get mean and CI
Store.obs.App1$Mean.Catch.App1=round(apply(Ag.sum, 1, mean,na.rm=T))
Store.obs.App1$Low95.Catch.App1=round(apply(Ag.sum, 1, function(x)quantile(x,probs=0.025,na.rm=T)))
Store.obs.App1$Up95.Catch.App1=round(apply(Ag.sum, 1, function(x)quantile(x,probs=0.975,na.rm=T)))
Obs.Catch.Missing.Fish.App1=Store.obs.App1


#EXPORT STUFF

#mean observed
write.table(Obs.Catch.Missing.Fish.App1,"Final.catch.obs.csv",sep=",",row.names = F)

#mean predicted for the unobserved
write.table(Pred.Catch.Missing.Fish.App1,"Final.catch.missing.csv",sep=",",row.names = F)




#plotting numbers vs period
RS.missing=reshape(Pred.Catch.Missing.Fish.App1, v.names = c("Low95.Catch.App1","Up95.Catch.App1","Mean.Catch.App1"), idvar = "Time",
                   timevar = "Fishing.region", direction = "wide")
RS.missing[2,2:4]=0

RS.obs=reshape(Obs.Catch.Missing.Fish.App1, v.names = c("Low95.Catch.App1","Up95.Catch.App1","Mean.Catch.App1"), idvar = "Time",
               timevar = "Fishing.region", direction = "wide")


tiff(file="Figure2AllFishers.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    #create tiff



#epsilon <- 0.02
fn.CI=function(x,y,z,jit,COL)
{
  segments(x,y,x,z,col=COL)
  segments(x,y,x,y,col=COL)
  segments(x,z,x,z,col=COL)
}

par(mai=c(1,1,.1,.1),oma=c(.1,.1,.1,.1),las=1,mgp=c(3,.5,0))
#mean missing
RS.missing$Time=as.numeric(as.character(RS.missing$Time))
RS.obs$Time=as.numeric(as.character(RS.obs$Time))
plot(RS.missing$Time,RS.missing$Mean.Catch.App1.SS1,ylim=c(0,750),xlab="Time period",
     ylab="Mean number caught (+/- 95% CI)",cex.lab=2,xaxt="n",cex=1.75,pch=15,
     xlim=c(0.9,3.2),cex.axis=1.5)
points(RS.missing$Time*(1+0.05),RS.missing$Mean.Catch.App1.SS2,col=1,cex=1.75,pch=17)
points(RS.missing$Time*(1-0.05),RS.missing$Mean.Catch.App1.WC,col=1,cex=1.75,pch=16)


#CI missing
fn.CI(RS.missing$Time,RS.missing$Low95.Catch.App1.SS1,RS.missing$Up95.Catch.App1.SS1,1,1)
fn.CI(RS.missing$Time*(1+0.05),RS.missing$Low95.Catch.App1.SS2,RS.missing$Up95.Catch.App1.SS2,1.03,1)
fn.CI(RS.missing$Time*(1-0.05),RS.missing$Low95.Catch.App1.WC,RS.missing$Up95.Catch.App1.WC,0.97,1)


#mean observed
points(RS.obs$Time*(1+0.01),RS.obs$Mean.Catch.App1.SS1,col=1,cex=1.75,pch=22,bg="white")
points(RS.obs$Time*(1+0.07),RS.obs$Mean.Catch.App1.SS2,col=1,cex=1.75,pch=24,bg="white")
points(RS.obs$Time*(1-0.06),RS.obs$Mean.Catch.App1.WC,col=1,cex=1.75,pch=21,bg="white")

#CI observed
fn.CI(RS.obs$Time*(1+0.01),RS.obs$Low95.Catch.App1.SS1,RS.obs$Up95.Catch.App1.SS1,1.02,1)
fn.CI(RS.obs$Time*(1+0.07),RS.obs$Low95.Catch.App1.SS2,RS.obs$Up95.Catch.App1.SS2,1.07,1)
fn.CI(RS.obs$Time*(1-0.06),RS.obs$Low95.Catch.App1.WC,RS.obs$Up95.Catch.App1.WC,0.96,1)


axis(1,at=c(1,2,3),labels=c("1988-1996","1997-2004","2005-2013"),tck=-0.02,cex.axis=1.50)

legend("topright",c("inferred, zone 1","inferred, zone 2","inferred, WC",
                    "observed, zone 1","observed, zone 2","observed, WC"),bty='n',pch=c(15,17,16,22,24,21),
       col=1,bg=c(1,1,1,"white","white","white"),cex=1.5,pt.cex=2)
dev.off()













################################################################################################

#PREVIOUS APPROACH (FIXED EFFECTS ONLY)
# for(i in 1:N.boots)
# {
#   dummy=Data.boot
#   fn.boot=function(a,b) round(runif(1,a,b))
#   dummy$Catch.boot=mapply(fn.boot,dummy$Min.per,dummy$Max.per)                                 
#   Dat.list[[i]]=dummy
#   GLM.list[[i]]=glm( Catch.boot ~ Fisher.ID+Time*Fishing.region+
#                      Effort*Fishing.region, data = dummy, family=poisson)
# 
# }

##################################################

#4. Predict unknown catch

###################################################


#4.1. Predict new data

  #Allocate FisherID value to unknown catch data set
Unk.Fisher.ID=Data.gn.pred$Fisher.ID


  #Different approaches for imputing Fisher ID

This.stuff=c("Fishing.region","Time","Effort")
ID=match(This.stuff,names(Data.gn.pred))

source.Dat=Data.gn.pred[,ID]
source.Dat=source.Dat[order(source.Dat$Fishing.region,source.Dat$Time),]  #sort new data to match FisherID order

    #Interviewed fishers
RegW.Per1=subset(Data.gn,Fishing.region=="WC"&Time==1)
RegW.Per2=subset(Data.gn,Fishing.region=="WC"&Time==2)
RegW.Per3=subset(Data.gn,Fishing.region=="WC"&Time==3)

RegSS1.Per1=subset(Data.gn,Fishing.region=="SS1"&Time==1)
RegSS1.Per2=subset(Data.gn,Fishing.region=="SS1"&Time==2)
RegSS1.Per3=subset(Data.gn,Fishing.region=="SS1"&Time==3)

RegSS2.Per1=subset(Data.gn,Fishing.region=="SS2"&Time==1)
RegSS2.Per2=subset(Data.gn,Fishing.region=="SS2"&Time==2)
RegSS2.Per3=subset(Data.gn,Fishing.region=="SS2"&Time==3)


    #Unobserved
Un.RegW.Per1=subset(Data.gn.pred,Fishing.region=="WC"&Time==1)
Un.RegW.Per2=subset(Data.gn.pred,Fishing.region=="WC"&Time==2)
Un.RegW.Per3=subset(Data.gn.pred,Fishing.region=="WC"&Time==3)

Un.RegSS1.Per1=subset(Data.gn.pred,Fishing.region=="SS1"&Time==1)
Un.RegSS1.Per2=subset(Data.gn.pred,Fishing.region=="SS1"&Time==2)
Un.RegSS1.Per3=subset(Data.gn.pred,Fishing.region=="SS1"&Time==3)

Un.RegSS2.Per1=subset(Data.gn.pred,Fishing.region=="SS2"&Time==1)
Un.RegSS2.Per2=subset(Data.gn.pred,Fishing.region=="SS2"&Time==2)
Un.RegSS2.Per3=subset(Data.gn.pred,Fishing.region=="SS2"&Time==3)


  #Create observed FisherID vectors by period and region
fn.levels=function(dat)id=(as.numeric(as.character(dat$Fisher.ID)))
Ves.W.1=fn.levels(RegW.Per1)
Ves.W.2=fn.levels(RegW.Per2)
Ves.W.3=fn.levels(RegW.Per3)
Ves.SS1.1=fn.levels(RegSS1.Per1)
Ves.SS1.2=fn.levels(RegSS1.Per2)
Ves.SS1.3=fn.levels(RegSS1.Per3)
Ves.SS2.1=fn.levels(RegSS2.Per1)
Ves.SS2.2=fn.levels(RegSS2.Per2)
Ves.SS2.3=fn.levels(RegSS2.Per3)


  #Create unobserved FisherID vectors by period and region
Un.Ves.W.1=fn.levels(Un.RegW.Per1)
Un.Ves.W.2=fn.levels(Un.RegW.Per2)
Un.Ves.W.3=fn.levels(Un.RegW.Per3)
Un.Ves.SS1.1=fn.levels(Un.RegSS1.Per1)
Un.Ves.SS1.2=fn.levels(Un.RegSS1.Per2)
Un.Ves.SS1.3=fn.levels(Un.RegSS1.Per3)
Un.Ves.SS2.1=fn.levels(Un.RegSS2.Per1)
Un.Ves.SS2.2=fn.levels(Un.RegSS2.Per2)
Un.Ves.SS2.3=fn.levels(Un.RegSS2.Per3)


  #Approach 1. Random selection within period and region

store.App.1=store.App.2=store.App.3=vector('list',length=N.boots)

fn.rand.sel=function(dat,dat1)
{
  N=length(dat)
  samp=sample(dat1,N,replace=T)
}

for (i in 1:N.boots)
{
  Dat.1=source.Dat
  New.FID.W.1=fn.rand.sel(Un.Ves.W.1,Ves.W.1)
  New.FID.W.2=fn.rand.sel(Un.Ves.W.2,Ves.W.2)
  New.FID.W.3=fn.rand.sel(Un.Ves.W.3,Ves.W.3)
  
  New.FID.SS1.1=fn.rand.sel(Un.Ves.SS1.1,Ves.SS1.1)
  New.FID.SS1.2=fn.rand.sel(Un.Ves.SS1.2,Ves.SS1.2)
  New.FID.SS1.3=fn.rand.sel(Un.Ves.SS1.3,Ves.SS1.3)
  
  New.FID.SS2.1=fn.rand.sel(Un.Ves.SS2.1,Ves.SS2.1)
  New.FID.SS2.2=fn.rand.sel(Un.Ves.SS2.2,Ves.SS2.2)
  New.FID.SS2.3=fn.rand.sel(Un.Ves.SS2.3,Ves.SS2.3)
  
  #convert to factor and attach new fisher id to data 
  ID=c(New.FID.SS1.1,New.FID.SS1.2,New.FID.SS1.3,New.FID.SS2.1,New.FID.SS2.2,New.FID.SS2.3,
       New.FID.W.1,New.FID.W.2,New.FID.W.3)
  ID=factor(ID,levels=sort(unique(ID)))
  Dat.1$Fisher.ID=ID
  
  #predict new data
  PRED.Min=predict(GLM.list[[i]],newdata=Dat.1,type="response",se.fit =T)
  Dat.1$Pred.Min.Catch=PRED.Min$fit
  store.App.1[[i]]=Dat.1
}



####################################################################
#   #Approaches 2 and 3
# dummy=0
# names(dummy)="1"
# Coefs.min=coef(slm3.min)[2:25]
# Coefs.min=c(dummy,Coefs.min)
# Coefs.max=coef(slm4)[2:25]
# Coefs.max=c(dummy,Coefs.max)
# 
# names(Coefs.min)=as.numeric(gsub("\\D", "",names(Coefs.min)))
# names(Coefs.max)=as.numeric(gsub("\\D", "",names(Coefs.max)))
# 
# fun.crap=function(DAT,WHAT)
# {
#   id.min=match(as.character(DAT), names(Coefs.min))
#   id.max=match(as.character(DAT), names(Coefs.max))
#   
#   if(WHAT=="min")
#   {
#     FisherID.min=names(sort(Coefs.min[id.min])[1])
#     FisherID.max=names(sort(Coefs.max[id.max])[1])
#   }
#   
#   if(WHAT=="max")
#   {
#     FisherID.min=names(rev(sort(Coefs.min[id.min]))[1])
#     FisherID.max=names(rev(sort(Coefs.max[id.max]))[1])
#   }
#   return(data.frame(FisherID.min=as.numeric(FisherID.min),FisherID.max=as.numeric(FisherID.max)))
# }
# 
# fn.min.max.sel=function(dat,dat1)
# {
#   N=length(dat)
#   Dat=data.frame(FID.min=rep(dat1[,1],N),FID.max=rep(dat1[,2],N))
#   return(Dat)
# }
# 
# 
#     #--Min within period and region---
# #1. extract coefficients
# MIN.FiID.W.1=fun.crap(RegW.Per1$Fisher.ID,"min")
# MIN.FiID.W.2=fun.crap(RegW.Per2$Fisher.ID,"min")
# MIN.FiID.W.3=fun.crap(RegW.Per3$Fisher.ID,"min")
# 
# MIN.FiID.SS1.1=fun.crap(RegSS1.Per1$Fisher.ID,"min")
# MIN.FiID.SS1.2=fun.crap(RegSS1.Per2$Fisher.ID,"min")
# MIN.FiID.SS1.3=fun.crap(RegSS1.Per3$Fisher.ID,"min")
# 
# MIN.FiID.SS2.1=fun.crap(RegSS2.Per1$Fisher.ID,"min")
# MIN.FiID.SS2.2=fun.crap(RegSS2.Per2$Fisher.ID,"min")
# MIN.FiID.SS2.3=fun.crap(RegSS2.Per3$Fisher.ID,"min")
# 
# 
# #2. replace FisherID
# New.FID.W.1=fn.min.max.sel(Un.Ves.W.1,MIN.FiID.W.1)
# New.FID.W.2=fn.min.max.sel(Un.Ves.W.2,MIN.FiID.W.2)
# New.FID.W.3=fn.min.max.sel(Un.Ves.W.3,MIN.FiID.W.3)
# 
# New.FID.SS1.1=fn.min.max.sel(Un.Ves.SS1.1,MIN.FiID.SS1.1)
# New.FID.SS1.2=fn.min.max.sel(Un.Ves.SS1.2,MIN.FiID.SS1.2)
# New.FID.SS1.3=fn.min.max.sel(Un.Ves.SS1.3,MIN.FiID.SS1.3)
# 
# 
# New.FID.SS2.1=fn.min.max.sel(Un.Ves.SS2.1,MIN.FiID.SS2.1)
# New.FID.SS2.2=fn.min.max.sel(Un.Ves.SS2.2,MIN.FiID.SS2.2)
# New.FID.SS2.3=fn.min.max.sel(Un.Ves.SS2.3,MIN.FiID.SS2.3)
# 
# 
#   #3. predict min catch new data
#   
#     #3.1. convert to factor and attach new fisher id to data 
# ID.min=c(New.FID.SS1.1[,1],New.FID.SS1.2[,1],New.FID.SS1.3[,1],
#            New.FID.SS2.1[,1],New.FID.SS2.2[,1],New.FID.SS2.3[,1],
#            New.FID.W.1[,1],New.FID.W.2[,1],New.FID.W.3[,1])
# ID.min=factor(ID.min,levels=sort(unique(ID.min)))
# Dat.2=source.Dat
# Dat.2$Fisher.ID=ID.min
# 
#     #3.2. prediction
# PRED.Min=predict(slm3.min,newdata=Dat.2,type="response",se.fit =T)
# Dat.2$Pred.Min.Catch=PRED.Min$fit
# Dat.2$Pred.Min.Catch.SE=PRED.Min$se.fit
# 
# 
#   #4. predict max catch new data
# 
#     #4.1. convert to factor and attach new fisher id to data 
# ID.max=c(New.FID.SS1.1[,2],New.FID.SS1.2[,2],New.FID.SS1.3[,2],
#          New.FID.SS2.1[,2],New.FID.SS2.2[,2],New.FID.SS2.3[,2],
#          New.FID.W.1[,2],New.FID.W.2[,2],New.FID.W.3[,2])
# ID.max=factor(ID.max,levels=sort(unique(ID.max)))
# Dat.2$Fisher.ID=ID.max
# 
#     #4.2. prediction
# PRED.Max=predict(slm4,newdata=Dat.2,type="response",se.fit =T)
# Dat.2$Pred.Max.Catch=PRED.Max$fit
# Dat.2$Pred.Max.Catch.SE=PRED.Max$se.fit
# 
# 
# 
# 
#   
#   #---Max within period and region----
#   
# #1. extract coefficients
# MAX.FiID.W.1=fun.crap(RegW.Per1$Fisher.ID,"max")
# MAX.FiID.W.2=fun.crap(RegW.Per2$Fisher.ID,"max")
# MAX.FiID.W.3=fun.crap(RegW.Per3$Fisher.ID,"max")
# 
# MAX.FiID.SS1.1=fun.crap(RegSS1.Per1$Fisher.ID,"max")
# MAX.FiID.SS1.2=fun.crap(RegSS1.Per2$Fisher.ID,"max")
# MAX.FiID.SS1.3=fun.crap(RegSS1.Per3$Fisher.ID,"max")
# 
# MAX.FiID.SS2.1=fun.crap(RegSS2.Per1$Fisher.ID,"max")
# MAX.FiID.SS2.2=fun.crap(RegSS2.Per2$Fisher.ID,"max")
# MAX.FiID.SS2.3=fun.crap(RegSS2.Per3$Fisher.ID,"max")
# 
# 
# #2. replace FisherID
# New.FID.W.1=fn.min.max.sel(Un.Ves.W.1,MAX.FiID.W.1)
# New.FID.W.2=fn.min.max.sel(Un.Ves.W.2,MAX.FiID.W.2)
# New.FID.W.3=fn.min.max.sel(Un.Ves.W.3,MAX.FiID.W.3)
# 
# New.FID.SS1.1=fn.min.max.sel(Un.Ves.SS1.1,MAX.FiID.SS1.1)
# New.FID.SS1.2=fn.min.max.sel(Un.Ves.SS1.2,MAX.FiID.SS1.2)
# New.FID.SS1.3=fn.min.max.sel(Un.Ves.SS1.3,MAX.FiID.SS1.3)
# 
# 
# New.FID.SS2.1=fn.min.max.sel(Un.Ves.SS2.1,MAX.FiID.SS2.1)
# New.FID.SS2.2=fn.min.max.sel(Un.Ves.SS2.2,MAX.FiID.SS2.2)
# New.FID.SS2.3=fn.min.max.sel(Un.Ves.SS2.3,MAX.FiID.SS2.3)
# 
# 
# #3. predict min catch new data
# 
# #3.1. convert to factor and attach new fisher id to data 
# ID.min=c(New.FID.SS1.1[,1],New.FID.SS1.2[,1],New.FID.SS1.3[,1],
#          New.FID.SS2.1[,1],New.FID.SS2.2[,1],New.FID.SS2.3[,1],
#          New.FID.W.1[,1],New.FID.W.2[,1],New.FID.W.3[,1])
# ID.min=factor(ID.min,levels=sort(unique(ID.min)))
# Dat.3=source.Dat
# Dat.3$Fisher.ID=ID.min
# 
# #3.2. prediction
# PRED.Min=predict(slm3.min,newdata=Dat.3,type="response",se.fit =T)
# Dat.3$Pred.Min.Catch=PRED.Min$fit
# Dat.3$Pred.Min.Catch.SE=PRED.Min$se.fit
# 
# 
# #4. predict max catch new data
# 
# #4.1. convert to factor and attach new fisher id to data 
# ID.max=c(New.FID.SS1.1[,2],New.FID.SS1.2[,2],New.FID.SS1.3[,2],
#          New.FID.SS2.1[,2],New.FID.SS2.2[,2],New.FID.SS2.3[,2],
#          New.FID.W.1[,2],New.FID.W.2[,2],New.FID.W.3[,2])
# ID.max=factor(ID.max,levels=sort(unique(ID.max)))
# Dat.3$Fisher.ID=ID.max
# 
# #4.2. prediction
# PRED.Max=predict(slm4,newdata=Dat.3,type="response",se.fit =T)
# Dat.3$Pred.Max.Catch=PRED.Max$fit
# Dat.3$Pred.Max.Catch.SE=PRED.Max$se.fit
# 
#   
#  #Compare predictions
# 
# 
#  plot(Dat.2$Pred.Min.Catch)
#  points(Dat.2$Pred.Max.Catch,col=2,pch=19)
# 
#  plot(Dat.3$Pred.Min.Catch)
#  points(Dat.3$Pred.Max.Catch,col=2,pch=19)
  
  
  
  
  #Approach 4. Steve's 



  
  
  
  


# ####### Random forest approach  ###############
# RF=randomForest(Min.catch.period ~ Fisher.ID+Time+Effort+Fishing.region, data = Data.gn,
#                 ntree=5000, importance=TRUE)
# 
# plot(RF)
# varImpPlot(RF)
# importance(RF) # importance of each predictor 
# 
# 
# #predictions
# Data.gn$Predicted.min.catch = predict(RF, Data.gn, type="response")
# table(Data.gn$Min.catch.period, Data.gn$Predicted.min.catch)
# prop.table(table(Data.gn$Min.catch.period, Data.gn$Predicted.min.catch),1)
# 
# plot(Data.gn$Min.catch.period,Data.gn$Predicted.min.catch,
#      ylim=c(0,max(Data.gn$Min.catch.period)),ylab="Predicted catch", xlab="Observed catch")
# lines(Data.gn$Min.catch.period,Data.gn$Min.catch.period,col=2,lwd=2)
# 
# 
# #Tree
# Tree=tree(Min.catch.period ~ Fisher.ID+Time+Effort+Fishing.region, data = Data.gn)
# 
# 
# #REPORT SECTION