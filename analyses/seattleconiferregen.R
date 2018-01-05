#Analyses for Seattle Confer regeneration manuscript
#Submitted to Urban Ecosystems
#July 2016
#Analyses by Ailene Ettinger
##LPseed trap data
library(lme4)
library(car)
#working directory: 
setwd("~/git/urbantreeregen")
##FIGURES
##First, seed trap figure, comparing urban forests to old growth forests at mt rainier
seeds<-read.csv("data/LPMoraSeedTrap.csv", header=T)
head(seeds)
LPseeds<-seeds[seeds$Site=="LP",]#just lincoln park seed trap data 
LPconseeds<-LPseeds[LPseeds$conifer=="y",]#just conifers at LP
LPconseeds$block=factor(LPconseeds$block)
LPconseedtot<-tapply(LPconseeds$ambientseedfall,LPconseeds$block,sum,na.rm = TRUE)#total number of conifer seeds per block at LP
LPconseedtotperm<-LPconseedtot/0.152###dimensions of LP seed traps=11" W x 21.37" L=.2794mX.542798m=0.1516578m2
tapply(LPconseeds$ambientseedfall,list(LPconseeds$species,LPconseeds$Year),sum,na.rm = TRUE)#total number of conifer seeds per block at LP
tapply(LPconseeds$ambientseedfall,LPconseeds$block,mean,na.rm = TRUE)
#Now mt. rainier old growth seed trap data from low elevation stands
TO04seeds<-seeds[seeds$Site=="TO04",]#just TO04 seed trap data 
TO11seeds<-seeds[seeds$Site=="TO11",]#just TO11 seed trap data 
TA01seeds<-seeds[seeds$Site=="TA01",]#just TA01 seed trap data 
MORAconseeds<-rbind(TO04seeds,TO11seeds,TA01seeds)#all MORA seeds
MORAconseeds$Siteblock <- paste(MORAconseeds$Site,MORAconseeds$block, sep = ".")
MORAconseedtot<-tapply(MORAconseeds$ambientseedfall,MORAconseeds$Siteblock,sum,na.rm = TRUE)#total number of conifer seeds per block at LP
MORAconseedtotperm<-MORAconseedtot/0.176###dimensions of MORA seed traps=0.176-m2 (from KRoiss&HilleRisLambers)
LPconseedmean<-mean(LPconseedtotperm)#average conifer seeds per block for LP
MORAconseedmean<-mean(MORAconseedtotperm)#average conifer seeds per block for MORA
LPconseedsd<-sd(LPconseedtotperm)#sd for LP seeds per block
MORAconseedsd<-sd(MORAconseedtotperm)#sd for MORA seeds per block
LPconseedn<-length(LPconseedtotperm)#n for LP blocks
MORAconseedn<-length(MORAconseedtotperm)#n for MORA blocks
LPconseedse<-LPconseedsd/sqrt(LPconseedn)#std error for LP mean seeds per block- use to make) error bars
MORAconseedse<-MORAconseedsd/sqrt(MORAconseedn)#std error for MORA mean seeds per block- use to make error bar
quartz(height=4,width=3)
bar<-barplot(c(LPconseedmean,MORAconseedmean), col=c("gray"), names.arg=c("Urban", "Old-growth"), ylab="Number of conifer seeds per m2",ylim=c(0,500), cex.lab=1.1, cex.axis=1.1, cex.names=1.1)
arrows(c(bar),c(LPconseedmean-LPconseedse,MORAconseedmean-MORAconseedse),c(bar),c(LPconseedmean+LPconseedse,MORAconseedmean+MORAconseedse),angle=90, code=0)
test=glmer(ambientseedfall~foresttype+(1|block),data=seeds,family=poisson)
Anova(test)
summary(test)
test.nb=glmer.nb(ambientseedfall~-1+foresttype+(1|block),data=seeds)
Anova(test.nb)# lower residual deviance and lower AIC, so use this model
summary(test.nb)
#Germination at Lincoln Park, Seattle 
data<-read.csv("data/LPgerm20112012.csv")
dim(data)#384 rows,  12 columns
data$year=as.factor(data$year)
data$block=as.factor(data$block)
data$seedsadd=NA
data[which(data$numseeds==0),]$seedsadd="n"
data[which(data$numseeds!=0),]$seedsadd="y"
abgrdat<-data[data$species=="ABGR",]
tshedat<-data[data$species=="TSHE",]
thpldat<-data[data$species=="THPL",]
seedadddata<-data[data$numseeds>0,]
controldata<-data[data$numseeds=="0",]
tapply(controldata$numgerm,list(controldata$year,controldata$species),sum,na.rm = T)
tapply(controldata$numgerm,list(controldata$block,controldata$plot,controldata$year),sum,na.rm = T)
tapply(seedadddata$numgerm,list(seedadddata$block,seedadddata$plot,seedadddata$year),sum,na.rm = T)

#try with data from both years, with year as a random effect
#try with all data- seeds added and not added so that one full model can be used
#Fit Poisson model for total germination (not surviving germinants at the end of the season)
#now fit a GLM with the following explanatory variables: seeds added, ivy presence,wood chips presence, and their interaction
abgrtest<-glmer(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=poisson,data=abgrdat) 
Anova(abgrtest,test="Chisq")#ivy only significant treatment
summary(abgrtest)#ivy presence marginally reduced total germination; seeds eadded had big effect
##fit ZIP for these models?
#Now for surviving germinants (at the end of the season), fit a GLM with the following explanatory variables: ivy presence,wood chips presence, and their interaction
y=cbind(abgrdat$survivinggerm,(abgrdat$numgerm-abgrdat$survivinggerm))
abgrsurvtest<-glmer(y~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=binomial,data=abgrdat,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
Anova(abgrsurvtest,test="Chisq")#ivy & woodchips significant
summary(abgrsurvtest)#significantly higher germination with wood chips; lower germination with ivy present
###
###Try ZIP model
library(pscl)
library(boot)
library(bbmle)
abgrtest.zip = zeroinfl(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy), data=abgrdat)#both poisson and binomial processes depend on seeds added treatment
summary(abgrtest.zip)
abgrtest.zip2 = zeroinfl(numgerm ~as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=abgrdat)#both poisson and binomial processes depend on seeds added treatment
summary(abgrtest.zip2)
abgrtest.zip3 = zeroinfl(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=abgrdat)#both poisson and binomial processes depend treatments
summary(abgrtest.zip3)
abgrtest.zip4 = zeroinfl(numgerm ~year+block+as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd)+year, data=abgrdat)#both poisson and binomial processes depend treatments
summary(abgrtest.zip4)
abgrtest.zip5 = zeroinfl(numgerm ~year+as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd)+year, data=abgrdat)#both poisson and binomial processes depend treatments
summary(abgrtest.zip5)
Anova(abgrtest.zip4)
hist(abgrdat$numgerm)
is.factor(abgrdat$block)
abgrtest<-glmer(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=poisson,data=abgrdat) 
Anova(abgrtest,test="Chisq")#ivy only significant treatment
summary(abgrtest)#ivy presence marginally reduced total germination; seeds eadded had big effect
abgrtest.glm<-glm(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy), family=poisson,data=abgrdat) 
abgrtest.glm2<-glm(numgerm ~block+year, family=poisson,data=abgrdat) 
summary(abgrtest.glm2)
Anova(abgrtest.glm2)
ICtab(abgrtest.zip, abgrtest,abgrtest.zip2,abgrtest.zip3,abgrtest.zip4,abgrtest.zip5,abgrtest.glm,abgrtest.glm2)#now fit a GLM with the following explanatory variables: seeds added, ivy presence,wood chips presence, and their interaction
y=cbind(abgrdat$survivinggerm,(abgrdat$numgerm-abgrdat$survivinggerm))
thplsurvtest<-glmer(y~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=binomial,data=thpldat, control=glmerControl(optCtrl=list(maxfun=100000))) 
confint(abgrtest.zip4)
###thpl
thpltest<-glmer(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=poisson,data=thpldat) 
Anova(thpltest,test="Chisq")#no significant treatmentd
summary(thpltest)#ivy presence marginally reduced total germination; seeds eadded had big effect
#Now for surviving germinants (at the end of the season)- can't fit these for tHPL
#now fit a GLM with the following explanatory variables: ivy presence,wood chips presence, and their interaction
y=cbind(thpldat$survivinggerm,(thpldat$numgerm-thpldat$survivinggerm))
thplsurvtest<-glmer(y~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=binomial,data=thpldat,control=glmerControl(optCtrl=list(maxfun=100000))) 

thpltest.zip3 = zeroinfl(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=thpldat)#both poisson and binomial processes depend treatments
summary(thpltest.zip3)

###tshe
tshetest<-glmer(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)+(1|block)+(1|year), family=poisson,data=tshedat) 
Anova(tshetest,test="Chisq")#ivy only significant treatment
summary(tshetest)#ivy presence marginally reduced total germination; seeds eadded had big effect
#Now for surviving germinants (at the end of the season)- can't fit these for tHPL
#now fit a GLM with the following explanatory variables: ivy presence,wood chips presence, and their interaction
y=cbind(tshedat$survivinggerm,(tshedat$numgerm-tshedat$survivinggerm))
tshesurvtest<-glm(y~seedsadd+woodchips*ivy,family=binomial,data=tshedat) 
Anova(tshesurvtest,test="Wald")#ivy & woodchips significant
summary(tshesurvtest)#significantly higher germination with wood chips; lower germination with ivy present
tshetest.zip = zeroinfl(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshetest.zip)
tshetest.zip2 = zeroinfl(numgerm ~as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshetest.zip2)
tshetest.zip3 = zeroinfl(numgerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshetest.zip3)
ICtab(tshetest.zip, tshetest,tshetest.zip2,tshetest.zip3)#now fit a GLM with the following explanatory variables: seeds added, ivy presence,wood chips presence, and their interaction
#now for surviving germs
tshesurvtest.zip = zeroinfl(survivinggerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshesurvtest.zip)
tshesurvtest.zip2 = zeroinfl(survivinggerm ~as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshesurvtest.zip2)
tshetest.zip3 = zeroinfl(survivinggerm ~as.factor(seedsadd)+as.factor(woodchips)*as.factor(ivy)|as.factor(seedsadd), data=tshedat)#both poisson and binomial processes depend treatments
summary(tshesurvtest.zip3)
ICtab(tshesurvtest.zip, tshesurvtest,tshesurvtest.zip2,tshesurvtest.zip3)#now fit a GLM with the following explanatory variables: seeds added, ivy presence,wood chips presence, and their interaction
##Zip doesn't work either, so i stuck with the glmms for now.

##Need to figure out models a bit more...
###Try doing N mixture model or ZIP model?

##summarize the data into a table, to make figures
##ABGR plots with seeds added:
addabgrdat=seedadddata[seedadddata$species=="ABGR",]
meanabgrgerm<-tapply(addabgrdat$numgerm,addabgrdat$treatment,mean, na.rm = TRUE)
meanabgrgermyear<-tapply(addabgrdat$numgerm,list(addabgrdat$treatment,addabgrdat$year),mean, na.rm = TRUE)
sdabgrgerm<-tapply(addabgrdat$numgerm,addabgrdat$treatment,sd)# to get the std deviation
sdabgrgermyear<-tapply(addabgrdat$numgerm,list(addabgrdat$treatment,addabgrdat$year),sd)# to get the std deviation

#and then
nabgrgerm<-tapply(addabgrdat$numgerm,addabgrdat$treatment,length)#the number of blocks for each treatment
nabgrgermyear<-tapply(addabgrdat$numgerm,list(addabgrdat$treatment,addabgrdat$year),length)# to get the std deviation

seabgrgerm<-sdabgrgerm/sqrt(nabgrgerm) #
seabgrgermyear<-sdabgrgermyear/sqrt(nabgrgermyear) #

#now, for controlplots with no seeds added:
abgrcontroldata<-controldata[controldata$species=="ABGR",]
meanabgrgermcontrol<-tapply(abgrcontroldata$numgerm,abgrcontroldata$treatment,mean, na.rm = TRUE)
meanabgrgermcontrolyear<-tapply(abgrcontroldata$numgerm,list(abgrcontroldata$treatment,abgrcontroldata$year),mean, na.rm = TRUE)

sdabgrgermcontrol<-tapply(abgrcontroldata$numgerm,abgrcontroldata$treatment,sd)# to get the std deviation
sdabgrgermcontrolyear<-tapply(abgrcontroldata$numgerm,list(abgrcontroldata$treatment,abgrcontroldata$year),sd, na.rm = TRUE)

#and then
nabgrgermcontrol<-tapply(abgrcontroldata$numgerm,abgrcontroldata$treatment,length)#the number of blocks for each treatment
nabgrgermcontrolyear<-tapply(abgrcontroldata$numgerm,list(abgrcontroldata$treatment,abgrcontroldata$year),length)

seabgrgermcontrol<-sdabgrgermcontrol/sqrt(nabgrgermcontrol) #
seabgrgermcontrolyear<-sdabgrgermcontrolyear/sqrt(nabgrgermcontrolyear) #

#now, for # of germinants surviving at the end of the season:
meanabgrsurvgerm<-tapply(addabgrdat$survivinggerm,addabgrdat$treatment,mean, na.rm = TRUE)
meanabgrsurvgermyear<-tapply(addabgrdat$survivinggerm,list(addabgrdat$treatment,addabgrdat$year),mean, na.rm = TRUE)
sdabgrsurvgerm<-tapply(addabgrdat$survivinggerm,addabgrdat$treatment,sd)# to get the std deviation
sdabgrsurvgermyear<-tapply(addabgrdat$survivinggerm,list(addabgrdat$treatment,addabgrdat$year),sd)# to get the std deviation

#and then
nabgrsurvgerm<-tapply(addabgrdat$survivinggerm,addabgrdat$treatment,length)#the number of blocks for each treatment
nabgrsurvgermyear<-tapply(addabgrdat$survivinggerm,list(addabgrdat$treatment,addabgrdat$year),length)# to get the std deviation

seabgrsurvgerm<-sdabgrsurvgerm/sqrt(nabgrsurvgerm) #
seabgrsurvgermyear<-sdabgrsurvgermyear/sqrt(nabgrsurvgermyear) #

#now, for controlplots with no seeds added:
meanabgrsurvgermcontrol<-tapply(abgrcontroldata$survivinggerm,abgrcontroldata$treatment,mean, na.rm = TRUE)
meanabgrsurvgermcontrolyear<-tapply(abgrcontroldata$survivinggerm,list(abgrcontroldata$treatment,abgrcontroldata$year),mean, na.rm = TRUE)

sdabgrsurvgermcontrol<-tapply(abgrcontroldata$survivinggerm,abgrcontroldata$treatment,sd)# to get the std deviation
sdabgrsurvgermcontrolyear<-tapply(abgrcontroldata$survivinggerm,list(abgrcontroldata$treatment,abgrcontroldata$year),sd, na.rm = TRUE)

#and then
nabgrsurvgermcontrol<-tapply(abgrcontroldata$survivinggerm,abgrcontroldata$treatment,length)#the number of blocks for each treatment
nabgrsurvgermcontrolyear<-tapply(abgrcontroldata$survivinggerm,list(abgrcontroldata$treatment,abgrcontroldata$year),length)

seabgrsurvgermcontrol<-sdabgrsurvgermcontrol/sqrt(nabgrsurvgermcontrol) #
seabgrsurvgermcontrolyear<-sdabgrsurvgermcontrolyear/sqrt(nabgrsurvgermcontrolyear) #

##THPL plots with seeds added:
addthpldat=seedadddata[seedadddata$species=="THPL",]
meanthplgerm<-tapply(addthpldat$numgerm,addthpldat$treatment,mean, na.rm = TRUE)
meanthplgermyear<-tapply(addthpldat$numgerm,list(addthpldat$treatment,addthpldat$year),mean, na.rm = TRUE)
sdthplgerm<-tapply(addthpldat$numgerm,addthpldat$treatment,sd)# to get the std deviation
sdthplgermyear<-tapply(addthpldat$numgerm,list(addthpldat$treatment,addthpldat$year),sd)# to get the std deviation

#and then
nthplgerm<-tapply(addthpldat$numgerm,addthpldat$treatment,length)#the number of blocks for each treatment
nthplgermyear<-tapply(addthpldat$numgerm,list(addthpldat$treatment,addthpldat$year),length)# to get the std deviation

sethplgerm<-sdthplgerm/sqrt(nthplgerm) #
sethplgermyear<-sdthplgermyear/sqrt(nthplgermyear) #

#now, for controlplots with no seeds added:
thplcontroldata<-controldata[controldata$species=="THPL",]
meanthplgermcontrol<-tapply(thplcontroldata$numgerm,thplcontroldata$treatment,mean, na.rm = TRUE)
meanthplgermcontrolyear<-tapply(thplcontroldata$numgerm,list(thplcontroldata$treatment,thplcontroldata$year),mean, na.rm = TRUE)

sdthplgermcontrol<-tapply(thplcontroldata$numgerm,thplcontroldata$treatment,sd)# to get the std deviation
sdthplgermcontrolyear<-tapply(thplcontroldata$numgerm,list(thplcontroldata$treatment,thplcontroldata$year),sd, na.rm = TRUE)

#and then
nthplgermcontrol<-tapply(thplcontroldata$numgerm,thplcontroldata$treatment,length)#the number of blocks for each treatment
nthplgermcontrolyear<-tapply(thplcontroldata$numgerm,list(thplcontroldata$treatment,thplcontroldata$year),length)

sethplgermcontrol<-sdthplgermcontrol/sqrt(nthplgermcontrol) #
sethplgermcontrolyear<-sdthplgermcontrolyear/sqrt(nthplgermcontrolyear) #

#now, for # of germinants surviving at the end of the season:
meanthplsurvgerm<-tapply(addthpldat$survivinggerm,addthpldat$treatment,mean, na.rm = TRUE)
meanthplsurvgermyear<-tapply(addthpldat$survivinggerm,list(addthpldat$treatment,addthpldat$year),mean, na.rm = TRUE)
sdthplsurvgerm<-tapply(addthpldat$survivinggerm,addthpldat$treatment,sd)# to get the std deviation
sdthplsurvgermyear<-tapply(addthpldat$survivinggerm,list(addthpldat$treatment,addthpldat$year),sd)# to get the std deviation

#and then
nthplsurvgerm<-tapply(addthpldat$survivinggerm,addthpldat$treatment,length)#the number of blocks for each treatment
nthplsurvgermyear<-tapply(addthpldat$survivinggerm,list(addthpldat$treatment,addthpldat$year),length)# to get the std deviation

sethplsurvgerm<-sdthplsurvgerm/sqrt(nthplsurvgerm) #
sethplsurvgermyear<-sdthplsurvgermyear/sqrt(nthplsurvgermyear) #

#now, for controlplots with no seeds added:
meanthplsurvgermcontrol<-tapply(thplcontroldata$survivinggerm,thplcontroldata$treatment,mean, na.rm = TRUE)
meanthplsurvgermcontrolyear<-tapply(thplcontroldata$survivinggerm,list(thplcontroldata$treatment,thplcontroldata$year),mean, na.rm = TRUE)

sdthplsurvgermcontrol<-tapply(thplcontroldata$survivinggerm,thplcontroldata$treatment,sd)# to get the std deviation
sdthplsurvgermcontrolyear<-tapply(thplcontroldata$survivinggerm,list(thplcontroldata$treatment,thplcontroldata$year),sd, na.rm = TRUE)

#and then
nthplsurvgermcontrol<-tapply(thplcontroldata$survivinggerm,thplcontroldata$treatment,length)#the number of blocks for each treatment
nthplsurvgermcontrolyear<-tapply(thplcontroldata$survivinggerm,list(thplcontroldata$treatment,thplcontroldata$year),length)

sethplsurvgermcontrol<-sdthplsurvgermcontrol/sqrt(nthplsurvgermcontrol) #
sethplsurvgermcontrolyear<-sdthplsurvgermcontrolyear/sqrt(nthplsurvgermcontrolyear) #
##TSHE plots with seeds added:
addtshedat=seedadddata[seedadddata$species=="TSHE",]
meantshegerm<-tapply(addtshedat$numgerm,addtshedat$treatment,mean, na.rm = TRUE)
meantshegermyear<-tapply(addtshedat$numgerm,list(addtshedat$treatment,addtshedat$year),mean, na.rm = TRUE)
sdtshegerm<-tapply(addtshedat$numgerm,addtshedat$treatment,sd)# to get the std deviation
sdtshegermyear<-tapply(addtshedat$numgerm,list(addtshedat$treatment,addtshedat$year),sd)# to get the std deviation

#and then
ntshegerm<-tapply(addtshedat$numgerm,addtshedat$treatment,length)#the number of blocks for each treatment
ntshegermyear<-tapply(addtshedat$numgerm,list(addtshedat$treatment,addtshedat$year),length)# to get the std deviation

setshegerm<-sdtshegerm/sqrt(ntshegerm) #
setshegermyear<-sdtshegermyear/sqrt(ntshegermyear) #

#now, for controlplots with no seeds added:
tshecontroldata<-controldata[controldata$species=="TSHE",]
meantshegermcontrol<-tapply(tshecontroldata$numgerm,tshecontroldata$treatment,mean, na.rm = TRUE)
meantshegermcontrolyear<-tapply(tshecontroldata$numgerm,list(tshecontroldata$treatment,tshecontroldata$year),mean, na.rm = TRUE)

sdtshegermcontrol<-tapply(tshecontroldata$numgerm,tshecontroldata$treatment,sd)# to get the std deviation
sdtshegermcontrolyear<-tapply(tshecontroldata$numgerm,list(tshecontroldata$treatment,tshecontroldata$year),sd, na.rm = TRUE)

#and then
ntshegermcontrol<-tapply(tshecontroldata$numgerm,tshecontroldata$treatment,length)#the number of blocks for each treatment
ntshegermcontrolyear<-tapply(tshecontroldata$numgerm,list(tshecontroldata$treatment,tshecontroldata$year),length)

setshegermcontrol<-sdtshegermcontrol/sqrt(ntshegermcontrol) #
setshegermcontrolyear<-sdtshegermcontrolyear/sqrt(ntshegermcontrolyear) #

#now, for # of germinants surviving at the end of the season:
meantshesurvgerm<-tapply(addtshedat$survivinggerm,addtshedat$treatment,mean, na.rm = TRUE)
meantshesurvgermyear<-tapply(addtshedat$survivinggerm,list(addtshedat$treatment,addtshedat$year),mean, na.rm = TRUE)
sdtshesurvgerm<-tapply(addtshedat$survivinggerm,addtshedat$treatment,sd)# to get the std deviation
sdtshesurvgermyear<-tapply(addtshedat$survivinggerm,list(addtshedat$treatment,addtshedat$year),sd)# to get the std deviation

#and then
ntshesurvgerm<-tapply(addtshedat$survivinggerm,addtshedat$treatment,length)#the number of blocks for each treatment
ntshesurvgermyear<-tapply(addtshedat$survivinggerm,list(addtshedat$treatment,addtshedat$year),length)# to get the std deviation

setshesurvgerm<-sdtshesurvgerm/sqrt(ntshesurvgerm) #
setshesurvgermyear<-sdtshesurvgermyear/sqrt(ntshesurvgermyear) #

#now, for controlplots with no seeds added:
meantshesurvgermcontrol<-tapply(tshecontroldata$survivinggerm,tshecontroldata$treatment,mean, na.rm = TRUE)
meantshesurvgermcontrolyear<-tapply(tshecontroldata$survivinggerm,list(tshecontroldata$treatment,tshecontroldata$year),mean, na.rm = TRUE)

sdtshesurvgermcontrol<-tapply(tshecontroldata$survivinggerm,tshecontroldata$treatment,sd)# to get the std deviation
sdtshesurvgermcontrolyear<-tapply(tshecontroldata$survivinggerm,list(tshecontroldata$treatment,tshecontroldata$year),sd, na.rm = TRUE)

#and then
ntshesurvgermcontrol<-tapply(tshecontroldata$survivinggerm,tshecontroldata$treatment,length)#the number of blocks for each treatment
ntshesurvgermcontrolyear<-tapply(tshecontroldata$survivinggerm,list(tshecontroldata$treatment,tshecontroldata$year),length)

setshesurvgermcontrol<-sdtshesurvgermcontrol/sqrt(ntshesurvgermcontrol) #
setshesurvgermcontrolyear<-sdtshesurvgermcontrolyear/sqrt(ntshesurvgermcontrolyear) #

##Figure 3, plot of total germination by treatment
quartz(height=7,width=6)
par(mfcol=c(3,1),mai=c(.6,.6,.1,.1), omi=c(.7,.1,.2,.2))
X=c(1,2,3,4)
X2=X+.1
X3=X-.1
#plot control (no seeds added) in white
plot(meanabgrgermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.1,10),type="p",bty="l",pch=21,cex=1.8,las=1,bg=c("white"), cex.main=1.3)
##add seeds added plots in gray

#error bars
for (i in 1:8){
  arrows(X3[i],meanabgrgermcontrol[i]-seabgrgermcontrol[i],X3[i],meanabgrgermcontrol[i]+seabgrgermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meanabgrgerm[i]-seabgrgerm[i],X2[i],meanabgrgerm[i]+seabgrgerm[i],length=.05,angle=90,code=0)}
points(meanabgrgermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meanabgrgerm~X2,bg="black",type="p",pch=21,cex=1.8)
mtext("Abies grandis",side=3,line=.8, adj=.1, font=3)
mtext("a)",side=3,line=.8, adj=-.05)
axis(2,at=c(0,2,4,6,8,10),las=1, cex.axis=1.3)
legend(.4,11, bty="n",legend=c("Seeds added", "No seeds added"), pch=21,pt.bg=c("black","white"),pt.cex=1.3)

#THPL
plot(meanthplgermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.006,.6),type="p",bty="l",pch=21,cex=1.8,las=1,bg="white", cex.main=1.3)
for (i in 1:8){
  arrows(X3[i],meanthplgermcontrol[i]-sethplgermcontrol[i],X3[i],meanthplgermcontrol[i]+sethplgermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meanthplgerm[i]-sethplgerm[i],X2[i],meanthplgerm[i]+sethplgerm[i],length=.05,angle=90,code=0)}
points(meanthplgermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meanthplgerm~X2,bg="black",type="p",pch=21,cex=1.8)
axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),las=1, cex.axis=1.3)
mtext("Thuja plicata",side=3,line=.8, adj=.1, font=3)
mtext("b)",side=3,line=.8, adj=-.05)
mtext("Mean number of total germinants", side=2, line=3,cex=1.2)
#TSHE
plot(meantshegermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.02,2.5),type="p",bty="l",pch=21,cex=1.8,las=1,bg="white", cex.main=1.3)
for (i in 1:8){
  arrows(X3[i],meantshegermcontrol[i]-setshegermcontrol[i],X3[i],meantshegermcontrol[i]+setshegermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meantshegerm[i]-setshegerm[i],X2[i],meantshegerm[i]+setshegerm[i],length=.05,angle=90,code=0)}
points(meantshegermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meantshegerm~X2,bg="black",type="p",pch=21,cex=1.8)
axis(2,at=c(0,.5,1,1.5,2,2.5),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=.8, adj=.1, font=3)
mtext("c)",side=3,line=.8, adj=-.05)
axis(1,at=c(1,2,3,4),labels=c("Hedera present","Hedera removed","Hedera present","Hedera removed"),tick = FALSE)
axis(1,at=c(1,2,3,4),labels=c("wood absent","wood absent","wood added","wood added"),tick = FALSE, line=1)
###Figure 4.  surviving germinants
quartz(height=7,width=6)
par(mfcol=c(3,1),mai=c(.6,.6,.1,.1), omi=c(.7,.1,.2,.2))
X=c(1,2,3,4)
#pl
#plot control (no seeds added) in white
plot(meanabgrsurvgermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.1,4),type="p",bty="l",pch=21,cex=1.8,las=1,bg=c("white"), cex.main=1.3)
##add seeds added plots in gray
#error bars
for (i in 1:8){
  arrows(X3[i],meanabgrsurvgermcontrol[i]-seabgrsurvgermcontrol[i],X3[i],meanabgrsurvgermcontrol[i]+seabgrsurvgermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meanabgrsurvgerm[i]-seabgrsurvgerm[i],X2[i],meanabgrsurvgerm[i]+seabgrsurvgerm[i],length=.05,angle=90,code=0)}
points(meanabgrsurvgermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meanabgrsurvgerm~X2,bg="black",type="p",pch=21,cex=1.8)
mtext("Abies grandis",side=3,line=.8, adj=.1, font=3)
mtext("a)",side=3,line=.8, adj=-.05)
axis(2,at=c(0,1,2,3,4),las=1, cex.axis=1.3)
legend(.4,4.3, bty="n",legend=c("Seeds added", "No seeds added"), pch=21,pt.bg=c("black","white"),pt.cex=1.3)

#THPL
plot(meanthplsurvgermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.006,.3),type="p",bty="l",pch=21,cex=1.8,las=1,bg="white", cex.main=1.3)
for (i in 1:8){
  arrows(X3[i],meanthplsurvgermcontrol[i]-sethplsurvgermcontrol[i],X3[i],meanthplsurvgermcontrol[i]+sethplsurvgermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meanthplsurvgerm[i]-sethplsurvgerm[i],X2[i],meanthplsurvgerm[i]+sethplsurvgerm[i],length=.05,angle=90,code=0)}
points(meanthplsurvgermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meanthplsurvgerm~X2,bg="black",type="p",pch=21,cex=1.8)
axis(2,at=c(0,.1,.2,.3),las=1, cex.axis=1.3)
mtext("Thuja plicata",side=3,line=.8, adj=.1, font=3)
mtext("b)",side=3,line=.8, adj=-.05)
mtext("Mean number of suviving germinants", side=2, line=3,cex=1.2)
#TSHE
plot(meantshesurvgermcontrol~X3,ylab="",xlab="",xlim=c(.5,4.5),col.axis="white",ylim=c(-.02,1),type="p",bty="l",pch=21,cex=1.8,las=1,bg="white", cex.main=1.3)
for (i in 1:8){
  arrows(X3[i],meantshesurvgermcontrol[i]-setshesurvgermcontrol[i],X3[i],meantshesurvgermcontrol[i]+setshesurvgermcontrol[i],length=.05,angle=90,code=0)}
for (i in 1:8){
  arrows(X2[i],meantshesurvgerm[i]-setshesurvgerm[i],X2[i],meantshesurvgerm[i]+setshesurvgerm[i],length=.05,angle=90,code=0)}
points(meantshesurvgermcontrol~X3,bg="white",type="p",pch=21,cex=1.8)
points(meantshesurvgerm~X2,bg="black",type="p",pch=21,cex=1.8)
axis(2,at=c(0,.2,.4,.6,.8,1),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=.8, adj=.1, font=3)
mtext("c)",side=3,line=.8, adj=-.05)
axis(1,at=c(1,2,3,4),labels=c("Hedera present","Hedera removed","Hedera present","Hedera removed"),tick = FALSE)
axis(1,at=c(1,2,3,4),labels=c("wood absent","wood absent","wood added","wood added"),tick = FALSE, line=1)

######Transplant survival
##TSHE
tshesurvgr=read.csv("data/TSHEtransplant.csv", header=T)
head(tshesurvgr)
dim(tshesurvgr)
unique(tshesurvgr$Status_2012Aug3)
tshesurvgr$Wood=as.factor(tshesurvgr$Wood)
tshesurvgr$Ivy=as.factor(tshesurvgr$Ivy)
tshesurvgr$Block=as.factor(tshesurvgr$Block)
tshesurvgr$Plot=as.factor(tshesurvgr$Plot)

tapply(tshesurvgr$Status_2012Aug3,list(tshesurvgr$Wood,tshesurvgr$Ivy,tshesurvgr$Plot),sum,na.rm=T)
tapply(tshesurvgr$Status_2012Aug3,list(tshesurvgr$Wood,tshesurvgr$Ivy,tshesurvgr$Plot),length)
sum(tapply(tshesurvgr$Status_2012Aug3,list(tshesurvgr$Wood,tshesurvgr$Ivy,tshesurvgr$Plot),length), na.rm=T)
###survival glm at end of study
survmod.tshe=glmer(Status_2012Aug3~Wood*Ivy + (1|Block), family=binomial,data=tshesurvgr)#get error message!
summary(survmod.tshe)#positive effect of wood, negative effect of ivy, but weaker
Anova(survmod.tshe)#both ivy and wood sig, but no interaction. 
#survmod.tshe=glm(Status_2012Aug3~Wood*Ivy, family=binomial,data=tshesurvgr) 
##growth at end of study
unique(tshesurvgr$Ht_cm_2011May10)
tshesurvgr$HI=tshesurvgr$Ht_cm_2012Aug3-tshesurvgr$Ht_cm_2011May10
himod.tshe=lmer(HI~as.factor(Wood)*as.factor(Ivy) + (1|Block),data=tshesurvgr)#get error message!
summary(himod.tshe)#positive effect of wood, negative effect of ivy
Anova(himod.tshe)#both wood sig, ivy & interaction not signiciant 
##THPL
thplsurvgr=read.csv("data/THPLtransplant.csv", header=T)
head(thplsurvgr)
dim(thplsurvgr)
unique(thplsurvgr$Ivy)
thplsurvgr$HI=thplsurvgr$Ht_cm_2012Aug3-thplsurvgr$Ht_cm_2011May26
thplsurvgr$Wood=as.factor(thplsurvgr$Wood)
thplsurvgr$Ivy=as.factor(thplsurvgr$Ivy)
thplsurvgr$Block=as.factor(thplsurvgr$Block)
thplsurvgr$Plot=as.factor(thplsurvgr$Plot)

tapply(thplsurvgr$Status_2012Aug3,list(thplsurvgr$Wood,thplsurvgr$Ivy,thplsurvgr$Block),sum,na.rm=T)

tapply(thplsurvgr$Status_2012Aug3,list(thplsurvgr$Wood,thplsurvgr$Ivy,thplsurvgr$Block),length)
sum(tapply(thplsurvgr$Status_2012Aug3,list(thplsurvgr$Wood,thplsurvgr$Ivy,thplsurvgr$Plot),length), na.rm=T)

###survival glm at end of study
survmod.thpl=glmer(Status_2012Aug3~Wood*Ivy + (1|Block), family=binomial,data=thplsurvgr,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))#get error message! doesn't converge- i think because perfect survival with wood!
summary(survmod.thpl)#positive effect of wood, negative effect of ivy, negative interaction
Anova(survmod.thpl)# ivy and wood and interaction sig
#survmod.tshe=glm(Status_2012Aug3~Wood*Ivy, family=binomial,data=tshesurvgr) 
thplsurvgr$Status_2012Aug3
##growth at end of study
thplsurvgr$HI=thplsurvgr$Ht_cm_2012Aug3-thplsurvgr$Ht_cm_2011Mar3
himod.thpl=lmer(HI~Wood*Ivy + (1|Block),data=thplsurvgr)#
summary(himod.thpl)#negative effect of wood, negative effect of ivy
Anova(himod.thpl)#both wood sig, ivy & interaction not signiciant#nothing significant
###Figure of survival/growth of transplants
tshesurvgr$trtmt=paste(tshesurvgr$Wood,tshesurvgr$Ivy, sep = ".")
tshesurvgr$trtmt2=NA
tshesurvgr[which(tshesurvgr$trtmt=="0.0"),]$trtmt2="2Ivy rem"
tshesurvgr[which(tshesurvgr$trtmt=="0.1"),]$trtmt2="1Control"
tshesurvgr[which(tshesurvgr$trtmt=="1.0"),]$trtmt2="4Wood,no ivy"
tshesurvgr[which(tshesurvgr$trtmt=="1.1"),]$trtmt2="3Wood+ivy"
proptshetrans<-tapply(tshesurvgr$Status_2012Aug3,list(tshesurvgr$trtmt2,tshesurvgr$Block),sum,na.rm = TRUE)/tapply(tshesurvgr$Status_2012Aug3,list(tshesurvgr$trtmt2,tshesurvgr$Block),length)
proptshetrans[3,4]=0#replace NA with 0
meantshetrans=rowMeans(proptshetrans,)
sdtshetrans=rbind(sd(proptshetrans[1,]),sd(proptshetrans[2,]),sd(proptshetrans[3,]),sd(proptshetrans[4,]))
setshetrans=sdtshetrans/sqrt(6)
thplsurvgr$trtmt=paste(thplsurvgr$Wood,thplsurvgr$Ivy, sep = ".")
thplsurvgr$trtmt2=NA
thplsurvgr[which(thplsurvgr$trtmt=="0.0"),]$trtmt2="2Ivy rem"
thplsurvgr[which(thplsurvgr$trtmt=="0.1"),]$trtmt2="1Control"
thplsurvgr[which(thplsurvgr$trtmt=="1.0"),]$trtmt2="4Wood, no ivy"
thplsurvgr[which(thplsurvgr$trtmt=="1.1"),]$trtmt2="3Wood+ivy"
propthpltrans<-tapply(thplsurvgr$Status_2012Aug3,list(thplsurvgr$trtmt2,thplsurvgr$Block),sum,na.rm = TRUE)/tapply(thplsurvgr$Status_2012Aug3,list(thplsurvgr$trtmt2,thplsurvgr$Block),length)
sdthpltrans=rbind(sd(propthpltrans[1,]),sd(propthpltrans[2,]),sd(propthpltrans[3,]),sd(propthpltrans[4,]))
sethpltrans=sdthpltrans/sqrt(6)
meanthpltrans=rowMeans(propthpltrans)

##Plot- interaction plot:
###Use triangles for no wood chips, circles for wood chips added, gray for ivy present, black for ivy absent
quartz(height=6,width=8)
par(mfcol=c(2,2),mai=c(.6,.6,.1,.1), omi=c(.7,.1,.2,.2))
#X=c(1,2,1,2)
#plot wood absent (triangles)
#THPL
plot(meanthpltrans~X,ylab="",xlab="",xlim=c(.5,2.5),col.axis="white",ylim=c(0,1),type="p",bty="l",pch=c(24,24,21,21),cex=1.8,las=1,bg="gray", cex.main=1.3)
axis(2,at=c(0,.2,.4,.6,.8,1),las=1, cex.axis=1.3)
mtext("Thuja plicata",side=3,line=.8, adj=0, font=3)
mtext("a)",side=3,line=.8, adj=-.05)

#error bars
for (i in 1:length(meanthpltrans)){
  arrows(X[i],meanthpltrans[i]-sethpltrans[i],X[i],meanthpltrans[i]+sethpltrans[i],length=.05,angle=90,code=0)}
points(meanthpltrans~X,pch=c(24,24,21,21),cex=1.8,las=1,bg="gray")
mtext("Survival",side=2,line=3, adj=.5, cex=1.2)
plot(meantshetrans~X,ylab="",xlab="",xlim=c(.5,2.5),col.axis="white",ylim=c(0,1),type="p",bty="l",pch=c(24,24,21,21),cex=1.8,las=1,bg=c("gray"), cex.main=1.3)
axis(2,at=c(0,.2,.4,.6,.8,1),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=.8, adj=0, font=3)
mtext("b)",side=3,line=.8, adj=-.05)
#error bars
for (i in 1:length(meantshetrans)){
  arrows(X[i],meantshetrans[i]-setshetrans[i],X[i],meantshetrans[i]+setshetrans[i],length=.05,angle=90,code=0)}
points(meantshetrans~X,pch=c(24,24,21,21),cex=1.8,las=1,bg="gray")
mtext("Survival",side=2,line=3, adj=.5, cex=1.2)
#mtext("Treatment",side=1,line=3.3, adj=.5)
axis(1,at=c(1,2,3,4),labels=c("Ivy","No ivy","Ivy,","No ivy,"),tick = FALSE)

##Growth
tshesurvgr2=na.omit(tshesurvgr)
growtshetrans2<-tapply(tshesurvgr2$HI,tshesurvgr2$trtmt2,mean)
setshetrans<-tapply(tshesurvgr2$HI,tshesurvgr2$trtmt2,sd)/sqrt(tapply(tshesurvgr2$HI,tshesurvgr2$trtmt2,length))
head(thplsurvgr)
thplsurvgr2=na.omit(thplsurvgr)
growthpltrans2<-tapply(thplsurvgr2$HI,thplsurvgr2$trtmt2,mean, na.rm=T)
sethpltrans<-tapply(thplsurvgr2$HI,thplsurvgr2$trtmt2,sd)/sqrt(tapply(thplsurvgr2$HI,thplsurvgr2$trtmt2,length))
#THPL
plot(growthpltrans2~X,ylab="",xlab="",xlim=c(.5,2.5),col.axis="white",ylim=c(0,5),type="p",bty="l",pch=c(24,24,21,21),cex=1.8,las=1,bg="gray", cex.main=1.3)
axis(2,at=c(0,1,2,3,4,5),las=1, cex.axis=1.3)
mtext("Thuja plicata",side=3,line=.8, adj=0, font=3)
mtext("c)",side=3,line=.8, adj=-.05)
#error bars
for (i in 1:length(growthpltrans2)){
  arrows(X[i],growthpltrans2[i]-sethpltrans[i],X[i],growthpltrans2[i]+sethpltrans[i],length=.05,angle=90,code=0)}
points(growthpltrans2~X,pch=c(24,24,21,21),cex=1.8,las=1,bg="gray")
legend(2,5.2, bty="n",legend=c("No Wood","Wood"), pch=c(24,21),pt.bg="gray", pt.cex=1.3)

#add to plot
plot(growtshetrans2~X,ylab="",xlab="",xlim=c(.5,2.5),col.axis="white",ylim=c(-.5,1),type="p",bty="l",pch=c(24,24,21,21),cex=1.8,las=1,bg="gray", cex.main=1.3)
axis(2,at=c(-.5,0,.5,1),las=1, cex.axis=1.3)
mtext("Tsuga heterophylla",side=3,line=.8, adj=0, font=3)
mtext("d)",side=3,line=.8, adj=-.05)
#error bars
for (i in 1:length(growtshetrans2)){
  arrows(X[i],growtshetrans2[i]-setshetrans[i],X[i],growtshetrans2[i]+setshetrans[i],length=.05,angle=90,code=0)}
points(growtshetrans2~X,pch=c(24,24,21,21),cex=1.8,las=1,bg="gray")
mtext("Height increment (cm)",side=2,line=3, adj=.5, cex=1.2)
axis(1,at=c(1,2,3,4),labels=c("Ivy","No ivy","Ivy,","No ivy,"),tick = FALSE)
#mtext("Treatment",side=1,line=3.3, adj=.5)

#get height differences for TSHE and THPL transplants
mean(tshesurvgr$Ht_cm_2011Feb14, na.rm=T)
mean(thplsurvgr$Ht_cm_2011Mar3, na.rm=T)


##Soil moisture and PAR data
setwd("~/Dropbox/SeattleConiferRegen/Manuscript/Analyses_ForManuscript/Ailene")
#Start with Soil Moisture (SM) Data
smdata<-read.csv("SMData.csv",header=T)

#take out rows that have no value for soil moisture:
smdata2<-smdata[is.na(smdata$SMAvg)==FALSE,]
#Looks at data:soil moisture over the season, by block
hist(smdata2$SMAvg)
smdata2$Block=as.factor(smdata2$Block)
smdata2$JulianDate2=as.factor(smdata2$JulianDate)
sm.mod<-lmer(SMAvg~IvyPresent*WCPresent+JulianDate(1|Block)+(1|JulianDate2), data=smdata2)
summary(sm.mod)
Anova(sm.mod, test="Chi")

sm.mod2<-lmer(SMAvg~IvyPresent*WCPresent+JulianDate+(1|Block)+(1|JulianDate2), data=smdata2)
summary(sm.mod2)
Anova(sm.mod2, test="Chi")
AIC(sm.mod,sm.mod2)
#AIC is lower for model without date included
#Plot soil moisture over time in each treatment
unique(smdata2$JulianDate)

#soil moisture over the season, by treatment
sm.mean=tapply(smdata2$SMAvg,list(smdata2$Trtmnt,smdata2$JulianDate),mean, na.rm=T)
sm.n=tapply(smdata2$SMAvg,list(smdata2$Trtmnt,smdata2$JulianDate),length)
sm.sd=tapply(smdata2$SMAvg,list(smdata2$Trtmnt,smdata2$JulianDate),sd, na.rm=T)
sm.se=sm.sd/sqrt(sm.n)
quartz(width=7,height=6)
par(mfrow=c(2,1),mai=c(.6,.6,.2,.6), omi=c(.7,.1,.2,.2))
#gray=no wood
#black=wood
#dashed=no ivy
#solid=ivy 
plot(colnames(sm.mean),sm.mean[1,],type="l", xlab="Date", ylab="Soil Moisture (%)", bty="l", ylim=c(0,25), lwd=2, col="gray", las=1)
polygon(c(as.numeric(colnames(sm.mean)),rev(as.numeric(colnames(sm.mean)))),c((sm.mean[1,]+sm.se[1,]),rev(sm.mean[1,]-sm.se[1,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(sm.mean)),rev(as.numeric(colnames(sm.mean)))),c((sm.mean[2,]+sm.se[2,]),rev(sm.mean[2,]-sm.se[2,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(sm.mean)),rev(as.numeric(colnames(sm.mean)))),c((sm.mean[3,]+sm.se[3,]),rev(sm.mean[3,]-sm.se[3,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(sm.mean)),rev(as.numeric(colnames(sm.mean)))),c((sm.mean[4,]+sm.se[4,]),rev(sm.mean[4,]-sm.se[4,])), col="lightgray", border="lightgray")
lines(colnames(sm.mean),sm.mean[2,],lwd=2, lty=2,col="darkgray")
lines(colnames(sm.mean),sm.mean[1,],lwd=2, lty=1,col="darkgray")
lines(colnames(sm.mean),sm.mean[3,],lwd=2, lty=1,col="black")
lines(colnames(sm.mean),sm.mean[4,],lwd=2, lty=2,col="black")
legend(11250,28, bty="n",legend=c("Ivy,Wood", "No Ivy, Wood","Ivy,No Wood","No Ivy, No Wood"), lty=c(1,2,1,2),col=c("black","black","darkgray","darkgray"), lwd=2, cex=.8)
###Now PAR
pardata<-read.csv("PARData.csv",header=T)
head(pardata)
#take out rows that have no value for soil moisture:
pardata2<-pardata[is.na(pardata$PARAvg)==FALSE,]
dim(pardata2)#there aren't any rows with missing data, so can use pardata
dim(pardata)
pardata2$Block=as.factor(pardata2$Block)
pardata2$JulianDate2=as.factor(pardata2$JulianDate)
par.mod<-lmer(PARAvg~IvyPresent*WCPresent+(1|Block)+(1|JulianDate2), data=pardata2)
summary(par.mod)
Anova(par.mod, test="Chi")

par.mod2<-lmer(PARAvg~IvyPresent*WCPresent+JulianDate+(1|Block)+(1|JulianDate2), data=pardata2)
summary(par.mod2)
Anova(par.mod2, test="Chi")
AIC(par.mod,par.mod2)

#AIC is lower for model without date included
#Plot PAR over time in each treatment
#soil moisture over the season, by treatment
par.mean=tapply(pardata2$PARAvg,list(pardata2$Trtmnt,pardata2$JulianDate),mean, na.rm=T)
par.n=tapply(pardata2$PARAvg,list(pardata2$Trtmnt,pardata2$JulianDate),length)
par.sd=tapply(pardata2$PARAvg,list(pardata2$Trtmnt,pardata2$JulianDate),sd, na.rm=T)
par.se=par.sd/sqrt(par.n)
##add par to figure
plot(colnames(par.mean),par.mean[1,],type="l", xlab="Date", ylab="Light (mmol m-2s-1)", bty="l", ylim=c(0,50), lwd=2, col="gray",las=1)
polygon(c(as.numeric(colnames(par.mean)),rev(as.numeric(colnames(par.mean)))),c((par.mean[1,]+par.se[1,]),rev(par.mean[1,]-par.se[1,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(par.mean)),rev(as.numeric(colnames(par.mean)))),c((par.mean[2,]+par.se[2,]),rev(par.mean[2,]-par.se[2,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(par.mean)),rev(as.numeric(colnames(par.mean)))),c((par.mean[3,]+par.se[3,]),rev(par.mean[3,]-par.se[3,])), col="lightgray", border="lightgray")
polygon(c(as.numeric(colnames(par.mean)),rev(as.numeric(colnames(par.mean)))),c((par.mean[4,]+par.se[4,]),rev(par.mean[4,]-par.se[4,])), col="lightgray", border="lightgray")
lines(colnames(par.mean),par.mean[2,],lwd=2, lty=2,col="darkgray")
lines(colnames(par.mean),par.mean[1,],lwd=2, lty=1,col="darkgray")
lines(colnames(par.mean),par.mean[3,],lwd=2, lty=1,col="black")
lines(colnames(par.mean),par.mean[4,],lwd=2, lty=2,col="black")

###Compare ivy, CWD cover in plots prior to experiment
setwd("~/Dropbox/SeattleConiferRegen/Analysis")
pretreat.dat=read.csv("PreTreatmentIvyCWDData.csv", header=T)
head(pretreat.dat)
summary(lm(Ivy_PercCov~as.factor(Block), data=pretreat.dat))#block 5 significant; lower than block 1
summary(lm(CWD_PercCov~as.factor(Block), data=pretreat.dat))#no differences
