## One compartment drug to soluble target binding model ####
rm(list=ls())
name='PKB1C_MC'

## load libraries & set background parameters #####
starttime=Sys.time()
library(doParallel) # note that this loads required library(foreach)
library(plyr)
library(ggplot2)
library(grid)
library(MASS)
cpu=4				# number of cpu compute cores to run in parallel
milli.nano = 1e6 # unit conversions
micro.nano = 1e3

## Trial design, dose event table ####
N     = 300	        # number of patients
MwD   = 150000      # molecular weight of drug 
MwL   = 190000      # molecular weight of ligand 
dose  = 24          # mg
ndoses= 4
tau   = 28          # days
tstart= -28         # days
firstdose = 0       # days
tend  = 7*tau       # days
Doses=data.frame(var   = rep("A1",ndoses),
                time   = seq(firstdose,(ndoses-1)*tau,by=tau),
                value  = rep(dose/MwD*milli.nano,ndoses),   # note dose unit conversion mg to nmoles
                method = rep("add",ndoses))
# Customisable dosing table
# Doses=data.frame(var   = c("A1" ,"A1" ,"A1" ,"A1" ,"A1" ,"A1" ,"A1" ,"A1" ),
#                 time   = c(0    ,14   ,28   ,42   ,56   ,70   ,84   ,98   ),
#                 value  = c(375  ,375  ,375  ,375  ,375  ,375  ,375  ,375  )/MwD*milli.nano,
#                 method = c("add","add","add","add","add","add","add","add"))

## Parameters (typical, population mean) from Xolair DTE2013 model (lowered KD) ####
tCLX= 0.200	# clearance drug
tCLE= 1.98	# clearance target
tCLC= 0.442	# clearance complexes [L/d]
tVX = 8.08	# volume [L]
tVC = 7.15	# volume for complexes [L]
tEPR= 655	  # target production rate [ug/d]
tKD = 1.53/10	# equilibrium binding dissociation constant [nM]
tka = 0.446	# absorption rate constant [1/d]

## Random effects parameters, variances, from Xolair DTE2013 model ####
vCLX = 0.117
  cv12 = 0.0750	  # covariance
vVX  = 0.0711
vCLC = 0.0348
  cv34 = -0.0189	# covariance
vEPR = 0.0627
vCLE = 0.0362
vVC  = 0.0422		
vka  = 0.330	
vKd  = 0.0465

## Create vectors of random patients' parameters ####
mu1=c(0,0)				          # means of zero for each
mu2=c(0,0)
si1=matrix(c(	vCLX,cv12,		# omega covariance matrix
		          cv12,vVX ),2,2)
si2=matrix(c(	vCLC,cv34,
		          cv34,vEPR),2,2)
x1=mvrnorm(n=N,mu1,si1)     # multivariate random normal distribution sampling
x2=mvrnorm(n=N,mu2,si2)     # can also include ,tol=1e-6,empirical=F,EISPACK=F

mvCLX = x1[,1]
mvVX  = x1[,2]
mvCLC = x2[,1]
mvEPR = x2[,2]
svCLE = rnorm(N,0,sqrt(vCLE))
svVC  = rnorm(N,0,sqrt(vVC))
svka  = rnorm(N,0,sqrt(vka))
svKd  = rnorm(N,0,sqrt(vKd))

parms=data.frame(            #collect patients' parameters to a dataframe
		CLX=rep(tCLX,N)*exp(mvCLX),
		CLE=rep(tCLE,N)*exp(svCLE),
		CLC=rep(tCLC,N)*exp(mvCLC),
		VX =rep(tVX ,N)*exp(mvVX ),
		VC =rep(tVC ,N)*exp(svVC ),
		EPR=rep(tEPR,N)*exp(mvEPR),
		KD =rep(tKD ,N)*exp(svKd ),
		ka =rep(tka ,N)*exp(svka ))

## DES: coupled differential equations and additional output variables function ########################
DES=function(t,y,parms){with(as.list(c(y,parms)),{
XT  =A2
ET  =A3
VE  =VX
COMP=((KD*VX*VE/VC+XT+ET)-((KD*VX*VE/VC+XT+ET)**2-4*XT*ET)**0.5)/2
XF  = XT-COMP                       # free drug amount
EF  = ET-COMP                       # free target amount
dA1 =-ka*A1                         # subcutaneous injection site amount
dA2 = ka*A1  -CLX/VX*XF-CLC/VC*COMP	# Total drug amount
dA3 = EPR/MwL*micro.nano-CLE/VE*EF-CLC/VC*COMP	# total target amount with unit conversion on input
addvar=c(
CTX = MwD/micro.nano*(XF/VX+COMP/VC),      # outputs in ng/mL or ug/L
CFX = MwD/micro.nano*(XF/VX),
CTE = MwL/micro.nano*(EF/VE+COMP/VC),
CFE = MwL/micro.nano*(EF/VE)     )
list(c(dA1,dA2,dA3),addvar)})}# output derivatives & additional variables

## Integrate using deSolve within parallel foreach loops ####
integ.starttime=Sys.time()
d=NULL
cluster=makeCluster(cpu)		  # Register a parallel backend
registerDoParallel(cluster)		# to run in parallel with foreach
d=foreach(id=seq(N),.combine=rbind) %dopar% {	# run N individuals in parallel and rbind into d dataframe
parm=unlist(parms[id,])		# get parameters for individual id from the parms data frame
VX =parms[id,'VX' ]			  # get individuals' parameters for setting target initial
EPR=parms[id,'EPR']
CLE=parms[id,'CLE']
init=c(A1=0,				# Initial conditions
	     A2=0,
	     A3=VX*EPR/CLE/MwL*micro.nano)    # target is production divided by elimination with unit conversion
library(deSolve)			# must have the library command in the foreach loop
data.frame(id=id,lsoda(func=DES,y=init,times=tstart:tend,parms=parm,events=list(data=Doses)))
} ; stopCluster(cluster)	# end of foreach loop and close the cluster
integ.time=paste('integ.time =',as.character.POSIXt(signif(Sys.time()-integ.starttime,3)))
# head(d) ; tail(d)

## Plot results ####
breaks=signif(10^(seq(-6,6,by=0.5)),1)
f1=ggplot(d)+
	geom_line(aes(time,CTX,group=id),colour="darkgreen",size=0.7,alpha=10/N) +	
	geom_line(aes(time,CFE,group=id),colour="red",size=0.7,alpha=10/N) +	
	geom_line(aes(time,CTE,group=id),colour="blue",size=0.7,alpha=10/N) +	
	geom_line(aes(time,CFX,group=id),colour="black",size=0.7,alpha=10/N) +	
	scale_y_log10("Concentration (ng/mL)",breaks=breaks) +
	scale_x_continuous("\nTime (days)", breaks=seq(tstart,tend,by=tau)) +
 	geom_segment(data=Doses,aes(x=time,xend=time,y=0,yend=value/30),
		colour="black",size=1,arrow=arrow(length=unit(0.1,"inches"),type="open",ends="last")) +
	theme_bw(base_size=11) +
	theme(legend.position="none")+
  coord_cartesian(ylim=c(min(d$CFE),max(d$CTX)))
f1
datetime=format(Sys.time(),"%Y-%m-%d-%H%M")
filename=paste(name,'_',datetime,'.pdf',sep='')
run.time=paste('run.time =',as.character.POSIXt(signif(Sys.time()-starttime,3)))  # integration time
pdf(file=filename,width=10,height=6)
plot.new()
grid_layout=grid.layout(100,100)
pushViewport(viewport(layout=grid_layout))
print(f1,vp=viewport(layout.pos.row=01:97,layout.pos.col=01:100)) 
grid.text(paste(getwd(),'/',filename,' ',integ.time,' ',run.time,sep=''),x=0.01,y=0.01,hjust=0,vjust=0,
	gp=gpar(fontsize=8,col="grey20"))
dev.off()



