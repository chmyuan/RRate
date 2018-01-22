# TODO: Simulation study to show the performance of reproducibility rate and false irreproducibility rate
# 
# Author: wjiangaa
###############################################################################
##rm(list = ls())##

startTime<-Sys.time()
setwd("D://Coding//What is the probability of replicating a statistically significant association in genome-wide association studies")

#install.packages('ROCR',dependencies=T)
library('ROCR')
#install.packages("xtable", dependencies=T)
library('xtable')
#install.packages("Hmisc", dependencies=T)
library("Hmisc")

##source('RRate.R')##
##source('genSimSmry.R')##

#Simulation Settings
n<-69033 #Sample Size
m<-2473441 #SNP number
kappa<-0.01 #prevalence
CCR<-1 #Control-to-case ratio
pi1<-0.05 #Causal SNP proportion
sigma02<-0.04 #log(OR)~N(0,sigma02) for nonnull
fL<-0.05 #Low MAF for noncausal SNP
fH<-0.5 #High MAF for noncausal SNP
outputBase<-'output'
simFilename1<-'T2D1.txt'
paramFilename1<-'param1.txt' #Ground truth parameters
simFilename2<-'T2D2.txt'
paramFilename2<-'param2.txt'
#alpha<-5e-5 #The significance level of the primary study
alpha<-5e-8
repRatio<-0.5 #The ratio between the sample size of replication study and primary study
alphaR<-5e-5 #The significance level of the replication study
runNum<-5 #Number of runs
#fixedRunNum<-10 #Number of runs when parameters were simulated and fixed
repRatioSeq<-seq(0.5,1.0,0.1)
binNum<-5 #Number of bins to partition [0,1]

zalpha2<-qnorm(1-alpha/2)
zalphaR2<-qnorm(1-alphaR/2)
pi0<-1-pi1

#loading Data
loadData<-function(filename, header=T){
  tmpData<-try(read.table(filename,header=header),silent=TRUE)
  if(class(tmpData)!='try-error'){
    return(tmpData)
  }else return(F)
}

calAUC<-function(x,y){ #Area under the curve
  return(sum(diff(x) * (y[-1] + y[-length(y)]) / 2))
}

dir.create(outputBase,showWarnings=F)
outputBase<-paste(outputBase,'RealData',sep='/')
dir.create(outputBase,showWarnings=F)

RMSE<-matrix(0,nrow=runNum,ncol=1)  #RR TNR 
avgRRresult<-matrix(0,nrow=runNum, ncol=5) #True avg RR, avg RR, CIlow, CIhigh, RP

corrResult<-matrix(0,nrow=runNum,ncol=4)

run<-0
while(TRUE){
  run<-run+1
  cat('run',run,'\n')
  #Generate simulation dataset
  outputDir<-paste(outputBase,'/run',run,sep='')
  dir.create(outputDir,showWarnings=F)
  simFile1<-paste(outputDir,simFilename1,sep='/')
  paramFile1<-paste(outputDir,paramFilename1,sep='/')
  

  #If CCR of replication study is different from CCR in primary study, then use function SEest for SE2 calculation
  
  #Obtain summary statistics of the simulation
  T2D1<-read.table('T2D1.txt',header=TRUE)
  pval<-T2D1$P_VALUE
  ORhat<-T2D1$OR
  SEhat<-(log(T2D1$OR_95U)-log(T2D1$OR_95L))/2*1.96
  SE2hat<-SEhat/sqrt(repRatio)
  MUhat<-log(ORhat)
  z<-MUhat/SEhat
  
  pi0hat<-snpNullS(pval, showEst=F, info=F, outputDir=outputDir)$pi0
  
  #install.packages("limma", dependencies=T)
  #	library("limma")
  #	pi0hat<-propTrueNull(pval,method='hist')
  #	if(pi0hat==1){
  #		run<-run-1
  #		next
  #	}
  
  sigIdx<-(pval<alpha)
  
  simFile2<-paste(outputDir,simFilename2,sep='/')
  paramFile2<-paste(outputDir,paramFilename2,sep='/')
  T2D2<-read.table('T2D2.txt',header=TRUE)
  SEhat2<-(log(T2D2$OR_95U)-log(T2D2$OR_95L))/2*1.96
  zR<-log(T2D2$OR)/SEhat2
  repIdx<-(sign(z[sigIdx])*zR[sigIdx]>zalphaR2)
  avgRRresult[run,5]<-sum(repIdx)/sum(sigIdx)
  
  ###### EB #####
  EBresult<-repRateEst(MUhat,SEhat, SE2hat,zalpha2,zalphaR2, output=T,dir=outputDir)
  sigPos<-match((1:m)[sigIdx],EBresult$idx)	
  paramEst<-cbind(EBresult$pi0, EBresult$sigma02)
  colnames(paramEst) <- c('pi0','sigma02')
  write.table(paramEst, paste(outputDir,'paramEst.txt',sep='/'), col.names=T,row.names=F,sep='\t',quote=F)
  
  dir.create(paste(outputDir,'GroundTruth',sep='/'),showWarnings=F)
  groundTruth<-repRate(MUhat[sigIdx],SEhat[sigIdx], SE2hat[sigIdx],sigma02,pi0,zalphaR2,output=T, dir=paste(outputDir,'GroundTruth',sep='/'))
  
  #scatter plot of RR vs RR_EB
  pdf(paste(outputDir,'/RR_EB.pdf',sep=''))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(groundTruth$RR,EBresult$RR[sigPos],cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
       main=expression(paste(hat(RR),' vs RR')),xlab='RR',ylab=expression(hat(RR)))
  abline(0,1)
  dev.off()
  
  RMSE[run,1]<-sqrt(sum((groundTruth$RR-EBresult$RR[sigPos])^2)/sum(sigIdx))
  
  avgRRresult[run,1]<-groundTruth$GRR
  avgRRresult[run,2]<-EBresult$GRR
  avgRRresult[run,3]<-EBresult$GRRlow
  avgRRresult[run,4]<-EBresult$GRRhigh
  
  #Average Power vs Sample size
  #	avgRR<-rep(0,length(repRatioSeq))
  avgRRhat<-rep(0,length(repRatioSeq))
  #	avgRRlow<-rep(0,length(repRatioSeq))
  #	avgRRhigh<-rep(0,length(repRatioSeq))
  RP<-rep(0,length(repRatioSeq))
  for(i in 1:length(repRatioSeq)){
    tmpSE2hat<-SEhat/sqrt(repRatioSeq[i])
    tmpResult<-repRate(MUhat[sigIdx],SEhat[sigIdx], tmpSE2hat[sigIdx],EBresult$sigma02,EBresult$pi0,zalphaR2, output=F)
    avgRRhat[i]<-tmpResult$GRR
    
    tmpSimFile2<-paste(outputDir,'/tmpSimData2_',i,'.txt',sep='')
    tmpParamFile2<-paste(outputDir,'/tmpParam2_',i,'.txt',sep='')
    #invisible(genSimSmry(tmpSimFile2,n=round(repRatioSeq[i]*n),m=m,kappa=kappa,CCR=CCR,pi1=pi1,paramFile=tmpParamFile2,f=f,OR=OR))#
    #tmpSmryStats2<-loadData(tmpSimFile2)#
	tmpSmryStats2<-sample(T2D1,size=round(repRatioSeq[i]*n),replace=F)
	tmpSE<-(log(tmpSmryStats2$OR_95U)-log(tmpSmryStats2$OR_95L))/2*1.96
    tmpZR<-log(tmpSmryStats2$OR)/tmpSE
    tmpRepIdx<-(sign(z[sigIdx])*tmpZR[sigIdx]>zalphaR2)
    RP[i]<-sum(tmpRepIdx)/sum(sigIdx)
  }
  
  ymin<-0.9*min(avgRRhat,RP)
  ymax<-min(1,1.1*max(avgRRhat,RP))
  pdf(paste(outputDir,'/GRR_SampleSize.pdf',sep=''))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(repRatioSeq,avgRRhat,ylim=c(ymin,ymax),type='b',pch=2,col='red',lty=1,lwd=2,
       main=expression(paste('GRR (GRP) vs ',n^(2)/n^(1))),
       ylab='GRR (GRP)',xlab=expression(n^(2)/n^(1)),cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  lines(repRatioSeq,RP,type='b',col='blue',lty=2,lwd=2)
  legend('bottomright',c('GRR','GRP'),
         col=c('red','blue'),pch=c(2,1),lty=c(2,1),lwd=2,cex=1.5)
  grid(lwd=2)
  dev.off()
  
  #Replication rate and Replication proportion
  sortedRR<-sort(EBresult$RR[sigPos])
  breaks<-sortedRR[ceiling(seq(1,length(sigPos),length.out=binNum+1))]
  #	binLabel<-cut(EBresult$RR[sigPos],breaks=breaks,labels=F,include.lowest=T)
  orderRR<-order(EBresult$RR[sigPos])
  binLabel<-rep(0,sum(sigIdx))
  binLabel[orderRR]<-ceiling((1:sum(sigIdx))/ceiling(sum(sigIdx)/binNum))
  binCount<-rep(0,binNum)
  binRep<-rep(0,binNum)
  sumRR<-rep(0,binNum)
  for(i in 1:sum(sigIdx)){
    binCount[binLabel[i]]<-binCount[binLabel[i]]+1
    sumRR[binLabel[i]]<-sumRR[binLabel[i]]+EBresult$RR[sigPos[i]]
    if(repIdx[i]==T) binRep[binLabel[i]]<-binRep[binLabel[i]]+1
  }
  nonZeroBin<-which(binCount>0)
  if(length(nonZeroBin)>2){
    #	x<-(breaks[nonZeroBin]+breaks[nonZeroBin+1])/2
    x<-sumRR[nonZeroBin]/binCount[nonZeroBin]
    y<-binRep[nonZeroBin]/binCount[nonZeroBin]
    rho<-cor(x,y)
    pdf(paste(outputDir,'/RP_RR.pdf',sep=''))
    lim<-c(min(x,y),max(x,y))
    par(mar=c(5.1,5.1,4.1,2.1))
    plot(x,y,cex=2,col='red',xlim=lim,ylim=lim,xlab='Replicatioin Rate (RR)',ylab='Replication Proportion (RP)'
         ,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,main=as.expression(bquote(paste('RP vs RR, ',rho,'=',.(format(rho,digits=3))))))
    abline(0,1)	
    dev.off()
  }
  
  #source('HLtest.R')##
  Hobj<-HLtest(repIdx,EBresult$RR[sigPos],g=binNum,dir=outputDir)
  HLtestContent<-cbind(Hobj$H,Hobj$pval_chi2,Hobj$pval_boot)
  colnames(HLtestContent) <- c("H", "p-value (chi2)", "p-value (bootstrap)")
  write.table(HLtestContent, paste(outputDir,'HLresult.txt',sep='/'), col.names=T,row.names=F,sep='\t',quote=F)
  
  HLtestTable<-xtable(HLtestContent,digits=3,caption='Hosmer-Lemeshow Test')
  align(HLtestTable) <- rep("c",4)
  print(HLtestTable,hline.after=0,latex.environments='center',
        type="latex", file=paste(outputDir,'HLresult.tex',sep='/'))
  print(HLtestTable,hline.after=0,latex.environments='center',
        type="html", file=paste(outputDir,'HLresult.html',sep='/'))
  
  #######Correlation Test##########
  corrResult[run,1]<-cor(repIdx,EBresult$RR[sigPos])
  corrResult[run,2]<-length(repIdx)-2
  corrResult[run,3]<-corrResult[run,1]*sqrt(corrResult[run,2])/sqrt(1-corrResult[run,1]^2)
  corrResult[run,4]<-1-pt(corrResult[run,3],df=corrResult[run,2])
  
  
  ########ROC/PR curve###########
  if(sum(repIdx)!=0 & sum(repIdx)!=sum(sigIdx)){
    #Precision-recall curve for RR
    pred<-prediction(EBresult$RR[sigPos],repIdx)
    perf<-performance(pred,"tpr","fpr")
    #		perf<-performance(pred,"prec","rec")
    #		perf@y.values[[1]][1]<-1
    AUC<-calAUC(perf@x.values[[1]],perf@y.values[[1]])
    
    pred2<-prediction(-pval[sigIdx],repIdx)
    perf2<-performance(pred2,"tpr","fpr")
    #		perf2<-performance(pred2,"prec","rec")
    #		perf2@y.values[[1]][1]<-1
    AUC2<-calAUC(perf2@x.values[[1]],perf2@y.values[[1]])	
    
    pdf(paste(outputDir,'/RR_ROC.pdf',sep=''))
    par(mar=c(5.1,5.1,4.1,2.1))
    ylim<-c(min(perf@y.values[[1]],perf2@y.values[[1]]),max(perf@y.values[[1]],perf2@y.values[[1]]))
    #		plot(perf,ylim=ylim,col='red',lty=1,lwd=2,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,
    #				main=as.expression(bquote(paste(hat(RR),', AUPRC=',.(format(AUC,digits=3))))))
    plot(perf,ylim=ylim,col='red',lty=1,lwd=2,main='ROC',cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    plot(perf2,col='blue',lty=2,lwd=2,add=T)
    abline(0,1)
    legend('bottomright',c(as.expression(bquote(paste('RR, AUC=',.(format(AUC,digits=3))))),
                           as.expression(bquote(paste(italic('p'),'-value, AUC=',.(format(AUC2,digits=3)))))), col=c('red','blue'),lty=c(1,2),cex=1.5,lwd=2)
    dev.off()
    
  }
  
  if(run==runNum) break
}

#RMSE output
avgRMSE<-colMeans(RMSE)

combineRMSE<-rbind(RMSE,avgRMSE)
rownames(combineRMSE) <- c(paste('run',1:runNum),'Avg.')
colnames(combineRMSE) <- c("RR")
RMSEtable<-xtable(combineRMSE,digits=3,caption='RMSE of RR')
align(RMSEtable) <- rep("c",2)
print(RMSEtable,hline.after=c(0,runNum),latex.environments='center',
      type="latex", file=paste(outputBase,'RMSE.tex',sep='/'))
print(RMSEtable,hline.after=c(0,runNum),latex.environments='center',
      type="html", file=paste(outputBase,'RMSE.html',sep='/'))

avgRRframe<-data.frame(avgRRresult[,1],avgRRresult[,2],
                       paste('(',format(avgRRresult[,3],digits=3),' , ',format(avgRRresult[,4],digits=3),')',sep=''), avgRRresult[,5])
rownames(avgRRframe) <- paste('run',1:runNum)
colnames(avgRRframe) <- c("Avg. RR", "Avg. $\\hat{RR}$", "95% CI","RP")
avgRRtable<-xtable(avgRRframe,digits=3,caption='Average RR estimation')
align(avgRRtable) <- rep("c",5)
print(avgRRtable,hline.after=0,latex.environments='center',
      type="latex", file=paste(outputBase,'avgRR.tex',sep='/'))
print(avgRRtable,hline.after=0,latex.environments='center',
      type="html", file=paste(outputBase,'avgRR.html',sep='/'))

corrFrame<-data.frame(corrResult)
rownames(corrFrame) <- paste('run',1:runNum)
colnames(corrFrame) <- c("r", "df", "t","p-value")
corrTable<-xtable(corrFrame,digits=3,caption='Correlation between replication status and RR')
align(corrTable) <- rep("c",5)
print(corrTable,hline.after=0,latex.environments='center',
      type="latex", file=paste(outputBase,'corr.tex',sep='/'))
print(corrTable,hline.after=0,latex.environments='center',
      type="html", file=paste(outputBase,'corr.html',sep='/'))

ymin<-min(avgRRresult[,1],avgRRresult[,3])
ymax<-max(avgRRresult[,1],avgRRresult[,4])
pdf(paste(outputBase,'/avgRR.pdf',sep=''))
plot(1:run,avgRRresult[,1],ylim=c(ymin-0.05,min(ymax+0.05,1)),type='b',pch=2,col='red',lty=2,lwd=1.5,
     main=expression(paste('Avg. RR ( Avg. ',hat(RR),' )')),ylab='Avg. RR',xlab='run #', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
errbar(1:run, avgRRresult[,2], avgRRresult[,4], avgRRresult[,3], add=T,type='b', 
       pch=1, lty=1, lwd=1.5, cap=.02,col='blue')
legend('bottomright',c('Avg. RR',expression(paste('Avg. ',hat(RR)))),
       col=c('red','blue'),pch=c(2,1),lty=c(2,1),lwd=1.5,cex=1.5)
grid(lwd=2)
dev.off()

cat('success!\n')
elapsedTime<-Sys.time()-startTime
print(elapsedTime)
save.image(paste(outputBase,'/.RData',sep=''))
