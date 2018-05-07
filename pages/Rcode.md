---
layout: page
---
<html lang="en-us">
<h2> Rcode

<p class="view"><a href="https://cran.r-project.org/web/packages/CpGassoc/index.html"> <h3>CpGassoc</h3> </a></p>
<p> R package wrote under the supervision of Dr. Conneely. Designed to test for association between methylation at CpG sites across the genome and a phenotype of interest, adjusting for any relevant covariates.</p>

<h3 id="EWASMEDIAT">Fast EWAS Mediation</h3>
<p>Code to run quick mediation analysis across a set of DNA methylation markers. Gives association between CpG and Y given A and Xm as well as association between A and CpG given Xm. Will give for each CpG site. Assuming linear models for all variables. Xm is matrix of covariates adjusting for (not including intercept), CpG is an Nx p matrix with rows corresponding to sample and columns corresponding to DNA methylation marker. </p> 


	fastMEDiat.func<-function(Y,Xm,A,CpG){
  		XX<-cbind(1,A,Xm)
  		numExplan<-ncol(XX)
  		XXproj<-solve(t(XX)%*%XX)%*%t(XX)
  		residY<-Y-XX%*%XXproj%*%Y
  		onDNA<-XXproj%*%CpG
  		residDNA<-CpG-XX%*%onDNA
  		par1<-colSums(residDNA^2)
  		SigEst<-par1/(length(Y)-numExplan)
	        VarEst<-SigEst*solve(t(XX)%*%XX)[2,2]
  		aest<-onDNA[2,]
  		Test.A.to.M<-aest/sqrt(VarEst)
  		Pval.A.to.M<-2*pt(abs(Test.A.to.M),df=length(Y)-numExplan,lower.tail=F)
		par2<-colSums(residDNA*c(residY))
  		Est2<-par2/par1
  		PredY<-residDNA*matrix(Est2,nrow=nrow(residDNA),ncol=ncol(residDNA),byrow=T)
  		newResid<-PredY-c(residY)
  		SSQ<-colSums(newResid^2)
  		SIGMA<-SSQ/(length(Y)-numExplan-1)
  		Tstat<-Est2*sqrt(par1/SIGMA)
  		Pval<-2*pt(abs(Tstat),lower.tail=F,df=length(Y)-numExplan-1)
  		return(cbind(c(Est2),c(Est2/Tstat),length(Y),numExplan,Pval,
               		aest,Test.A.to.M,Pval.A.to.M))
		}
