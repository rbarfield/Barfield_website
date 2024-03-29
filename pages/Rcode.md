---
layout: page
---

<h2>

Rcode

<p class="view">

<a href="https://cran.r-project.org/web/packages/CpGassoc/index.html">

<h3>

CpGassoc

</h3>

</a>

</p>

<p>

R package written under the supervision of Dr. Conneely. Designed to test for association between methylation at CpG sites across the genome and a phenotype of interest, adjusting for any relevant covariates.

</p>

<h3 id="DNAmpop">

DNA methylation population stratification

</h3>

<p>

Code and information on DNA methylation population stratification work done from Barfield et al 2014. Can be <a href="http://genetics.emory.edu/research/conneely/annotation-r-code.html">found here</a>. Paper can be found <a href="https://www.ncbi.nlm.nih.gov/pubmed/24478250">here</a>

</p>

<h3 id="EWASMEDIAT">

Fast EWAS Mediation

</h3>

<p>

Code to run quick mediation analysis across a set of DNA methylation markers. Gives association between CpG and Y given A and Xm as well as association between A and CpG given Xm. Will give for each CpG site. Assuming linear models for all variables. Xm is matrix of covariates adjusting for (not including intercept), CpG is an Nx p matrix with rows corresponding to sample and columns corresponding to DNA methylation marker.

</p>

<figure>

<figcaption>

Fast EWAS Mediation Function for Normal Outcome

</figcaption>

```{=html}
<pre>
    <code>
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
    </code>
  </pre>
```
</figure>

<h3 id="LDA MR Code">

LDA MR-Egger Code

</h3>

<p>

Code to run LDA MR-Egger Regression as discussed in Barfield et al (2018). X is the vector of joint eQTL effects. Y is the vector of joint GWAS effects.W is the inverse of the covariance of the joint GWAS effects (i.e. var(Y)).

</p>

<figure>

<figcaption>

LDA MR-Egger Code

</figcaption>

```{=html}
<pre>
    <code>
    LDA.MREgger<-function(X,Y,W){
     bX<-cbind(1,X)
     bread<-solve(crossprod(bX,W)%*%bX)
     theEsts<-bread%*%crossprod(bX,W%*%Y)
     theresid<-c(Y-theEsts[1]-X*theEsts[2])
     Sig.Est<-c(crossprod(theresid,W%*%theresid))/(length(X)-2)
     finresults<- cbind(theEsts,diag(bread)*Sig.Est)
     TestStat<-theEsts/sqrt(finresults[,2])
     Pvals<-2*pt(abs(TestStat),df = nrow(bX)-2,lower.tail = F)
     return(cbind(finresults,TestStat,Pvals))
     }
    </code>
  </pre>
```
</figure>
