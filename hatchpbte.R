#R-code to estimate phos using maximum likelihood theory. Statistical properties of the estimator
#are estimated using maximum likelihood theory and bootstrapping. This code treats the general case
#of inputs from several hatcheries with potentially different visible marking fractions and different
#parentage-based tag fractions.
#AUTHOR: Richard A. Hinrichsen
#CONTACT: rich@hinrichsenenvironmental.com

#Variables and parameters used in the analysis
#inputs
#Nsamp = total sample size
#x1 = number of marked spawners in sample of size Nsamp
#n = total subsample sizes (number tested for PBT)
#n1 = number of marked spawners tested for PBT
#y = number of marked spawners tested that were PBT (hatchery-specific)
#z = number of unmarked spawners tested that were PBT (hatchery-specific)
#lambda = marking rate (lambda) (hatchery-specific)
#ppbt=fraction of marked fish that are also parentage-based tagged (hatchery-specific)
#BOOT = FALSE for theoretical results, TRUE for Boostrap results
#NBOOT = number of Bootstrap simulations (needed if BOOT=TRUE)
#
#
#Select intermediate variables
#nhatch=number of hatcheries supplying spawners in the wild
#I = Observed Fisher Information Matrix
#x2 = number of unmarked spawners in sample of size Nsamp
#n2 = number of unmarked spawners tested for PBT

#Results
#phosi= maximum likelihood estimates of the proportions of hatchery-origin spawners
#SE.phosi= standard errors of the estimates of the proportions of hatchery-origin spawners
#phos = maximum likelihood estimate of the proportion of hatchery-origin spawners
#SE.phos = standard error (SE) of the estimate of phos
#CV.phos = Coefficient of variation of the estimate of phos
#BIAS.phos = relative bias of the estimate of phos (NA if BOOT=FALSE)

#top level function
phos.pbte.main<-function(Nsamp=1000,x1=50,n=200,n1=50,
y=c(16,15,8,8),
z=c(0,1,14,14),
lambda=c(1,.95,.5,.5),
ppbt=c(.95,.95,.95,.95),
BOOT=FALSE,NBOOT=1000){
x2<-Nsamp-x1
n2<-n-n1
check.input(x1=x1,x2=x2,n1=n1,n2=n2,y=y,z=z,lambda=lambda,ppbt=ppbt,BOOT=BOOT,NBOOT=NBOOT)
if(!BOOT){res<-phos.pbte.estimates(x1=x1,x2=x2,n1=n1,n2=n2,y=y,z=z,lambda=lambda,ppbt=ppbt,suppress=FALSE)}
if(BOOT){res<-phos.pbte.estimates2(NBOOT=NBOOT,x1=x1,x2=x2,n1=n1,n2=n2,y=y,z=z,lambda=lambda,ppbt=ppbt)}
return(res)
}

#check input to phos.pbte.main,phos.pbte.estimates,phos.pbte.estimates2
check.input<-function(x1,x2,n1,n2,y,z,lambda,ppbt,BOOT,NBOOT)
{
Nsamp<-x1+x2
if(!is.logical(BOOT))stop("BOOT must be TRUE or FALSE")
if(floor(Nsamp)!=Nsamp){stop("Nsamp must be a positive integer")}
if(Nsamp<=0){stop("Nsamp must be a positive integer")}
if(BOOT){
if(floor(NBOOT)!=NBOOT){stop("NBOOT must be a positive integer")}
if(NBOOT<=0){stop("NBOOT must be a positive integer")}
}
#check dimension of inputs
k2<-length(lambda);k3<-length(ppbt);k4<-length(y);k5<-length(z)
mytest<-abs(k2-k3)+abs(k2-k4)+abs(k2-k5)
if(mytest>0) stop("dimensions of lambda, ppbt, y, and z must match")
#check that subsample size is less than sample size
n<-n1+n2
if(n>Nsamp)stop(paste("n=n1+n2 must be less than or equal to Nsamp=x1+x2=",Nsamp))
#check sample and subsample subsample
if(x1<0)stop("x1 must be nonnegative")
if(x2<0)stop("x2=Nsamp-x1 must be nonnegative")
if(n1<0)stop("n1 must be nonnegative")
if(n2<0)stop("n2=n-n1 must be nonnegative")
if(n1>x1)stop(paste("n1 must not exceed the number of VM spawners in sample=",x1))
if(n2>x2)stop(paste("n2 must not exceed the number of ~VM spawners in the sample=",x2=Nsamp-x1))
#check that all ppbts are between zero and 1.0
if(sum(ppbt<0))stop("ppbts must all be greater than or equal to zero")
if(sum(ppbt>1))stop("ppbts must all be less than or equal to 1.0")
#check that all lambdas are between zero and 1.0
if(sum(lambda<0))stop("lambdas must all be greater than or equal to zero")
if(sum(lambda>1))stop("lambdas must all be less than or equal to one")
#check for consistency between y,x1,z,and x2,and lambda
if(sum(y<0))stop("ys must all be greater than or equal to zero")
if(sum(z<0))stop("zs must all be greater than or equal to zero")
if(sum(y)>n1)stop("ys must not sum to greater than n1")
if(sum(z)>n2)stop("zs must not sum to greater than n2=n-n1")
iii<-(y>0)&(lambda==0)
if(sum(iii))stop("y cannot be > 0 if no releases are visibly marked")
iii<-(z>0)&(lambda==1)
if(sum(iii))stop("z cannot be > 0 if all releases are visibly marked")
if((sum(lambda)==0)&(x1>0))stop("x1 cannot be > 0 if no hatchery releases are visibly marked")
if((sum(lambda*(1-ppbt))==0)&(sum(y)!=n1)){
stop("sum(y) should equal n1 since all ppbt=1 for all hatcheries with lambda>0")}
#When a group of hatcheries has zero expected PBT recoveries, they must use constant VM fractions
iii<-n1*lambda*ppbt+n2*(1-lambda)*ppbt==0
onelambda2<-sum(mean(lambda[iii])==lambda[iii])==sum(iii)
onelambda2<-onelambda2&(mean(lambda[iii])>0)
if((!onelambda2)&sum(iii))stop("Expected tag recoveries must not be zero when marking fractions differ")
return(NULL)
}

#Maximum Likelihood Theory results
phos.pbte.estimates<-function(x1,x2,n1,n2,y,z,lambda,ppbt,suppress){
nhatch<-length(lambda)
Nsamp<-x1+x2

#An important case for combining cells occurs when the E(y+z)=0
#in this case, constant lambdas over these cells saves the estimation.
#This also takes care of the case of a single lambda for all hatcheries
iii<-n1*lambda*ppbt+n2*(1-lambda)*ppbt==0
onelambda2<-sum(lambda[iii]==mean(lambda[iii]))==sum(iii)
if((sum(iii)>1)&onelambda2){
if(!suppress)warning("collapsing cells with expected tag recoveries of zero into single cell since lambda is constant")
lambda1.new<-mean(lambda[iii])
ppbt1.new<-0.0
lambda.new<-c(lambda1.new,lambda[!iii])
ppbt.new<-c(ppbt1.new,ppbt[!iii])
y.new<-c(0,y[!iii])
z.new<-c(0,z[!iii])
nhatch.new<-length(lambda.new)
res<-phos.pbte.estimates(x1=x1,x2=x2,n1=n1,n2=n2,y=y.new,z=z.new,lambda=lambda.new,ppbt=ppbt.new,suppress=suppress)
phosi<-rep(NA,nhatch)
SE.phosi<-rep(NA,nhatch)
phosi[!iii]<-res$phosi[2:nhatch.new]
SE.phosi[!iii]<-res$SE.phosi[2:nhatch.new]

myres<-list(BOOT=FALSE,
NBOOT=NA,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi,
SE.phosi=SE.phosi,
phos=res$phos,
SE.phos=res$SE.phos,
CV.phos=res$CV.phos,
BIAS.phos=NA)
return(myres)
}

#get initial estimate of phosi
phosi.init<-init(x1,x2,n1,n2,y,z,lambda,ppbt,suppress=suppress)

myres<-get.estimates2(phosi.init,x1,x2,n1,n2,y,z,lambda,ppbt,suppress=suppress)
phos<-myres$phos
phos.var<-myres$phos.var
phosi<-myres$phosi
SE.phosi<-sqrt(myres$phosi.var)
SE.phos<-sqrt(phos.var)
CV.phos<-SE.phos/phos
SE.phos<-as.numeric(SE.phos)
CV.phos<-as.numeric(CV.phos)

myres<-list(BOOT=FALSE,
NBOOT=NA,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi,
SE.phosi=SE.phosi,
phos=phos,
SE.phos=SE.phos,
CV.phos=CV.phos,
BIAS.phos=NA)
return(myres)
}

#Fisher Information Matrix (general case)
getI<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
theta1<-n1/Ex1
theta2<-n2/Ex2

v<-lambda
I<- v%*%t(v)*Nsamp*(1-theta1)/sum(v*phosi)
I<-I +v%*%t(v)*Nsamp*(1-theta2)/(1-sum(v*phosi))

v<-(1-ppbt)*lambda
I<-I+v%*%t(v)*Nsamp*theta1/sum(v*phosi)

v<-lambda*(1-ppbt)+ppbt
I<-I+v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))

#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*(theta1*lambda+theta2*(1-lambda))/phosi
diag(I)<-mydiag
return(I)
}

#Fisher Information matrix used when sum((1-ppbt)*lambda) is zero
getI2<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex1<-Nsamp*sum(lambda*phosi)
Ex2<-Nsamp-Ex1
theta1<-n1/Ex1
theta2<-n2/Ex2

v<-lambda
I<- v%*%t(v)*Nsamp*(1-theta1)/sum(v*phosi)
I<-I +v%*%t(v)*Nsamp*(1-theta2)/(1-sum(v*phosi))

v<-ppbt
I<-I+v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))

#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*(theta1*lambda+theta2*(1-lambda))/phosi
diag(I)<-mydiag
return(I)
}

#Fisher Information matrix (used when all lambdas are zero)
getI3<-function(Nsamp,n1,n2,lambda,ppbt,phosi){
Ex2<-Nsamp
theta2<-n2/Ex2
v<-ppbt
I<-v%*%t(v)*Nsamp*theta2/(1-sum(v*phosi))
#fix diagonal
mydiag<-diag(I)+Nsamp*ppbt*theta2/phosi
diag(I)<-mydiag
return(I)
}

#Bootstrap estimates of standard error and bias
#consider cases where some of the phosi are missing
phos.pbte.estimates2<-function(NBOOT,x1,x2,n1,n2,y,z,lambda,ppbt){
nhatch<-length(lambda)
Nsamp<-x1+x2
#get MLE using theoretical results
res<-phos.pbte.estimates(x1,x2,n1,n2,y,z,lambda,ppbt,suppress=FALSE)
phosi<-res$phosi
phosi.orig<-res$phosi
phos<-res$phos
phos.sim<-rep(NA,NBOOT)
phosi.sim<-matrix(NA,nrow=NBOOT,ncol=nhatch)
if(!is.na(phos)){
for(ii in 1:NBOOT){
iii<-is.na(phosi)
phosi[iii]<-rep(phos-sum(phosi,na.rm=T),sum(iii))/sum(iii)
mysim<-pbtsim1(phosi,Nsamp,n1,n2,res$lambda,res$ppbt)
my.n1<-n1
my.n2<-n2
if(mysim$x1<n1)my.n1<-mysim$x1
if(mysim$x2<n2)my.n2<-mysim$x2
res<-phos.pbte.estimates(x1=mysim$x1,x2=mysim$x2,n1=my.n1,n2=my.n2,
y=mysim$y,z=mysim$z,lambda=lambda,ppbt=ppbt,suppress=TRUE)
phos.sim[ii]<-res$phos
phosi.sim[ii,]<-res$phosi

}}

SE.phosi<-apply(phosi.sim,c(2),var,na.rm=T)
SE.phosi<-sqrt(SE.phosi)
SE.phos<-sqrt(var(phos.sim,na.rm=T))
CV.phos<-SE.phos/phos
mymean<-mean(phos.sim,na.rm=T)
BIAS.phos<-(mymean-phos)/phos

myres<-list(BOOT=TRUE,
NBOOT=NBOOT,
Nsamp=Nsamp,
x1=x1,
x2=x2,
n=n1+n2,
n1=n1,
n2=n2,
y=y,
z=z,
lambda=lambda,
ppbt=ppbt,
phosi=phosi.orig,
SE.phosi=SE.phosi,
phos=phos,
SE.phos=SE.phos,
CV.phos=CV.phos,
BIAS.phos=BIAS.phos)
return(myres)
}

#simulate data when phosi is available
pbtsim1<-function(phosi,Nsamp,n1,n2,lambda,ppbt){
m<-length(phosi)
#first get binomial sample of fish marked and unmarked
P<-sum(phosi*lambda)
x1<-rbinom(n=1,size=Nsamp,prob=P)
x2<-Nsamp-x1
#next use multinomial distribution to simulate
#how many fish are pbt and how many are not pbt
py<-ppbt*phosi*lambda/P
pz<-ppbt*phosi*(1-lambda)/(1-P)
if(P>0){y<-rmultinom(n=1,size=min(n1,x1),prob=c(py,1-sum(py)))}
if(P==0){y<-matrix(0,ncol=1,nrow=length(phosi))}
if(P<1){z<-rmultinom(n=1,size=min(n2,x2),prob=c(pz,1-sum(pz)))}
if(P==1){z<-matrix(0,ncol=1,nrow=length(phosi))}

return(list(Nsamp=Nsamp,x1=x1,x2=x2,n1=n1,n2=n2,y=y[1:m,1],z=z[1:m,1]))
}

#Use R.A. Fisher"s scoring algorithm to estimate phosi
#phosi represents the intial guess on input
get.estimates2<-function(phosi.init,x1,x2,n1,n2,y,z,lambda,ppbt,suppress){
Nsamp<-x1+x2
NTRIAL<-100
tolx<-1.e-5
nhatch<-length(phosi.init)
nhatch.orig<-nhatch
w<-y+z
nas<-rep(NA,nhatch)
phosi<-nas
phosi.var<-nas
phos<-NA
phos.var<-NA

#in a rare case init can return zeroes even though the true estimate is not zero
#which would defeat Fisher"s scoring method
jjj<-phosi.init==0
phosi.init[jjj]<-phosi.init[jjj]+.00001
#find zeroes that occur when lambda=0 and n1*lambda*ppbt+n2*(1-lambda)*ppbt>0
jjj<-(n1*lambda*ppbt+n2*(1-lambda)*ppbt>0)&(lambda==0)&(w==0)
if(sum(jjj)>=1){
y<-y[!jjj]
z<-z[!jjj]
ppbt<-ppbt[!jjj]
lambda<-lambda[!jjj]
phosi.init<-phosi.init[!jjj]
nhatch<-length(y)}

#check for special cases
#check to see if all pbts are one (special case)
pbttest<-sum(ppbt==1)==nhatch
#modify slightly
pbttest<-sum((ppbt-1)*lambda)==0
#check to see if all lambdas are zero (special case)
lambdatest<-sum(lambda==0)==nhatch

#get the right Fisher Information function and dlike function
if(lambdatest){
my.getI<-getI3
dlike<-dlike3
}
else{

if(pbttest){
my.getI<-getI2
dlike<-dlike2
}
else{
my.getI<-getI
dlike<-dlike1
}
}

#use R.A. Fisher"s scoring algorithm to find where the partial derivatives of the
#log-likelihood are zero. Use Fisher Information matrix in Newton"s method to approx -Hessian
phosi<-phosi.init
errf<-0.0
alpha<-0.9
for(ii in 1:NTRIAL){
I<-my.getI(Nsamp,n1,n2,lambda,ppbt,phosi)
df<-dlike(phosi,x1,x2,n1,n2,y,z,lambda,ppbt)
size<-prod(dim(I))
if(size==0){
if(!suppress)warning("dimension of I is 0 x 0")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
if(is.na(rcond(I))){
if(!suppress)warning("condition number of I is NA")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
if(rcond(I)<1.e-15){
if(!suppress)warning("computationally singular information matrix")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
delx<-solve(I)%*%df
phosi<-phosi+delx*(1-alpha)
phosi<-abs(phosi)
errx<-sum(abs(delx))/sum(abs(phosi))
alpha<-alpha*alpha
if(errx<=tolx)break
}
if(ii==NTRIAL){
if(!suppress)warning("maximum number of iterations was reached")
return(list(phosi=nas,phosi.var=nas,phos=NA,phos.var=NA))}
phos<-sum(phosi)
e<-rep(1,nhatch)
myvar<-solve(I)
phos.var<-t(e)%*%myvar%*%e
phosi.var<-diag(myvar)
full.phosi.var<-rep(0,nhatch.orig)
full.phosi<-rep(0,nhatch.orig)
#recall that jjj represents hatcheries with estimates of phosi=0
#reduction occurs only when sum(jjj>=1)
if(sum(jjj)>=1){
full.phosi.var[!jjj]<-phosi.var
full.phosi[!jjj]<-phosi
}
if(sum(jjj)<1){
full.phosi.var<-phosi.var
full.phosi<-phosi

}
return(list(phosi=full.phosi,phosi.var=full.phosi.var,phos=phos,phos.var=phos.var))
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used in most general case
dlike1<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum1<-sum(lambda*phosi)
sum2<-sum((1-ppbt)*lambda*phosi)
sum3<-sum(phosi*ppbt)
res<-lambda*(x1-n1)/sum1-lambda*(Nsamp-x1-n2)/(1-sum1)+y/phosi
res<-res+(n1-sum(y))*(1-ppbt)*lambda/sum2+z/phosi
res<-res-(n2-sum(z))*(lambda*(1-ppbt)+ppbt)/(1-sum2-sum3)
return(res)
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used in the case where sum(lambda*(1-ppbt))===0
dlike2<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum1<-sum(lambda*phosi)
sum3<-sum(phosi*ppbt)
res<-lambda*(x1-n1)/sum1-lambda*(Nsamp-x1-n2)/(1-sum1)+y/phosi
res<-res+z/phosi
res<-res-(n2-sum(z))*ppbt/(1-sum3)
return(res)
}

#get gradient of the log likelihood function
#phosi is the current best estimate of phosi
#used when lambda=0 at all hatcheries
dlike3<-function(phosi,x1,x2,n1,n2,y,z,lambda,ppbt){
#estimate phosi
Nsamp<-x1+x2
sum3<-sum(phosi*ppbt)
res<-z/phosi
res<-res-(n2-sum(z))*ppbt/(1-sum3)
return(res)
}

#get initial estimates of phosi
#by equating x1,y,and z to their expectations
#this yields 2*nhatch +1 equations with nhatch unknowns
#which is, in general, overdetermined.
init<-function(x1,x2,n1,n2,y,z,lambda,ppbt,suppress=FALSE){
Nsamp<-x1+x2
nhatch<-length(lambda)
A1<-lambda
B1<-x1/Nsamp
A2<--n1*diag(lambda*ppbt,ncol=nhatch,nrow=nhatch)
LAMBDAMAT<-matrix(lambda,ncol=nhatch,nrow=nhatch)
LAMBDAMAT<-t(LAMBDAMAT)
YDIAG<-diag(y,nrow=nhatch,ncol=nhatch)
ZDIAG<-diag(z,nrow=nhatch,ncol=nhatch)
A2<-A2+YDIAG%*%LAMBDAMAT
B2<-rep(0,nhatch)
A3<-n2*diag((1-lambda)*ppbt,nrow=nhatch,ncol=nhatch)
A3<-A3+ZDIAG%*%LAMBDAMAT
B3<-z
A<-rbind(A1,A2,A3)
B<-c(B1,B2,B3)
if(rcond(t(A)%*%A)<1.e-15){
if(!suppress)warning("matrix t(A)%A in init() is computationally singular")}
phosi<-solve(t(A)%*%A)%*%t(A)%*%B
return(abs(phosi))
}