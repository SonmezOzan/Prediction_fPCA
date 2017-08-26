###############################################
# Function that computes the matrix X from the 
# matrix Y (pxn) of principal components

X<-function(Y){
p=length(Y[,1])
n=length(Y[1,])
X0=NULL
for(i in 1:(n-1)){
X0=rbind(X0,diag(1,p)%x%t(Y[,i]))
}
X0
}

###############################################
# Function that computes the matrix Z from the
# matrix Y (pxn) of principal components

Z=function(Y){
n=length(Y[1,])
Z0=NULL
for(i in 2:n){
Z0=cbind(Z0,t(Y[,i]))
}
t(Z0)
}


###############################################
# Function that computes \hat\beta from Y by 
# least squares

hat.beta.ls<-function(Y){
((solve(t(X(Y))%*%X(Y)))%*%t(X(Y)))%*%Z(Y)
}

###############################################
# Function that computes \hat B from Y

hat.B.ls<-function(Y,skip=TRUE){
n=length(Y[1,])
bb=hat.beta.ls(Y)
p=sqrt(length(bb))
t(matrix(bb,ncol=p))
}


###############################################
# Function that computes \Delta from Y

Delta<-function(Y){
n=length(Y[1,])
B=hat.B.ls(Y)
Delta0=NULL
Z0=NULL
for(i in 2:n){
Z0=cbind(Z0,Y[,i]-B%*%Y[,i-1])
}
Z0
}

###############################################
# Computes cross covariances of data matrix D (pxn)
# at given lag 

cro.cov<-function(D,lag){
n=length(D[1,])
mean=1/n*colSums(D)
D0<-D-mean
D1<-D0[,(1+lag):n]
D2<-D0[,1:(n-lag)]
D1%*%t(D2)/n
}

###############################################
# Computes matrix (nxn) having all zeros but
# in k-th upper offdiagonal, where it has ones 

off.diag<-function(n,k){
k<-n-k
d0<-diag(1,k)
for(i in k:(n-1)){
d0<-rbind(cbind(rep(0,i),d0),rep(0,i+1))
}
d0
}




###########################################################
#
# Prediction with regressors
# We define function that predicts Yn+1 using the
# elements n-m+1,....n of Y and a vector of regressors R. For estimtion
# of the necessary covariances we use the entire sample

###############################################
# Computes cross covariances of data matrix Y (pxn) and R (rxn)
# at given lag 

cro.cov.2<-function(Y,R,lag){
n=length(Y[1,])
meanY=1/n*colSums(Y)
meanR=1/n*colSums(R)
Y0<-Y-meanY
R0<-R-meanR
Y1<-Y0[,(1+lag):n]
R1<-R0[,1:(n-lag)]
Y1%*%t(R1)/n
}


##############################################
fcX<-function(Y,R,p){

r=length(R[,1])
d=length(Y[,1])
n=length(Y[1,])
# the left hand side of the pred equations excluding covariates for the moment
lhs.Y=matrix(0,d*p,d*p)
# the left hand side of the pred equations 
lhs=matrix(0,d*p+r,d*p+r)
# the right hand side of the pred equations
rhs=matrix(0,d,d*p+r)


for(i in 1:p){
rhs[,((i-1)*d+1):(i*d)]=cro.cov(Y,i)
}
rhs[,(d*p+1):(d*p+r)]=cro.cov.2(Y,R,1)

lhs.Y=lhs.Y+diag(p)%x%cro.cov(Y,0)
if(p>1){
	for(i in 1:(p-1)){
	lhs.Y=lhs.Y+off.diag(p,i)%x%cro.cov(Y,i)+t(off.diag(p,i))%x%t(cro.cov(Y,i))
	}
}
lhs[1:(d*p),1:(d*p)]=lhs.Y
for(i in 1:p){
lhs[((i-1)*d+1):(i*d),(d*p+1):(d*p+r)]<-cro.cov.2(Y,R,i-1)
lhs[(d*p+1):(d*p+r),((i-1)*d+1):(i*d)]<-cro.cov.2(R,Y,i-1)
}
lhs[(d*p+1):(d*p+r),(d*p+1):(d*p+r)]<-cro.cov.2(R,R,0)


PHI=rhs%*%solve(lhs)

PHI
}


