library(fda)

# FUNCTION PREDICTS FAR WITH KOKOSZKA & REIMHERR METHODOLOGY
# Input: 
# z = functional time series, 
# Pmax = maximum order, 
# nharmy = maximum number of harmonics to determine dy 
# nharmx = maximum number of harmonics to determine dx
# sig = Pmax-dim vector with significance levels for the individual tests
# ExpVar = % of variance to be explained

method.KR <- function(z, Pmax=NULL, nharmx=NULL, nharmy=NULL, sig=NULL, ExpVar=NULL){

#####################################################################################
# DEFAULT OPTIONS
#####################################################################################

# if no Pmax is specified use upper limit of 3 tests
if(is.null(Pmax)){ Pmax = 3 }

# if number of harmonics for x and y are not specified use 60 and 30 respectively
if(is.null(nharmx)){ nharmx = 30 }
if(is.null(nharmy)){ nharmy = 10 }

# if significance levels are not specified use 5%, 3% and 2% for the first, second and third test, respectively
if(is.null(sig)){ sig = c(.05,.03,.02) } # to achieve approximately 10% level for the three tests

# if not otherwise specified, use 80% of variation explained to select number of FPC to be used
if(is.null(ExpVar)){ ExpVar = .8 } # % of variance to be explained

#####################################################################################
# PREPROCESSING OF FUNCTIONAL TIME SERIES
#####################################################################################

N = 100 # number of functions used for order selection
y = z[1:N] # first 100 functions

H = dim(z$coefs)[2]-N # number of functions used for prediction
yp = z[(N+1):(N+H)] # Remaining 10% of functions for prediction

K = 48 # number of evenly spaced sample points per curve
t = (0:(K-1))/(K-1) # points at which functional time series y is evaluated

zd = eval.fd(t,z)
yd = eval.fd(t,y) # values of discretized y at grid t

nb = 21 # number of basis functions for the data
fbf = create.fourier.basis(rangeval=c(0,K)/K, nbasis=nb) # basis for data

#####################################################################################
# FPCA FOR Y
#####################################################################################

ydfun = Data2fd(y=yd, argvals=t, fbf) # functions generated from discretized y
ydfun = center.fd(ydfun)
y.fd = pca.fd(ydfun,nharm=nharmy)

dy = min(which(cumsum(y.fd$varprop)>ExpVar)) # number of harmonics for y
dy = max(dy,2) # use at least dy = 2 because dy = 1 is not interpreted as a matrix(?)
#per = sum(y.fd$varprop[1:dy])

y.fd = pca.fd(ydfun,nharm=dy) # repeat FPCA with dy harmonics
Y = y.fd$scores
Y.harm = y.fd$harmonics
per = sum(y.fd$varprop)

#Y = y.fd$scores[,1:dy]
#Y.harm = y.fd$harmonics[1:dy]

#####################################################################################
# DEFINITION OF OUTPUT MATRIX
#####################################################################################

P = Pmax # Upper bound for order p to be tested

OUT = matrix(0,nrow=P+1,ncol=6)
colnames(OUT)<-c("Test Stat","P-Value","mean(MSE)", "med(MSE)","sd(MSE)","rej")
rownames(OUT)<-c("p=0","p=1","p=2","p=3")

#####################################################################################
# PREDICTION FOR p=0
#####################################################################################

MSE0 = c()
for (j in 1:H){
	MSE0 = c(MSE0,inprod(yp[j],yp[j]))
	}

OUT[1,1:5] = c(NA,NA,mean(MSE0),median(MSE0),sd(MSE0))
OUT[1,6] = 1

#####################################################################################
# LOOP FOR REPEATED TESTS 
#####################################################################################

for(p in 1:P){

#####################################################################################
# FPCA FOR X (THE STACKED Y)
#####################################################################################
	
	tp = (0:(p*K-1))/(p*K-1)
	
	xd = yd[,p:(N-1)]
	if(p > 1){
		for (j in 2:p){xd = rbind(xd,yd[,(p-j+1):(N-j)])} # stringing the p observations together
	}

	xdfun = Data2fd(y=xd, argvals=tp, fbf)
	xdfun = center.fd(xdfun)
	dx = dy*p # number of harmonics for x
	#dx = min(which(cumsum(pca.fd(xdfun,nharm=nharmx)$varprop)>ExpVar))
	
	X.fd = pca.fd(xdfun,dx)
	X = X.fd$scores
	X.harm = X.fd$harmonics
	
	Y1 = Y[(p+1):N,] # delete first p elements of Y to make dimensions match
	
#####################################################################################
# SOLVING THE FUNCTIONAL LINEAR MODEL
#####################################################################################

	lm1 = lm(Y1~X)
	Ceps = (1/(N-p-dx-1))*(t(lm1$res)%*%lm1$res)

#####################################################################################
# COMPUTATION OF TEST STATISTIC AND P-VALUE
#####################################################################################
	
	M = 400 # number of points to approximate integral of harmonics over [(p-1)/(p),1] interval
	t.tmp = seq(from=(p-1)/(p), to = 1, length.out = (M+1))
	t.tmp = t.tmp[2:(M+1)]
	v.tmp = eval.fd(t.tmp,X.fd$harmonics)
	V = t(t(v.tmp)%*%v.tmp)/(M*p)
	eg.v = eigen(V)
	dv = length(which(eg.v$values>0.9))
	alpha = eg.v$vectors[,(1:dv)]
	lambda = diag(t(X)%*%X)/N
	Psi = diag((N*lambda)^(-1))%*%(t(X)%*%Y1)
	I_dy = diag(rep(1,times=dy))
	av_tmp = as.vector(t(alpha)%*%Psi)
	
	Tst.Stat = (N-p)*av_tmp%*%solve( (I_dy%x%t(alpha))%*%(Ceps%x%diag(1/lambda))%*%(I_dy%x%alpha) )%*%av_tmp
	
	PValue = 1-pchisq(Tst.Stat,df=(dv*dy))
	
#####################################################################################
# PREDICTION FOR p 
#####################################################################################

	xdp = zd[,N:(N+H-1)] 
	if(p > 1){
		for (j in 2:p){ xdp = rbind(xdp,zd[,(N+1-j):(N+H-j)]) }
	}
		
	xdpfun = Data2fd(y=xdp, argvals=tp, fbf)
	xdpfun = center.fd(xdpfun)
		
	Xp = inprod(xdpfun,X.harm)
	yhat.v =  Xp %*% Psi
		
	MSE = c()
	for (j in 1:H){
		yhat = 0*Y.harm[1]
		for (i in 1:dy){ yhat = yhat + yhat.v[j,i]*Y.harm[i] }
		MSE = c(MSE,inprod(yp[j]-yhat,yp[j]-yhat))
				
		OUT[p+1,1] = Tst.Stat
		OUT[p+1,2] = PValue
		OUT[p+1,3] = mean(MSE)
		OUT[p+1,4] = median(MSE)
		OUT[p+1,5] = sd(MSE)
		if(PValue<sig[p]) {OUT[p+1,6] = 1} 
	}
}

#####################################################################################
# PREPARE OUTPUT 
#####################################################################################

OUT1 = matrix(NA,nrow=1,ncol=9) # more columns needed if Pmax is larger than 3
# colnames(OUT1)<-c("p.hat","dy","mean(MSE)","med(MSE)","sd(MSE)","P-Value","P-Value","P-Value","varprop")

OUT1[1,2] = dy
OUT1[1,9] = per

if(min(OUT[,6]) == 0){	
OUT1[1,1] = which.min(OUT[,6])-2
OUT1[1,3:5] = OUT[OUT1[1,1]+1,3:5]
OUT1[1,6:(6+OUT1[1,1])] = t(OUT[2:(2+OUT1[1,1]),2]) 
}
else {
OUT1[1,1] = Pmax 
OUT1[1,3:5] = OUT[OUT1[1,1]+1,3:5]
OUT1[1,6:(6+OUT1[1,1]-1)] = t(OUT[2:(2+OUT1[1,1]-1),2]) 
}

OUT1 # Return 
}