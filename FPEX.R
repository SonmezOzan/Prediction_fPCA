library(vars)
library(fda)


#########################
# Computation of FPE
#########################

# Given the residuals of the fitted model, this computes the FPE as
# (n+pd+r)/(n-(pd+r))*tr(\hat\Sigma), where \hat\Sigma is the covariance
# estimated from residuls. Each residual vector defines one row in 
# the res matrix. Here r is the dimension of the exogeneous variables

FPEX.trace<-function(res,p=2,r=2){
res=t(t(res))
d=length(res[1,])
n=length(res[,1])
if(d==1)
{
out=(p*d+n)/(n-p*d)*var(res)
}
else
{
out=(p*d+r+n)/(n-p*d-r)*sum(diag(cov(res)))
}
out
}


method.FPEX<-function(FTS,EXO,D=10,Pmax){

##########################################################
# Find the optimal order p and dimenion d in an FARX model
##########################################################

# determine length of FTS
m = dim(FTS[[1]])[2]

# determine the dimension of the exogeneous variables
# EXO is an nxr vector, we assume it to be centered
r = dim(EXO)[2]

# perform the PCA on 90% of FTS data
m = 100

# center the data

FTS=center.fd(FTS)

# we compute the FPC scores of the first m observations, with D harmonics
pca.FTS=pca.fd(FTS[1:m],nharm=D)

# compute the total variance
vartot=sum(pca.FTS$values)

# the matrix below will contain the different (FPE + vartot - var.explain) values
# for p in 0:Pmax and d in 1:D
values=matrix(0,D,(Pmax+1))

for(d in 1:D)
	{
	scores=pca.FTS$scores[,1:d]
	var.explain=sum(pca.FTS$values[1:d])
		for(p in 0:Pmax)
		{
			if(d==1)
			{
				res=arima(scores, xreg=EXO, order = c(p, 0, 0))$residuals
			}
			else
			{
#!!!!!!!!!!!!!!!!!! 
# p=0 has to be still treated
				if(p==0)
				{
					mean=t(matrix(rep(colMeans(scores),m),d))
					res=scores-mean
				}
				else
				{
					# we have to give colnames to the scores, 
					# otherwise we get warnings from vars resid function
					colnames(scores)<-as.character(seq(1:d))
					res=resid(VAR(scores,exogen=EXO,p=p,type="const"))
				}
			}
			values[d,p+1]=FPEX.trace(res=res,p=p,r=r)+vartot-var.explain
		}
	}

# compute the estimates hat.p and hat.d for optimal order p and dimension d
hat.p=(which.min(values)-1)%/%D
hat.d=which.min(values)%%D
if(hat.d==0) {hat.d=hat.d+D}

out=c(hat.p,hat.d)
out
}


