#library(vars)
#library(fda)


#########################
# Computation of FPE
#########################

# Given the residuals of the fitted model, this computes the FPE as
# (n+pd)/(n-pd)*tr(\hat\Sigma), where \hat\Sigma is the covariance
# estimated from residuls. Each residual vector defines one row in 
# the res matrix 

FPE.trace<-function(res,p=2){
res=t(t(res))
d=length(res[1,])
n=length(res[,1])
if(d==1)
{
out=(p*d+n)/(n-p*d)*var(res)
}
else
{
out=(p*d+n)/(n-p*d)*sum(diag(cov(res)))
}
out
}


method.FPE<-function(FTS,D=21,Pmax){

#######################################
# Our procedure applied to given FAR(2)
#######################################

# determine length of FTS
n = dim(FTS[[1]])[2]

# perform the PCA on 90% of FTS data
m = 100

# center the data

# we compute the mean from the first m observations
  means<-rowMeans(FTS[[1]][,1:m])

# and subtract it from all observations
  FTS[[1]][,1:n]<-(FTS[1:n][[1]]-means)  	

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
				res=arima(scores, order = c(p, 0, 0))$residuals
			}
			else
			{
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
					res=resid(VAR(scores,p=p,type="const"))
				}
			}
			values[d,p+1]=FPE.trace(res=res,p=p)+vartot-var.explain
		}
	}

# compute the estimates hat.p and hat.d for optimal order p and dimension d
hat.p=(which.min(values)-1)%/%D
hat.d=which.min(values)%%D
if(hat.d==0) hat.d=hat.d+D


# compute the 1 step ahead prdiction of the VAR for the coefficients
MSE=c()

for(k in (m+1):n)
	{
	if(hat.d==1)
	{
		VAR.pre=predict(arima(pca.FTS$scores[,1],order=c(hat.p,0,0)),n.ahead=1)$pred[1]
		yhat=VAR.pre
	}
	else
	{
		vec.ts=pca.FTS$scores[,1:hat.d]
		if(p==0)
			{
			yhat=colMeans(vec.ts)
			}
		else
			{
			# we need to give colnames to the scores 
			# to avoid warnings from vars predict function below
			colnames(vec.ts)<-as.character(seq(1:hat.d))
			VAR.pre=predict(VAR(vec.ts,p=hat.p),n.ahead=1,type="const")$fcst
			yhat=c()
			for(i in 1:hat.d)
				{
				yhat=c(yhat,VAR.pre[[i]][1])
				}
			}
	}
	FAR.pre=pca.FTS$harmonics[1]*0
	for(i in 1:hat.d)
		{
		FAR.pre=FAR.pre+yhat[i]*pca.FTS$harmonics[i]
		}

	MSE=c(MSE,inprod(FTS[k]-FAR.pre,FTS[k]-FAR.pre))

	pca.FTS$scores=rbind(pca.FTS$scores,inprod(FTS[k],pca.FTS$harmonics))

}

out = t(c(hat.p,hat.d,values[hat.d,hat.p+1],mean(MSE),median(MSE),sd(MSE),sum(pca.FTS$varprop[1:hat.d])))

out
}

