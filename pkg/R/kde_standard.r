#(c) Martin Rittner 2014

skde<-function(data,bandwidth,n,from,to,kernel="gaussian",...){
	#returns kde estimates at n equally spaced points between from and to,
	#calculated with density() from package stats, with given data, bandwidth
	#and kernel

	if(missing(bandwidth))bandwidth<-optimal_bw(data,n=2^10)
	if(missing(n))n<-2^12
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)
	tf<-approxfun(density(data,bw=bandwidth,kernel=kernel))
	ret<-tf(ev)
	return(matrix(c(ev,ret),nrow=2,byrow=TRUE))
}

#########################################
#following Jann, 2007
#this works and allows manual definition of kernel functions, but is slow...
jkde<-function(data,bandwidth,n,from,to,kernel="gaussian",...){
	data<-data[!is.na(data)]
	if(missing(bandwidth))bandwidth<-optimal_bw(data,n=2^10)/2
	if(missing(n))n<-2^11
	kernel<-match.arg(kernel,c("gaussian","epanechnikov","triangular","epan2"))
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)
	ret<-fk(data,ev,bandwidth,kernel)
	return(matrix(c(ret$x,ret$y),nrow=2,byrow=TRUE))
}

#adaptive kde - after Jann, 2007, following Abramson, 1982
# FIXME: automatic bandwidth calculation
akde<-function(data,bandwidth,n,from,to,kernel="gaussian",...){
	if(missing(bandwidth))bandwidth<-optimal_bw(data,n=2^10)
	if(missing(n))n<-2^12
	kernel<-match.arg(kernel,c("gaussian","epanechnikov","triangular"))
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)
	tf<-approxfun(x=ev,y=fak(data,ev,bandwidth,kernel)/length(data))
	ret<-tf(ev)
	return(matrix(c(ev,ret),nrow=2,byrow=TRUE))
}

# helper functions from JANN paper:
# (1)
fk<-function(data,vals,h,kernel){
	#n<-length(data)
	#zs<-outer(vals,data,FUN="-")/h
	# not the prettiest, but the fastest: do NOT pre-calculate anything, let R do
	# the optimisation on the fly
	return(approx(x=vals,xout=vals,
		y=sapply(vals,function(x){return(
			sum(sapply(data,function(Xi){return(
				Kfun((x-Xi)/h,kernel=kernel)/h
			)}))/length(data)
		)})
	)
	)
}

# (2)
fak<-function(data,vals,h,kernel){
	initial<-fk(data,vals,h,kernel)
	G<-exp(mean(log(c(initial))))
	tf<-approxfun(x=vals,y=initial/length(data))
	h2<-approxfun(x=data,y=h*sapply(data,function(Xi){return(h*sqrt(G/tf(Xi)))}))
	return(sapply(vals,function(x){
		return(sum(sapply(data,function(Xi){
			return(Kfun((x-Xi)/h2(Xi),kernel=kernel)/h2(Xi))
		})))
	}))
}

# (3)
lambda<-function(Xi,h,G,data,kernel){
	return(sqrt(G/fk(data,Xi,h,kernel)))
}

# after Table 1
Kfun<-function(z,kernel="gaussian"){
	# 	d<-switch(EXPR=c(kernel),
	# 		gaussian=dnorm(z),
	# 		epanechnikov=
	# 			if(abs(z)<sqrt(5)){
	# 				3/4*(1-z^2/5)/sqrt(5)
	# 			}else{0},
	# 		triangular=
	# 			if(abs(z)<1){
	# 				1-abs(z)
	# 			}else{0},
	# 		NULL
	# 	)
	#not as pretty, but >80% faster:
	d<-rep(0,length(z))
	if(kernel=="gaussian"){
		d<-dnorm(z)
	}else if(kernel=="epanechnikov"){
		ids<-abs(z)<sqrt(5)
		d[ids]<-0.75*(1-z[ids]^2/5)/sqrt(5)
	}else if(kernel=="triangular"){
		ids<-abs(z)<1
		d[ids]<-1-abs(z[ids])
	}else if(kernel=="epan2"){
		ids<-abs(z)<1
		d[ids]<-0.75*(1-z[ids]^2)
	}else if(kernel=="rectangular"){
		ids<-abs(z)<1
		d[ids]<-0.5
	}
	return(d)
}
