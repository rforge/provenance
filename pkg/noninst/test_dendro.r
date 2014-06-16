#############################################
## stacked KDE-Plots, minimal stacking height:

require(provenance)
data<-read.csv("/media/ldata/martin/LOESS/data/zrn_UPb/raw_data_Tb/Tb_all.csv",header=TRUE,stringsAsFactors=FALSE)
data<-as.list(data)
for(i in 1:length(data)){data[[i]]<-data[[i]][!is.na(data[[i]])]}

# ##for two series:
# ya<-t(pkde(data$Tb02,MIN=0))
# yb<-t(pkde(data$Tb09,MIN=0))
#
# afa<-approxfun(x=ya,rule=2)
# afb<-approxfun(x=yb,rule=2)
#
# xvals<-c(0:3500)
# dcomp<-data.frame(x=xvals,ya=afa(xvals),yb=afb(xvals))
#
# distance<-max((dcomp$ya-dcomp$yb))
# offset<-0.2*max(dcomp[2:length(dcomp)])
#
# ggplot(data=dcomp)+geom_line(aes(x=x,y=ya),colour="blue")+geom_line(aes(x=x,y=yb+distance+offset),colour="red")

# ## for all:
# nx<-2^12
# xvals<-seq(0,max(sapply(data,max)),length.out=nx)
#
# redata<-data.frame(x=xvals)
# for(i in seq(along.with=data)){
# 	askde<-t(pkde(data[[i]],MIN=0,bandwidth=20))
# 	#askde[,2]<-length(data[[i]])*askde[,2]/sum(askde[,2])
# 	askde[,2]<-askde[,2]/max(askde[,2])
# 	asfun<-approxfun(x=askde,rule=2)
# 	redata[names(data)[i]]<-asfun(xvals)
# }
#
# stdata<-redata
# labpos<-data.frame(sample=names(stdata)[2:length(stdata)],x=max(stdata$x),y=0)
# offset<-0.3*max(stdata[2:length(stdata)])
# for(i in (length(stdata)-1):2){
# 	distance<-max(stdata[[i+1]]-stdata[[i]])
# 	stdata[[i]]<-stdata[[i]]+distance+offset
# 	labpos[labpos$sample==names(stdata)[i],]$y<-min(stdata[[i]])
# }
#
# require(reshape2)
# meltdata<-melt(stdata,id.vars="x",variable.name="sample",value.name="d")
# ggplot(data=meltdata,aes(x=x,y=d))+geom_line(aes(group=sample))+geom_text(data=labpos,aes(x=x,y=y,label=sample),hjust=0,vjust=0)

#reorder to maximise stacking:
diss = dissimilarity(data)
closest<-hclust(as.dist(diss),method="ward.D")
# stdata<-redata[c("x",closest$labels[closest$order])]
# labpos<-data.frame(sample=names(stdata)[2:length(stdata)],x=min(stdata$x),y=0)
# offset<-0.3*max(stdata[2:length(stdata)])
# for(i in (length(stdata)-1):2){
# 	distance<-max(stdata[[i+1]]-stdata[[i]])
# 	stdata[[i]]<-stdata[[i]]+distance+offset
# 	labpos[labpos$sample==names(stdata)[i],]$y<-min(stdata[[i]])
# }
# meltdata<-melt(stdata,id.vars="x",variable.name="sample",value.name="d")
# ggplot(data=meltdata,aes(x=x,y=d))+geom_line(aes(group=sample))+geom_text(data=labpos,aes(x=x,y=y,label=sample),hjust=1.1,vjust=0)

# require(ggdendro)

## generate segments of tree based on pre-calculated coordinates for leaves
## e.g. to match age spectra...

pos<-2	# braching direction: 1 top-down, 2 right-to-left, 3 bottom-up, 4 left-to-right
htype<-1	# 1 uniform height steps, 2 calculated by hclust w/base 0, 3 same as 2, w/base = height from hclust

# clustree<-function(atree,orientation=c(1:4),type=c(1,2),positions=NULL){
# 	#takes a hclust structure and returns rendered coordinates to draw the tree,
# 	#according to orientation and type parameters
# 	#positions takes y-positions in the order of drawing layout (atree$order), i.e. ascending values
#
# 	if(class(atree)!="hclust")stop("atree must be of class hclust")
# 	stepval<-switch(orientation,1,1,-1,-1)
# 	angle<-switch(orientation,270,0,90,0)
# 	if(length(positions)==0){
# 		positions<-match(c(1:length(atree$labels)),atree$order)
# 	}else if(length(positions)==length(atree$labels)){
# 		#nothing to do!?
# 	}else{
# 		stop("invalid positions")
# 	}
# 	# TODO: if positions not given, optionally compute from height/distances in hclust
#
# 	subtree<-function(mergetable,curpos,branches){
# 		newheight<-curpos+1
# 		stree<-data.frame()
# 		markers<-NULL
# 		for(i in c(1:dim(mergetable)[1])){
# 			if(!all(is.na(mergetable[i,])) && mergetable[i,1]<0 && mergetable[i,2]<0){
# 				if(type==1){
# 					x<-c(branches[-mergetable[i,1],"height"],branches[-mergetable[i,2],"height"],newheight)
# 					xend<-c(newheight,newheight,newheight)
# 				}else{
# 					x<-c(branches[-mergetable[i,1],"height"],branches[-mergetable[i,2],"height"],atree$height[i])
# 					xend<-rep(atree$height[i],3)
# 				}
# 				y<-c(branches[-mergetable[i,1],"position"],branches[-mergetable[i,2],"position"],branches[-mergetable[i,2],"position"])
# 				yend<-c(branches[-mergetable[i,1],"position"],branches[-mergetable[i,2],"position"],branches[-mergetable[i,1],"position"])
# 				markers<-c(markers,i)
# 				stree<-rbind(stree,data.frame(x=x,y=y,xend=xend,yend=yend))
# 			}
# 		}
# 		branches[-mergetable[markers,1],"position"]<-(branches[-mergetable[markers,1],"position"]+branches[-mergetable[markers,2],"position"])/2
# 		if(type==1){
# 			branches[-mergetable[markers,1],"height"]<-newheight
# 		}else{
# 			branches[-mergetable[markers,1],"height"]<-atree$height[markers]
# 		}
# 		branches[-mergetable[markers,2],2:3]<-NA
# 		#mergetable[mergetable %in% markers]<-mergetable[markers,1]
# 		mergetable[match(markers,mergetable)]<-mergetable[markers,1]
# 		mergetable[markers,]<-NA
# 		if(!all(is.na(mergetable)))stree<-rbind(stree,subtree(mergetable,newheight,branches))
#
# 		return(stree)
# 	}
#
# 	ret<-list()
# 	branches<-data.frame(id=c(1:length(positions)),height=as.double(rep(0,length.out=length(positions))),position=as.double(positions))
# 	ret$segments<-subtree(matrix(as.double(atree$merge),ncol=2),0.0,branches)
# 	ret$labels<-data.frame(x=-0.2,y=positions,labels=atree$labels,hj=1,vj=0.5,an=angle)
# 	#TODO: set labels' x-value based on tree height
# 	#TODO: scale height to stepval? use heights from atree?
# 	if(orientation==1){
# 		ret$segments<-transform(ret$segments,x=y,xend=yend,y=x,yend=xend)
# 		ret$labels<-transform(ret$labels,x=y,y=x,vj=0,hj=0.5)
# 	}else if(orientation==3){
# 		xmax<-max(c(ret$segments$x,ret$segments$xend))
# 		ret$segments<-transform(ret$segments,x=y,xend=yend,y=xmax-x,yend=xmax-xend)
# 		ret$labels<-transform(ret$labels,x=y,y=xmax-x,vj=0)
# 	}else if(orientation==4){
# 		xmax<-max(c(ret$segments$x,ret$segments$xend))
# 		ret$segments<-transform(ret$segments,x=xmax-x,xend=xmax-xend)
# 		ret$labels<-transform(ret$labels,x=xmax-x,hj=0)
# 	}
#
# 	return(ret)
# }

#testtree<-clustree(closest,pos,htype)
#ggplot()+geom_segment(data=testtree$segments,aes(x=x,y=y,xend=xend,yend=yend))+geom_text(data=testtree$labels,aes(x=x,y=y,label=labels,hjust=hj,vjust=vj,angle=an))

# testclust<-closest
# testclust$merge<-matrix(c(-7,-6,-5,1,-2,4,-4,-1,-3,2,3,5),ncol=2)
# testclust$order<-c(7,4,6,1,2,5,3)
# testclust$labels<-LETTERS[1:7]
# testclust$height<-c(1,1,1,2,2,3)
# testclust$method="complete"

#testtree<-clustree(testclust,2,htype)


tpos<-c(1,3:37,40)
ids<-match(c(1:length(closest$labels)),closest$order)
tpos<-tpos[ids]
testtree<-clustree(closest,4,2,tpos)
