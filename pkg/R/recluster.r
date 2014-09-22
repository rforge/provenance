
clustree<-function(atree,orientation=c(1:4),type=c(1,2),positions=NULL,classes=1){
	#takes a hclust structure and returns rendered coordinates to draw the tree,
	#according to orientation and type parameters. Also returns leaf labels.
	#generates segments of tree based on pre-calculated coordinates for leaves and in any orientation,
	#unlike plot.hclust() or package ggdendro.
	#positions takes y-positions in the order of drawing layout (atree$order), i.e. ascending values
	#orientations ... braching direction: 1 top-down, 2 right-to-left, 3 bottom-up, 4 left-to-right
	#type ... 1 uniform height steps, 2 calculated by hclust w/base 0, 3 same as 2, w/base = height from hclust
	#classes ... cut tree into classes sub-branches, return as column in $segments and $leaves

	#TODO: type==3

	if(class(atree)!="hclust")stop("atree must be of class hclust")
	stepval<-switch(orientation,1,1,-1,-1)
	angle<-switch(orientation,270,0,90,0)
	if(length(positions)==0){
		positions<-match(c(1:length(atree$labels)),atree$order)
	}else if(length(positions)==length(atree$labels)){
		#nothing to do!?
	}else{
		stop("invalid positions")
	}
	# TODO: if positions not given, optionally compute from height/distances in hclust
	# TODO: if a class/branch has only one member, plot in this colour

	#internal recursive function to parse the tree
	subtree<-function(mergetable,curpos,branches){
		newheight<-curpos+1
		stree<-data.frame()
		markers<-NULL
		for(i in c(1:dim(mergetable)[1])){
			if(!all(is.na(mergetable[i,])) && mergetable[i,1]<0 && mergetable[i,2]<0){
				if(type==1){
					x<-c(branches[-mergetable[i,1],"height"],branches[-mergetable[i,2],"height"],newheight)
					xend<-c(newheight,newheight,newheight)
				}else{
					x<-c(branches[-mergetable[i,1],"height"],branches[-mergetable[i,2],"height"],atree$height[i])
					xend<-rep(atree$height[i],3)
				}
				y<-c(branches[-mergetable[i,1],"position"],branches[-mergetable[i,2],"position"],branches[-mergetable[i,2],"position"])
				yend<-c(branches[-mergetable[i,1],"position"],branches[-mergetable[i,2],"position"],branches[-mergetable[i,1],"position"])
				markers<-c(markers,i)
				branch<-ifelse(branches[-mergetable[i,1],"branch"]==branches[-mergetable[i,2],"branch"],branches[-mergetable[i,1],"branch"],0)
				stree<-rbind(stree,data.frame(x=x,y=y,xend=xend,yend=yend,branch=branch))
			}
		}
		branches[-mergetable[markers,1],"position"]<-(branches[-mergetable[markers,1],"position"]+branches[-mergetable[markers,2],"position"])/2
		if(type==1){
			branches[-mergetable[markers,1],"height"]<-newheight
		}else{
			branches[-mergetable[markers,1],"height"]<-atree$height[markers]
		}
		branches[-mergetable[markers,2],2:3]<-NA
		mergetable[match(markers,mergetable)]<-mergetable[markers,1]
		mergetable[markers,]<-NA
		if(!all(is.na(mergetable)))stree<-rbind(stree,subtree(mergetable,newheight,branches))

		return(stree)
	}

	ret<-list()
	bclass<-cutree(atree,k=classes)
	branches<-data.frame(id=c(1:length(positions)),height=as.double(rep(0,length.out=length(positions))),position=as.double(positions),branch=bclass)
	ret$segments<-subtree(matrix(as.double(atree$merge),ncol=2),0.0,branches)
	ret$labels<-data.frame(x=-0.2,y=positions,labels=atree$labels,hj=1,vj=0.5,an=angle,branch=bclass)
	#TODO: set labels' x-value based on tree height
	#TODO: scale height to stepval? use heights from atree?
	if(orientation==1){
		ret$segments<-transform(ret$segments,x=y,xend=yend,y=x,yend=xend)
		ret$labels<-transform(ret$labels,x=y,y=x,vj=0,hj=0.5)
	}else if(orientation==3){
		xmax<-max(c(ret$segments$x,ret$segments$xend))
		ret$segments<-transform(ret$segments,x=y,xend=yend,y=xmax-x,yend=xmax-xend)
		ret$labels<-transform(ret$labels,x=y,y=xmax-x,vj=0)
	}else if(orientation==4){
		xmax<-max(c(ret$segments$x,ret$segments$xend))
		ret$segments<-transform(ret$segments,x=xmax-x,xend=xmax-xend)
		ret$labels<-transform(ret$labels,x=xmax-x,hj=0)
	}

	return(ret)
}
