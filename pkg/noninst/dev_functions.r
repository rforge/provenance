#additional functions and function overrides for testing
require(provenance)
require(CDFt)

#########################################
## from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
#########################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
										 ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
																			layout.pos.col = matchidx$col))
		}
	}
}

#########################################


# define a function to calculate the dissimilarity matrix
dissimilarity <- function(data,metric="K-S") {
	# TODO: make method for calculation chooseable to allow for Sircombe & Hazelton, 2004 - and others?
	# TODO: calculate based on log-transformed ages?
	metric<-match.arg(metric,c("K-S","C-v-M"))
	n = length(data)
	# instantiate the dissimilarity matrix as an empty data frame
	diss = as.data.frame(mat.or.vec(n,n))
	rownames(diss) = names(data)
	for (i in 1:n){   # loop through all possible pairs of samples
		for (j in 1:n){  # calculate the dissimilarity
			if(metric=="K-S"){
				diss[i,j]<-KolmogorovSmirnov(data[[i]],data[[j]])
			}else if(metric=="C-v-M")
				diss[i,j]<-CramerVonMisesTwoSamples(data[[i]],data[[j]])
		}
	}
	return(diss)
}


dissSH<-function(data,dataerrors){
	#calculate dissimilarity matrix according to Sircombe & Hazelton

	smpls<-names(data)
	#smpls<-unique(smpls)
	cs<-sapply(smpls,FUN=function(x){return(sig2.con(data[[x]],dataerrors[[x]]))})
	#calculated bandwidths:
	Run.c2<-max(cs)

	#calculate (dis?)similarity matrix:
	nsmpl<-length(smpls)
	diss<-matrix(0,nsmpl,nsmpl)
	for(i in 1:nsmpl){
		for(j in 1:nsmpl){
			diss[i,j]<-dXY(data[[smpls[i]]],dataerrors[[smpls[i]]],data[[smpls[j]]],dataerrors[[smpls[j]]],Run.c2)*1000
		}
	}
	row.names(diss)<-smpls
	return(as.dist(diss))
}

clrTransform<-function(x){
	# performs a centered log-ratio transform according to Aitchinson, 1986, on
	# compositional data
	# this is a hack aiming to emulate the results given by clr() of the not
	# furtner maintained package compositions

	clr<-function (x, ...)
	{
		W <- oneOrDataset(x)
		nmv <- is.NMV(W)
		LOG <- unclass(log(ifelse(nmv, W, 1)))
		erg <- ifelse(nmv, LOG - gsi.rsum(LOG)/gsi.rsum(nmv), 0)
		rmult(gsi.simshape(erg, x), orig = x)
	}

	cenLR<-function (x)
	{
		if (dim(x)[2] == 1) {
			res <- list(x.clr = x, gm = rep(1, dim(x)[1]))
		}
		else {
			geometricmean <- function(x) {
				if (any(na.omit(x == 0)))
					0
				else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
			}
			gm <- apply(x, 1, geometricmean)
			x.clr <- log(x/gm)
			res <- list(x.clr = x.clr, gm = gm)
		}
		class(res) <- "clr"
		return(res)
	}


}

#########################################
## new plotKDE()

plotKDE2<-function(data,title,limits=c(0,3000),breaks=NA_real_,bandwidth=NA,
	fcolour=NA,splitat=NA,plotonly=names(data),classes=NA,
	periods=FALSE,hist=FALSE,markers="none",order="auto",
	logx=FALSE,method="botev",...){
	# !!! NOT FINISHED !!! - but useable... (periods and order don't work yet)

	# data      ... input ages, either as a data.frame() with samples in separate
	#                columns for creating stacked KDEs, or for a single KDE plot,
	#                a numeric vector or a data.frame() with one column (or only
	#                one selected, see "plot" below)
	# title     ... optional string giving the title for single KDE plot, or
	#                vector of strings, same length as (selected columns of) data
	# limits    ... optional vector of length 2 giving the age range to plot, or
	#                vector of length 4 giving the age ranges for the two split
	#                plots (see "splitat" below). Set to -1 for automatic min/max.
	# breaks    ... x-axis breaks, optional
	# bandwidth ... optional bandwidth overriding automatic bandwidth calculation.
	#                By default, half bandwidth is used for younger ages in split
	#                plot; vector of length 2 to override split plots separately.
	#                Set to -1 (or FALSE) to plot with individual, "optimal" band-
	#                widths.
	# fcolour   ... optional fill colour, either single colour to be applied to
	#                all KDEs, or vector of colour definitions of same length as
	#                classes.
	# splitat   ... optional age to split plot at. Half-plots will occupy same
	#                amount of space. Limits of length 4 will override splitat.
	# plotonly  ... optional string vector of column names to select from data.
	#                Other columns will be ignored. If given, title and classes
	#                must be of same length as resulting selection.
	# periods   ... optional - plot major geological boundaries for guidance?
	# classes   ... optinal vector of same length as (selected) data, giving a
	#                classification for each data column, determining fill colour.
	#                If unspecified, colours are chosen randomly, unless given in
	#                fcolour. Set classes=-1 to plot all samples in the same
	#                colour (setable by fcolour).
	# hist      ... optionally underlie a histogram? Maximum count will be scaled
	#                to max height of KDE, binwidth is calculated optimal or given
	#                kernel bandwidth.
	# markers   ... optional, either "dash" or "circle" - plot markers at
	#                x-position (age) of data.
	# order     ... optional, have multiple KDEs stacked in "auto"matically in
	#                alphabetical order of names(data) or classes (default), in
	#                the order as entered (order="asis") or set manual (named
	#                numeric vector of same length as unique classes/names).
	# logx      ... optional, should x-axis (ages) be plottet in log-scale?
	# method    ... optional, method to be used for KDE calculation. Possible
	#                values: "botev" for use of Botev et al., 2010, "R" for use of
	#                "standard" density() function of package stats.

	# TODO: how to make sure which class gets which colour? Choose order of classes!?
	#       add handling of length(fcolour)==length(classes) versus length(fcolour)==length(unique(classes))
	#       how to fix assertain order of classes and colours in latter case?
	# TODO: Maybe put classes and colour together in one df/key-value list/...
	# TODO: maybe make ratio of space that left/right half occupy chooseable (is that possible?)?
	# TODO: plot period boundaries?
	# CHECK: should "n=xxx" in title only count ages lying within x-limits?
	# TODO: make"n=xxx" and/or whole title optional, make change of font possible?
	# TODO: make legend optional (despite classification?) - achievable with theme()
	# TODO: implement optional KDE calculation according to Sircombe & Hazelton 2004 - and others?
	# TODO: chooseable, additional algorithms for bw calculation
	# FIXME: if bandwidth=-1, hist=TRUE will result in an error!
	# FIXME: histograms are same for all plots!
	# TODO: take data.frame, list, or geoDataSet as input, format accordingly in here
	#       (check geoDataSet for errors, use for spot indicators? (rectangles/gaussians?))
	# TODO: add custom theme. Either remove horizontal lines only, or try lighter background, or b&w
	# FIXME: better handling of "title"(s)
	# FIXME: use coord_cartesian(xlim=,ylim=) for zooming, instead of limits on scales)
	# CHECK: if custom range is selected, is bandwidth calculated on data with, or all data? Use all? Or not?
	# TODO: better adaption of size of titles, dashes and circles, maybe add transparend rectangles, gaussians
	# TODO: tidy up bandwidth calculation in split KDEs, with potential for histograms in mind...

	# density values smaller than this value will be cur out in final plot
	# prevents very long, drawn-out tails of effectively zero density
  #cutoffy<-1e-5
  cutoffy<-0
  nplot<-2^11
	ncalc<-2^10

	#data frame containing period boundaries and names
	tsdf<-data.frame(stringsAsFactors=FALSE,
		bdry=c(0.0,2.59,23.0,66.0,145,201.3,252.2,298.9,358.9,419.2,443.8,485.4,541.0,2500,4000),
		name=c("Quaternary","Neogene","Paleogene","Cretaceous","Jurassic","Triassic",
			"Permian","Carboniferous","Devonian","Silurian","Ordovician","Cambrian","Proterozoic","Archean","Hadean"),
		short=c("Q","N","Pg","K","J","Tr","P","Ca","D","S","O","Cm","Pz","A","pA"))

	# CHECK: maybe try-catch unknown method and use default?
	method<-match.arg(method,c("botev","R","standard","adaptive"))
	#kdefun<-switch(method,botev=pkde,R=skde,jann=jkde,abramson=akde)

	#check and preformat data
	if(length(data)==0)stop("no data")
	if(!is.list(data)){
		data<-list(age=data)
		plotonly<-c("age")
	}
	plotonly<-names(data) %in% plotonly
	if(length(plotonly)==0)stop("data column selection yielded no data")
	data<-data[plotonly]

	#convert to list and cut out NA, NaN, Inf and -Inf values
	data<-as.list(data)
	for(i in 1:length(data)){
		data[[i]]<-as.numeric(data[[i]])
		rejects<-is.na(data[[i]])|!is.finite(data[[i]])
		data[[i]]<-data[[i]][!rejects]
		if(length(data[[i]])==0){
			warning(sprintf("%s contained no numeric data - removed",names(data)[i]))
			data[[i]]<-NULL
		}
	}
	if(length(data)==0)stop("no data left after data column selection and NA removal")

	#set/check limits, adapt if splitat is given:
	if(length(limits)==2){
		if(!is.na(splitat)){
			if((splitat>limits[1])&&(splitat<limits[2])){
				limits<-c(limits[1],splitat,splitat,limits[2])
			}else{
				warning("splitat outside age limits - ignored")
			}
		}
	}else if(length(limits)==1 && limits==-1){
		limits<-c(min(unlist(data)),max(unlist(data)))
	}else if(length(limits)!=4){
		warning("invalid limits parameter - using default")
		limits<-c(0,3000)
	}

	# Check on breaks...
	if(any(!is.numeric(breaks)))warning("non-numeric break values - using default")
	breaks<-as.numeric(breaks)
	if(logx){
		# breaks and data transformation for logarithmic scale
		for(i in c(1:length(data))){
			data[[i]]<-log10(data[[i]])
		}
		# TODO: should not use data, but minimal positive value of KDE-x-values after calculation
		if(any(limits<=0))limits[limits<=0]<-min(unlist(data))
		if(any(limits<=0))limits[limits<=0]<- 1 # use minimum provided age, unless erroneous data - use 1 Ma
		if(anyNA(breaks)){
			if(length(limits)==2){
				breaks=floor(trans_breaks("log10", function(x) 10^x,n=8)(limits))
			}else{
				breaks=floor(c(trans_breaks("log10", function(x) 10^x,n=5)(limits[c(1,2)]),
					trans_breaks("log10", function(x) 10^x,n=5)(limits[c(3,4)])))
			}
		}else{
			if(any(breaks<=0))breaks<-breaks[breaks>0]
			if(min(breaks)<limits[1])limits[1]<-min(breaks)
			if(max(breaks)>limits[length(limits)])limits[length(limits)]<-max(breaks)
		}
		limits<-log10(limits)
		labels<-breaks
		#cheat a bit: can't actually plot at 0 in log scale!
		breaks[breaks<=0]<-1e-10
	}else{
		# set/adjust breaks for lin scale
		if(anyNA(breaks)){
			if(length(limits)==2){
				breaks<-pretty_breaks(n=8)(limits)
			}else{
				breaks<-c(pretty_breaks(n=5)(limits[c(1,2)]),pretty_breaks(n=5)(limits[c(3,4)]))
			}
		}else{
			if(min(breaks)<limits[1])limits[1]<-min(breaks)
			if(max(breaks)>limits[length(limits)])limits[length(limits)]<-max(breaks)
		}
		labels<-breaks
	}


	##### new: prepare data.frames plotpars and classpars, that keep all the information together...
	#at this point, data is already cropped to data[plotonly]...
	plotpars<-data.frame(series=names(data))
	plotpars$ndata<-sapply(data,length)

	#titles and plotting coordinates
	plotpars$labx<-limits[length(limits)] 	#coodinates for plotting labels...
	plotpars$laby<-0.98	#*max(data2$d,na.rm=TRUE)
	plotpars$section<-ifelse(length(limits)!=4,1,2)
	#check/generate title(s)
	if((names(data)=="age")||is.null(names(data))){
		if(missing(title)||is.na(title)){
			plotpars$title<-""
		}else{
			plotpars$title<-title
		}
	}else	if(length(data)>=1){
		if(missing(title)||is.na(title)||(length(title)!=length(data))){
			plotpars$title<-plotpars$series
		}else{
			plotpars$title<-""
		}
	}else{
		stop("invalid data")
	}

	#bandwidth(s):
	bw1<-NA
	#choose best bandwidth as 25%-quantile of all calculated bandwidths, or as input:
	if(length(bandwidth)==1 && is.na(bandwidth)){
		plotpars$bw<-sapply(data,optimal_bw,n=ncalc)
		bw1<-quantile(plotpars$bw,probs=c(0.25))
	}else if(length(bandwidth)==1 && (bandwidth==FALSE||bandwidth==-1)){
		plotpars$bw<-sapply(data,optimal_bw,n=ncalc)
		bw1<-FALSE
	}else if(length(bandwidth)==1 && is.numeric(bandwidth)){
		plotpars$bw<-rep(bandwidth,length(data))
		bw1<-bandwidth
	}else if(length(bandwidth)==2 && is.numeric(bandwidth)){
		plotpars$bw<-rep(bandwidth,length(data))
		bw1<-bandwidth[2]
	}else{
		stop("type or length of bandwith unsuitable")
	}

	#adding classes:
	if(!is.na(classes) && (length(classes)==length(data))){
		plotpars$class<-"n/a"
		for(i in 1:length(data)){plotpars$class[plotpars$series==names(data)[i]]<-classes[i]}
	}else if(!is.na(classes) && (classes==-1)){
		plotpars$class<-"sample"
		classes<-"sample"
	}else{
		if(!is.na(classes))warning("invalid classes - ignoring")
		plotpars$class<-plotpars$series
		classes<-plotpars$series
	}
	# TODO: plotpars$position here

	print(plotpars)

	classpars<-data.frame(class=unique(classes))
	classpars$order<-c(1:length(classpars$class))
	if(order=="auto"){
		classpars$order<-order(classpars$class)
		
	}else if(order=="asis"){
		#nothing to do?
	}else if(is.numeric(order) && length(order)==length(classpars$class) && !anyDuplicated(order)){
		if(all(classpars$class %in% names(order))){
			for(i in 1:dim(classpars)[1])classpars$order[i]<-order[names(order)==classpars$class[i]]
		}
	}else{
		warning("invalid 'order' parameter - using default")
	}


	print(classpars)

	#calculate KDEs, put in data. frame
	data2<-data.frame(x=NULL,d=NULL,smpl=NULL,section=NULL,stringsAsFactors=FALSE)	#collected kdes go in here
	lbs<-data.frame(smpl=NULL,n=NULL,stringsAsFactors=FALSE)	#labels for plotting
	for(i in 1:length(data)){	#loop over elements of data
		spl<-names(data)[i]	#current column name = sample name
		lbs<-rbind(lbs,data.frame(smpl=spl,n=length(data[[i]])))
		if(bw1==FALSE){
			ckde<-gkde(data=data[[i]],n=nplot,method=method,...)
		}else{
			# TODO: does not work:
			#if(logx){bw1<-log10(bw1)}
			ckde<-gkde(data=data[[i]],n=nplot,bandwidth=bw1,method=method,...)
		}
		dd<-data.frame(t(ckde))	#reformat the output...

		#calculate second kde for ages within limits[3:4], with optional different bandwidth
		#cut output to range(s) within limits
		if(length(limits)==4){
			dd<-dd[((dd[[1]]>=limits[3])&(dd[[1]]<=limits[4])),]
			dd$section<-2
			if(length(bandwidth)==2){
				bw2<-bandwidth[1]
			}else{
				if(bw1==FALSE){
					bw2<-FALSE
				}else{
					bw2<-bw1/2
				}
			}
			if(bw2==FALSE){
				dd2<-data.frame(t(gkde(data=data[[i]],n=nplot,method=method,...)))
			}else{
				#if(logx){bw2<-log10(bw2)}
				dd2<-data.frame(t(gkde(data=data[[i]],n=nplot,bandwidth=bw2,method=method,...)))
			}
			dd2$section<-1
			dd<-rbind(dd,dd2[((dd2[[1]]>=limits[1])&(dd2[[1]]<=limits[2])),])
		}else{
			dd<-dd[(dd[[1]]>=limits[1])&(dd[[1]]<=limits[2]),]
			dd$section<-1
		}
		names(dd)<-c("x","d","section")
		dd$smpl<-spl
		dd$d<-dd$d/max(dd$d,na.rm=TRUE)	#normalise density values to max==1
		data2<-rbind(data2,dd)	#collate new kde into data2
	}
	#cut out values effectively 0:
	data2<-data2[data2$d>=cutoffy,]
	data2<-data2[!is.na(data2$d),]

# 	#adding new column "type" for classification, based on classes parameter:
# 	if(!is.na(classes) && (length(classes)==length(data))){
# 		data2$type<-"n/a"
# 		for(i in 1:length(data)){data2$type[data2$smpl==names(data)[i]]<-classes[i]}
# 	}else if(!is.na(classes) && (classes==-1)){
# 		data2$type<-"sample"
# 	}else{
# 		if(!is.na(classes))warning("invalid classes - ignoring")
# 		data2$type<-data2$smpl
# 	}

	#adding new column "type" for classification, based on classes
	data2$type<-"n/a"
	for(i in 1:dim(plotpars)[1]){data2$type[data2$smpl==plotpars$series[i]]<-plotpars$class[i]}

	#sort data so that the same "type"s plot together:
	data2<-within(data2,i<-order(data2$type,data2$smpl,data2$x))
	data2<-data2[data2$i,]
	data2$smpl<-factor(data2$smpl,levels=unique(data2$smpl))

	#check/set fcolour:
	ntypes<-length(unique(data2$type))
	if(length(fcolour)==0||all(is.na(fcolour))){
		if(ntypes==1){
			fcolour<-"#4444FF66"
		}else{
			fcolour<-hsv(seq(0.08,1,length.out=ntypes),s=0.5,v=0.95,alpha=0.8)
		}
	}else if(length(fcolour)==1){
		if(ntypes==1){
			#nothing to do?
		}else{
			fcolour<-rep(fcolour,ntypes)
		}
	}else if(length(fcolour)==length(classes)){
		#match the right colours to the already sorted classes - there should be a better way to do this!
		#ucols<-unique(fcolour)
		#names(ucols)<-unique(classes)
		#fcolour<-ucols[match(classes,names(ucols))]
	}else{
		warning("fcolour invalid - using default")
		fcolour<-hsv(seq(0.08,1,length.out=ntypes),s=0.5,v=0.95,alpha=0.8)
	}

# 	#populate data frame for labels
# 	lbs$x<-limits[length(limits)] 	#coodinates for plotting labels...
# 	lbs$y<-0.98*max(data2$d,na.rm=TRUE)
# 	lbs$section<-ifelse(length(limits)!=4,1,2)
# 	#check/generate title(s)
# 	if((names(data)=="age")||is.null(names(data))){
# 		if(missing(title)||is.na(title)){
# 			lbs$title<-paste0("n=",sapply(data,length))
# 		}else{
# 			lbs$title<-paste0(title,", n=",sapply(data,length))
# 		}
# 	}else	if(length(data)>=1){
# 		if(missing(title)||is.na(title)||(length(title)!=length(data))){
# 			lbs$title<-paste0(names(data),", n=",sapply(data,length))
# 		}else{
# 			lbs$title<-paste0(title,", n=",sapply(data,length))
# 		}
# 	}else{
# 		stop("invalid data")
# 	}

	#molten input data for markers, also used for histograms:
	dm<-melt(data,stringsAsFactors=FALSE)
	names(dm)[grep("value",names(dm),invert=TRUE)]<-"smpl"
	names(dm)[grep("value",names(dm))]<-"age"
	dm$section<-1
	dm$bw<-bw1
	if(length(limits)==4){
		dm$section[((dm$age>=limits[3])&(dm$age<=limits[4]))]<-2
		bw2<-bw1/2
		if(length(bandwidth)==2)bw2<-bandwidth[1]
		dm$bw[dm$section==1]<-bw2
		dm<-dm[(dm$age>=limits[1]&dm$age<=limits[2])|(dm$age>=limits[3]&dm$age<=limits[4]),]
	}else{
		dm<-dm[dm$age>=limits[1]&dm$age<=limits[2],]
	}

	#create plot:
	lw<-rel(0.3)
	g<-ggplot()
	#density:
	g<-g+geom_density(data=data2,aes(x=x,y=d,fill=type),stat="identity",size=lw)
	#histogram:
	if(hist){
		# 		h<-hist(x=dm$age,plot=FALSE,breaks=(limits[2]-limits[1])/bw1)
		# 		hdf<-data.frame(x=h$mids,y=h$count/max(h$count,na.rm=TRUE))
		# 		g<-g+geom_histogram(data=hdf,aes(x=x,y=y),stat="identity",fill=NA,colour="black",binwidth=bw1)
		g<-g+stat_bin(data=dm,aes(x=age,y=..ncount..),binwidth=bw,fill=NA,colour="grey40",size=lw)
	}
	#	g<-g+geom_vline(xintercept=splitat)
	#data markers:
	if(markers=="dash"){
		g<-g+geom_segment(data=dm,aes(x=age,xend=age,y=-0.02,yend=-0.06))
	}else if(markers=="circle"){
		g<-g+geom_point(data=dm,aes(x=age,y=-0.05),colour="#00000022",size=rel(2.5))
	}
	#title:
	g<-g+geom_text(data=plotpars,aes(x=labx,y=laby,label=sprintf("%s, n=%d",title,ndata),group=class),hjust=1,vjust=1,size=rel(4.2))
	#x-axis label, tick marks, age range...:
	if(logx){
		# 		scale_x_continuous(name="Ma",limits=limits,breaks=breaks,labels=breaks)+
		g<-g+scale_x_continuous(name="Ma",breaks=log10(breaks),labels=labels)
	}else{
		g<-g+scale_x_continuous(name="Ma",breaks=breaks,labels=labels)
	}
	g<-g+
		#colour scale:
		scale_fill_manual(values=fcolour)+
		#we need no legend title:
		guides(fill=guide_legend(title=NULL))+
		#stack plots by sample name:
		facet_grid(smpl ~ section,scales="free_x")+
		#blank out y-axis:
		theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
		#remove facet strips, horizontal grid lines, make background very light grey
		theme(strip.text=element_blank(),strip.background=element_blank(),
			panel.grid.major.y=element_blank(),panel.grid.minor.y=element_blank(),
			panel.background=element_rect(fill="#F4F4F4"))
	#if only one colour - need no legend for colours:
	if((length(unique(data2$type))==1)||(length(unique(fcolour))==1))g<-g+theme(legend.position="none")

	return(g)
}


#########################################

# define a function to plot the MDS coordinates and connect the nearest neighbours, ggplot2-style
plotMDS <- function(mds,diss,col="",sym="",nearest=TRUE,labels=TRUE,symbols=TRUE,fcolour=NA,
										stretch=FALSE,axes=FALSE,expand=FALSE) {
	#quick and dirty, makes big assumptions about what mds looks like if it's not a two-column matrix.
	#mds can be: data.frame(x,y,any,further,columns,...)
	#
	# mds     ... data.frame() or matrix with at least 2 columns representing x and y
	#              coordinates calculated by cmdscale() or isoMDS()
	# diss    ... (dis-)similarities of samples, only required if nearest==TRUE
	# col     ... string giving the name of column in mds to be used for colour scale
	# sym     ... string giving the name of column in mds to be used for symbol scale
	# nearest ... boolean - plot lines connecting nearest neighbours?
	# labels  ... boolean - plot data labels (taken from row.names of mds)?
	# symbols ... boolean - plot data points (useful if plotting only labels)?
	# fcolour ... symbol fill colours - not very functional yet...
	# stretch ... boolean - strech plot to full extend, or use fixed x/y scale ratio
	# axes    ... boolean, plot axis labels?
	# expand  ... boolean, expand plot area to be square (only valid if stretch==FALSE)
	
	# TODO: let labels optionally be taken from a column in mds
	# TODO: if length(col) too great for brewer, switch to alternative colour scale
	# TODO: more symbols for symbol scale
	# TODO: check diss if nearest==TRUE
	
	require(ggplot2)
	require(scales)
	defcol<-"lightblue"
	ssize<-9
	tsize<-3
	ddf<-as.data.frame(mds)
	if(!all(c("x","y") %in% names(ddf)))names(ddf)[1:2]<-c("x","y")
	# create a new (empty) plot
	#	plot(mds[,1],mds[,2],type='n')
	p<-ggplot()
	# draw lines between closest neighbours
	#	for (j in 1:nrow(mds)) {
	#		lines(c(mds[j,1],x1[j]),c(mds[j,2],y1[j]),lty=1)
	#		lines(c(mds[j,1],x2[j]),c(mds[j,2],y2[j]),lty=2)
	#	}
	if(nearest){
		# indices of nearest and second nearest neighbours:
		i = t(apply(as.matrix(diss),1,function(x) order(x))[2:3,])
		# plot coordinates for the lines
		x1 = as.vector(ddf[i[,1],"x"])
		y1 = as.vector(ddf[i[,1],"y"])
		x2 = as.vector(ddf[i[,2],"x"])
		y2 = as.vector(ddf[i[,2],"y"])
		p<-p+geom_segment(data=ddf,x=ddf[["x"]],y=ddf[["y"]],xend=x1,yend=y1,linetype="solid")
		p<-p+geom_segment(data=ddf,x=ddf[["x"]],y=ddf[["y"]],xend=x2,yend=y2,linetype="dashed")
	}
	# plot the configuration as labeled circles
	#	points(mds[,1],mds[,2],pch=21,cex=2.5,col='red',bg='white')
	# TODO: a lot of checks on col and sym, sensible automatic assumptions
	if(symbols){
		if(length(ddf)>2){
			if(col==""){
				if(sym==""){
					p<-p+geom_point(data=ddf,aes_string(x="x",y="y"),fill=defcol,colour="black",size=ssize,shape=21)
				}else{
					p<-p+geom_point(data=ddf,aes_string(x="x",y="y",shape=sym),fill=defcol,colour="black",size=ssize)
				}
			}else{
				if(sym==""){
					p<-p+geom_point(data=ddf,aes_string(x="x",y="y",fill=col),colour="black",size=ssize,shape=21)
				}else{
					p<-p+geom_point(data=ddf,aes_string(x="x",y="y",fill=col,shape=sym),colour="black",size=ssize)
				}
			}
		}else{
			c2=defcol
			if(col!=""){
				if(length(col)>1){
					warning("more than one colour value specified - using default")
					c2<-defcol
				}else{
					# TODO: check if col IS a colour value
					c2<-col
				}
			}
			p<-p+geom_point(data=ddf,aes(x=x,y=y),fill=c2,colour="black",size=ssize,shape=21)
		}
	}
	#	text(mds[,1],mds[,2],row.names(mds))
	#add labels:
	if(labels)p<-p+geom_text(data=ddf,aes(x=x,y=y),label=row.names(ddf),size=tsize)
	
	#apply colour scale, symbol scale, format plot:
	
	cols<-c(rgb(t(col2rgb(defcol))/255))  #default colour
	#print(cols)
	#browser()
	if(col!=""){
		if(length(fcolour)==0||is.na(fcolour)){
			cols<-brewer_pal(type="div",palette=6)(length(unique(ddf[[col]])))
		}else if(length(fcolour)==1){
			#nothing to do
		}else if(length(fcolour)==length(unique(ddf[[col]]))){
			ucols<-fcolour
			names(ucols)<-unique(ddf[[col]])
			cols<-ucols[match(sort(unique(ddf[[col]])),names(ucols))]
			# FIXME: allow reordering here...
		}else{
			warning("invalid fill colour parameter - using default")
			#...to be precise: "leaving" default....
		}
	}
	#print(cols)
	if(length(cols)>1)p<-p+scale_fill_manual(values=cols)
	if(sym!=""){
		bks<-unique(ddf[[sym]])	# TODO: warn if too many, define new ones?
		p<-p+scale_shape_manual(values=rep(c(21,22,24,23),length.out=10)[1:length(bks)])
	}
	p<-p+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
		theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank())+
		theme(panel.background=element_rect(fill=NA,colour="black"),plot.background=element_rect(fill=NA,colour=NA))
	if(!stretch){
		p<-p+coord_equal()
		if(expand){
			xrange<-range(ddf$x)
			yrange<-range(ddf$y)
			maxrange<-c(min(xrange[1],yrange[1]),max(xrange[2],yrange[2]))
			p<-p+xlim(maxrange)+ylim(maxrange)
		}
	}
	if(col!=""){
		p<-p+guides(fill=guide_legend(title=NULL,order=2,override.aes=list(fill=cols,shape=21)),shape=guide_legend(title=NULL,order=1))
	}
	if(!axes){
		p<-p+theme(axis.ticks=element_blank(), axis.text=element_blank())
	}
	return(p)
}


#########################################
#following Jann, 2007
#this works and allows manual definition of kernel functions, but is slow...
jkde<-function(data,bandwidth,n,from,to,kernel="gaussian",...){
	data<-data[!is.na(data)]
	if(missing(bandwidth))bandwidth<-optimal_bw(data,n=2^10)
	if(missing(n))n<-2^11
	kernel<-match.arg(kernel,c("gaussian","epanechnikov","triangular","epan2"))
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)
	ret<-fk(data,ev,bandwidth,kernel)
	return(matrix(c(ev,ret),nrow=2,byrow=TRUE))
}

#adaptive kde - after Jann, 2007, following Abramson, 1982
# FIXME: automatic bandwidth calculation
akde<-function(data,bandwidth,n,from,to,kernel="gaussian",...){
	if(missing(bandwidth))bandwidth<-optimal_bw(data,n=2^10)
	if(missing(n))n<-2^11
	kernel<-match.arg(kernel,c("gaussian","epanechnikov","triangular","epan2"))
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)
	#tf<-approxfun(x=ev,y=fak(data,ev,bandwidth,kernel)/length(data))
	#ret<-tf(ev)
	ret<-fak(data,ev,bandwidth,kernel)
	return(matrix(c(ev,ret),nrow=2,byrow=TRUE))
}

# helper functions from JANN paper:
# (1)
fk<-function(data,vals,h,kernel){
	#n<-length(data)
	#zs<-outer(vals,data,FUN="-")/h
	# not the prettiest, but the fastest: do NOT pre-calculate anything, let R do
	# the optimisation on the fly
	return(sapply(vals,function(x){return(
			sum(sapply(data,function(Xi){return(
				Kfun((x-Xi)/h,kernel=kernel)/h
			)}))/length(data)
		)})
	)
}

fak<-function(data,vals,h,kernel){
	initial<-fk(data,data,h,kernel)
	G<-exp(mean(log(c(initial))))
	tf<-approx(x=data,y=initial,xout=data)
	h2<-approxfun(x=data,y=h*sapply(seq(data),function(Xi){return(h*sqrt(G/tf$y[Xi]))}))
#	h2<-approxfun(x=data,y=h*(h*sqrt(G/initial)))
	return(sapply(vals,function(x){
		return(sum(sapply(data,function(Xi){
			h2xi<-h2(Xi)
			return(Kfun((x-Xi)/h2xi,kernel=kernel)/h2xi)
		})))
	}))
}

lambda<-function(Xi,h,G,data,kernel){
	return(sqrt(G/fk(data,Xi,h,kernel)))
}

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

#general KDE function, calls specific bandwith calculations according to "method"
#following Jann, 2007, Abramson, 1982
gkde<-function(data,bandwidth=optimal_bw(data,n=2^10),n=2^11,from,to,kernel="gaussian",method="botev",...){
	data<-data[!is.na(data)]
	if(length(data)==0)stop("no valid data")
	kernel<-match.arg(kernel,c("gaussian","epanechnikov","triangular","epan2"))
	method<-match.arg(method,c("botev","R","standard","adaptive"))
	minimum=min(data)
	maximum=max(data)
	Range=maximum-minimum
	if(missing(from))from<-minimum-Range/10
	if(missing(to))to<-maximum+Range/10
	ev<-seq(from=from,to=to,length.out=n)

	ret<-NULL
	if(method=="botev"){
		ret<-pkde(data=data,bandwidth=bandwidth,n=n)[2,]
	}else if(method=="R"){
		ret<-density(x=data,bw=bandwidth,kernel=kernel,from=from,to=to,n=n)$y
	}else if(method=="standard"){
		ret<-fk(data,ev,bandwidth,kernel)
	}else if(method=="adaptive"){
		ret<-fak(data,ev,bandwidth,kernel)
	}
	return(matrix(c(ev,ret),nrow=2,byrow=TRUE))
}

#########################################

plotDendrogram<-function(data,orientation=4,type=2,positions=c(1:length(data)),metric=c("K-S","C-v-M"),method="ward.D",classes=1){
  #plots a dendrogram for the given data, according to specified metric
  #wrapper for clustree(), see there for further parameters
  #classes ... cuts tree into classes branches, colours leaves accordingly
  
  #TODO: improve colouring of branches to include all lines of one branch
  #      maybe implement in clustree rather, and return as extra column in $segments
  
  diss<-dissimilarity(data,metric=metric)
  closest<-hclust(as.dist(diss),method=method)
  
  ids<-match(c(1:length(closest$labels)),closest$order)
  positions<-positions[ids]
  tree<-clustree(closest,orientation,type,positions,classes)
  
  g<-ggplot()+
    geom_segment(data=tree$segments,aes(x=x,y=y,xend=xend,yend=yend,colour=factor(branch)))+
    geom_text(data=tree$labels,aes(x=x,y=y,label=labels,hjust=hj,vjust=vj,angle=an,colour=factor(branch)),size=rel(4))+
    scale_colour_manual(values=c("#000000",brewer_pal(type="qual",palette="Paired")(length(unique(tree$labels$branch)))),guide="none")
  return(g)
}

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
  ret$labels<-data.frame(x=0,y=positions,labels=atree$labels,hj=1.2,vj=0.5,an=angle,branch=bclass)
  #TODO: set labels' x-value based on tree height
  #TODO: scale height to stepval? use heights from atree?
  if(orientation==1){
    ret$segments<-transform(ret$segments,x=y,xend=yend,y=x,yend=xend)
    ret$labels<-transform(ret$labels,x=y,y=x,hj=-0.2)
  }else if(orientation==3){
    xmax<-max(c(ret$segments$x,ret$segments$xend))
    ret$segments<-transform(ret$segments,x=y,xend=yend,y=xmax-x,yend=xmax-xend)
    ret$labels<-transform(ret$labels,x=y,y=xmax-x,hj=-0.2)
  }else if(orientation==4){
    xmax<-max(c(ret$segments$x,ret$segments$xend))
    ret$segments<-transform(ret$segments,x=xmax-x,xend=xmax-xend)
    ret$labels<-transform(ret$labels,x=xmax-x,hj=-0.2)
  }
  
  return(ret)
}
