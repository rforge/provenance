# main KDE plotting function
# (c) Martin Rittner 2013-2014
# based on:
# Pieter Vermeesch, 2012
# "On the visualisation of detrital age distributions"
# Chemical Geology 312–313 (2012) p.190–194

plotKDE<-function(data,title,limits=c(0,3000),breaks=NA_real_,bandwidth=NA,
									fcolour=NA,splitat=NA,plotonly=names(data),classes=NA,
									periods=FALSE,hist=FALSE,markers="none",order=TRUE,
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
	# order     ... optional, have multiple KDEs stacked in alphabetical order of
	#                names(data) (default), or in order as entered (order=FALSE)?
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
	cutoffy<-1e-5

	#data frame containing period boundaries and names
	tsdf<-data.frame(stringsAsFactors=FALSE,
			 bdry=c(0.0,2.59,23.0,66.0,145,201.3,252.2,298.9,358.9,419.2,443.8,485.4,541.0,2500,4000),
			 name=c("Quaternary","Neogene","Paleogene","Cretaceous","Jurassic","Triassic",
				 			 "Permian","Carboniferous","Devonian","Silurian","Ordovician","Cambrian","Proterozoic","Archean","Hadean"),
			 short=c("Q","N","Pg","K","J","Tr","P","Ca","D","S","O","Cm","Pz","A","pA"))

	# CHECK: maybe try-catch unknown method and use default?
	method<-match.arg(method,c("botev","R","jann","abramson"))
	kdefun<-switch(method,botev=pkde,R=skde,jann=jkde,abramson=akde)

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
		if(any(limits<=0))limits[limits<=0]<-min(unlist(data))
		if(any(limits<=0))limits[limits<=0]<- 0.1 # use minimum provided age, unless erroneous data - use 1 Ma
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
	}
	labels<-breaks

	#choose best bandwidth as 25%-quantile of all calculated bandwidths, or as input:
	bw1<-NA
	bws<-NA
	if(length(bandwidth)==1 && is.na(bandwidth)){
		bws<-sapply(data,optimal_bw,n=2^10)
		bw1<-quantile(bws,probs=c(0.25))
	}else if(length(bandwidth)==1 && (bandwidth==FALSE||bandwidth==-1)){
		bws<-sapply(data,optimal_bw,n=2^10)
		bw1<-FALSE
	}else if(length(bandwidth)==1 && is.numeric(bandwidth)){
		bws<-rep(bandwidth,length(data))
		bw1<-bandwidth
	}else if(length(bandwidth)==2 && is.numeric(bandwidth)){
		bw1<-bandwidth[2]
	}else{
		stop("type or length of bandwith unsuitable")
	}

	#calculate KDEs, put in data. frame
	data2<-data.frame(x=NULL,d=NULL,smpl=NULL,section=NULL,stringsAsFactors=FALSE)	#collected kdes go in here
	lbs<-data.frame(smpl=NULL,n=NULL,stringsAsFactors=FALSE)	#labels for plotting
	for(i in 1:length(data)){	#loop over elements of data
		spl<-names(data)[i]	#current column name = sample name
		lbs<-rbind(lbs,data.frame(smpl=spl,n=length(data[[i]])))
		if(bw1==FALSE){
			ckde<-kdefun(data=data[[i]],n=2^12,...)
		}else{
			# TODO: does not work:
			#if(logx){bw1<-log10(bw1)}
			ckde<-kdefun(data=data[[i]],n=2^12,bandwidth=bw1,...)
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
				dd2<-data.frame(t(kdefun(data=data[[i]],n=2^12,...)))
			}else{
				#if(logx){bw2<-log10(bw2)}
				dd2<-data.frame(t(kdefun(data=data[[i]],n=2^12,bandwidth=bw2,...)))
			}
			dd2$section<-1
			dd<-rbind(dd,dd2[((dd2[[1]]>=limits[1])&(dd2[[1]]<=limits[2])),])
		}else{
			dd<-dd[(dd[[1]]>=limits[1])&(dd[[1]]<=limits[2]),]
			dd$section<-1
		}
		names(dd)<-c("x","d","section")
		dd$smpl<-spl
		dd$d<-dd$d/max(dd$d,na.rm=TRUE)	#normalise density values
		data2<-rbind(data2,dd)	#collate new kde into data2
	}

	#cut out values effectively 0:
	data2<-data2[data2$d>cutoffy,]
	data2<-data2[!is.na(data2$d),]

	#adding new column "type" for classification, based on classes parameter:
	if(!is.na(classes) && (length(classes)==length(data))){
		data2$type<-"n/a"
		for(i in 1:length(data)){data2$type[data2$smpl==names(data)[i]]<-classes[i]}
	}else if(!is.na(classes) && (classes==-1)){
		data2$type<-"sample"
	}else{
		if(!is.na(classes))warning("invalid classes - ignoring")
		data2$type<-data2$smpl
	}

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
	  ucols<-unique(fcolour)
	  names(ucols)<-unique(classes)
	  fcolour<-ucols[match(classes,names(ucols))]
	}else{
		warning("fcolour invalid - using default")
		fcolour<-hsv(seq(0.08,1,length.out=ntypes),s=0.5,v=0.95,alpha=0.8)
	}

	#populate data frame for labels
	lbs$x<-limits[length(limits)] 	#coodinates for plotting labels...
	lbs$y<-0.98*max(data2$d,na.rm=TRUE)
	lbs$section<-ifelse(length(limits)!=4,1,2)
	#check/generate title(s)
	if((names(data)=="age")||is.null(names(data))){
		if(missing(title)||is.na(title)){
			lbs$title<-paste0("n=",sapply(data,length))
		}else{
			lbs$title<-paste0(title,", n=",sapply(data,length))
		}
	}else	if(length(data)>=1){
		if(missing(title)||is.na(title)||(length(title)!=length(data))){
			lbs$title<-paste0(names(data),", n=",sapply(data,length))
		}else{
			lbs$title<-paste0(title,", n=",sapply(data,length))
		}
	}else{
		stop("invalid data")
	}

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
	g<-g+geom_text(data=lbs,aes(x=x,y=y,label=title),hjust=1,vjust=1,size=rel(4.2))
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
