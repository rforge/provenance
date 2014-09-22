##### To load Kaja's zrn data, follow this:
#you need to install the "provenance" package, Menu Tools/Install Packages..., choose "Package Archive File".
#after opening this file, go through it step-by-step, using Ctrl-Enter

##( - Giulia, ignore this and see further down...)

## path to the data - change according to where you saved it
basepath<-"/home/martin/Documents/LA-Lab-related/FOR_/for_Kaja/zrn_data"

## find all ..._preferred.csv files in this folder, load them into a data object called 'data':
# list all .._preferred.csv files
fls<-list.files(path=basepath,pattern=".*_preferred.csv",full.names=TRUE)
# prepare data object
data<-list()
# load the individual data files
for(f in fls){
	# extract sample names from filename
	spl<-sub("_preferred.csv","",basename(f))
	# load data
	cdat<-read.table(file=f,header=TRUE,sep=",",stringsAsFactors=FALSE,skip=6)
	# skip line without a 'preferred age'
	cdat<-cdat[!is.na(cdat$preferred.age),]
	# copy the remaining ages into the 'data' object
	data[[spl]]<-cdat$preferred.age
}


##if you have all data in one Excel table, columns=samples, as in Giulia's case:
##you might have to install the "gdata" package, see Tools/Install Packages... - choose "CRAN"
##un-comment these lines below by marking them and press Shift-Ctrl-c:

# require(gdata)
# inpath="...wherever_your_file_is/filename.xls"
#
# #read xls file - this assumes data to be in the first sheet, with headers (sample names) in the first line
# data<-read.xls(inpath,stringsAsFactors=FALSE)
# #the columns aren't all the same length - cut out empty values
# data<-as.list(alldata)
# for(i in 1:length(data)){data[[i]]<-data[[i]][!is.na(data[[i]])]}


##### from here on, it should work the same for both of you:
##### plotting KDEs:
# we need the 'provenance' package
require(provenance)

plotKDE(data,bandwidth=18)
plotKDE(data,logx=TRUE)
plotKDE(data,hist=TRUE)
plotKDE(data,markers="circle")
plotKDE(data,markers="circle",splitat=600)
plotKDE(data,markers="circle",splitat=600,bandwidth=15)

#to choose output file format, just change extension. width and height are in inches.
#For raster files, you can also give the resolution.
ggsave(paste(basepath,"test.pdf",sep="/"),width=10,height=8)
ggsave(paste(basepath,"test.svg",sep="/"),width=10,height=8)
ggsave(paste(basepath,"test.png",sep="/"),width=10,height=8,dpi=300)

##### plotting MDSs:
#not much to see with only few data points...
require(provenance)
require(MASS)

#calculate dissimilarity matrix
diss = dissimilarity(data)
#calculate MDS coordinate, re-shape into a data.frame
mds<-isoMDS(as.matrix(diss))$points
mds<-as.data.frame(mds)
names(mds)<-c("x","y")
#add sample names to data.frame
mds$smpl<-row.names(mds)

plotMDS(mds,diss,col="smpl",nearest=TRUE)
