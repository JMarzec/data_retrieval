################################################################################
#
#   File name: geo.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk ), code modified from Rossalind Cutts
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Getting data from GEO using GEOquery package
#
#	Command line use: R --file=./geo.R --args "GSE45016" "raw" "Y" "/scratch/jack/data/PhD/raw/Affy_U133Plus2" - this will get the raw files if they exist
#
#	First arg: geo id
#	Second arg: whether data is raw or matrix format
#	Third arg: tidy data arg = Y will try to remove extra names from GEO patient ids
#	Fourth arg: the directory for the data
#
################################################################################

#===============================================================================
#    Functions
#===============================================================================

##### Returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

##### Returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

##### Returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

##### Get raw files from GEO
getRawFiles <-function(GEO, DataDir) {

	getGEOSuppFiles(GEO, makeDirectory = FALSE, baseDir = DataDir)

	untar(list.files(pattern="*.tar"))
	lapply(list.files(pattern="*.gz"),gunzip)
}

##### Get the study design info and normalised expression data
getGEODetails <- function(geo_eset) {

	print(geo_eset)

	##### Get probes annotation
	names(fData(geo_eset))
	geneSymbol=as.character(fData(geo_eset)[,grep("Symbol",names(fData(geo_eset)))])

	##### show(pData(phenoData(geo_eset[[1]]))[1:5,c(1,6,8)])

	Name=as.character(pData(phenoData(geo_eset))[,c(1)])

	if(TidyNames == "Y") {
		Name=gsub("-","_",Name)
		Name=gsub(" ","_",Name)
	}

	##### Get raw data file name
	FileName=gsub(".gz","",sapply(strsplit(as.character(pData(phenoData(eset))[,names(pData(phenoData(eset)))=="supplementary_file"]),"/"), function(elt) elt[length(elt)]))

	##### Get sample source name
	Characteristics=sapply(pData(phenoData(eset))[,names(pData(phenoData(eset)))=="source_name_ch1"],trim)

	##### Get sample characteristics
	for (i in 1:length(names(pData(phenoData(eset)))) ) {

		if ( length(grep("characteristics",names(pData(phenoData(eset)))[i])) != 0 ) {

			Characteristics=data.frame(Characteristics,pData(phenoData(eset))[i])

		}
	}
	names(Characteristics)[1]="Source_name"

	phenodata=data.frame(Name,FileName,Characteristics)
	print(phenodata)
	sampleNames(geo_eset)=gsub(" ","_",as.character(pData(phenoData(geo_eset))[,c(1)]))

	exprsData=exprs(geo_eset)
	head(exprsData)

	##### Get log2 expression values
	#exprsData=log(exprsData,base=2)

	return(list(phenodata,exprsData,geneSymbol))
}

#===============================================================================
#    Main
#===============================================================================

library("GEOquery")

Args <- commandArgs();

GEO=Args[4]
FileType=Args[5]
TidyNames=Args[6]
DataDir=paste(Args[7], GEO, sep="/")

print(GEO)
print(FileType)
print(TidyNames)
print(DataDir)

##### Set/create a directory for data download
if (file.exists(DataDir)){
	setwd(DataDir)
} else {
	dir.create(DataDir, recursive=TRUE);
	setwd(DataDir)
}

##### Report used parameters to a file
write(Args, file = "parameters.txt", append = FALSE, sep="\t")

if(FileType=="raw") {
	getRawFiles(GEO, DataDir)
}

gse_eset <- getGEO(GEO,GSEMatrix=TRUE, destdir=DataDir, AnnotGPL=TRUE)
print(gse_eset)

if(length(gse_eset) == 1) {

	##### Get the first platform here
	eset=gse_eset[[1]]
	geo_data=getGEODetails(eset)
	phenodata=geo_data[[1]]
	exprsData=geo_data[[2]]

	write.table(phenodata,file=paste0(GEO, "_target.txt"),sep="\t",row.names=FALSE,quote=FALSE)
	write.table(exprsData,paste0(GEO, "_normalised.txt"),sep="\t", quote=FALSE, row.names=TRUE,col.names=NA)

} else {
	##### Get info for all platforms
	for (i in names(gse_eset)) {
		print(i)

		eset=gse_eset[[i]]
		gse_val=unlist(strsplit(i,"_"))[1]
		geo_data=getGEODetails(eset)
		phenodata=geo_data[[1]]
		exprsData=geo_data[[2]]

		name=paste0("target_",gse_val,".txt")
		write.table(phenodata,file=name,sep="\t",row.names=FALSE,quote=FALSE)

		sampleNames(eset)=gsub(" ","_",as.character(pData(phenoData(eset))[,c(1)]))
		name=paste0("normalised_",gse_val,".txt")
		write.table(exprsData,file=name,sep="\t", quote=FALSE, row.names=TRUE,col.names=NA)
	}
}

##### Remove the compressed raw data file
if(FileType=="raw") {
    rawFile2rm=paste("rm", paste(GEO,"RAW.tar", sep="_"), sep=" ")
system(rawFile2rm)
}
