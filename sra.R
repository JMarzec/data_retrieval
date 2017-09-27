################################################################################
#
#   File name: sra.R
#
#   Authors: Jacek Marzec  ( j.marzec@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Getting data from SRA using SRAdb package
#
#	Command line use: R --file=./sra.R --args "E-MTAB-567" "Y" "/scratch/jack/data/PhD/raw/Illum_HiSeq2000_RNAseq" - this will get the raw files if they exist
#
#	First arg: accession number
#	Second arg: tidy data arg = Y will try to remove extra names from GEO patient ids
#	Third arg: the directory for the data
#
################################################################################

#===============================================================================
#    Functions
#===============================================================================

##### Returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

##### Get the study design of study recoreded in GEO
getGEODetails <- function(geo_eset) {

	print(geo_eset)

	Name=as.character(pData(phenoData(geo_eset))[,c(1)])

	if(TidyNames == "Y") {
		Name=gsub("-","_",Name)
		Name=gsub(" ","_",Name)
	}

	##### Get sample source name
	Characteristics=sapply(pData(phenoData(eset))[,names(pData(phenoData(eset)))=="source_name_ch1"],trim)

	##### Get sample characteristics
	for (i in 1:length(names(pData(phenoData(eset)))) ) {

		if ( length(grep("characteristics",names(pData(phenoData(eset)))[i])) != 0 ) {

			Characteristics=data.frame(Characteristics,pData(phenoData(eset))[i])
		}
	}
	names(Characteristics)[1]="Source_name"

	phenodata=data.frame(Name,Characteristics)
	print(phenodata)
	sampleNames(geo_eset)=gsub(" ","_",as.character(pData(phenoData(geo_eset))[,c(1)]))

	return(phenodata)
}

#===============================================================================
#    Main
#===============================================================================

library("SRAdb")

Args <- commandArgs();

Acc=Args[4]
TidyNames=Args[5]
DataDir=paste(Args[6], Acc, sep="/")

print(Acc)
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

##### Download the most recent SRAdb database file
sqlfile = getSRAdbFile()

##### Connect the SRAdb database to R
sra_con = dbConnect(SQLite(), sqlfile)

#####  Search SRA database with specified query
rs = getSRA(search_terms = Acc, out_types = c("run", "study"), sra_con = sra_con)

cat(paste("Study accession number:", rs$study[1], "\n", sep=" "))

#####  Retrieve data file information
SRAfiles <- listSRAfile(rs$study[1], sra_con = sra_con, fileType = "fastq", srcType = "ftp")

cat(paste(nrow(SRAfiles), "fastq files", "for", nrow(rs), "samples\n", sep=" "))

#####  Check the number of fastq files
if ( nrow(SRAfiles)/2 < nrow(rs) ) {

	cat(" !!! WARNING: Some fastq files may be missing !!!")
}

##### Get fastq files info
FASTQinfo <- getFASTQinfo(in_acc = rs$study[1], srcType = "ftp")

missingFiles=NULL
existingFiles=NULL

##### Get fastq files from SRA
for (i in 1:nrow(FASTQinfo)) {

	##### Check if there are any files missing
	if ( FASTQinfo$file.name[i]=="Not available" ) {

		cat(" ##################################################\n", "#\n")
		cat(" #  !!! fastq file", "for", FASTQinfo$run[i], "is missing !!!\n #\n")
		cat(" ##################################################", "\n\n")

		missingFiles <- rbind(missingFiles, FASTQinfo[i,])

	##### Check if the files for current run were already downloaded
	} else if ( file.exists(paste(DataDir,FASTQinfo$file.name[i],sep="/")) ) {

		cat(" ##################################################\n", "#\n")
		cat(" #  ! fastq file", "for run ", FASTQinfo$run[i], "already exist !\n #\n")
		cat(" ##################################################", "\n\n")

	} else {
		cat(" Downloading fastq files for run:", FASTQinfo$run[i], "\n")
		getFASTQfile(in_acc = FASTQinfo$run[i], destDir = getwd(),srcType = "ftp")

		##### Record samples for which fastq files were downloaded
		existingFiles <- c(existingFiles, FASTQinfo$run[i])
	}
}

write.table(data.frame(missingFiles),file=paste(DataDir,"/missing_fastq_files_", Acc,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

##### Get experiment info for records in GEO
if ( length(grep("^GSE", Acc)) == 1 ) {

	library("GEOquery")
	gse_eset <- getGEO(Acc,GSEMatrix=TRUE)

	##### Get the first platform here
	eset=gse_eset[[1]]
	phenodata=getGEODetails(eset)

	write.table(phenodata,file=paste("phenodata_",Acc,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
	write.table(FASTQinfo,file=paste(Acc,".sdrf.txt", sep=""),sep="\t",row.names=FALSE,quote=FALSE)

##### Get experiment info for records in ArrayExpress
} else {

	library("ArrayExpress")
	rawset = ArrayExpress(Acc, path = DataDir, save = FALSE)
}

##### Remove the SRAdb database file
SRAdb2rm=paste("rm", sqlfile, sep=" ")
system( SRAdb2rm )
