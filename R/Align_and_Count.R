setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_ResultsCode/")
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
######################################

## built an index (hash table) for read mapping
if(!file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index.files")){
    buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
               reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta")}
##buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome_reference_index",
##            reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome.fa")

#list our _dereplicated_ files (see tally in bash notes). 
FilesR1 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R1_001.fastq.unique.gz$", full.names=TRUE)
FilesR2 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R2_001.fastq.unique.gz$", full.names=TRUE)

## 2919Single_S12(FilesR1[15]) was excluded, keep crashing (memory pb?), same pb with 2809Si [9] + 2812Do [11]
FilesR1 <- FilesR1[-c(9,11,15,19)] # exclude also "undetermined"
FilesR2 <- FilesR2[-c(9,11,15,19)]

## Create the alignments 
Mismatches <- c(10,20,30,50,70)
## sapply(Mismatches, function (MM){
##         Rsubread::align(index="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
##                       readfile1=FilesR1,
##                       readfile2=FilesR2,
##                       type="dna",
##                       maxMismatches=MM,
##                       indels=10,
##                       nthreads=length(FilesR1),
##                       output_file=paste(FilesR1, "_", MM, ".BAM", sep="")
##                       )})


### Create the featurecounts objects and save them as tab-delimited file which includes identifier, length and counts for each bait in each library

## bams_10<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
##                    pattern="10.BAM$", full.names=TRUE)
## FC_10<- featureCounts(bams_10,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
##                      isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
##                      allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
##                      isPairedEnd=TRUE,reportReads=TRUE)
##write.table(x=data.frame(FC_10$annotation[,c("GeneID","Length")],FC_10$counts,stringsAsFactors=FALSE),file="counts_10.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("BaitID", "LengthBait", paste("Lib",substr(names(FC_10$stat[-1]),82,87), sep="_")))
## bams_20<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
##                   pattern="20.BAM$", full.names=TRUE)
## FC_20<- featureCounts(bams_20,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
##                     isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
##                     allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
##                     isPairedEnd=TRUE,reportReads=TRUE)
#write.table(x=data.frame(FC_20$annotation[,c("GeneID","Length")],FC_20$counts,stringsAsFactors=FALSE),file="counts_20.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("BaitID", "LengthBait", paste("Lib",substr(names(FC_20$stat[-1]),82,87), sep="_")))
## bams_30<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
##                   pattern="30.BAM$", full.names=TRUE)
## FC_30<- featureCounts(bams_30,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
##                     isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
##                     allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
##                     isPairedEnd=TRUE,reportReads=TRUE)
#write.table(x=data.frame(FC_30$annotation[,c("GeneID","Length")],FC_30$counts,stringsAsFactors=FALSE),file="counts_30.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("BaitID", "LengthBait", paste("Lib",substr(names(FC_30$stat[-1]),82,87), sep="_")))
## bams_50<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
##                   pattern="50.BAM$", full.names=TRUE)
## FC_50<- featureCounts(bams_50,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
##                     isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
##                     allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
##                     isPairedEnd=TRUE,reportReads=TRUE)
#write.table(x=data.frame(FC_50$annotation[,c("GeneID","Length")],FC_50$counts,stringsAsFactors=FALSE),file="counts_50.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("BaitID", "LengthBait", paste("Lib",substr(names(FC_50$stat[-1]),82,87), sep="_")))
## bams_70<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
##                   pattern="70.BAM$", full.names=TRUE)
## FC_70<- featureCounts(bams_70,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
##                     isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
##                     allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
##                     isPairedEnd=TRUE,reportReads=TRUE)
#write.table(x=data.frame(FC_70$annotation[,c("GeneID","Length")],FC_70$counts,stringsAsFactors=FALSE),file="counts_70.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("BaitID", "LengthBait", paste("Lib",substr(names(FC_70$stat[-1]),82,87), sep="_")))

## Load the tab delimited files, results from featureCounts
counts10 <- read.table("counts_10.txt", header = TRUE);counts10 <- counts10[-ncol(counts10)]
counts20 <- read.table("counts_20.txt", header = TRUE);counts20 <- counts20[-ncol(counts20)]
counts30 <- read.table("counts_30.txt", header = TRUE);counts30 <- counts30[-ncol(counts30)]
counts50 <- read.table("counts_50.txt", header = TRUE);counts50 <- counts50[-ncol(counts50)]
counts70 <- read.table("counts_70.txt", header = TRUE);counts70 <- counts70[-ncol(counts70)]

## Calculate proportion assigned to baits
propmapped
## To be continued.....


## For each library, get the number of efficient baits
coverage <- 10
l <- ncol(counts10)

## How many baits capture sequences?
myvec1 <- vector(mode="numeric", length=0);myvec2 <- vector(mode="numeric", length=0);myvec3 <- vector(mode="numeric", length=0);myvec4 <- vector(mode="numeric", length=0);myvec5 <- vector(mode="numeric", length=0)
for (i in (3:l)) {
    myvec1 <- c(myvec1,sum(counts10[i]>coverage))
    myvec2 <- c(myvec2,sum(counts20[i]>coverage))
    myvec3 <- c(myvec3,sum(counts30[i]>coverage))
    myvec4 <- c(myvec4,sum(counts50[i]>coverage))
    myvec5 <- c(myvec5,sum(counts70[i]>coverage))
}
mytab <- data.frame(names(counts10)[3:l],myvec1,myvec2,myvec3,myvec4,myvec5)
names(mytab) <- c("Libraries", "MM=10", "MM=20", "MM=30", "MM=50", "MM=70")
pdf("Fig_NumberEfficientBaitsCov10.pdf", height=11, width=8.5, title ="Baits with a minimum coverage of 10")
grid.table(mytab)
dev.off()
## same with coverage = 1 saved in "Fig_NumberEfficientBaitsCov1.pdf"

## Which baits? Selection of the efficient baits
## To be continued.....
GoodBaits <- counts10$BaitID[counts10[3]>coverage]

## Plot a heatmap of the "cool" baits
NamesB <- counts10[,1]
df <- counts10[-c(1,2)]
mymat <- as.matrix(df)
rownames(mymat) <- NamesB



## 1 -> baits with at least a sum of coverage of 500 for all the libraries

pheatmap(mymat[rowSums(mymat)>500,])


help(pheatmap)
table(rowSums(counts10[3:l])>500)
#countsDF[grepl("apico", rownames(countsDF)),]

pheatmap(log10(countsDF[rowSums(countsDF)>500,]+0.1))
