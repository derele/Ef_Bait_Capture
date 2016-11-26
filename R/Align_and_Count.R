setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_ResultsCode/")
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library("Biostrings")
######################################

## I. built an index (hash table) for read mapping 
if(!file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index.files")){
    buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
               reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta")}
if(!file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome_reference_index.files")){
    buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome_reference_index",
            reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome.fa")}

## II. List our _dereplicated_ files (see tally in bash notes). 
FilesR1 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R1_001.fastq.unique.gz$", full.names=TRUE)
FilesR2 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_
ed",
                      pattern="R2_001.fastq.unique.gz$", full.names=TRUE)
## erase "Undetermine" library
FilesR1 <- FilesR1[-19];FilesR2 <- FilesR2[-19]
## and files crashing R : 2808Do,2809Do,2812Double,2919Single, 2809Si, 2812Si
FilesR1 <- FilesR1[-c(5,8,11,15,9,12)]
FilesR2 <- FilesR2[-c(5,8,11,15,9,12)]

## III. Create the alignments 
## Alignements were created for Mismatches <- c(0,1,3,10,20,30,50,70)
Mismatches <- NA
sapply(Mismatches, function (MM){
        Rsubread::align(index="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
                      readfile1=FilesR1,
                      readfile2=FilesR2,
                      type="dna",
                      maxMismatches=MM,
                      indels=10,
                      nthreads=length(FilesR1),
                      output_file=paste(FilesR1, "_", MM, ".BAM", sep="")
                      )})

## IV. Create the featureCounts objects for a given "maxmismatches"
MM <- NA
thebams <- list.files(path="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
                      pattern=paste("_",MM,"onlyGenome.BAM$",sep=""), full.names=TRUE)
theFC <- featureCounts(thebams, annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                       isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE,
                       GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)

#################################
## FeactureCount with a bam from bed from psl from BLAT

setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/PslDnax/")
# import psl output of blat
mytab <- read.table("2TRAnna_S16_R1_001.unique.dnax.psl", skip=5)
names(mytab) <- c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert","tNumInsert", "tBaseInsert",
                  "strand", "qName","qSize", "qStart","qEnd", "tName", "tSize", "tStart","tEnd", "blockCount",
                  "blockSizes", "qStarts", "tStarts")
head(mytab)
####




featureCounts("../All_alignment_Blat/2TRAnna_S16_R1_001.unique.dnax.bam",
              annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
              isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE,
              GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)



#################################

## V. FeatureCounts with a GTF that contains genomeoffbaits + mito + api + baits
## V.1. Prepare the good GTF

## load the Genome
fastaFile <- readDNAStringSet("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)
## load the baits
baits <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf")

## get the contigs names
Contigs <- as.character(df$seq_name)

## contigs NOT used for the baits design :
ContigsOff <- setdiff(as.character(df$seq_name),levels(baits$V1))

## get the length of these contigs :
lengthContig <- nchar(as.character(df$sequence))

## in a data frame :
datcontig <- data.frame(Contigs, lengthContig)
datcontig$Contigs <- as.character(datcontig$Contigs)

## Phase 1: exclude the genes used for baits :
datcontigOFF <- datcontig[!(datcontig$Contigs %in% levels(baits$V1)),]

Offbaits <- data.frame(datcontigOFF$Contigs,"rtracklayer", "sequence_feature",0,datcontigOFF$lengthContig,".","+",".","gene_id",datcontigOFF$Contigs,";","bait_id",paste("Off_target_",seq(1:nrow(datcontigOFF)),sep=""))

names(Offbaits) <- names(baits)

Fullbaits <- rbind(baits,Offbaits)

## Phase 2: add the piece without contigs, from contigs where baits where designed

### Function to "Fill the Gaps" of the gff file

mydata <- data.frame(NA,NA,NA)
if (baits$V4[1]!=0) {
    mydata <- data.frame(baits$V1[1],0, baits$V4[1])
}
for (i in 2:nrow(baits)) {
    if (baits$V1[i]==baits$V1[i-1]) {
        if (baits$V4[i] != baits$V5[i-1]+1) {
            new <- data.frame(baits$V1[i],baits$V5[i-1]+1,baits$V4[i]-1)
            names(new) <- names(mydata)
            mydata <- rbind(mydata,new)
        }
    }
    else {
        if (baits$V4[i]!=0) {
            new <- data.frame(baits$V1[i],0,baits$V4[i]-1)
            names(new) <- names(mydata)
            mydata <- rbind(mydata,new)
        }
    }
}

vec <- paste("Off_target_",
             seq(from=(nrow(datcontigOFF)+1),to=(nrow(datcontigOFF) + nrow(mydata)) ),  sep="")

Offbaits2 <- data.frame(mydata$baits.V1.1,"rtracklayer", "sequence_feature",mydata$X0, mydata$baits.V4.1., ".", "+", ".", "gene_id", mydata$baits.V1.1.,";","bait_id", vec)

names(Offbaits2) <- names(Fullbaits)

Fullbaitsfinal <- rbind(Fullbaits,Offbaits2)

# Check for NAs?
sapply(Fullbaitsfinal, function(x) sum(is.na(x)))

## Export my GTF file
write.table(x=Fullbaitsfinal, file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf", sep="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
########################################

## Execute the featureCount with this new GTF
for (MM in c(0,1,3,10,30,40,50,70)){

MM <- 0
    thebams <- list.files(path="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
                          pattern=paste("_",MM,".BAM$",sep=""), full.names=TRUE)
    theFC <- featureCounts(thebams, annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf", isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)
    theFC2 <- featureCounts(thebams, annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf", isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)
    




    theFC2$stat$X.SAN.Alices_sandpit.sequencing_data_dereplicated.All_alignments_Rsubread_align..2672Single_S17_R1_001.fastq.unique.gz_0onlyGenome.BAM
    theFC$stat$X.SAN.Alices_sandpit.sequencing_data_dereplicated.All_alignments_Rsubread_align..2672Single_S17_R1_001.fastq.unique.gz_0onlyGenome.BAM
    



## Make a tab summarizing
Mymat <- cbind(stat0,stat1,stat10,stat30)
Mymat <- Mymat[ , order(names(Mymat))]
pdf("Mymat.pdf", height=10, width=50)
grid.table(Mymat)
dev.off()


#####################################


### For each library, get the number of efficient baits
coverage <- 10
l <- ncol(counts0)

## How many baits capture sequences?
myvec1 <- vector(mode="numeric", length=0);myvec2 <- vector(mode="numeric", length=0);myvec3 <- vector(mode="numeric", length=0);myvec4 <- vector(mode="numeric", length=0);myvec5 <- vector(mode="numeric", length=0);myvec6 <- vector(mode="numeric", length=0)
for (i in (3:l)) {
    myvec1 <- c(myvec1,sum(counts10[i]>coverage))
    myvec2 <- c(myvec2,sum(counts20[i]>coverage))
    myvec3 <- c(myvec3,sum(counts30[i]>coverage))
    myvec4 <- c(myvec4,sum(counts50[i]>coverage))
    myvec5 <- c(myvec5,sum(counts70[i]>coverage))
    myvec6 <- c(myvec5,sum(counts0[i]>coverage))
}
mytab <- data.frame(names(counts10)[3:l],myvec1,myvec2,myvec3,myvec4,myvec5)
names(mytab) <- c("Libraries", "MM=10", "MM=20", "MM=30", "MM=50", "MM=70")
pdf("Fig_NumberEfficientBaitsCov0.pdf", height=11, width=8.5, title ="Baits with a minimum coverage of 10")
grid.table(mytab)
dev.off()
## same with coverage = 1 saved in "Fig_NumberEfficientBaitsCov1.pdf"

## Which baits? Selection of the efficient baits
## To be continued.....
GoodBaits <- counts0$BaitID[counts0[3]>coverage]
head(GoodBaits)

## Plot a heatmap of the "cool" baits
NamesB <- counts0[,1]
df <- counts0[-c(1,2)]
mymat <- as.matrix(df)
rownames(mymat) <- NamesB

## Baits with at least a sum of coverage of 500 for all the libraries
pdf("Fig_heatmap500.pdf",paper="a4")

pheatmap(log10(mymat[rowSums(mymat)>500,]+0.1))

dev.off()
