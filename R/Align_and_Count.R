#########################################################################
## Align with Rsubread & read the blat alignments done before,         ##
## then count the efficient baits and store them, with R:FeatureCounts.##
#########################################################################
library(Rsubread)
library(reshape)
library(ggplot2)
library(ggthemes)
library(pheatmap)
library(gridExtra)
library("Biostrings")
########################

## Import the total number of sequences per libraries
Numberseq <- read.csv("/SAN/Alices_sandpit/sequencing_data_dereplicated/NumberOfSequences.csv", sep=" ", header=FALSE)

##Remove the libraries with not enough sequences :
Numberseq <- Numberseq[-which(Numberseq$V1 < 500000),]

## Rename human friendly :
Numberseq$V2 <- apply(Numberseq[2], 1,function(x){return(sub("_00.*.","",x))})
Numberseq$V3 <- apply(Numberseq[2], 1, function(x){return(sub("_S.*.","",x))})
names(Numberseq) <- c("totalseq","libraries","libshort")
## Same but with the mean of the R1 and R2 reads
Numberseq_paired <- aggregate(Numberseq$totalseq ~ Numberseq$libshort, FUN = mean )
names(Numberseq_paired) <- c("libshort", "meanpairedseq")

########### Create a function calculating the count of features ###########
MyfeatureCounts <- function(SamFiles, paired_or_not) {
    return(featureCounts(SamFiles,
                         annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_IRanges_IMP_120_mitoapi.gtf",
                         isGTFAnnotationFile=TRUE,
                         GTF.featureType= "sequence_feature",
                         useMetaFeatures=FALSE,
                         GTF.attrType="bait_id",
                         allowMultiOverlap=TRUE,
                         largestOverlap=TRUE, ## if multiple overlap, biggest region counted
                         isPairedEnd=paired_or_not,
                         countMultiMappingReads=FALSE)
           )
    }

######################
## I. FeactureCount with a sam from BLAT alignment 80% similitude nucleotide level
## --> Is it trustable? The manual says 95%, we set up manually 80%
setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/Blat_Dna_Dnax/GOOD/")
SamsR1 <- list.files(pattern="R1_001.fastq.unique.fasta.fa_blatDna.sam$", full.names=TRUE)
SamsR2 <- list.files(pattern="R2_001.fastq.unique.fasta.fa_blatDna.sam$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
F1_R1 <- MyfeatureCounts(SamsR1, FALSE)
F1_R2 <- MyfeatureCounts(SamsR2, FALSE)

######################
## II. FeactureCount with a bam from Rsubread:align alignment 80% similitude nucleotide level (30 MM max allowed, 10 indels allowed)

## II.1. built an index (hash table) for read mapping 
if(!file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index.files")){
    buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
               reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta")}
if(!file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome_reference_index.files")){
    buildindex(basename="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome_reference_index",
            reference="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome.fa")}

## II.2. List our _dereplicated_ files (see tally in bash notes). 
FilesR1 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R1_001.fastq.unique.gz$", full.names=TRUE)
FilesR2 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated",
                      pattern="R2_001.fastq.unique.gz$", full.names=TRUE)

## erase "Undetermine" library
FilesR1 <- FilesR1[-grep("Undete", FilesR1)]; FilesR2 <- FilesR2[-grep("Undete", FilesR2)]

## and files crashing R : 2808Do,2809Do,2812Double,2919Single, 2809Si, 2812Si
## FilesR1crashed <- FilesR1[grep(paste(setdiff(substr(FilesR1,50,55),substr(thebams,80,85)),
##                                     collapse= "|"), FilesR1)]
##FilesR2crashed <- FilesR2[grep(paste(setdiff(substr(FilesR2,50,55),substr(thebams,80,85)),
##                                     collapse= "|"), FilesR1)]

## II.4. Create the alignments 
## Alignments were created for Mismatches <- c(0,1,3,10,20,30,50,70)
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

## II.5. Create the featureCounts objects for a given "maxmismatches"
MM <- 30
thebams <- list.files(path="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align", pattern="_30.BAM$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
F2 <- MyfeatureCounts(thebams, FALSE)

######################
## III. FeactureCount with a sam from BLAT X alignment 80% similitude protein level
setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/Blat_Dna_Dnax/GOOD/")
SamsXR1 <- list.files(pattern="R1_001.fastq.unique.fasta.fa_blatDnax.sam$", full.names=TRUE)
SamsXR2 <- list.files(pattern="R2_001.fastq.unique.fasta.fa_blatDnax.sam$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
F3_R1 <- MyfeatureCounts(SamsXR1, FALSE)
F3_R2 <- MyfeatureCounts(SamsXR2, FALSE)

########### Create a tab with the coverage of the mapped reads positions ###########
## Usage : mytabcov(df) with df an output of R Feature Count

mytabcov <- function(df){
    MyDF <- data.frame(matrix(unlist(df$count), nrow=nrow(df$count), byrow=F),stringsAsFactors=FALSE)
    ## Add the libraries names
    names(MyDF) <- colnames(df$count)
    ## Add the regions names
    MyDF$baits <- make.names(rownames(df$count), unique=TRUE)
    ## Add the length of the regions
    MyDF$Length <- df$annotation$Length
    ## Melt all libraries
    MyDF <- melt(MyDF, id=c("baits","Length"))
    ## Calculate the coverage per region
    MyDF$cov <- MyDF$value/MyDF$Length
    substr(MyDF$baits, nchar(MyDF$baits), nchar(MyDF$baits))==0
    ## Subset by group : mito, api, off target, on target
    MyDF$Group <- ifelse(grepl("mito", MyDF$baits, ignore.case = T), "Mitochondria",
                  ifelse(grepl("apico", MyDF$baits, ignore.case = T),  "Apicoplast",
                  ifelse(substr(MyDF$baits, nchar(MyDF$baits), nchar(MyDF$baits))==0,  "On target genome",
                  ifelse(substr(MyDF$baits, nchar(MyDF$baits), nchar(MyDF$baits))==1,  "One bait away",
                  ifelse(substr(MyDF$baits, nchar(MyDF$baits), nchar(MyDF$baits))==2,  "Two baits away", "Three baits away")))))
    ## Mean coverage of the nuclear genome
    MyDFn <- MyDF[(!grepl("Mitochondria",MyDF$Group,)),]
    MyDFn <- MyDFn[(!grepl("Apicoplast",MyDFn$Group,)),]
    MeanNucl <- aggregate(MyDFn$cov ~ MyDFn$variable, FUN=mean); MeanNucl$Group <- "MeanNuclGen"
    MyDF2 <- aggregate(MyDF$cov ~ MyDF$variable + MyDF$Group, FUN=mean)
    names(MyDF2) <- c("Libraries","Group", "Coverage")
    names(MeanNucl)<- c("Libraries", "Coverage","Group")
    MyDF2 <- rbind(MyDF2, MeanNucl)
    return(MyDF2)
}

## Blat : M1
## NB : we use the average value between R1 and R2 alignments
M1.1 <- mytabcov(F1_R1)
M1.1$Libraries <-substr(sub("_S.*.","",M1.1$Libraries),3,50)
M1.2 <- mytabcov(F1_R2)
M1.2$Libraries <-substr(sub("_S.*.","",M1.2$Libraries),3,50)

M1 <- rbind(M1.1,M1.2)
M1 <- aggregate(M1$Coverage ~ M1$Libraries + M1$Group, FUN=mean)
M1$Alignment <- "Blat DNA"
names(M1) <- c("Libraries", "Group", "Coverage", "Alignment")

## Rsubread:align : M2
M2 <- mytabcov(F2)
M2$Libraries <- sub(".*align.","",sub("_S.*.","",M2$Libraries))
M2$Alignment <- "Rsubread"
names(M2) <- c("Libraries", "Group", "Coverage", "Alignment")

## Blatx : M3
M3.1 <- mytabcov(F3_R1)
M3.1$Libraries <-substr(sub("_S.*.","",M3.1$Libraries),3,50)
M3.2 <- mytabcov(F3_R2)
M3.2$Libraries <-substr(sub("_S.*.","",M3.2$Libraries),3,50)

M3 <- rbind(M3.1,M3.2)
M3 <- aggregate(M3$Coverage ~ M3$Libraries + M3$Group, FUN=mean)
M3$Alignment <- "Blat DNAX"
names(M3) <- c("Libraries", "Group", "Coverage", "Alignment")

##### ALL :
Mtot <- rbind(M1, M2, M3)

## Coverage per Mb
Mtot$Coverage <- Mtot$Coverage *10000
Mtot$Group <- as.factor(Mtot$Group)
Mtot$Group <- factor(Mtot$Group, levels = levels(Mtot$Group)[c(2,1,3:5,7,6)])

myplotAllCOV <- ggplot(Mtot, aes(x=Libraries, y=Coverage, fill=Group))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c( "black","darkgreen","darkblue","#990000", "#CC0000", "#FF0000", "#FF3366"))+
    ylab("Coverage per 10 000 bases for each type of region")+
    scale_y_log10()+
    facet_grid(Alignment ~ .)
myplotAllCOV

## In a table, to compare with coverage at target region
Mtot <- Mtot[order(Mtot$Libraries, Mtot$Alignment),]
Mtot$Coverage <- round(Mtot$Coverage)

Mwide <- reshape(Mtot, timevar = "Group",
        idvar = c("Libraries", "Alignment"),
        direction = "wide")

for (i in c(3,5:9)) {
    Mwide[i] <- Mwide[i]/Mwide[4]
    }

M <- melt(Mwide, id=c("Libraries","Alignment"))

M$value <- as.numeric(as.character(M$value))
M$value <- round(M$value,1)

M <- M[order(M$Libraries, M$Alignment),]

TableFunc <- function(Tab){
    Tab <- Tab[-2]
    Tab <- reshape(Tab, timevar = "Libraries",
                   idvar = "variable",
                   direction = "wide")
    rownames(Tab) <- substring(Tab$variable,10,100)
    colnames(Tab) <- substring(colnames(Tab),7,100)
    Tab <- Tab[-1] # remove useless column
    Tab <- Tab[-2,] # remove mean
    return(Tab)}

## For Blat DNA
TableFunc(M[M$Alignment=="Blat DNA",])

## For Blat DNAX
TableFunc(M[M$Alignment=="Blat DNAX",])

## For Rubread
TableFunc(M[M$Alignment=="Rsubread",])

########################################## 
## Then, which baits have which coverage?#
##################### ####################

### For each library, get the number of "efficient" baits

### Q°1 : how many different baits captured some reads?
library(stringr)

## Which baits? Selection of the efficient baits
FindBait <- function(FC, isblat){
    ## Keep only the baits
    X <- FC$counts
    X <- X[-(grepl("mito", rownames(C1.1), ignore.case = T)),]
    X <- X[-(grepl("mito", rownames(C1.1), ignore.case = T)),]
    X <- X[substr(rownames(X), nchar(rownames(X)), nchar(rownames(X)))==0,]
    if (isblat==TRUE) {
        colnames(X) <- substr(colnames(X),3,8)
    } else {
        colnames(X) <- substr(colnames(X),81,86)
    }
    ## How many baits have a coverage of min 10? Calculate on each libraries?
    names <- colnames(X)
    mybiglist <- list()
    for (i in 1:ncol(X)){
        vec <- names(X[X[,i]>10,i]) # baits with a coverage > 10 for lib i are stored
        name <- names[i]
        tmp <- list(Bait=vec)
        mybiglist[[name]] <- tmp
        print(names[i])
        print(table(X[,i]>10))
    }
    ##->let's select the libraries with more than 1000 Good Baits
    GL <- which(sapply(mybiglist, sapply,length)>1000) ## list the "good" libraries 
    ## Which baits overlap?
    mysmalllist <- mybiglist[GL]
    print(rapply(mysmalllist, length)) # number of baits for each library
    myTab <- table(unlist(mysmalllist))# NUMBER OF LIBRARIES CAPTURED BY EACH BAIT
    print(length(myTab[myTab==length(mysmalllist)]))# Number of baits in common to ALL libraries
    ## -> Let's select baits good in AT LEAST 2 good libraries
    ## Which are:
    BAITS <- names(myTab[myTab>=2])
}

BaitsC1.1 <- FindBait(F1_R1, TRUE)
BaitsC1.2 <- FindBait(F1_R2, TRUE)
BaitsC2 <- FindBait(F2, FALSE)
BaitsC3.1 <- FindBait(F3_R1, TRUE)
BaitsC3.2 <- FindBait(F3_R2, TRUE)

## Baits obtained from ALL aligners
BaitsTot <-c(BaitsC1.1,BaitsC1.2, BaitsC2, BaitsC3.1, BaitsC3.2)

TabTot <- table(BaitsTot)

## Questions: are the selected baits the same in all alignments?
length(TabTot) == length(TabTot[TabTot==5])

### -> How many baits do I select?
length(TabTot[TabTot>=5])
## Stored in :
BaitsTot <- names(TabTot[TabTot==5])
          
########## Which contigs contain the baits chosen? #########
a <- sapply(BaitsTot,function(x) strsplit(x,"_"))
a <- data.frame(table(sapply(a,function(x) x[2])))

ggplot(data=a, aes(x=Freq))+
    geom_histogram(aes(fill = ..count..), binwidth=1)+
    scale_fill_gradient("Count", low = "green", high = "red")+
    theme_minimal()+
    labs(title = "Histogram of the contigs containing chosen baits")+
    theme(axis.text = element_text(size=20),
          title = element_text(size=20))+
    annotate("text", x = 100, y = 20, size = 10, label = paste("N = ",nrow(a)))

#############################
## Heatmaps of "good" baits##
#############################

## "Good" baits stored in BaitsTot
adf <- as.data.frame(C1.2.baits)
subdf <- adf[rownames(adf) %in% BaitsTot,] # select only the good baits
pheatmap(log10(as.matrix(subdf)+0.1),show_rownames=F)


############## Are there baits with a really high coverage? ##############
## Distribution of coverage among baits : 
adf <- as.data.frame(C1.2.baits)
subdf <- adf[rownames(adf) %in% BaitsTot,] # select only the good baits
subdf$Baits <-rownames(subdf) # prepare data frame
subdf <- melt(subdf, id="Baits")
ggplot(data= subdf, aes(x=variable, y=value, fill=variable))+
    geom_jitter(height = 0, alpha=.05)+
    theme_minimal()+
    geom_violin(scale="width", alpha=.8)+
    labs(title = "Number of reads aligned to targeted positions per library")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
          axis.text.y = element_text(size=10),
          axis.title = element_blank(),
          title = element_text(size=10),
          legend.position= "none")+
    scale_y_continuous(breaks=seq(0,900,50))

######## Mean & median of coverage of the selected baits?
med <- aggregate(value~variable, data=subdf, median)
mea <- round(aggregate(value~variable, data=subdf, mean)[2],0)
Meda <- cbind(med, mea)
names(Meda) <- c("Libraries", "coverage median", "coverage mean")
Meda
