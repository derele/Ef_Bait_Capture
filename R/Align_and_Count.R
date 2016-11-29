library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library("Biostrings")
######################################

## Import the total number of sequences per libraries
Numberseq <- read.csv("/SAN/Alices_sandpit/sequencing_data_dereplicated/NumberOfSequences.csv", sep=" ", header=FALSE)
## Rename human friendly :
myfunc <- function(x){
    return(sub("_00.*.","",x))
}
Numberseq$V2 <- apply(Numberseq[2], 1,myfunc)
myfunc2 <- function(x){
    return(sub("_S.*.","",x))
}
Numberseq$V3 <- apply(Numberseq[2], 1,myfunc2)
names(Numberseq) <- c("totalseq","libraries","libshort")
## Same but with the mean of the R1 and R2 reads
Numberseq_paired <- aggregate(Numberseq$totalseq ~ Numberseq$libshort, FUN = mean )
names(Numberseq_paired) <- c("libshort", "meanpairedseq")
## Remove the number of sequences too low
Number_paired <- Numberseq_paired[-which(Numberseq_paired$meanpairedseq < 500000),]

## FUNCTIONS ##
########### Create a function calculating the count of features ###########
## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
MyfeatureCounts <- function(SamFiles, paired_or_not, multiple_overlap_or_not) {
    return(featureCounts(SamFiles,
                         annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf",
                         isGTFAnnotationFile=TRUE,
                         GTF.featureType= "sequence_feature",
                         useMetaFeatures=FALSE,
                         GTF.attrType="bait_id",
                         allowMultiOverlap=multiple_overlap_or_not,
                         isPairedEnd=paired_or_not,
                         countMultiMappingReads=FALSE)
           )
    }

########### Create a tab with the info of the mapped reads positions ###########
## Usage : makemytab(df) with df being a FC$count object and "alignment" name 
makemytab <- function(df,alignment){
    mito <- df[grepl("mito_only_1", rownames(df)),]
    api <- df[grepl("apico_only_1", rownames(df)),]
    ## Remove apico and mito lines
    dfnucl <- df[-c(nrow(df),nrow(df)-1),]
    ## Count the sum of the off target sequences
    Efal_genome_Off_Target <- colSums(dfnucl[grepl("Off", rownames(dfnucl)),])
    Efal_genome_On_Target <- colSums(dfnucl[!grepl("Off", rownames(dfnucl)),])
    ## Create tab
    BigTab <- data.frame(Efal_genome_On_Target,Efal_genome_Off_Target,mito, api)
    BigTab$TotalAligned <- rowSums(BigTab)
    BigTab$not_aligned <- NA
    ## Rename human friendly style
    if (alignment == "blat") {
        rownames(BigTab) <- substr(sub("_00.*.","", rownames(BigTab)),3,20)
    } else if (alignment == "Rsubread") {
        rownames(BigTab) <- sub(".*_align.","", sub("_S.*.","", rownames(BigTab)))
    }
    BigTab$libraries <- rownames(BigTab)
    return(BigTab)
}

######################################

## I. FeactureCount with a sam from BLAT alignment 80% similitude nucleotide level
## --> Is it trustable? The manual says 95%, we set up manually 80%
setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/Blat_Dna_Dnax/GOOD/")
SamsR1 <- list.files(pattern="R1_001.fastq.unique.fasta.fa_blatDna.sam$", full.names=TRUE)
SamsR2 <- list.files(pattern="R2_001.fastq.unique.fasta.fa_blatDna.sam$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
F1_R1 <- MyfeatureCounts(SamsR1, FALSE, FALSE)
F1_R2 <- MyfeatureCounts(SamsR2, FALSE, FALSE)

## Usage : makemytab(df) with df being a FC$count object
Tab1 <- rbind(makemytab(F1_R1$count, alignment="blat"),makemytab(F1_R2$count, alignment="blat"))
rn <- rownames(Tab1)
Tab1 <- Tab1[order(rn), ]
Tab1 <- merge(Numberseq, Tab1, by = intersect(names(Numberseq), names(Tab1)))
Tab1$not_aligned <- Tab1$totalseq - Tab1$TotalAligned
write.table(x=Tab, file="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_Dna.csv", row.names = FALSE, quote=FALSE)

### Some visualisation :
Tabplot1 <- melt(Tab1, id=c("libraries","TotalAligned","libshort","totalseq"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_Dna.pdf")
myplot1 <- ggplot(Tabplot1, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    ggtitle("Blat with alignment on nucleotide level, 80% identity")+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()
myplot1
dev.off()
## NB : if Multipleoverlap not allowed, it's not comparable (one read can touch On and Off targets)

#################################
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
FilesR2 <- list.files(path = "/SAN/Alices_sandpit/sequencing_data_
ed",
                      pattern="R2_001.fastq.unique.gz$", full.names=TRUE)
## erase "Undetermine" library
FilesR1 <- FilesR1[-19];FilesR2 <- FilesR2[-19]
## and files crashing R : 2808Do,2809Do,2812Double,2919Single, 2809Si, 2812Si
FilesR1 <- FilesR1[-c(5,8,11,15,9,12)]
FilesR2 <- FilesR2[-c(5,8,11,15,9,12)]

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
thebams <- list.files(path="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align",
                      pattern="_30.BAM$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
## F2 <- MyfeatureCounts(thebams, FALSE, FALSE)

## Usage : makemytab(df) with df being a FC$count object and "alignment" name 
Tab2 <- makemytab(F2$count, "Rsubread")
rn <- rownames(Tab2)
Tab2 <- Tab2[order(rn), ]
### Extra step !!
Tab2$libshort <- Tab2$libraries
## Paired!
Tab2 <- merge(Numberseq_paired, Tab2, by = intersect(names(Numberseq_paired), names(Tab2)))
## Paired!
Tab2$not_aligned <- Tab2$meanpairedseq - Tab2$TotalAligned
write.table(x=Tab, file="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Rsubread_align.csv", row.names = FALSE, quote=FALSE)

### Some visualisation :
Tabplot2 <- melt(Tab2, id=c("libraries","TotalAligned","libshort","meanpairedseq"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Rsubread_align.pdf")
myplot2 <- ggplot(Tabplot2, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    ggtitle("Rsubread:align, MM=30 (20%), MI=10")+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()
myplot2
dev.off()
## NB : if Multipleoverlap not allowed, it's not comparable (one read can touch On and Off targets)

#################################
## COMPARISONS :
## First, stuck together the reads for R1 and R2 when applicable
Tabplot1bis <- aggregate(Tabplot1$value ~ Tabplot1$libshort + Tabplot1$variable, FUN=sum)
names(Tabplot1bis) <- c("libraries", "variable", "value")

Tabplot2bis <- Tabplot2[c("libraries", "variable", "value")]

Tabplot1bis$aligner <- "Blat_DNA"
Tabplot2bis$aligner <- "Rsubread_align"
TabAll <- rbind(Tabplot1bis, Tabplot2bis)

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/All_aligners.pdf", width=20)
myplotAll <- ggplot(TabAll, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()+
    facet_grid(aligner ~ .)
myplotAll
dev.off()







##################### Old : align by Rsubread:align for different MM --> not checked

## Execute the featureCount with this new GTF
for (MM in c(0,1,3,10,30,40,50,70)){

MM <- 0
    thebams <- list.files(path="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/",
                          pattern=paste("_",MM,".BAM$",sep=""), full.names=TRUE)
    theFC <- featureCounts(thebams, annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf", isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)
    theFC2 <- featureCounts(thebams, annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf", isGTFAnnotationFile=TRUE, GTF.featureType= "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id", allowMultiOverlap=TRUE,isPairedEnd=TRUE,reportReads=TRUE)

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
