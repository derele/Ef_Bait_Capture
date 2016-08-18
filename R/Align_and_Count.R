setwd("~/Eimeria_Wild_coding/Bait_capture/")
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
###################### IDEAS ###################### 
#Hierarchical clustering on the readcounts per bait.  V
#A heatmap just on this raw data?  V
#Pairs plot to identify correlating (good) libraries.  V
#Table of pairwise correlation coefficients.  V
#And a truly random 120nt baits gff would be awesome ;-)
#Just random no intron or etc selection needed.
#My guess is the hierarchical clustering will show a cluster of "working" baits. And the correlations will show "working" libraries.
#column scaling of the libraries might give the best clustering.
#is there reads to non baits region in the genomes?
################################################### 

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
FilesR1 <- FilesR1[-c(9,11,15)]
FilesR2 <- FilesR2[-c(9,11,15)]

## Create the alignments ## SOME BUGS with R, ran outside ## /home/alice/subread-1.5.0-p3-source/bin/subread-align -T 16 -i Efal_mtapi_reference_index -r 2812Single_S10_R1_001.fastq.unique.gz -R 2812Single_S10_R2_001.fastq.unique.gz -I 10 -M 3 -t 1 -o 2812Single_S10_R1_001.fastq.unique.gz_3.BAM
Mismatches <- c(10,20,30,50,70)
if(file.exists("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/2672Single_S17_R1_001.fastq.unique.gz_BAM")){
    sapply(Mismatches, function (MM){
        Rsubread::align(index="/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi_reference_index",
                      readfile1=FilesR1,
                      readfile2=FilesR2,
                      type="dna",
                      maxMismatches=MM,
                      indels=10,
                      nthreads=length(FilesR1),
                      output_file=paste(FilesR1, "_", MM, ".BAM", sep="")
                      )
})}

## Create the featurecounts objects
bams_10<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/",
                  pattern="10.BAM$", full.names=TRUE)
FC_10<- featureCounts(bams_10,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                    isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
                    allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
                    isPairedEnd=TRUE,reportReads=TRUE,nthreads=20)
bams_20<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/",
                  pattern="20.BAM$", full.names=TRUE)
FC_20<- featureCounts(bams_20,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                    isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
                    allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
                    isPairedEnd=TRUE,reportReads=TRUE,nthreads=20)
bams_30<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/",
                  pattern="30.BAM$", full.names=TRUE)
FC_30<- featureCounts(bams_30,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                    isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
                    allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
                    isPairedEnd=TRUE,reportReads=TRUE,nthreads=20)
bams_50<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/",
                  pattern="50.BAM$", full.names=TRUE)
FC_50<- featureCounts(bams_50,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                    isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
                    allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
                    isPairedEnd=TRUE,reportReads=TRUE,nthreads=20)
bams_70<- list.files(path = "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments/",
                  pattern="70.BAM$", full.names=TRUE)
FC_70<- featureCounts(bams_70,annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf",
                    isGTFAnnotationFile=TRUE, GTF.featureType = "sequence_feature", useMetaFeatures=FALSE, GTF.attrType="bait_id",
                    allowMultiOverlap=TRUE,## important to allow reads to map multiple baits
                    isPairedEnd=TRUE,reportReads=TRUE,nthreads=20)

## First look into the data
## Define mismatches max
MM <- 10
FC <-get(paste("FC_", MM, sep=""))
############
countsDF <- FC$counts
colnames(countsDF) <- paste("Lib_",substr(colnames(countsDF),51,56),sep="") # shorten names 
libraries <- 1:10 ; coverage <- 1:100 ##
tabcountcov <- matrix(, nrow = 10, ncol = 100)
for (l in libraries){
    for (c in coverage){
        tabcountcov[l,c] <- (table(rowSums(countsDF > c) > l))[2]}}
colnames(tabcountcov) <- as.character(1:100)
df <- melt(tabcountcov[,c(1,2,5,10,15,20,30,50,100)]); names(df) <- c("nLib","cov","nBaits")
ggplot(data=df, aes(x=nLib, y=nBaits, group=cov)) + geom_line(aes(color=factor(cov))) + scale_color_discrete(name="Divers levels\nof coverage") + ggtitle(paste(as.character(MM),"mismatches maximum allowed in the alignment"))

##Plot a heatmap
table(rowSums(countsDF)>500)
rowSums(countsDF)>500

countsDF[grepl("apico", rownames(countsDF)),]

pheatmap(log10(countsDF[rowSums(countsDF)>500,]+0.1))
##pheatmap(log10(countsDF[rowSums(countsDF)>50,]+0.1))
##pheatmap(log10(countsDF[rowSums(countsDF)>10,]+0.1))

## Plot correlation heatmap + pairwise coef [adapted from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization]
create_heatmap_correl <- function(mydata){
    cormat <- round(cor(mydata),2)# I.Compute the correlation matrix
    get_upper_tri <- function(cormat){# II.Get upper triangle of the correlation matrix
        cormat[lower.tri(cormat)]<- NA; return(cormat)}
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- na.omit(melt(upper_tri))# II.Melt the correlation matrix
    ggplot(data = melted_cormat, aes(X2,X1, fill = value))+ # IV.Plot
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "blue", high = "red",
                             mid = "white",midpoint = 0,
                             limit = c(-1,1), space = "Lab",
                             name="Pearson\nCorrelation") +
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
        coord_fixed()+ geom_text(aes(X2, X1, label = value),
                                 color = "black", size = 4) +  #Add correlation coefficients on the heatmap
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.ticks = element_blank(),
              legend.justification = c(1, 0),
              legend.position = c(0.6, 0.7),
              legend.direction = "horizontal")+
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5)) }

create_heatmap_correl(countsDF[!grepl("apico|mito", rownames(countsDF)),])

#create_heatmap_correl(countsDF)
