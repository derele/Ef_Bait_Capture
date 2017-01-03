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

###############
## FUNCTIONS ##
###############

########### Create a function calculating the count of features ###########
## Usage : MyfeatureCounts(SamFiles, paired_or_not) 
MyfeatureCounts <- function(SamFiles, paired_or_not) {
    return(featureCounts(SamFiles,
                         annot.ext="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_complete.gtf",
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

########### Create a tab with the info of the mapped reads positions ###########
## Usage : makemytab(df) with df being a FC$count object and "alignment" name 
makemytab <- function(df,alignment){
    mito <- df[grepl("mito_only_1", rownames(df)),]
    api <- df[grepl("apico_only_1", rownames(df)),]
    ## Remove apico and mito lines
    dfnucl <- df[-c(nrow(df),nrow(df)-1),]
    ## Count the sum of the off target sequences
    ## Efal_genome_Off_Target <- colSums(dfnucl[grepl("Off", rownames(dfnucl)),])
    Efal_genome_0_to_100 <- colSums(dfnucl[grepl("0_to_100", rownames(dfnucl)),])
    Efal_genome_100_to_300<- colSums(dfnucl[grepl("100_to_300", rownames(dfnucl)),])
    Efal_genome_300_to_500<- colSums(dfnucl[grepl("300_to_500", rownames(dfnucl)),])
    Efal_genome_more_than_500 <- colSums(dfnucl[grepl("more_than_500", rownames(dfnucl)),])
    Efal_genome_distant <- colSums(dfnucl[grepl("distant", rownames(dfnucl)),])
        Efal_genome_On_Target <- colSums(dfnucl[!grepl("Off", rownames(dfnucl)),])
    ## Create tab
    BigTab <- data.frame(Efal_genome_On_Target,
                         Efal_genome_0_to_100,
                         Efal_genome_100_to_300,
                         Efal_genome_300_to_500,
                         Efal_genome_more_than_500,
                         Efal_genome_distant,mito, api)
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


########### Create a tab with the coverage of the mapped reads positions ###########
## Usage : mytabcov(df) with df an output of R Feature Count

mytabcov <- function(df,whichpart){
    MyDF <- data.frame(matrix(unlist(df$count), nrow=83690, byrow=F),stringsAsFactors=FALSE)
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
    ## Subset by group : mito, api, off target, on target
    MyDF$Group <- ifelse(grepl("mito", MyDF$baits, ignore.case = T), "Mitochondria",
                  ifelse(grepl("apico", MyDF$baits, ignore.case = T),  "Apicoplast",
                  ifelse(grepl("0_to_100", MyDF$baits, ignore.case = T),  "Genome 0 to 100bp from bait",
                  ifelse(grepl("100_to_300",MyDF$baits, ignore.case = T),  "Genome 100 to 300bp from bait",
                  ifelse(grepl("300_to_500",MyDF$baits, ignore.case = T),  "Genome 300 to 500bp from bait",
                  ifelse(grepl("more_than_300",MyDF$baits, ignore.case = T),  "Genome more than 500bp from bait",
                  ifelse(grepl("distant",MyDF$baits, ignore.case = T),  "Genome on a contig with no bait", "Genome ON target")))))))
    ## For the second part, a category is added
    if (whichpart==2){
        MyDF[MyDF$baits %in% BaitsTot,]$Group <- "Selected_baits"
    }
    MyDF2 <- aggregate(MyDF$cov ~ MyDF$variable + MyDF$Group, FUN=mean)
    names(MyDF2) <- c("Libraries","Group", "Coverage")
    return(MyDF2)
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

## Usage : makemytab(df) with df being a FC$count object
Tab1 <- rbind(makemytab(F1_R1$count, alignment="blat"),makemytab(F1_R2$count, alignment="blat"))
rn <- rownames(Tab1)
Tab1 <- Tab1[order(rn), ]
Tab1 <- merge(Numberseq, Tab1, by = intersect(names(Numberseq), names(Tab1)))
Tab1$not_aligned <- Tab1$totalseq - Tab1$TotalAligned
write.table(x=Tab1, file="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_Dna.csv", row.names = FALSE, quote=FALSE)

### Some visualisation :
Tabplot1 <- melt(Tab1, id=c("libraries","TotalAligned","libshort","totalseq"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_Dna.pdf")
myplot1 <- ggplot(Tabplot1, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    ggtitle("Blat with alignment on nucleotide level, 80% identity")+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    #scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()
myplot1
dev.off()

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
FilesR1 <- FilesR1[-19]; FilesR2 <- FilesR2[-19]
## and files crashing R : 2808Do,2809Do,2812Double,2919Single, 2809Si, 2812Si
##FilesR1 <- FilesR1[-c(5,8,11,15,9,12)]
##FilesR2 <- FilesR2[-c(5,8,11,15,9,12)]

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
write.table(x=Tab2, file="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Rsubread_align.csv", row.names = FALSE, quote=FALSE)

### Some visualisation :
Tabplot2 <- melt(Tab2, id=c("libraries","TotalAligned","libshort","meanpairedseq"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Rsubread_align.pdf")
myplot2 <- ggplot(Tabplot2, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    ggtitle("Rsubread:align, MM=30 (20%), MI=10")+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#    scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()
myplot2
dev.off()
## NB : if Multipleoverlap not allowed, it's not comparable (one read can touch On and Off targets)

######################
## III. FeactureCount with a sam from BLAT X alignment 80% similitude protein level
setwd("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/Blat_Dna_Dnax/GOOD/")
SamsXR1 <- list.files(pattern="R1_001.fastq.unique.fasta.fa_blatDnax.sam$", full.names=TRUE)
SamsXR2 <- list.files(pattern="R2_001.fastq.unique.fasta.fa_blatDnax.sam$", full.names=TRUE)

## Usage : MyfeatureCounts(SamFiles, paired_or_not, multiple_overlap_or_not) 
F3_R1 <- MyfeatureCounts(SamsXR1, FALSE)
F3_R2 <- MyfeatureCounts(SamsXR2, FALSE)

## Usage : makemytab(df) with df being a FC$count object
Tab3 <- rbind(makemytab(F3_R1$count, alignment="blat"),makemytab(F3_R2$count, alignment="blat"))
rn <- rownames(Tab1)
Tab3 <- Tab1[order(rn), ]
Tab3 <- merge(Numberseq, Tab3, by = intersect(names(Numberseq), names(Tab3)))
Tab3$not_aligned <- Tab3$totalseq - Tab3$TotalAligned
write.table(x=Tab3, file="/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_DnaX.csv", row.names = FALSE, quote=FALSE)

### Some visualisation :
Tabplot3 <- melt(Tab3, id=c("libraries","TotalAligned","libshort","totalseq"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Blat_DnaX.pdf")
myplot3 <- ggplot(Tabplot3, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    ggtitle("Blat with alignment on protein, 80% identity")+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#    scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()
myplot3
dev.off()
## NB : if Multipleoverlap not allowed, it's not comparable (one read can touch On and Off targets)

######################
## COMPARISONS :
## First, stick together the reads for R1 and R2 when applicable
Tabplot1bis <- aggregate(Tabplot1$value ~ Tabplot1$libshort + Tabplot1$variable, FUN=sum)
names(Tabplot1bis) <- c("libraries", "variable", "value")

Tabplot2bis <- Tabplot2[c("libraries", "variable", "value")]

Tabplot3bis <- aggregate(Tabplot3$value ~ Tabplot3$libshort + Tabplot3$variable, FUN=sum)
names(Tabplot3bis) <- c("libraries", "variable", "value")

Tabplot1bis$aligner <- "Blat_DNA"
Tabplot2bis$aligner <- "Rsubread_align"
Tabplot3bis$aligner <- "Blat_DNAX"
TabAll <- rbind(Tabplot1bis, Tabplot2bis, Tabplot3bis)

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/All_aligners.pdf", width=20)
myplotAll <- ggplot(TabAll, aes(x=libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
 #   scale_fill_manual(values=c("red","orange","darkblue","darkgreen","grey"))+
    ylab("Number of reads per category")+
    scale_y_log10()+
    facet_grid(aligner ~ .)
myplotAll
dev.off()

##################################
#######Now with coverages########

## Blat : M1
## NB : we use the average value between R1 and R2 alignments
M1.1 <- mytabcov(F1_R1,"first")
M1.1$Libraries <-substr(sub("_S.*.","",M1.1$Libraries),3,50)
M1.2 <- mytabcov(F1_R2,"first")
M1.2$Libraries <-substr(sub("_S.*.","",M1.2$Libraries),3,50)

M1 <- rbind(M1.1,M1.2)
M1 <- aggregate(M1$Coverage ~ M1$Libraries + M1$Group, FUN=mean)
M1$Alignment <- "Blat DNA"
names(M1) <- c("Libraries", "Group", "Coverage", "Alignment")

## Rsubread:align : M2
M2 <- mytabcov(F2,"first")
M2$Libraries <- sub(".*align.","",sub("_S.*.","",M2$Libraries))
M2$Alignment <- "Rsubread"
names(M2) <- c("Libraries", "Group", "Coverage", "Alignment")

## Blatx : M3
M3.1 <- mytabcov(F3_R1,"first")
M3.1$Libraries <-substr(sub("_S.*.","",M3.1$Libraries),3,50)
M3.2 <- mytabcov(F3_R2,"first")
M3.2$Libraries <-substr(sub("_S.*.","",M3.2$Libraries),3,50)

M3 <- rbind(M3.1,M3.2)
M3 <- aggregate(M3$Coverage ~ M3$Libraries + M3$Group, FUN=mean)
M3$Alignment <- "Blat DNAX"
names(M3) <- c("Libraries", "Group", "Coverage", "Alignment")

##### ALL :
Mtot <- rbind(M1, M2, M3)

## Coverage per Mb
Mtot$Coverage <- Mtot$Coverage * 1000000

Mtot$Group <- factor(Mtot$Group, levels = c("Genome ON target", "Genome 0 to 100bp from bait", "Genome 100 to 300bp from bait", "Genome 300 to 500bp from bait", "Genome on a contig with no bait", "Mitochondria", "Apicoplast"))

pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/All_aligners_coverage.pdf", width=20)
myplotAllCOV <- ggplot(Mtot, aes(x=Libraries, y=Coverage, fill=Group))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values=c("#990000", "#CC0000", "#FF0000", "#FF3366", "#FF6699", "darkgreen","darkblue"))+
    ylab("Coverage per Mb for each type of region")+
    scale_y_log10()+
    facet_grid(Alignment ~ .)
myplotAllCOV
dev.off()

## In a table, to compare with coverage at target region
Mtot <- Mtot[order(Mtot$Libraries, Mtot$Alignment),]
Mtot$Coverage <- round(Mtot$Coverage)

Mwide <- reshape(Mtot, timevar = "Group",
        idvar = c("Libraries", "Alignment"),
        direction = "wide")

for (i in 3:9) {
    Mwide[i] <- Mwide[i]/Mwide[8]
    }

M <- melt(Mwide, id=c("Libraries","Alignment"))

M$variable <- factor(M$variable, levels = c("Coverage.Genome ON target", "Coverage.Genome 0 to 100bp from bait", "Coverage.Genome 100 to 300bp from bait", "Coverage.Genome 300 to 500bp from bait", "Coverage.Genome on a contig with no bait", "Coverage.Mitochondria", "Coverage.Apicoplast"))

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
    Tab <- Tab[-1]
    return(Tab)}

## For Blat DNA
pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/BlatDna_coverage_ratio_table.pdf",height=4, width=20)
Tab <- TableFunc(M[M$Alignment=="Blat DNA",])
grid.table(Tab)
dev.off()

## For Blat DNAX
pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/BlatDnaX_coverage_ratio_table.pdf",height=4, width=20)
Tab <- TableFunc(M[M$Alignment=="Blat DNAX",])
grid.table(Tab)
dev.off()

## For Rubread
pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Rsubread_coverage_ratio_table.pdf",height=4, width=20)
Tab <- TableFunc(M[M$Alignment=="Rsubread",])
grid.table(Tab)
dev.off()

############## Percentage enrichment by fold ##############







########################################## 
## Then, which baits have which coverage?#
##################### ####################

### For each library, get the number of "efficient" baits

### Q°1 : how many different baits captured some reads?
library(stringr)

## To erase?
## Ex 1 : F1_R1
#DfBaits <- as.data.frame(as.table(F1_R1$counts))
#names(DfBaits) <- c("region","libraries","coverage")
#DfBaits$libraries <- substr(DfBaits$libraries,3,8)
# Define a type of region (apico, mito, EfaB=bait, OT=Off Target)
#DfBaits$type <- str_split_fixed(DfBaits$region, "_",3)[,1]


## For later heatmaps of the baits vs a heatmap
C1.1 <- F1_R1$counts
colnames(C1.1) <- substr(colnames(C1.1),3,8)
C1.1.baits <- C1.1[grep("EfaB_", rownames(C1.1)), ]

C1.2 <- F1_R2$counts
colnames(C1.2) <- substr(colnames(C1.2),3,8)
C1.2.baits <- C1.2[grep("EfaB_", rownames(C1.2)), ]

C2 <- F2$counts
colnames(C2) <- substr(colnames(C2),81,86)##
C2.baits <- C2[grep("EfaB_", rownames(C2)), ]

C3.1 <- F3_R1$counts
colnames(C3.1) <- substr(colnames(C3.1),3,8)
C3.1.baits <- C3.1[grep("EfaB_", rownames(C3.1)), ]

C3.2 <- F3_R2$counts
colnames(C3.2) <- substr(colnames(C3.2),3,8)
C3.2.baits <- C3.2[grep("EfaB_", rownames(C3.2)), ]

## Which baits? Selection of the efficient baits
FindBait <- function(counts){
    ## How many baits have a coverage of min 10? Calculate on each libraries?
    names <- colnames(counts)
    mybiglist <- list()
    for (i in 1:ncol(counts)){
        vec <- names(counts[counts[,i]>10,i]) # baits with a coverage > 10 for lib i are stored
        name <- names[i]
        tmp <- list(Bait=vec)
        mybiglist[[name]] <- tmp
        print(names[i])
        print(table(counts[,i]>10))
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

BaitsC1.1 <- FindBait(C1.1.baits)

BaitsC1.2 <- FindBait(C1.2.baits)

BaitsC2 <- FindBait(C1.2.baits)

BaitsC3.1 <- FindBait(C3.1.baits)

BaitsC3.2 <- FindBait(C3.2.baits)

## Baits obtained from ALL aligners
BaitsTot <-c(BaitsC1.1,BaitsC1.2, BaitsC2, BaitsC3.1, BaitsC3.2)

TabTot <- table(BaitsTot)

## Questions: are the selected baits the same in all alignments?
length(TabTot) == length(TabTot[TabTot==5])

### -> How many baits do I select?
length(TabTot[TabTot>=5])
## Stored in :
BaitsTot <- names(TabTot[TabTot==5])

#############################
## Heatmaps of "good" baits##
#############################

## "Good" baits stored in BaitsTot
pdf("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Heatmap_GB_blatR1.pdf")

adf <- as.data.frame(C1.1.baits)
subdf <- adf[rownames(adf) %in% BaitsTot,] # select only the good baits
pheatmap(log10(as.matrix(subdf)+0.1),show_rownames=F, main="Targeted positions, blat R1")

dev.off()

############## Are there baits with a really high coverage? ##############
## Distribution of coverage among baits : 
png("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/Coverage_among_baits_Rsubreadalign.png")
adf <- as.data.frame(C2.baits)
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
dev.off()



########## DE NOVO METAGENOME ASSEMBLY ##########
T1 <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_tblastn_Efalci.txt", header=T)

T2 <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Efalci.txt", header=T)

# T3 <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu.txt", header=T)

T <- cbind(T1[-2], T2$percentmatch)
names(T) <- c("Libraries", "N scaffolds", "tblastn E. falci","blastn E. falci")

T <- melt(T, id=c("Libraries","N scaffolds"))

png("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Results_80_plot.png", width=800, height=800)
ggplot(T, aes(x=Libraries, y=value, fill=variable))+
    geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7)+
    theme_minimal()+
    labs(fill = "")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=30),
          axis.text.y = element_text(size=30),
          axis.title = element_text(size=30),
          legend.text= element_text(size=30))+
    ylab("Percentage of scaffolds > 80% similarities to reference")
dev.off() 














