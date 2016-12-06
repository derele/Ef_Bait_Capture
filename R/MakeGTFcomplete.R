## Complete the GTF
#################################
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library("Biostrings")
library(stringr)


## PART I. Make a GTF that contains genomeoffbaits + mito + api + baits :

## load the Genome
fastaFile <- readDNAStringSet("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_genome.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

## load the baits
baitstot <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf")

## !! Remove api and mito, back at the end!!
baits <- baitstot[-c(nrow(baitstot),nrow(baitstot)-1),]

## get the contigs names
Contigs <- as.character(df$seq_name)

## get the length of these contigs :
lengthContig <- nchar(as.character(df$sequence))

## in a data frame :
datcontig <- data.frame(Contigs, lengthContig)
datcontig$Contigs <- as.character(datcontig$Contigs)

## Phase 1: exclude the contigs used for baits :
datcontigOFF <- datcontig[!(datcontig$Contigs %in% levels(baits$V1)),]

Offbaits <- data.frame(datcontigOFF$Contigs,"rtracklayer", "sequence_feature",1,datcontigOFF$lengthContig,".","+",".","gene_id",datcontigOFF$Contigs,";","bait_id",paste("Off_target_",seq(1:nrow(datcontigOFF)),sep=""))

names(Offbaits) <- names(baits)

Fullbaits <- rbind(baits,Offbaits)

## Phase 2: add the piece without contigs, from contigs where baits where designed

## Reorder "baits"

## extract the contig number in 1 column
baits2 <- cbind(baits, str_split_fixed(baits$V1, "_", 2))
colnames(baits2)[15] <- "num"
baits2$num <- as.numeric(as.character(baits2$num))

## and order by contigs, and within contigs
baits2 <- with(baits2, baits2[order(num,V4),])
## then discard useless columns
baits2 <- baits2[-c(14,15)]

### Function to "Fill the Gaps" of the gff file
mydata <- data.frame(NA,NA,NA)

## initialisation
if (baits2$V4[1]!=1) {
    mydata <- data.frame(baits2$V1[1],1, baits2$V4[1]-1)
}
#then fill the gaps
for (i in 2:nrow(baits2)) {
    if (baits2$V1[i]==baits2$V1[i-1]) {
        if (baits2$V4[i] != baits2$V5[i-1]+1) {
            new <- data.frame(baits2$V1[i],baits2$V5[i-1]+1,baits2$V4[i]-1)
            names(new) <- names(mydata)
            mydata <- rbind(mydata,new)
        }
    }
    else {
        if (baits2$V4[i]!=1) {
            new <- data.frame(baits2$V1[i],1,baits2$V4[i]-1)
            names(new) <- names(mydata)
            mydata <- rbind(mydata,new)
        }
    }
}

# Add "Off_target" to the number of the contigs that are off target
vec <- paste("Off_target_",
             seq(from=(nrow(datcontigOFF)+1),to=(nrow(datcontigOFF) + nrow(mydata)) ),  sep="")

## Bring them in a data frame
Offbaits2 <- data.frame(mydata$baits2.V1.1, "rtracklayer",
                        "sequence_feature", mydata$X1,
                        mydata$baits2.V4.1., ".", "+", ".",
                        "gene_id", mydata$baits2.V1.1.,";",
                        "bait_id", vec)

names(Offbaits2) <- names(Fullbaits)

Fullbaitsfinal <- rbind(Fullbaits,Offbaits2,
                        baitstot[nrow(baitstot),], baitstot[nrow(baitstot)-1,])
## NB : add the mito and api that I delete earlier for counting reasons

# Check for NAs?
sapply(Fullbaitsfinal, function(x) sum(is.na(x)))

## Finalize my GTF file format

head(Fullbaitsfinal)
Fullbaitsfinal$V10 <- paste0("\"",Fullbaitsfinal$V10,"\"")
Fullbaitsfinal$V13 <- paste0("\"",Fullbaitsfinal$V13,"\"")

Fullbaitsfinal$last <- with(Fullbaitsfinal, paste(V9,V10, sep=" "))
Fullbaitsfinal$last2 <- with(Fullbaitsfinal, paste(V11,V12,V13, sep=" "))
Fullbaitsfinal$LAST <- with(Fullbaitsfinal, paste(last,last2, sep=""))

Fullbaitsfinal <- Fullbaitsfinal[-c(9,10,11,12,13,14,15)]
head(Fullbaitsfinal)

## Export my GTF file
write.table(x=Fullbaitsfinal, file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf", sep="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)


#################################################
## PART II.  Same but with Off-target in chunks :
#################################################
## Complete the GTF
#################################
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library("Biostrings")
library(stringr)

## Read in the previously written new GTF
GTFfile <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget.gtf")

## Remove 2 last lines (apicoplast + mito), to add at the end!!
GTFfile <- GTFfile[-c(nrow(GTFfile),nrow(GTFfile)-1),]

## Reorder the file to get it contig after contig, bait by bait
### Add a column with the contig numbers (split V1)
GTFfile$contig <- (str_split_fixed(GTFfile$V1, "_", 2))[,2]

### and reorder by contig, and within contig
head(GTFfile)
GTFordered <- GTFfile[with(GTFfile, order(as.numeric(contig), V4)),]




##########
essai <- head(GTFordered, 50)

essai$category <- NA # Categories will be:
## "lonely", "<100", "100-300", "300-500", ">500"

## initialization:
if (str_split_fixed(essai$V13, "_", 3)[1,1]=="Off")  # Is it an "Off"?
    {
        if (essai[1,1]!=essai[2,1])
        {
            essai$category[1] <- "lonely"
        }
}    
## full loop:
for (i in 2:nrow(essai)) {
    if (str_split_fixed(essai$V13, "_", 3)[i,1]=="Off")  # Is it an "Off"?
    {
        if (essai[i,1]!=essai[i-1,1] & essai[i,1]!=essai[i+1,1]) # Is it a seq without baits?
        {
            essai$category[i] <- "lonely"
        }
        else
        {
            if (essai[i,1]==essai[i-1,1] & essai[i,1]==essai[i+1,1]) # Is it a seq between 2 baits?
            {
                essai$category[i] <- "between"
            }
            else
            {
                if (essai[i,1]==essai[i-1,1]) # Is it a seq AFTER a bait?
                {
                    essai$category[i] <- "after"
                }
                else
                {
                    if (essai[i,1]==essai[i+1,1]) # Is it a seq BEFORE a bait?
                    {
                        essai$category[i] <- "before"
                    }
                }
            }
        }
    }
    else
    {
        essai$category[i] <- "BAIT"
    }
}
            

## What is the distance to bait?

### Calculate the length of each sequence
essai$length <- essai$V5 - essai$V4 +1

### Then assess the distance to bait of each piece. Split if necessary.
essai$distance_to_bait <- NA
for (i in 1:nrow(essai)) {
    if (essai$category[i]=="lonely") # If the sequence is not on a contig with a bait, it's a "distant" sequence
    {
        essai$distance_to_bait[i] <- "distant"
    }
    if (essai$
}



tail(essai)
head(essai)
essai
