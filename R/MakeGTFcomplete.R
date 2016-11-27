## Complete the GTF
#################################
library(Rsubread)
library(reshape)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library("Biostrings")
library(stringr)

## Make a GTF that contains genomeoffbaits + mito + api + baits :

## load the Genome + api + mito
fastaFile <- readDNAStringSet("/SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa")
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

## load the baits
baits <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_feature_counted.gtf")

## get the contigs names
Contigs <- as.character(df$seq_name)

## contigs NOT used for the baits design :
#ContigsOff <- setdiff(as.character(df$seq_name),levels(baits$V1))
### ???###

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

## !! Not api and mito, back later!!
baits2 <- baits[-c(nrow(baits),nrow(baits)-1),]

## extract the contig number in 1 column
baits2 <- cbind(baits2, str_split_fixed(baits2$V1, "_", 2))
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
                        baits[nrow(baits),], baits[nrow(baits)-1,])
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
