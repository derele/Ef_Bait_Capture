#####################
## Complete the GTF##
#####################
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
GTF <- GTFfile[-c(nrow(GTFfile),nrow(GTFfile)-1),]

## Reorder the file to get it contig after contig, bait by bait
### Add a column with the contig numbers (split V1)
GTF$contig <- (str_split_fixed(GTF$V1, "_", 2))[,2]

## Add a column with the subcontig sequence number
GTF$subcontig <- (str_split_fixed(GTF$V13, "_", 3))[,3]

### and reorder by contig, and within contig
GTFordered <- GTF[with(GTF, order(as.numeric(contig), V4)),]

########## LET'S GO!!!!
myDF <- GTFordered

## Add a column that will be filled with categories
myDF$category <- NA # Categories will be:
## "lonely", "<100", "100-300", "300-500", ">500"

### Calculate the length of each sequence
myDF$length <- myDF$V5 - myDF$V4 +1

## initialization:
myDF$category[1] <- "lonely"

## full loop:
for (i in 2:(nrow(myDF)-1)) {
    if (str_split_fixed(myDF$V13, "_", 3)[i,1]=="Off")  # Is it an "Off"?
    {if (myDF[i,1]!=myDF[i-1,1] & myDF[i,1]!=myDF[i+1,1]) # Is it a seq without baits?
     {myDF$category[i] <- "lonely"}
     else
        {if (myDF[i,1]==myDF[i-1,1] & myDF[i,1]==myDF[i+1,1]) # Is it a seq between 2 baits?
         {myDF$category[i] <- "between"}
         else
            {if (myDF[i,1]==myDF[i-1,1]) # Is it a seq AFTER a bait?
             {myDF$category[i] <- "after"}
             else
                {if (myDF[i,1]==myDF[i+1,1]) # Is it a seq BEFORE a bait?
                 {myDF$category[i] <- "before"}
                }
            }
        }
    }
    else
    {myDF$category[i] <- "BAIT"}
}

## Check last line, add manualy
myDF$category[nrow(myDF)] <- "lonely"

## Check lines with no distance
myDF[myDF$category %in% NA,]

## Save
write.table(x=myDF, file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_complete_TEMP.gtf")


myDF <- read.table(file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_complete_TEMP.gtf")

############ GOOD SO FAR

# table(myDF[myDF$V5-myDF$V4 <= 0,][17])
#### NB : PROBLEMS FOR SEQUENCES OF LENGTH :
# 201  202  203  601  602  603  604  605 1001 1002 1003 1004 1006 1007
# 15   18   11    3    2    7    5    6    4    7    5    3    4    5

### Then assess the distance to bait of each piece. Split if necessary.
myDF$distance_to_bait <- NA

for (i in 1:nrow(myDF)) {
    if (myDF$category[i]=="lonely") { # If the sequence is not on a contig with a bait, it's a "distant" sequence
        myDF$distance_to_bait[i] <- "distant"
    } else {
        if (myDF$category[i]=="BAIT") { # If the sequence is a bait, then... distance to bait is nulle
            myDF$distance_to_bait[i] <- "in_bait"
        } else {
            if (myDF$category[i]=="between") { # if the sequence is between 2 baits
                if (myDF$length[i]<=203) { # If the sequence is smaller than 203bp
                    myDF$distance_to_bait[i]<-"0_to_100" # Any base is distant to max 100bp to a bait
                } else {
                    if (myDF$length[i]<=605) { # If the sequence is smaller than 605bp, we cut it
                        j <- nrow(myDF)+1 # We fill at the end of the file
                        myDF[j,] <- myDF[i,] # New lines are created to split the sequence
                        myDF[j+1,] <- myDF[i,]
                        myDF[j+2,] <- myDF[i,]
                        myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                        myDF[j,18] <- "0_to_100"
                        myDF$subcontig[j] <- 1
                        myDF[j+2,4] <- myDF[j+2,5]-100 # Last piece is 100 long, close to 1 bait
                        myDF[j+2,18] <- "0_to_100"
                        myDF$subcontig[j+2] <- 3
                        myDF[j+1,4] <- myDF[j,5]+1 # Last piece is further from bait (100<d<300)
                        myDF[j+1,5] <- myDF[j+2,4]-1
                        myDF[j+1,18] <- "100_to_300"
                        myDF$subcontig[j+1] <- 2
                        myDF[i,18] <- "To remove"
                    } else {
                        if (myDF$length[i]<=1007) { # If the sequence is smaller than 1007bp, we cut it
                            j <- nrow(myDF)+1 # We fill at the end of the file
                            myDF[j,] <- myDF[i,] # Init.
                            myDF[j+1,] <- myDF[i,]
                            myDF[j+2,] <- myDF[i,]
                            myDF[j+3,] <- myDF[i,]
                            myDF[j+4,] <- myDF[i,]
                            myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                            myDF[j,18] <- "0_to_100"
                            myDF$subcontig[j] <- 1
                            myDF[j+4,4] <- myDF[j+4,5]-100 # Last piece is 100 long, close to 1 bait
                            myDF[j+4,18] <- "0_to_100"
                            myDF$subcontig[j+4] <- 5
                            myDF[j+1,4] <- myDF[j,5]+1
                            myDF[j+1,5] <- myDF[j+1,4]+200 # 2nd piece is 200 long, 100to300 from bait
                            myDF[j+1,18] <- "100_to_300"
                            myDF$subcontig[j+1] <- 2
                            myDF[j+3,5] <- myDF[j+4,4]-1
                            myDF[j+3,4] <- myDF[j+3,5]-200 # 4nd piece is 200 long, 100to300 from bait
                            myDF[j+3,18] <- "100_to_300"
                            myDF$subcontig[j+3] <- 4
                            myDF[j+2,4] <- myDF[j+1,5]+1
                            myDF[j+2,5] <- myDF[j+3,4]-1 # Middle piece is further away [300 to 500bp]
                            myDF[j+2,18] <- "300_to_500"
                            myDF$subcontig[j+2] <- 3
                            myDF[i,18] <- "To remove"
                        } else {
                            j <- nrow(myDF)+1 # We fill at the end of the file
                            myDF[j,] <- myDF[i,] # Init.
                            myDF[j+1,] <- myDF[i,]
                            myDF[j+2,] <- myDF[i,]
                            myDF[j+3,] <- myDF[i,]
                            myDF[j+4,] <- myDF[i,]
                            myDF[j+5,] <- myDF[i,]
                            myDF[j+6,] <- myDF[i,]
                            myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                            myDF[j,18] <- "0_to_100"
                            myDF$subcontig[j] <- 1
                            myDF[j+6,4] <- myDF[j+6,5]-100 # Last piece is 100 long, close to 1 bait
                            myDF[j+6,18] <- "0_to_100"
                            myDF$subcontig[j+6] <- 7
                            myDF[j+1,4] <- myDF[j,5]+1
                            myDF[j+1,5] <- myDF[j+1,4]+200 # 2nd piece is 200 long, 100to300 from bait
                            myDF[j+1,18] <- "100_to_300"
                            myDF$subcontig[j+1] <- 2
                            myDF[j+5,5] <- myDF[j+6,4]-1
                            myDF[j+5,4] <- myDF[j+5,5]-200 # 6st piece is 200 long, 100to300 from bait
                            myDF[j+5,18] <- "100_to_300"
                            myDF$subcontig[j+5] <- 6
                            myDF[j+2,4] <- myDF[j+1,5]+1
                            myDF[j+2,5] <- myDF[j+2,4]+200 # 3rd piece is 200 long, 300to500 from bait
                            myDF[j+2,18] <- "300_to_500"
                            myDF$subcontig[j+2] <- 3
                            myDF[j+4,5] <- myDF[j+5,4]-1
                            myDF[j+4,4] <- myDF[j+4,5]-200 # 5th piece is 200 long, 300to500 from bait
                            myDF[j+4,18] <- "300_to_500"
                            myDF$subcontig[j+4] <-5 
                            myDF[j+3,4] <- myDF[j+2,5]+1 # Middle piece is further away [more than 500]
                            myDF[j+3,5] <- myDF[j+4,4]-1
                            myDF[j+3,18] <- "more_than_500"
                            myDF$subcontig[j+3] <- 4
                            myDF[i,18] <- "To remove"
                        }
                    }
                }
            } else {
                if (myDF$category[i]=="before") { # if the sequence is before a bait
                    if (myDF$length[i]<=100) { # If the sequence is smaller than 100bp
                        myDF$distance_to_bait[i]<-"0_to_100" # Any base is distant to max 100bp to a bait
                    } else {
                        if (myDF$length[i]<=300) { # If the sequence is smaller than 300bp, we cut it
                            j <- nrow(myDF)+1 # We fill at the end of the file
                            myDF[j,] <- myDF[i,] # New lines are created to split the sequence
                            myDF[j+1,] <- myDF[i,]
                            myDF[j,4] <- myDF[j,5]-100 # First piece is 100 long, close to 1 bait
                            myDF[j,18] <- "0_to_100"
                            myDF$subcontig[j] <- 2
                            myDF[j+1,5] <- myDF[j,4]-1 # Last piece is further from bait (100<d<300)
                            myDF[j+1,18] <- "100_to_300"
                            myDF$subcontig[j+1] <- 1
                            myDF[i,18] <- "To remove"
                        } else {
                            if (myDF$length[i]<=500) { # If the sequence is smaller than 500bp, we cut it
                                j <- nrow(myDF)+1 # We fill at the end of the file
                                myDF[j,] <- myDF[i,] # Init.
                                myDF[j+1,] <- myDF[i,]
                                myDF[j+2,] <- myDF[i,]
                                myDF[j,4] <- myDF[j,5]-100 # First piece is 100 long, close to 1 bait
                                myDF[j,18] <- "0_to_100"
                                myDF$subcontig[j] <- 3
                                myDF[j+1,5] <- myDF[j,4]-1 # 2nd piece is 200 long, 100to300 from bait
                                myDF[j+1,4] <- myDF[j+1,5]-200
                                myDF[j+1,18] <- "100_to_300"
                                myDF$subcontig[j+1] <- 2
                                myDF[j+2,5] <- myDF[j+1,4]-1 # 3rd piece is 200 long, further away [300 to 500bp]
                                myDF[j+2,18] <- "300_to_500"
                                myDF$subcontig[j+2] <- 1
                                myDF[i,18] <- "To remove"
                            } else {
                                j <- nrow(myDF)+1 # We fill at the end of the file
                                myDF[j,] <- myDF[i,] # Init.
                                myDF[j+1,] <- myDF[i,]
                                myDF[j+2,] <- myDF[i,]
                                myDF[j+3,] <- myDF[i,]
                                myDF[j,4] <- myDF[j,5]-100 # First piece is 100 long, close to 1 bait
                                myDF[j,18] <- "0_to_100"
                                myDF$subcontig[j] <- 4
                                myDF[j+1,5] <- myDF[j,4]-1 # 2nd piece is 200 long, 100to300 from bait
                                myDF[j+1,4] <- myDF[j+1,5]-200
                                myDF[j+1,18] <- "100_to_300"
                                myDF$subcontig[j+1] <- 3
                                myDF[j+2,5] <- myDF[j+1,4]-1 # 3rd piece is 200 long, further away [300 to 500bp]
                                myDF[j+2,4] <- myDF[j+2,5]-200
                                myDF[j+2,18] <- "300_to_500"
                                myDF$subcontig[j+2] <- 2
                                myDF[j+3,5] <- myDF[j+2,4]-1 # 4th piece is further than 500bp from bait
                                myDF[j+3,18] <- "more_than_500"
                                myDF$subcontig[j+3] <- 1
                                myDF[i,18] <- "To remove"
                            }
                        }
                    }
                } else {
                    if (myDF$category[i]=="after") { #if the sequence is after a bait
                        if (myDF$length[i]<=100) { # If the sequence is smaller than 100bp
                            myDF$distance_to_bait[i]<-"0_to_100" # Any base is distant to max 100bp to a bait
                        } else {
                            if (myDF$length[i]<=300) { # If the sequence is smaller than 300bp, we cut it
                                j <- nrow(myDF)+1 # We fill at the end of the file
                                myDF[j,] <- myDF[i,] # New lines are created to split the sequence
                                myDF[j+1,] <- myDF[i,]
                                myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                                myDF[j,18] <- "0_to_100"
                                myDF$subcontig[j] <- 1
                                myDF[j+1,4] <- myDF[j,5]+1 # Last piece is further from bait (100<d<300)
                                myDF[j+1,18] <- "100_to_300"
                                myDF$subcontig[j+1] <- 2
                                myDF[i,18] <- "To remove"
                            } else {
                                if (myDF$length[i]<=500) { # If the sequence is smaller than 500bp, we cut it
                                    j <- nrow(myDF)+1 # We fill at the end of the file
                                    myDF[j,] <- myDF[i,] # New lines are created to split the sequence
                                    myDF[j+1,] <- myDF[i,]
                                    myDF[j+2,] <- myDF[i,]
                                    myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                                    myDF[j,18] <- "0_to_100"
                                    myDF$subcontig[j] <- 1
                                    myDF[j+1,4] <- myDF[j,5]+1 # 2nd piece is further from bait (100<d<300)
                                    myDF[j+1,5] <- myDF[j+1,4]+200
                                    myDF[j+1,18] <- "100_to_300"
                                    myDF$subcontig[j+1] <- 2
                                    myDF[j+2,4] <- myDF[j+1,5]+1 #3rd piece is between 300 and 500 from bait
                                    myDF[j+2,18] <- "300_to_500"
                                    myDF$subcontig[j+2] <- 3
                                    myDF[i,18] <- "To remove"
                                } else {
                                    j <- nrow(myDF)+1 # We fill at the end of the file
                                    myDF[j,] <- myDF[i,] # New lines are created to split the sequence
                                    myDF[j+1,] <- myDF[i,]
                                    myDF[j+2,] <- myDF[i,]
                                    myDF[j+3,] <- myDF[i,]
                                    myDF[j,5] <- myDF[j,4]+100 # First piece is 100 long, close to 1 bait
                                    myDF[j,18] <- "0_to_100"
                                    myDF$subcontig[j] <- 1
                                    myDF[j+1,4] <- myDF[j,5]+1 # 2nd piece is further from bait (100<d<300)
                                    myDF[j+1,5] <- myDF[j+1,4]+200
                                    myDF[j+1,18] <- "100_to_300"
                                    myDF$subcontig[j+1] <- 2
                                    myDF[j+2,4] <- myDF[j+1,5]+1 #3rd piece is between 300 and 500 from bait
                                    myDF[j+2,5] <- myDF[j+2,4]+200                                     
                                    myDF[j+2,18] <- "300_to_500"
                                    myDF$subcontig[j+2] <- 3
                                    myDF[j+3,4] <- myDF[j+2,4]+1 # 4th piece is further away (more than 500)
                                    myDF[j+3,18] <- "more_than_500"
                                    myDF$subcontig[j+3] <- 4
                                    myDF[i,18] <- "To remove"
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

### CHECK
table(myDF[myDF$V5-myDF$V4 <= 0,][17])

# Save
write.table(x=myDF, file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_complete_TEMP2.gtf")

## Remove the "to remove"
myDF <- myDF[!grepl("remove", myDF$distance_to_bait),]

## Pass the "distance to bait" in name of bait
myDF$V13 <- paste(myDF$V13, myDF$distance_to_bait, sep="_")

## Re-order by contig, and within contig
myDF <- myDF[with(myDF, order(as.numeric(contig), V4)),]

## Erase bad columns
myDFfinal <- myDF[-c(14:18)]

## Add back the last 2 lines (apicoplast + mitocondria)

## Prepare my GTF file for export
myDFfinal <- rbind(myDFfinal, GTFfile[c(nrow(GTFfile),nrow(GTFfile)-1),])

##myDFfinal <- read.table("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_COMPLETE.gtf", sep="\t")

# Check for NAs?
sapply(myDFfinal, function(x) sum(is.na(x)))

## Finalize my GTF file format

head(myDFfinal)

## Troubleshooting : check when V5 - V4 < 100 (should be never)

myDFfinal[(myDFfinal$V5 - myDFfinal$V4) < 0,]

diff <- myDFfinal$V5 - myDFfinal$V4
table(diff[diff<0])

head(myDFfinal)
myDFfinal$V10 <- paste0("\"",myDFfinal$V10,"\"")
myDFfinal$V13 <- paste0("\"",myDFfinal$V13,"\"")

myDFfinal$last <- with(myDFfinal, paste(V9,V10, sep=" "))
myDFfinal$last2 <- with(myDFfinal, paste(V11,V12,V13, sep=" "))
myDFfinal$LAST <- with(myDFfinal, paste(last,last2, sep=""))

myDFfinal <- myDFfinal[-c(9,10,11,12,13,14,15)]

myDFfinal$LAST <- gsub("Off_target", "OT", myDFfinal$LAST)

## Export
write.table(x=myDFfinal, file="/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_AND_offtarget_complete.gtf", sep="\t", col.names=FALSE, row.names=FALSE, quote = FALSE)
