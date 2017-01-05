#################################################
########## DE NOVO METAGENOME ASSEMBLY ##########
#################################################
## January 2017##
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(ggplot2)
#################################################

## Import the blast results
T_blastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6")

T_tblastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf.b6")

## Import the libraries names
lib <- read.csv("/SAN/Alices_sandpit/sequencing_data_dereplicated/NumberOfSequences.csv", sep=" ", header=FALSE)
lib <- unique(substr(sub("_S.*.","",lib$V2),1,6))
lib <- lib[-19] # remove "undetermine"

## Percentage similarities wanted:
X <- 80

## 1. With blastn
names(T_blastn) <- c("Query label","Target label","Percent identity","Alignment length",
                    "Number of mismatches","Number of gap opens", "Start position in query",
                    "End position in query","Start position in target",
                    "End position in target","E-value","Bit score")

## Select only the >X% similarities
T_blastn <- T_blastn[T_blastn$'Percent identity' >= X, ]

## Initialisation
All <- vector()
## Loop over libraries

## loop over a list lapply lib


for (i in lib){
    ## By library (to automatise)
    T <- T_blastn[grep(i, T_blastn$'Target label'), ]
    ## Create GRanges object :
    gr <- GRanges(seqnames = T$'Query label',
                  ranges = IRanges(T$'Start position in query',
                                   T$'End position in query'))
    ## Reduce align the ranges and merge overlapping ranges to produce a simplified set:
    red <- reduce(gr)
    ## Mesure the total length of the alignment:
    len <- sum(width(red))
    All <- c(All,len)
}
DF <- data.frame(lib,All)
DF$blast <- "blastn"

## 2. With tblastn
names(T_tblastn) <- c("Query label","Target label","Percent identity","Alignment length",
                    "Number of mismatches","Number of gap opens", "Start position in query",
                    "End position in query","Start position in target",
                    "End position in target","E-value","Bit score")

## Select only the >X% similarities
T_tblastn <- T_tblastn[T_tblastn$'Percent identity' >= X, ]

## Initialisation
All <- vector()
## Loop over libraries
for (i in lib){
    ## By library (to automatise)
    T <- T_tblastn[grep(i, T_tblastn$'Target label'), ]
    ## Create GRanges object :
    gr <- GRanges(seqnames = T$'Query label',
                  ranges = IRanges(T$'Start position in query',
                                   T$'End position in query'))
    ## Reduce align the ranges and merge overlapping ranges to produce a simplified set:
    red <- reduce(gr)
    ## Mesure the total length of the alignment:
    len <- sum(width(red))
    All <- c(All,len)
}
DF2 <- data.frame(lib,All)
DF2$blast <- "tblastn"


## Plot
DF3 <- rbind(DF,DF2) 

ggplot(DF3, aes(x=lib, y=All, fill=blast))+
    geom_bar(stat="identity", position="dodge")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))+
    labs(title = "Length of the alignment showing >X% identity")


DF2
