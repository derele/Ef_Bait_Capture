#################################################
########## DE NOVO METAGENOME ASSEMBLY ##########
#################################################
## January 2017##
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(ggplot2)
library(reshape)
#################################################

## Import the blast results
blastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6")

tblastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf.b6")

## Import the libraries names
lib <- read.csv("/SAN/Alices_sandpit/sequencing_data_dereplicated/NumberOfSequences.csv",
                sep=" ", header=FALSE)
lib <- unique(substr(sub("_S.*.","",lib$V2),1,6))
lib <- lib[!lib%in%"Undete"] # remove "undetermine"

## Percentage similarities wanted:
X <- 80

## 1. With blastn
bl.names <- c("Query label","Target label","Percent identity","Alignment length",
              "Number of mismatches","Number of gap opens", "Start position in query",
              "End position in query","Start position in target",
              "End position in target","E-value","Bit score")

names(blastn) <- bl.names
names(tblastn) <- bl.names

## Select only the >X% similarities
blastn <- blastn[blastn$'Percent identity' >= X, ]
tblastn <- tblastn[tblastn$'Percent identity' >= X, ]

## Initialisation
cov.sums <- lapply(lib, function (x) {
    bn <- blastn[grep(x, blastn$'Target label'), ]
    bt <- tblastn[grep(x, tblastn$'Target label'), ]
    get.range.sum <- function (bl.tab){
        ## By library (to automatise) Create GRanges object :
        gr <- GRanges(seqnames = bl.tab[, 'Query label'],
                      ranges = IRanges(bl.tab[, 'Start position in query'],
                                       bl.tab[, 'End position in query']))
        ## Reduce align the ranges and merge overlapping ranges to produce
        ## a simplified set:
        red <- reduce(gr)
        ## Mesure the total length of the alignment:
        sum(width(red))
    }
    list(BN=get.range.sum(bn), BT=get.range.sum(bt))
})

names(cov.sums) <- lib

cov.sums <- melt(cov.sums)

ggplot(cov.sums, aes(x=L1, y=value, fill=L2))+
    geom_bar(stat="identity", position="dodge")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))+
    labs(title = "Length of the alignment showing >X% identity")


