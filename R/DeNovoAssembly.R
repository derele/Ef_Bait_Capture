#################################################
########## DE NOVO METAGENOME ASSEMBLY ##########
#################################################
## January 2017##
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(ggplot2)
library(reshape)
library(rtracklayer)
library(Biostrings)
#################################################

## Import the libraries names
lib <- read.csv("/SAN/Alices_sandpit/sequencing_data_dereplicated/NumberOfSequences.csv",
                sep=" ", header=FALSE)
lib <- unique(substr(sub("_S.*.","",lib$V2),1,6))
lib <- lib[!lib%in%"Undete"] # remove "undetermine"

## Import the blast results
blastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6")
#tblastn <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf.b6")
bl.names <- c("Query label","Target label","Percent identity","Alignment length",
              "Number of mismatches","Number of gap opens", "Start position in query",
              "End position in query","Start position in target",
              "End position in target","E-value","Bit score")
names(blastn) <- bl.names
names(tblastn) <- bl.names

## Choose now the best hits (keep the first line when duplicates):
blastn <- blastn [!duplicated(B[,c("Query label","Target label")]),]

## Percentage similarities wanted:
X <- 80

## Select only the >X% similarities
blastn <- blastn[blastn$'Percent identity' >= X, ]
tblastn <- tblastn[tblastn$'Percent identity' >= X, ]

## Apply for all libraries
cov.sums <- lapply(lib, function (x) {
    bn <- blastn[grep(x, blastn$'Target label'), ]
    bt <- tblastn[grep(x, tblastn$'Target label'), ]
    get.range.sum <- function (bl.tab){
        ## By library, create GRanges object :
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

## Proportion of the "baits parts" covered:
# gtf <- import.gff("/SAN/Alices_sandpit/MYbaits_Eimeria_IRanges_IMP.gtf") pb to find?
gtf <- import.gff("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_Alice.gtf")

## Remove first apicoplast and mitochondria
gtf <- gtf[!grepl("mito|apico", gtf$bait_id),]

## Apply for all libraries
cov.baits <- lapply(lib, function (x) {
    bn <- blastn[grep(x, blastn$'Target label'), ]
    ## By library, create GRanges object :
    gr <- GRanges(seqnames = bn[, 'Query label'],
                  ranges = IRanges(bn[, 'Start position in query'], bn[, 'End position in query']),
                  strand = Rle(strand(rep("+",nrow(bn)))))
    ## Reduce align the ranges and merge overlapping ranges to produce
    ## a simplified set:
    red <- reduce(gr)
    ## Extracts the elements in the query that overlap at least one
    ## element in the subject between gtf and red
    ol <- intersect(gtf,red)
    ## Proportion of the alignment along the baits by the total length of the baits    
    b.b <- sum(width(reduce(ol)))/sum(width(gtf))*100
    list(b.b)
})

names(cov.baits) <- lib

cov.baits <- melt(cov.baits)

ggplot(cov.baits, aes(x=L1, y=value))+
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))



