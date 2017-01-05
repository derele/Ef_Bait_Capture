## we actually use IRanges mostly but rtracklayer imports this and has
## the additional nice function import.gff
library(rtracklayer)
library(Biostrings)

## for a good intro see
## http://genomicsclass.github.io/book/pages/iranges_granges.html

gff <- import.gff("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_Alice.gtf")
## seqinfo: 242 sequences from an unspecified genome; no seqlengths
## you need to fix this to get the ranges up to the end of contigs
## best just readFastq and add the width (both bioc jargon ;-) to the
## gff
EfG <- readDNAStringSet("/SAN/Alices_sandpit/Efal_genome.fa")

seqlengths(seqinfo(gff)) <- width(EfG[seqnames(seqinfo(gff))])
isCircular(seqinfo(gff)) <- rep(FALSE, length(seqinfo(gff)))
genome(seqinfo(gff)) <- "E.falciformis"

add.gff.flanks <- function (gff, len, add){
    flanking.gff.S <- flank(gff, len, start=TRUE)
    flanking.gff.S <- trim(flanking.gff.S)
    flanking.gff.E <- flank(gff, len, start=FALSE)
    flanking.gff.E <- trim(flanking.gff.E)
    flanking.gff <- unique(c(flanking.gff.S, flanking.gff.E))
    ## avoid queries giving hits
    OL <- findOverlaps(flanking.gff, gff)
    elementMetadata(flanking.gff)$bait_id <-
                                    paste(elementMetadata(flanking.gff)$bait_id,
                                          add, sep="")
    c(gff, flanking.gff[-from(OL)])
}

len.list <- c(120, 200, 400)

gff.imputed <- gff
elementMetadata(gff.imputed)$bait_id <- paste(elementMetadata(gff.imputed)$bait_id,
                                              "none_0", sep="_")

for(i in seq_along(len.list)){
    gff.imputed <- add.gff.flanks(gff.imputed, len.list[i],
                                  paste(len.list[i], i, sep="_")) 
}

table(gsub(".*_(\\d)", "\\1", elementMetadata(gff.imputed)$bait_id))


export.gff(gff.imputed, "/SAN/Alices_sandpit/MYbaits_Eimeria_IRanges_IMP.gtf")