## we actually use IRanges mostly but rtracklayer imports this and has
## the additional nice function import.gff
library(rtracklayer)
library(Biostrings)

## for a good intro see
## http://genomicsclass.github.io/book/pages/iranges_granges.html

gff <- import.gff("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120_Alice_nuclgenonly.gtf")

## seqinfo: 244 sequences from an unspecified genome; no seqlengths
## you need to fix this to get the ranges up to the end of contigs
## best just readFastq and add the width (both bioc jargon ;-) to the
## gff
EfG <- readDNAStringSet("/SAN/Alices_sandpit/Efal_genome.fa")

seqlengths(seqinfo(gff)) <- width(EfG[seqnames(seqinfo(gff))])
isCircular(seqinfo(gff)) <- rep(FALSE, length(seqinfo(gff)))
genome(seqinfo(gff)) <- "E.falciformis"

## Add "none_0" to the baits themselves
elementMetadata(gff)$bait_id <- paste(elementMetadata(gff)$bait_id,
                                              "none_0", sep="_")

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

len.list <- c(120, 120, 120)

gff.imputed <- gff

## And a suffix for the flanking regions
for(i in seq_along(len.list)){
    gff.imputed <- add.gff.flanks(gff.imputed, len.list[i],
                                  paste(len.list[i], i, sep="_")) 
}

gff.imputed <- gff.imputed[!width(ranges(gff.imputed))==0,]

table(gsub(".*_(\\d)", "\\1", elementMetadata(gff.imputed)$bait_id))

## Add a term when there are flanks both sides
gff.imputed[duplicated(gff.imputed$bait_id),]$bait_id <- paste("other", gff.imputed[duplicated(gff.imputed$bait_id),]$bait_id, sep="_")

# export.gff(gff.imputed, "/SAN/Alices_sandpit/MYbaits_Eimeria_IRanges_IMP.gtf")

export.gff(gff.imputed, "/SAN/Alices_sandpit/MYbaits_Eimeria_IRanges_IMP_120.gtf")
