## we actually use IRanges mostly but rtracklayer imports this and has
## the additional nice function import.gff
library(rtracklayer)


## for a good intro see  http://genomicsclass.github.io/book/pages/iranges_granges.html

gff <- import.gff("/SAN/Alices_sandpit/MYbaits_Eimeria_V1.single120.gff")
## seqinfo: 242 sequences from an unspecified genome; no seqlengths
## you need to fix this to get the ranges up to the end of contigs
## best just readFastq and add the width (both bioc jargon ;-) to the
## gff

## so now these three lines do the magic:
flanking.gff <- flank(gff, 120)
OL <- findOverlaps(flanking.gff, gff)
## avoid queries giving hits
gff.120 <- flanking.gff[-from(OL)]


## Now the programming: construct a function from this which takes a
## vector of distances from baits and constructs a list of gffs

## Then you can append someting like the distance to each name,
## flatten the list and write to a file

## Then find on or two of the very many other ways to do this with
## IRanges 
