## In bash

## The last database entries for each mapping is reported and stored in positions
positions=$(for file in $(ls /SAN/Alices_sandpit/sequencing_data_dereplicated/*.fasta.fa_dnax.psl); do echo $file; tail -n2 $file | cut -f 10; done | sed -n '/psl/{n;p}')

printf '%s\n' "${positions[@]}"

for pos in  ${positions[@]}; do echo $pos; done

for pos in  ${positions[@]}; do grep -aob $pos /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa ; done > /SAN/Alices_sandpit/sequencing_data_dereplicated/positions

## find /SAN/Alices_sandpit/sequencing_data_dereplicated/*unique.fasta.fa > librariesList

## for fname in *unique.gz; do zcat "$fname" | echo $((`wc -l`/4)); done > nbseqperlib


##then finish in R
############################################

## After extracting the position in the genome of the processed contigs, the libraries sizes,
## and the libraries list, calculate the progress of BLAT

positions <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/positions", sep=":")
libList <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/librariesList", sep=c("/","_"))
libList <- substr(libList$V5,1,6)
nbsequences <- read.table("/SAN/Alices_sandpit/sequencing_data_dereplicated/nbseqperlib")
table <- data.frame(libList,positions$V1,nbsequences$V1)
table$progress <- table$positions.V1/43698993
table$progress_reads <- table$progress * table$nbsequences.V1
## How many reads have been processed?
sum(table$progress_reads)
## How many reads to go?
sum(table$nbsequences.V1) - sum(table$progress_reads)
## What is the percentage of reads that have been processed?
PC <-round(sum(table$progress_reads)/sum(table$nbsequences.V1)*100,2)
PC

## Plot that:
library(ggplot2)
library(ggtheme)
table4plot <- table
table4plot
library(reshape)
table4plot <- melt(table4plot, id=c("libList","positions.V1","progress"))
ggplot(table4plot, aes(x=libList, y=value, fill=variable)) +
    geom_bar(position='dodge', stat='identity') +
    annotate("text", label=paste(PC),x = 10, y = 15000000, size = 8, colour = "red")+
    theme_minimal()

## Plot just the number of sequences per library
tablemini <- table[c(1,3)]
pdf("/home/alice/Figures/Numberreadsperlib.pdf")
ggplot(tablemini, aes(x=libList, y=nbsequences.V1)) +
    geom_bar(position='dodge', stat='identity', fill="black")
    theme_minimal()
dev.off()
