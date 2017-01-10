## SNPs calling from chosen baits
## "/SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_All/BaitsChosen.txt"

## Extract the bam files ONLY for the chosen baits :

## Convert the GTF/GFF file to a BED format
cd /SAN/Alices_sandpit/

tail -6100 MYbaits_Eimeria_Chosen_Baits_6100.gtf 1 > temp.gtf
gff2bed < temp.gtf >  MYbaits_Eimeria_Chosen_Baits_6100.bed
rm temp.gtf

## Limit my alignments to the "chosen baits" region of the genome:
cd /SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Rsubread_align/

for file in *30.BAM; do samtools view -hL /SAN/Alices_sandpit/MYbaits_Eimeria_Chosen_Baits_6100.bed $file > $file.inregion.sam; done

## Convert sam to bam
for file in *inregion.sam; do samtools view -b -S -o $file.bam $file; done

## Sort and index
for file in *sam.bam; do samtools sort $file -o $file.sorted.bam; done

for file in *sorted.bam; do samtools index $file; done

## Identifying Genomic Variants
## Using samtools mpileup and then bcftools and generating a vcf file.
for file in *sorted.bam; do samtools mpileup -g -f /SAN/Alices_sandpit/Efal_genome.fa $file > $file.bcf; done

for file in *bcf; do bcftools call -c -v $file > $file.vcf; done


## First easy pre analyses
ls *.vcf | while read file
do
    echo $file  >> SNPcounting.txt
    grep -v "^#" $file | wc -l >> SNPcounting2.txt
done

paste SNPcounting.txt SNPcounting2.txt > SNPcountingtemp.txt
mv SNPcountingtemp.txt SNPcounting.txt
rm SNPcounting2.txt
