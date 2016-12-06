## Transform pls > bed > sam/bam:
## Usage : sh PsltoBamAlice.sh file.psl

## NB : translation of "reconstructed_apicoplast" into "EfaB_7" between bed and sam. Why? No clue. Be careful!!

echo "Welcome, I am ready to change your psl file to a sam/bam one"
FILE=$1
echo "Here is the head of your file.psl"
head $FILE

## Keep best hit in the .psl NB TIME CONSUMING!!
echo ""
echo Step 1/4 : keep the best hit with trinity...
perl /tools/trinityrnaseq-Trinity-v2.3.2/util/misc/blat_util/blat_top_hit_extractor.pl $FILE > "`basename "$FILE" .psl`_trinity.psl"
FILETRY="`basename "$FILE" .psl`_trinity.psl"

# pslToBed input.psl input.bed
echo ""
echo Step 2/4 : pslToBed...
/home/alice/pslToBed "$FILETRY" "`basename "$FILETRY" _trinity.psl`.bed"
FILEBED="`basename "$FILETRY" _trinity.psl`.bed"

## make a .genome file [DONE]
## cat /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Efal_mtapi.genome #ok

# bedtools bedtobam -bed12 -i input.bed -g hg38.chrom.sizes > output.bam
echo ""
echo Step 3/4 : bedtools bedtobam...
bedtools bedtobam -bed12 -i "$FILEBED" -g /SAN/Alices_sandpit/sequencing_data_dereplicated/All_alignments_Blat/Efal_mtapi.genome > "`basename "$FILEBED" .bed`.bam"
FILEBAM="`basename "$FILEBED" .bed`.bam"

# samtools view output.bam > output.sam
echo ""
echo Step 4/4 : bam to sam with samtools view...
samtools view "$FILEBAM" > "`basename "$FILEBAM" .bam`.sam"
