#####################
##GENERAL bash code##
#####################
awk '{print $3}' MYbaits_Eimeria_V1.single120_feature_counted.gtf | head   #check a column

## Prepare GFF file
# add a counter to the end of the line of a gff file, and save it in a new file
awk '{if($1==c) {i++} else {i=1} ;c=$1; print $0 "_" i}' MYbaits_Eimeria_V1.single120.gff > MYbaits_Eimeria_V1.single120_feature_counted_alice.gff
most MYbaits_Eimeria_V1.single120_feature_counted.gff   # visualize
mv filetemp filefinal   #change name

# remove files by part of name:
find | grep prpm | xargs rm -f
#1- Find the files in the current directory beginning by Alignments.
find -name "Alignments.*"
#or
ls | grep "Alignments*"
#2- Put all the files with a commun name in a folder
mv -t ./FeatureCounts_results  `ls|grep "Alignments.*"` #good but also moved just "Alignment" file
#Copy file feature count
cp Alignments.Alignment_2672Si.featureCounts summaryFeatureCount_2672Si
#and modify it (keep only the lines with Assigned)
grep "Assigned" Alignments.Alignment_pattern.featureCounts >> summaryFeatureCount_pattern

#fastx_quality_stats -i 2672Single_S17_R1_001.fastq.gz -o 2672Single_S17_R1_001.fastq.gz_stats.txt
              $ fastq_quality_boxplot_graph.sh -i bc54_stats.txt -o bc54_quality.png -t "My Library"
                      $ fastx_nucleotide_distribution_graph.sh -i bc54_stats.txt -o bc54_nuc.png -t "My Library"

python --version # check the version of python

# Change names of extensions
for file in *.bam.sam; do mv "$file" "`basename "$file" .bam.sam`.sam"; done  

# Find rank of Efab_
grep ">" /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa | grep -n "25398";

# number of characters in fasta file :
wc -m < /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa

# Number of sequences par fasta file :
for i in *unique.fasta.fa; do wc -l $i >> NumberOfSequences; done #NB must be divided by 2!
awk '$1/=2' < NumberOfSequences > NumberOfSequences.csv

##########################
##Transform a GFF to GTF##
##########################
# Same, apart from the 9th field:
# Suppress 9th column !!!! TAB SEP FIELDS
awk -v OFS='\t' '{$9=""; print $0}' MYbaits_Eimeria_V1.single120_feature_counted.gff > temporaryAlice.gtf
#Remove the first 3 lines and change the file
awk '{if(NR>3)print $0 "gene_id " "\"" $1 "\"" "; " "transcript_id " "\"" $1 "\"" ";"}' temporaryAlice.gtf > MYbaits_Eimeria_V1.single120_alicechange.gtf
##########################

####################
## part 1 : ALIGN ##
####################

# step 1 : Dereplicate with Tally
tally -i out1.gz -j out2.gz -o out1.unique.gz -p out2.unique.gz --pair-by-offset
# Find the dereplicated files in the current directory
ls | grep "fastq.unique."
# Put all the files with a commun name in a folder
mv -t Dereplicated  `ls|grep "fastq.unique."`

## At one point : check quality more formaly

# step 2 : Align and Count

## --> FIRST ALIGNER TESTED : Rsubread:align (cf R file /home/alice/Ef_Bait_Capture/R/Align_and_Count.R)

## --> SECOND ALIGNER TESTED : Blat (with dna and dnax) & store in /Blat_Dna and /Blat_Dnax

cd /SAN/Alices_sandpit/sequencing_data_dereplicated/

time find What_I_want | parallel blat Efal_mtapi.fasta.fa {} -t=dna -q=dna -minIdentity=80 -dots=10000 All_alignments_Blat/Blat_Dna_Dnax/{}_blatDna.psl

#or

time find What_I_want | parallel blat Efal_mtapi.fasta.fa {} -t=dnax -q=dnax -minIdentity=80 -dots=10000 All_alignments_Blat/Blat_Dna_Dnax/{}_blatDnax.psl

########
## Transform pls > bed > sam/bam:
## Usage : sh /home/alice/Ef_Bait_Capture/Bashprograms/PsltoBamAlice.sh
########
for file in what_I_Want; do sh /home/alice/Ef_Bait_Capture/Bashprograms/PsltoBamAlice.sh $file; done

## example
for file in *psl; do sh /home/alice/Ef_Bait_Capture/Bashprograms/PsltoBamAlice.sh $file; done

## NB error in the sams files : EfaB_7 replace reconstructed_apicoplast!!
for i in  *.sam; do sed -i.bak 's/\<EfaB_7\>/reconstructed_apicoplast/g' $i; done

## --> ALIGNER TO TEST : Bowtie, Bwa at first
## Compare the TIME needed to align (i) and the RESULTS (ii)

####################
## part 2 : COUNT ##
####################

## NB : let's do it via R, and create a table directly from it... Here are some comparison

# Count the proportion of reads mapped to the reference genome
cat file.sam | cut –f3 | sort | uniq | wc -l > Nm # number of reads in the alignment
# Count sequences in the fasq.unique.gz files:
#for fname in *unique.gz; do zcat "$fname" | echo $((`wc -l`/4)); done
zcat file.unique.gz | echo $((`wc -l`/4)) > Nt
# Ratio
Nm/Nt

##########################################
## part 3 : de novo metagenome assembly ##
##########################################

# step 1 : create the assembly with SPAdes
/tools/SPAdes-3.7.1-Linux/bin/spades.py --meta -o Results_SPAdes_f9Ann/ -1 f9Anna_S15_R1_001.fastq.unique.gz -2 f9Anna_S15_R2_001.fastq.unique.gz 
## ===== Assembling finished.
##  * Corrected reads are in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/corrected/
##   * Assembled contigs are in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/contigs.fasta
##    * Assembled scaffolds are in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/scaffolds.fasta
##     * Assembly graph is in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/assembly_graph.fastg
##      * Paths in the assembly graph corresponding to the contigs are in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/contigs.paths
##       * Paths in the assembly graph corresponding to the scaffolds are in /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/scaffolds.paths
##======= SPAdes pipeline finished.

# step 2 : move all the scaffolds.fasta in one unique directory
scp De_Novo_Assembly/SPAdes_metagenomic_assembly/Results_SPAdes_f9Ann/scaffolds.fasta  De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/scaf_f9Ann.fa
 
# step 3 : add an identifier from the library to all lines with >
awk '/>/ {$0=$0 "_2672Si"}1' scaf_2672Si.fa > scaf_2672Si.fasta
#etc. NB for later : to automatise

# step 4 : create one big file with all the scaffolds
cat *fasta > allscaf.fa

# step 5 : create a blast database with the assembly
makeblastdb -in SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa -out Blast/blastdb_all19scaf -parse_seqids -dbtype 'nucl'

## Building a new DB, current time: 09/07/2016 15:33:53
## New DB name:   /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf
## New DB title:  SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa
## Sequence type: Nucleotide
## Keep MBits: T
## Maximum file size: 1000000000B
## Adding sequences from FASTA; added 18903594 sequences in 841.087 seconds.

###################
# step 6 : tblastn with Efalciprot as query and the metagenome assemblies as reference
tblastn -query /SAN/db/blastdb/Eimeria_falciformis/proteins.fa -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf.b6 -outfmt 6

###################
# step 7 : blastn with /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa as query and the metagenome assemblies as reference
blastn -query /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6 -outfmt 6
