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
## 11 = BLAST archive format (ASN.1)
tblastn -query /SAN/db/blastdb/Eimeria_falciformis/proteins.fa -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf -outfmt 11

## Try with BLAST 6 format
tblastn -query /SAN/db/blastdb/Eimeria_falciformis/proteins.fa -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf.b6 -outfmt 6

## Or transform it in outfmt 6 with blastformatter
blast_formatter -archive /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf -outfmt 6 -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf_transfo.b6

# With 80% similarities, how many lines?
## lib names stored in lib.txt

# grep for each library the number of hits
while read p; do
    awk '$3>80' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/tblastn/blast_19scaf_transfo.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result.txt

# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result2.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result.txt

# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result3.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/result3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_tblastn_Efalci.txt

# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_tblastn_Efalci.txt


###################
# step 7 : blastn with /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa as query and the metagenome assemblies as reference
blastn -query /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6 -outfmt 6

# With 80% similarities, how many lines?
## lib names stored in lib.txt

# grep for each library the number of hits
while read p; do
    awk '$3>80' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt

# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt

# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Efalci.txt

# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Efalci.txt


###################
# step 8 : blastn with mus musculus genome (Reference genome: Mus musculus (assembly GRCm38.p5)) as query and the metagenome assemblies as reference
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_genomic.fna.gz"
zcat GCF_000001635.25_GRCm38.p5_genomic.fna.gz > GCF_000001635.25_GRCm38.p5_genomic.fna

blastn -query /SAN/Alices_sandpit/GCF_000001635.25_GRCm38.p5_genomic.fna -db /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastdb_all19scaf -out /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf_musmusculus.b6 -outfmt 6

# Error: NCBI C++ Exception:
# T0 "/home/amu/src/ncbi-blast+/ncbi-blast+-git/c++/src/corelib/ncbiobj.cpp", line 977: Critical: ncbi::CObject::ThrowNullPointerException() - Attempt to access NULL pointer.
# Stack trace:
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CStackTraceImpl::CStackTraceImpl() offset=0x5B
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CStackTrace::CStackTrace(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&) offset=0x1F
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CException::x_GetStackTrace() offset=0x94
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CException::SetSeverity(ncbi::EDiagSev) offset=0x64
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CObject::ThrowNullPointerException() offset=0xAF
# /usr/lib/ncbi-blast+/libxblast.so ???:0 ncbi::blast::CBlastTracebackSearch::Run() offset=0xCA2
# /usr/lib/ncbi-blast+/libxblast.so ???:0 ncbi::blast::CLocalBlast::Run() offset=0x1422
# blastn ???:0 CBlastnApp::Run() offset=0x15FD
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CNcbiApplication::x_TryMain(ncbi::EAppDiagStream, char const*, int*, bool*) offset=0x133
# /usr/lib/ncbi-blast+/libxncbi.so ???:0 ncbi::CNcbiApplication::AppMain(int, char const* const*, char const* const*, ncbi::EAppDiagStream, char const*, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&) offset=0x6AD
# blastn ???:0 main offset=0x66
# /lib/x86_64-linux-gnu/libc.so.6 ???:0 __libc_start_main offset=0xF1
# blastn ???:0 _start offset=0x29

# With 80% similarities, how many lines?
## lib names stored in lib.txt

# grep for each library the number of hits
while read p; do
    awk '$3>80' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf_musmusculus.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt

# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt

# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu.txt

# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu.txt


# With 80% similarities, how many lines?
## lib names stored in lib.txt

# grep for each library the number of hits
while read p; do
    awk '$3>80' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf_musmusculus.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt

# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt

# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt

# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu.txt

# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu.txt

######## With 90% similarities?
# grep for each library the number of hits
while read p; do
    awk '$3>90' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf_musmusculus.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt
# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt
# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu_90.txt
# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu_90.txt



# grep for each library the number of hits
while read p; do
    awk '$3>90' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt
# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt
# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_90.txt
# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_90.txt


######## With 99% similarities?
# grep for each library the number of hits
while read p; do
    awk '$3>99' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf_musmusculus.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt
# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt
# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu_99.txt
# Add headers
sed -i 1i"lib match>80 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_Mmuscu_99.txt


# grep for each library the number of hits
while read p; do
    awk '$3>99' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/blastn/blast_19scaf.b6 | grep "$p" | wc -l
    done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt
# grep pour each library the number of sequences in the metagenome assembly
while read p; do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/SPAdes_metagenomic_assembly/All_scaffold_fasta_alignments/allscaf.fa | grep ">" | grep "$p" | wc -l
done < /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/lib.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp2.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt
# Divide
awk '$3 != 0 { print $2/$3*100 }' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt
# All in one
paste -d ' ' /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Temp.txt /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/temp3.txt > /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_99.txt
# Add headers
sed -i 1i"lib match>99 tot percentmatch" /SAN/Alices_sandpit/sequencing_data_dereplicated/De_Novo_Assembly/Blast/Result_blastn_99.txt
