# Finding similar pairs of sequences

# Write a shell script that uses BLAST to identify pairs of sequences in a protein FASTA file with high similarity. You should be able to supply the name of the FASTA file and$

# 1.  turn the FASTA file into a BLAST database using makeblastdb
FILENAME=$1
CUTOFF=$2
echo ""
echo blastp is running, please be patient...

makeblastdb -in $FILENAME -out mydatabase -dbtype 'prot'

#2. do an all-vs-all BLAST search by running blastp with the database you have built and the FASTA file as the query
# get blastp to produce the output in tabular format (outfmt 6)
blastp -query $FILENAME -db mydatabase -outfmt 6 | column -t > $FILENAME.out
echo ""

echo The sequences resulting of blastp are stored in $FILENAME.out
echo There are $(cat $FILENAME.out | wc -l) pairs of sequences there
echo ""

#3. use the tools from the previous session to identify rows where the identity score is greater than the cutoff
#  extract just the sequence names from those rows
awk -v cutoff=$CUTOFF '$3 >= cutoff' $FILENAME.out | tr -s " " | cut -d " " -f 1,2 > $FILENAME.out2

echo The sequences resulting of blastp with an idscore greater than $CUTOFF are stored in $FILENAME.out2
echo There are $(cat $FILENAME.out2 | wc -l) pairs of sequences there
echo ""

# Figure out all the steps on the command line first before you start to put them into a shell script. You'll probably have to consult the help for makeblastdb and blastp.

# Exclude the sequences aligned to themselves
awk '$2 != $1 { print $0 }' $FILENAME.out2 > $FILENAME.outfinal

echo The final file, excluding the sequences matching to themselves, is stored in $FILENAME.outfinal
echo There are $(cat $FILENAME.outfinal | wc -l) pairs of sequences there





