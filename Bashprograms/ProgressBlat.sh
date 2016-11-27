#!/bin/bash

# The last database entries each mapping is reporting :
for file in $(ls /SAN/Alices_sandpit/sequencing_data_dereplicated/*.fasta.fa_dnax.psl); do echo $file; tail -n2 $file | cut -f 10; done | sed -n '/psl/{n;p}' > position

# and
IFS=$'\n'       # make newlines the only separator
set -f          # disable globbing
for i in $(cat $1); do
    cat /SAN/Alices_sandpit/sequencing_data_dereplicated/Efal_mtapi.fasta.fa | sed -e '/{$i}/,$d' | wc -c;
done
