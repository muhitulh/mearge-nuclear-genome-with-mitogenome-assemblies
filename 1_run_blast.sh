#!/bin/bash

# Iterate through each file that matches the pattern *_ncl_assembly.fasta
for ncl_file in *_ncl_assembly.fasta; do
  # Extract the base name of the file without the _ncl_assembly.fasta extension
  db_name=$(basename "$ncl_file" _ncl_assembly.fasta)

  # Create a blast database using the ncl_file
  makeblastdb -in "$ncl_file" -out "$db_name" -dbtype 'nucl' -hash_index

  # Define the mito_file name as "${db_name}_mito_assembly.fasta"
  mito_file="${db_name}_mito_assembly.fasta"

  # Run blastn to perform a nucleotide-nucleotide search
  # Query: mito_file
  # Database: db_name
  # Output format: 6 (tabular)
  # Output file: "${db_name}_blast_hits.tsv"
  # E-value threshold: 1e-10
  # Number of threads: 16
  blastn -query "$mito_file" -task blastn -db "$db_name" -outfmt 6 -out "${db_name}_blast_hits.tsv" -evalue 1e-10 -num_threads 16

  # Sort the blast hits based on the 12th column in descending order
  # Save the sorted results to "${db_name}_blast_hits_sorted.tsv"
  sort -n -r -k 12 "${db_name}_blast_hits.tsv" > "${db_name}_blast_hits_sorted.tsv"
done
# here -n: this option tells sort to perform a numerical sort, -r: sort the input in reverse order, -k: column 12

mkdir database_file
mv *.nhd *.nhi *.nhr *.nin *.nog *.nsd *.nsi *.nhd *.nsq *_hits.tsv database_file/
