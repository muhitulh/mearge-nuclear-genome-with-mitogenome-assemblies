# blast-and-mearge-nuclear-genome-with-mitogenome

- create environment and install packages

```
conda create -n merging
conda activate merging

conda install -c bioconda blast
conda install -c bioconda bedtools
conda install -c bioconda samtools
conda install -c bioconda bioawk
```

- To run the blast command on a list of FASTA files (assuming both nuclear and mito assembly in same folder): `1_run_blast.sh`
```
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





#!/bin/bash

for ncl_file in *_ncl_assembly.fasta; do
  db_name=$(basename "$ncl_file" _ncl_assembly.fasta)
  makeblastdb -in "$ncl_file" -out "$db_name" -dbtype 'nucl' -hash_index

  mito_file="${db_name}_mito_assembly.fasta"
  blastn -query "$mito_file" -task blastn -db "$db_name" -outfmt 6 -out "${db_name}_blast_hits.tsv" -evalue 1e-10 -num_threads 16

  sort -n -r -k 12 "${db_name}_blast_hits.tsv" > "${db_name}_blast_hits_sorted.tsv"
done

mkdir database_file
mv *.nhd *.nhi *.nhr *.nin *.nog *.nsd *.nsi *.nhd *.nsq *_hits.tsv database_file/
```
![image](https://github.com/muhitulh/blast-and-mearge-nuclear-genome-with-mitogenome/assets/67751990/e4ec3119-d270-40df-9e76-69c3a8cc52dc)

The column headers are as follows:

Column 1: qseqid, query sequence id. If you open you query file, you will find the id after ‘>’ sign

Column 2: sseqid, subject i.e., reference sequence id

Column 3: pident, percentage of identical matches

Column 4: length, alignment length

Column 5: mismatch, number of mismatches

Column 6: gapopen, number of gap openings

Column 7: qstart, start of alignment in query

Column 8: qend, end of alignment in query

Column 9: sstart, start of alignment in subject

Column 10: send, end of alignment in subject

Column 11: evalue, expect value

Column 12: bitscore, bit score



- calculate length all together: `2_run_calculate_length.py`
```
import glob

# loop through all filtered blast hits files
for blast_hits_file in glob.glob('*_blast_hits_filtered.tsv'):
    # create a dictionary to store start and end positions for each ID
    id_dict = {}

    # read the file and process each line
    with open(blast_hits_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            id = fields[1]
            start = int(fields[6])
            end = int(fields[7])

            # check if the ID is already in the dictionary
            if id in id_dict:
                overlap = False
                # check if the current interval overlaps with any existing intervals
                for interval in id_dict[id]:
                    if start <= interval[1] and end >= interval[0]:
                        # update the existing interval with the overlapped coordinates
                        interval[0] = min(start, interval[0])
                        interval[1] = max(end, interval[1])
                        overlap = True
                        break
                # if there is no overlap, add the new interval to the list of intervals for this ID
                if not overlap:
                    id_dict[id].append([start, end])
            else:
                # add the first interval for this ID
                id_dict[id] = [[start, end]]

    # calculate the total length of each ID and write to output file
    output_file = blast_hits_file.replace('_blast_hits_filtered.tsv', '_lengths.tsv')
    with open(output_file, 'w') as out:
        for id in id_dict:
            total_length = 0
            for interval in id_dict[id]:
                total_length += interval[1] - interval[0]
            out.write(f"{id}\t{total_length}\n")
```
![image](https://github.com/muhitulh/blast-and-mearge-nuclear-genome-with-mitogenome/assets/67751990/e8803772-7d76-49cd-a4e1-1e7d700c3c90)


- index spade assembly fasta file: `3_run_index_ncl_assembly.sh`
```
#!/bin/bash
for fasta_file in *_ncl_assembly.fasta; do samtools faidx "$fasta_file"; done
```
![image](https://github.com/muhitulh/blast-and-mearge-nuclear-genome-with-mitogenome/assets/67751990/64d41489-7e58-418e-96db-959ae026e65e)
here, col 1 = contig ID, col2 = length (bases long) of contigs, col3 = the third column indicates the first base of "contig1" starts in the FASTA file.

- `4_run_exclude_overlapped.py`

```
import glob

for lengths_file in glob.glob('*_lengths.tsv'):
    contigs_file = lengths_file.replace('_lengths.tsv', '_ncl_assembly.fasta.fai')
    exclude_file = lengths_file.replace('_lengths.tsv', '_exclude.txt')
    with open(lengths_file, 'r') as f1, open(contigs_file, 'r') as f2, open(exclude_file, 'w') as f3:
        lengths = {}
        for line in f2:
            contig_id, length = line.split()[0:2]
            lengths[contig_id] = int(length)
        for line in f1:
            contig_id, length = line.split()[0:2]
            length = int(length)
            coverage = length / lengths[contig_id] * 100
            if coverage > 90:  # Adjust the threshold as needed
                f3.write(contig_id + '\n')
```



- `5_run_extract_ncl_fasta_after_excluding.sh`
```
for f in *_ncl_assembly.fasta; do
    exclude="${f%_ncl_assembly.fasta}_exclude.txt"
    output="${f%_ncl_assembly.fasta}_ncl_assembly_filtered.fasta"
    awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$f" | grep -vFf "$exclude" - | tr "\t" "\n" > "$output"
done
```


- `6_run_reverse_sort_200.sh`
reverse size sort filt.fasta, renumber, append isolate ID to prefix e.g. SN15_scf001
```
for f in *_ncl_assembly_filtered.fasta; do
    prefix="${f%_ncl_assembly_filtered.fasta}"
    bioawk -c fastx -v prefix="$prefix" '{if(length($seq)>=200) print ">" prefix "_scaffold-" (++i)" "length($seq)"\n"$seq }' < "$f" > "${prefix}_ncl_assembly_filtered_sorted.fasta"
done
```
![image](https://github.com/muhitulh/blast-and-mearge-nuclear-genome-with-mitogenome/assets/67751990/6d5a0e81-6071-4ec4-b505-ef99807170ca)


- `7_run_reverse_mito_sort.sh`
- mitoz.fasta id will need changing too.. e.g. SN15_mtDNA_1/2
```
for f in *_mito_assembly.fasta; do
    prefix="${f%_mito_assembly.fasta}"
    bioawk -c fastx -v prefix="$prefix" '{ print ">" prefix "_mito-" (++i)" "length($seq)"\n"$seq }' < "$f" > "${prefix}_mito_assembly_sorted.fasta"
done
```

- `8_run_final_fasta.sh`

```
for f in *_ncl_assembly_filtered_sorted.fasta; do
  prefix="${f%_ncl_assembly_filtered_sorted.fasta}"
	mito_file="${f%_ncl_assembly_filtered_sorted.fasta}_mito_assembly_sorted.fasta"
cat "$f" "$mito_file" > "${prefix}_final_assembly.fasta"
done
```
![image](https://github.com/muhitulh/blast-and-mearge-nuclear-genome-with-mitogenome/assets/67751990/465c88b4-c7f4-4633-a6f5-5eb9493b263f)
