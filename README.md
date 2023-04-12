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

- To run the blast command on a list of FASTA files (assuming both nuclear and mito assembly in same folder): '1_run_blast.sh'
```
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

- index spade assembly fasta file: ‘3_run_index_ncl_assembly.sh’
```
#!/bin/bash
for fasta_file in *_ncl_assembly.fasta; do samtools faidx "$fasta_file"; done
```


- '4_run_exclude_overlapped.py'

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

