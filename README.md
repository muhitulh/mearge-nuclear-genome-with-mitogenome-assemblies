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

- To run the blast command on a list of FASTA files (assuming both nuclear and mito assembly in same folder)
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

