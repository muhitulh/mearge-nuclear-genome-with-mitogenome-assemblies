
#!/bin/bash
for fasta_file in *_ncl_assembly.fasta; do samtools faidx "$fasta_file"; done
