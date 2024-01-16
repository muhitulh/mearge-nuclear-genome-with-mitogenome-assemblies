for f in *_ncl_assembly_filtered_sorted.fasta; do
  prefix="${f%_ncl_assembly_filtered_sorted.fasta}"
  mito_file="${f%_ncl_assembly_filtered_sorted.fasta}_mito_assembly_sorted.fasta"
  cat "$f" "$mito_file" > "${prefix}_final_merged_mito_ncl_assembly.fasta"
done

###create final folder and mv all in that folder
mkdir complete_assembly
mv *final_merged_mito_ncl_assembly.fasta ./complete_assembly

##remove intermadiate files
rm *_filtered.tsv
rm *_filtered.fasta
rm *_sorted.tsv
rm *_lengths.tsv
rm *_sorted.fasta
rm *_exclude.txt
rm *_assembly.fasta.fai



