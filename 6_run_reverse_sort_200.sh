for f in *_ncl_assembly_filtered.fasta; do
 prefix="${f%_ncl_assembly_filtered.fasta}"
 bioawk -c fastx -v prefix="$prefix" '{if(length($seq)>=200) print ">" prefix "_scaffold-" (++i)" "length($seq)"\n"$seq }' < "$f" > "${prefix}_ncl_assembly_filtered_sorted.fasta"
done
