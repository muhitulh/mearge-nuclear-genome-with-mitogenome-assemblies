for f in *_mito_assembly.fasta; do
    prefix="${f%_mito_assembly.fasta}"
    bioawk -c fastx -v prefix="$prefix" '{ print ">" prefix "_mito-" (++i)" "length($seq)"\n"$seq }' < "$f" > "${prefix}_mito_assembly_sorted.fasta"
done
