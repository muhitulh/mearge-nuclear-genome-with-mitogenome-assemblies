for f in *_ncl_assembly.fasta; do
    exclude="${f%_ncl_assembly.fasta}_exclude.txt"
    output="${f%_ncl_assembly.fasta}_ncl_assembly_filtered.fasta"
    awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' "$f" | grep -vFf "$exclude" - | tr "\t" "\n" > "$output"
done
