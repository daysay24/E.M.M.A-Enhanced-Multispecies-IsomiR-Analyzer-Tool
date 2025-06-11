miranda ... ... >"${f%.*}".out &"
# Filter Miranda Output: With my ~200 miRNAs and my 75k UTRs, I got ~1.5 million putative interactions. I filter based on the highest scores (sort), and arbitrarily cut off at highest scoring 50k interactions
cat *out |grep ">>" |sort -k5,5nr |awk 'NR<50000' |sed 's/>>//g' |cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" |tr "," "\t") - >MirandaOutput.tab






