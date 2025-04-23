import pandas as pd

# get the coordinates of UTR coordinates of isomiRs 
awk '$3=="three_prime_UTR"' asu.gff3 > UTRs.gff3

# get the UTR sequence from coordinate
# bedtools doesn't use the Parent= name, so make it a tab output, paste with gff3 names, and convert it to fasta.  Make sure you set -s to get the strand information correct.  I also cut UTRs that were 10bp or shorter.  
bedtools getfasta -fi s.jap.fa -bed UTRs.gff3 -s -tab | paste <(cut -f 9 UTRs.gff3|sed 's/Parent=transcript://g' ) - |awk 'length($3)>10' |cut -f 1,3 |sed 's/\t/#/1' |sed 's/^/>/g' |tr "#" "\n " >NamedUTRs.fasta

# split into smaller files 
fasta-splitter --n-parts 36 mature.fasta

# Run miranda in a loop, otherwise run times are much longer
for f in *part*; 
do echo "miranda "$f" NamedUTRs.fasta >"${f%.*}".out &";
done >runMiranda.sh
sh runMiranda.sh &
# Filter Miranda Output: With my ~200 miRNAs and my 75k UTRs, I got ~1.5 million putative interactions. I filter based on the highest scores (sort), and arbitrarily cut off at highest scoring 50k interactions
cat *out |grep ">>" |sort -k5,5nr |awk 'NR<50000' |sed 's/>>//g' |cat <(echo "Seq1,Seq2,Tot Score,Tot Energy,Max Score,Max Energy,Strand,Len1,Len2,Positions" |tr "," "\t") - >MirandaOutput.tab






