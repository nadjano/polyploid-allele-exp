
cd *tsv

grep -v '*' *_seperate_* | cut -f1 | sort | uniq | wc -l 
grep -v '*' *_competetive_* | cut -f1 | sort | uniq | wc -l 

cat  ${mm_parameters}_diff.tsv  | cut -f1 | sort | uniq | wc -l 

cat  *read_diff.tsv  | cut -f1 | sort | uniq | wc -l