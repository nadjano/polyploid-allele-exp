grep -v '*' Orangutan_seperate_P.tsv | cut -f1 | sort | uniq | wc -l 


grep -v '*' Orangutan_competetive_P.tsv | cut -f1 | sort | uniq | wc -l 

grep -v '*' Orangutan_competetive_N200.tsv | cut -f1 | sort | uniq | wc -l 

grep -v '*' Orangutan_seperate_N200.tsv | cut -f1 | sort | uniq | wc -l 

grep -v '*' Atlantic_seperate_N200.tsv | cut -f1 | sort | uniq | wc -l 


cat  RIL_updated_P_diff.tsv  | cut -f1 | sort | uniq | wc -l