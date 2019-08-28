#!/bin/bash

#site number
j=1

#parameter setting
samplistA="SRR493366 SRR493367 SRR493368"
idxname="mouse"
bs=100
total=100
lg=115

#create directoty
mkdir ../anti_bisulfite
mkdir ../results
mkdir ../results_methy
mkdir ../m5c_results

#edit the control file anti_bisulfite_single_batch.ctl
sed '/length/s/lvalue/'"$lg"'/g' anti_bisulfite_single_batch.ctl > anti_bisulfite_single_batch_ls.ctl

#site differential analysis
echo "Start ..."
echo "Start sort file ..."
date

for name in ${samplistA}; do
  sed '1d' ../inputmapping/${name}"_pe.txt" | sort -k 3 > ../inputmapping/${name}"_pe_sort.txt"
done

echo "End sort file."
date

while [ $j -le $total ]; do
  echo "Start the $j site ..."
  echo "Start expression level estimating ..."
  date
  
  for name in ${samplistA}; do
    mv ../inputmapping/${name}"_pe_sort.txt" ../inputmapping/demo_pe.txt
    ./anti_bisulfite_single_batch anti_bisulfite_single_batch_ls.ctl ../inputmapping/site_info.txt $j
    ./selmethy anti_bisulfite_1.fq anti_bisulfite_2.fq
    mkdir ../anti_bisulfite/$name
    mv anti_bisulfite_1.fq ../anti_bisulfite/$name/
    mv anti_bisulfite_2.fq ../anti_bisulfite/$name/
    mv methy_1.fq ../anti_bisulfite/$name/
    mv methy_2.fq ../anti_bisulfite/$name/
    mv ../inputmapping/demo_pe.txt ../inputmapping/${name}"_pe_sort.txt"
    rm ../inputmapping/trans_anno_out
    rm out
    rm methylation_summary_sam
    mkdir ../results/$name
    ./kallisto quant -i ../kallisto_index/${idxname}"_transcripts.idx" -o ../results/$name/ -b $bs ../anti_bisulfite/$name/anti_bisulfite_1.fq ../anti_bisulfite/$name/anti_bisulfite_2.fq
    ./kallisto h5dump -o ../results/$name/bootstrap/ ../results/$name/abundance.h5
    mkdir ../results_methy/$name
    ./kallisto quant -i ../kallisto_index/${idxname}"_transcripts.idx" -o ../results_methy/$name/ -b $bs ../anti_bisulfite/$name/methy_1.fq ../anti_bisulfite/$name/methy_2.fq
    ./kallisto h5dump -o ../results_methy/$name/bootstrap/ ../results_methy/$name/abundance.h5    
  done
  
  echo "End expression level estimating."
  date
  
  ./m5c.sh
  ./m5c_filter m5c_A.tsv
  mv filter_out m5c_out
  rm m5c_A.tsv
  ./combination_diff m5c_out ../inputmapping/site_info.txt $j
  mv diff_out_single.tsv m5c_out_single
  ./combination_diff summary_A.tsv ../inputmapping/site_info.txt $j
  mv diff_out_single.tsv summary_A_single.tsv
  ./combination_diff summary_methy_A.tsv ../inputmapping/site_info.txt $j
  mv diff_out_single.tsv summary_methy_A_single.tsv
  if [[ $j == 1 ]]; then
    cat m5c_out_single >> ../m5c_results/m5c_out_single_all
    cat summary_A_single.tsv >> ../m5c_results/summary_A_single_all.tsv
    cat summary_methy_A_single.tsv >> ../m5c_results/summary_methy_A_single_all.tsv
  else
    sed '1d' m5c_out_single >> ../m5c_results/m5c_out_single_all
    sed '1d' summary_A_single.tsv >> ../m5c_results/summary_A_single_all.tsv
    sed '1d' summary_methy_A_single.tsv >> ../m5c_results/summary_methy_A_single_all.tsv
  fi
  rm m5c_out_single
  rm summary_A_single.tsv
  rm summary_methy_A_single.tsv
  rm m5c_out
  rm summary_A.tsv
  rm summary_methy_A.tsv
  
  echo "Start deleting ..."
  date
  
  for name in ${samplistA}; do
    rm -rf ../anti_bisulfite/$name
    rm -rf ../results/$name
    rm -rf ../results_methy/$name
  done
 
  echo "End deleting."
  echo "End the $j site."
  date  
  
  let j=j+1
done

#delete directory and file
rm -rf ../anti_bisulfite
rm -rf ../results
rm -rf ../results_methy
rm anti_bisulfite_single_batch_ls.ctl
for name in ${samplistA}; do
  rm ../inputmapping/${name}"_pe_sort.txt"
done

echo "End."