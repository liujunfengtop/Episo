#!/bin/bash

#parameter setting
samplistA="SRR493366 SRR493367 SRR493368"
idxname="mouse"
bs=100

#create directoty
mkdir ../anti_bisulfite
mkdir ../results
mkdir ../results_methy
mkdir ../m5c_results

#differential analysis
echo "Start ..."
echo "Start expression level estimating ..."
date

for name in ${samplistA}; do
  mv ../inputmapping/${name}"_pe.txt" ../inputmapping/demo_pe.txt
  ./anti_bisulfite anti_bisulfite.ctl
  ./selmethy anti_bisulfite_1.fq anti_bisulfite_2.fq
  mkdir ../anti_bisulfite/$name
  mv anti_bisulfite_1.fq ../anti_bisulfite/$name/
  mv anti_bisulfite_2.fq ../anti_bisulfite/$name/
  mv methy_1.fq ../anti_bisulfite/$name/
  mv methy_2.fq ../anti_bisulfite/$name/
  mv ../inputmapping/demo_pe.txt ../inputmapping/${name}"_pe.txt"
  rm out
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
mv m5c_out ../m5c_results/
mv summary_A.tsv ../m5c_results/
mv summary_methy_A.tsv ../m5c_results/

#delete directory
rm -rf ../anti_bisulfite
rm -rf ../results
rm -rf ../results_methy

echo "End."
