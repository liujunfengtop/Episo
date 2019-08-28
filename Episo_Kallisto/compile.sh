#!/bin/bash

#create directoty
mkdir ../m5c_command

#compile C code
echo "Start compile ..."
date

gcc -o ../m5c_command/compare-paired compare-paired.c -lm
gcc -o ../m5c_command/anti_bisulfite anti_bisulfite.c
gcc -o ../m5c_command/anti_bisulfite_single_batch anti_bisulfite_single_batch.c -lm
gcc -o ../m5c_command/combination_diff combination_diff.c
gcc -o ../m5c_command/compute_m5c compute_m5c.c -lm
gcc -o ../m5c_command/sel_compare sel_compare.c
gcc -o ../m5c_command/m5c_filter m5c_filter.c
gcc -o ../m5c_command/selmethy selmethy.c -lm
gcc -o ../m5c_command/sum_counts sum_counts.c -lm
cp anti_bisulfite.ctl ../m5c_command
cp anti_bisulfite_single_batch.ctl ../m5c_command
cp m5c.sh ../m5c_command/
cp m5c_whole.sh ../m5c_command/
cp m5c_single.sh ../m5c_command/
cp kallisto ../m5c_command/
cp bismark_genome_preparation ../m5c_command/
cp bismark-liu ../m5c_command/
chmod u=rwx,g=rx,o=x ../m5c_command/m5c.sh
chmod u=rwx,g=rx,o=x ../m5c_command/m5c_whole.sh
chmod u=rwx,g=rx,o=x ../m5c_command/m5c_single.sh
chmod u=rwx,g=rx,o=x ../m5c_command/kallisto
chmod u=rwx,g=rx,o=x ../m5c_command/bismark_genome_preparation
chmod u=rwx,g=rx,o=x ../m5c_command/bismark-liu

echo "End compile."
date
