#!/bin/bash
echo  "python run_evaluation_on_assembly.py \
--in assembly.fasta \
--cluster_dir output/ \
--inactive" | \
qsub \
-pe multislot 8 \
-cwd \
-P fair_share \
-N ass_eval \
-l vf=5G \
-l idle=1 \
-l arch=lx-amd64 \
-o genome_eval.out \
-e genome_eval.err
