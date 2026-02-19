#!/bin/bash
#/work/apps/Bandage # export all graph nodes to FASTA (positive only)
/data/dana/kraken2_to_bandage.sh all_positive_graph_nodes.fasta /data/scratch/refdbs/krakendb/standard_16gb > colors_standard.csv
/data/dana/kraken2_to_bandage.sh all_positive_graph_nodes.fasta /data/scratch/refdbs/krakendb/16S_SILVA138_k2db > colors_16S.csv
/data/dana/kraken2_to_bandage.sh all_positive_graph_nodes.fasta /data/scratch/refdbs/krakendb/pluspfp_08gb > colors_plus8.csv

