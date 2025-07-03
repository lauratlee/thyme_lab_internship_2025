# runs a blast query in the command line based off of zebrafish data
# example usage (to be run from /outputs/blast/): python ../../scripts/run_blast.py ../../gpcr_pocket_dir/AVPR2/Class_F/avpr2aa_pocket.fasta 2.0

import os, sys

query_file = sys.argv[1]
threshold = sys.argv[2]

# get zebrafish gene name
name = os.path.splitext(os.path.basename(query_file))[0]
gene = name.split('_')[0]

out_file = f"{gene}_{threshold}.txt"

os.system(f"blastp -query {query_file} -db {threshold}/zf_db_{threshold} -out {threshold}/{out_file} -outfmt 6 -max_target_seqs 10 -evalue 10")


# Example blast command that was manually input
# blastp -query ../gpcr_pocket_dir/AVPR2/Class_F/avpr2aa_pocket.fasta -db blast/zf_db_2.0 -out avpr2aa_2.0.txt -outfmt 6 -max_target_seqs 10 -evalue 10
