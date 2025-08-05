#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 800
#BSUB -R "rusage[mem=100000]"
#BSUB -e err_log.err
#BSUB -o rosetta.out
/pi/summer.thyme-umw/2024_intern_lab_space/rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @flags.txt #edit flags name
