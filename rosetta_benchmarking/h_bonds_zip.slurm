#!/bin/bash
#BSUB -n 1
#BSUB -R "span[hosts=1]"
#BSUB -W 800
#BSUB -R "rusage[mem=100000]"
tar -czvf rosetta_system_dir_close_res.tar.gz system_dir_h_bonds/
