#!/usr/bin/env pymol
load /home/laura.lee5-umw/thyme_lab_internship_2025/tm-align_alignments/gpcr_pocket_dir/S1PR5/Class_F/s1pr5a.pdb_all, format=pdb
hide all
show stick
color blue, chain A
color red, chain B
set ray_shadow, 0
set stick_radius, 0.3
set sphere_scale, 0.25
show stick, not polymer
show sphere, not polymer
bg_color white
set transparency=0.2
zoom polymer

