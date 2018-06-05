#!/bin/bash
vmd << EOF 

mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_init.pdb
mol addfile ../analysis_data/min_eq_nucl.dcd waitfor all

display projection orthographic
display rendermode GLSL

mol modstyle 0 top NewCartoon
mol modselect 0 top "nucleic"
mol modcolor 0 top ColorID 15

mol color ColorID 3
mol selection "segname CHA CHE"
mol representation NewCartoon
mol addrep top

mol color ColorID 7
mol selection "segname CHB CHF"
mol representation NewCartoon
mol addrep top

mol color ColorID 9
mol selection "segname CHC CHG"
mol representation NewCartoon
mol addrep top

mol color ColorID 10
mol selection "segname CHD CHH"
mol representation NewCartoon
mol addrep top

#mol addrep top
# set sel_all [atomselect top "all and noh"]	

# set frame0_helix [atomselect top "alpha_helix and noh" frame 0]	 
# set sel_helix [atomselect top "alpha_helix and noh"]	

# set frame0_helix_CA [atomselect top "alpha_helix and name CA" frame 0]	 
# set sel_helix_CA [atomselect top "alpha_helix and name CA"]	

# set frame0_helix_sidech [atomselect top "sidechain and alpha_helix and noh" frame 0]	 
# set sel_helix_sidech [atomselect top "sidechain and alpha_helix and noh"]	

# set frame0_dna [atomselect top "nucleic and noh" frame 0]	 
# set sel_dna [atomselect top "nucleic and noh"]	

# set frame0_dna_bb [atomselect top "backbone and nucleic and noh" frame 0]	 
# set sel_dna_bb [atomselect top "backbone and nucleic and noh"]

wait 360000000000000
EOF