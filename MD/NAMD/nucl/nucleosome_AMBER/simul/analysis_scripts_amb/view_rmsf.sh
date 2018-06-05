#!/bin/bash

cat << EOF > /tmp/pymol_scr.pml

load ../analysis_data/nucl_rmsf.pdb
spectrum b, blue_white_red, minimum=0.4, maximum=2
as cartoon
#cartoon putty

EOF

pymolm -u /tmp/pymol_scr.pml