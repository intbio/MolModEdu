#!/bin/bash

cat << EOF > /tmp/pymol_scr.pml

load ../analysis_data/nucl_bfactors.pdb
spectrum b, blue_white_red, minimum=1, maximum=5
as cartoon
#cartoon putty

EOF

pymolm -u /tmp/pymol_scr.pml