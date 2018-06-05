#!/bin/bash
cat << EOF > movie$1.tcl
set first $1
set last $2
EOF
cat movie.tcl >> movie$1.tcl
vmd-1.9.1 -dispdev text  -e movie$1.tcl
rm movie$1.tcl

