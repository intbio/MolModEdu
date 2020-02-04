# split states from NMR ensemble and save it as the individual pdb
dir=$(pwd)
project='NpXynWT_apo_340K_MD_multi'
snapshots=${dir}/${project}.pdb
output=${dir}/${project}

mkdir ${output}

for pdb in ${snapshots}; do
 pdb_name2=$(basename "$pdb")
 pdb_name3="${pdb_name2/.pdb}"
 pdb_name="${pdb_name3/snapshots.}"
 grep -n 'MODEL\|ENDMDL' ${pdb} | cut -d: -f 1 | \awk '{if(NR%2) printf "sed -n %d,",$1+1; else printf "%dp '${pdb}' > '${output}'/'${pdb_name}'_%04d.pdb\n", $1-1,NR/2;}' |  bash -s
done
