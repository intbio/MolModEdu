#!/bin/bash
#
# this script downloads pdb defined in pdb_array in the current folder using wget

home=$(pwd)
mkdir ${home}/pdb_xyl

#pdb_array=('2c1f' '2vg9' '4ixl' '4ixl' '2f6b' '2f6b' '2dck' '2nqy' '1qh7' '1h4g' '2nqy' '1qh7' '1qh6' '2dcj' '2dcj' '1qh6' '1h4g' '1h4h' '1igo' '1igo' '1h4h' '1h4h' '1h4h' '1f5j' '1f5j' '3wp5' '3wp4' '3wp3' '3wp6' '3mf6' '2dfc' '4s2d' '4s2h' '3mf9' '3wp3' '4s2g' '5hxv' '5hxv' '4s2f' '5hxv' '3mfa' '5hxv' '1pvx' '1m4w' '1ref' '2vuj' '1red' '3b5l' '5ej3' '5hxv' '3aks' '1yna' '5ej3' '5hxv' '3akr' '1xyp' '1enx' '1xyp' '3akp' '3akt' '3akp' '1ree' '5hxv' '1xyo' '1enx' '3lgr' '3akt' '1xnk' '2vgd' '5hxv' '2dfb' '4xpv' '1h1a' '5k7p' '1ree' '1te1' '2jic' '4hko' '4xqw' '5gyg' '1xnd' '5gyb' '5gy8' '1red' '2d97' '5gya' '1hix' '5gv1' '1xnk' '5gya' '5gyi' '2d98' '5gyg' '5gy9' '1ref' '5gye' '5gyc' '5gyh' '4hkw' '4hkl' '2vul' '5gy9' '3akq' '1h1a' '3mfc' '3zse' '5gyf' '1hix' '4hk8' '5jrm' '5jrn' '4xqd' '4xq4' '4xqd' '4hk9'); 
pdb_array=('2B45' '2B46'); 


for pdb in "${pdb_array[@]}"; do
  pdb_title=$(basename "$pdb")
  #echo "${pdb_title} is downloaded"
  wget "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=${pdb}" -O ${home}/pdb_xyl/${pdb_title}.pdb
done
  
