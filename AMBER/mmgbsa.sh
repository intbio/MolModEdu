# this is template of shell script that uses amber-tools for MMGBSA analysis
# it demonstrates how to calculate protein-ligand affinity (free energy of binding) for amber trajectories
# see the tutorial for details and all possible options: http://ambermd.org/tutorials/advanced/tutorial3/

# set workdir
analysis="/mydir/Analysis"
# system to be analysed
sim=${analysis}/NpXynCYS_com_300K
# indicate name of the trajectory
simulation="NpXynCYS_com_300K"

# source amber
source /opt/gbamod/application/amber16/amber.sh

# create separate topologies for protein ligand and complexes
ante-MMPBSA.py -p ${sim}/protein.prmtop -c ${sim}/protein.prmcrd -c complex.prmtop -r receptor.prmtop -l ligand.prmtop -s :WAT:Cl-:NA+:K+ -n :ROH:4XB:0XB > ante_${simulation}.log
# run mmgbsa using 20 CPUS to calculate free energy of binding
# since we decompose results to obtain individual contribution of each amino acids into the resulted dG
# so we don't calculate entropy for dG in this example to speed up the calculation
mpirun -np 20 MMPBSA.py.MPI -O -i mmgbsa_params.in -o mmgbsa_${simulation}.dat -do DECOMP_${simulation}.dat -sp ${sim}/protein.prmtop -cp complex.prmtop -rp receptor.prmtop -y ${sim}/${simulation}.nc  -lp ligand.prmtop > progress_${simulation}.log
