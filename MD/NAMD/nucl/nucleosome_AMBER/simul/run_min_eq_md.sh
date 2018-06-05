#!/bin/bash
#PBS -N min_eq_md_1kx5nt_ge
#PBS -k oe
#PBS -m be
#

#Comment:
# this file is aimed at making a full minimization/relaxation/md protocol in one script
#
#
#

#Set Runtype here
RUNTYPE="PBS_BIOWULF" #"PBS_BIOWULF" or "LOCAL" or "PBS_BIOWULF_IPATH" 

if [ $RUNTYPE == "LOCAL" ]
# this is how you run it qsub -v np=64 -l nodes=8:ib pbsnamd.sh
then

RUN_NAMD="`which namd2` +p8"

fi


if [ $RUNTYPE == "PBS_BIOWULF" ]
# this is how you run it: qsub -v np=384 -l nodes=48:ib run_min_eq_md.sh
then
NAMD_VER=2.9
NETWORK=ib
TRANS=ibverbs

PATH=/usr/local/NAMD/$NAMD_VER/$NETWORK/$TRANS:$PATH

cd $PBS_O_WORKDIR

# Create host file (required)
make-namd-nodelist

RUN_NAMD="charmrun ++nodelist /home/panch/namd.$PBS_JOBID ++p $np ++remote-shell ssh `which namd2`"

fi





#main commands
echo "Running profile:" $RUNTYPE
echo "Namd command is:" $RUN_NAMD
#minimize and equilibrate
mkdir -p output
cd output

$RUN_NAMD ../input/min_equil.conf > min_eq.log

#run md
$RUN_NAMD ../input/md.conf > md.log


#rm *.BAK

if [ RUNTYPE == "PBS_BIOWULF" ]
then
rm -f ~/namd.$PBS_JOBID
fi



