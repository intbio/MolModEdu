#!/bin/sh
# this should be adapted to your supercomputer

# !! NB you have to change these pathways according to your system !!
#SYSTEM="16_GPCR_MARTINI"
CASE="test"
HOME=$(pwd)
WORKDIR=${HOME}


# !! NB you have to change these variables according to your version of Gromacs, MPI and number of used CPUs for MD run
PREP="gmx grompp"

PROG="gmx_mpi mdrun"


# We are performing equilibration and production run of the pre-minimized system;
# this is example sequence of the workflow, which can be customized;
cd $WORKDIR
#$PREP -f ${WORKDIR}/mdp/equil.mdp -c ${WORKDIR}/system_start.gro -r ${WORKDIR}/system_start.gro -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil1_${CASE} -maxwarn 1
#$PROG -v -deffnm equil1_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/equil_long.mdp -c ${WORKDIR}/equil1_${CASE} -r ${WORKDIR}/system_start.gro -t ${WORKDIR}/equil1_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil2_${CASE} -maxwarn 1
$PROG -v -deffnm equil2_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/equil_long2.mdp -c ${WORKDIR}/equil2_${CASE} -r ${WORKDIR}/system_start.gro -t ${WORKDIR}/equil2_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil3_${CASE} -maxwarn 2
$PROG -v -deffnm equil3_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/equil_long3.mdp -c ${WORKDIR}/equil3_${CASE} -r ${WORKDIR}/system_start.gro -t ${WORKDIR}/equil3_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil4_${CASE} -maxwarn 2
$PROG -v -deffnm equil4_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/equil_long3.mdp -c ${WORKDIR}/equil4_${CASE} -r ${WORKDIR}/system_start.gro -t ${WORKDIR}/equil4_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil5_${CASE} -maxwarn 2
$PROG -v -deffnm equil5_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/equil_long3.mdp -c ${WORKDIR}/equil5_${CASE} -r ${WORKDIR}/system_start.gro -t ${WORKDIR}/equil5_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/equil6_${CASE} -maxwarn 2
$PROG -v -deffnm equil6_${CASE} -rdd 1.8

$PREP -f ${WORKDIR}/mdp/dynamic.mdp -c ${WORKDIR}/equil6_${CASE} -t ${WORKDIR}/equil6_${CASE} -p ${WORKDIR}/topol.top -n ${WORKDIR}/index.ndx -o ${WORKDIR}/md_${CASE}_prod_5ms -maxwarn 2
$PROG -v -deffnm md_${CASE}_prod_5ms -rdd 1.8



exit
