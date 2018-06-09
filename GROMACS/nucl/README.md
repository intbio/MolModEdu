# Tutorial on simulating molecular dynamics of a nucleosome core particle in Gromacs
## Directory structure
[docs](docs) - documentation for the current tutorial  
[prep](prep) - scripts to prepare the structure for simulations  
[simul](simul) - directory to perform simulations  
[analysis](analysis) - scripts to perform trajectory analysis  

## Contents
1. [**Introduction**](#Introduction) 
   + [Objectives](#Objectives)
   + [Prerequisites](#Prerequisites)
     +[]()





<a name="Introduction"/>

## 1. Introduction
This tutorial introduces students to molecular dynamics simulations method using GROMACS by the example of simulating the nucleosome core particle. 

<a name="Objectives"/>

### Objectives

- To gain an understanding of Molecular Dynamics Simulations method;
- To learn how to prepare a molecular system for MD simulations from a PDB file;
- To understand how to choose ForceField and correct stimulation conditions;
- To learn how to run MD simulations using GROMACS on a parallel computer cluster;
- To learn how to visualize and analyze MD trajectories;
- To provide a reusable set of scripts and examples that students may reuse to simulate the system of interest.

### Prerequisites
Basic knowledge of the following areas is required:
- Unix-like Operating System (e.g. Linux), recommended tutorial [here](http://swcarpentry.github.io/shell-novice/);
- Terminal / Command line interface and Bash scripting, recommended tutorial [here](http://swcarpentry.github.io/shell-novice/); 
- Basic Physical Chemistry and Biochemistry;
- PDB file structure, recommended tutorial [here](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction)

### Hardware/Software requirements for the tutorial
- Access to a Linux machine, ideally Ubuntu;
- Gromacs v. , [intallations instructions]();
- VMD v. , [intallations instructions]();
- Chimera v. , [intallations instructions]();
- Anaconda Python, [intallations instructions]();
...

### Reference materials and further reading

- Bash CheatSheet [1](https://gist.github.com/LeCoupa/122b12050f5fb267e75f) or [2](https://devhints.io/bash)
- [GROMACS manual](http://ftp.gromacs.org/pub/manual/manual-5.0.4.pdf)
- Quick and easy tutorial [Lysozyme in Water](http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/lysozyme/index.html)
- [Understanding Molecula Dynamics, Frenkel and Smith, 2002](https://www.sciencedirect.com/science/book/9780122673511)


## 2. System description and simulation strategy choice
### Nucleosome core particle and its PDB structures

### Dealing with flexible histone tails

### Ionic conditions and choosing simulation box size

### Force field choice

Ion parameters chosen as [(Yoo & Aksiementiev, JPC, 2012)](https://pubs.acs.org/doi/abs/10.1021/jz201501a)

## 3. Installing software

## 4. Obtaining force field files

## 5. Pereparing system for simulation

