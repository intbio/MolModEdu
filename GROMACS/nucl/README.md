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

#### Manuals and CheatCheets
- Bash CheatSheet [[1]](https://gist.github.com/LeCoupa/122b12050f5fb267e75f) or [[2]](https://devhints.io/bash)
- [GROMACS manual](http://ftp.gromacs.org/pub/manual/manual-5.0.4.pdf)


#### Textbooks 

- [Computer Simulation of Liquids](https://books.google.ru/books/about/Computer_Simulation_of_Liquids.html?id=O32VXB9e5P4C&redir_esc=y), Allen and Tildesley, 1989
- [Understanding Molecula Dynamics](https://www.sciencedirect.com/science/book/9780122673511), Frenkel and Smith, 2002
- [Molecular Driving Forces](https://books.google.ru/books/about/Molecular_Driving_Forces.html?id=hdeODhjp1bUC&redir_esc=y), Dill, 2003 
- [Intermolecular and Surface forces](https://www.sciencedirect.com/science/book/9780123751829), Israelachvili, 1985

#### Other useful tutorials

- Quick and easy tutorial [Lysozyme in Water](http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/lysozyme/index.html)

## 2. System description and simulation strategy choice
### Nucleosome core particle and its PDB structures
write a few words about NCP
Nucleosome structural studies davy 2011

будем использова 1kx5, 1,9 А
структуры с каким разрешением стоит брать для МД - бокс - ссылку на другой файл в докс - как оценить качество пдб структуры
при каком разрешении, какие начинатся проблемы (ссылка на пдб 101) x-ray resolution
[Protein crystallography for non-cristallographers](https://www.ncbi.nlm.nih.gov/pubmed/18034855), Wlodawer A., 2008

### Understanding your PDB structure 

- checking occupancy and b-factor
открыть в химере и сделать скрины, how to examine a pdb structure
chains number, 
there are dna, protein, water and ions - присутствуют структуры
наличие связанной воды - удалять или не удалять, как затем достраивать атомы водорода
разные ионы - что с ними делать
исследование occupancy - раскраскаа по ней и б-фактору( какие части со слишком большим факторам - их положению возможно нельз доверять, + формула для него

### Dealing with flexible histone tails

### Ionic conditions and choosing simulation box size

### Force field choice

Ion parameters chosen as [(Yoo & Aksiementiev, JPC, 2012)](https://pubs.acs.org/doi/abs/10.1021/jz201501a)

## 3. Installing software

## 4. Obtaining force field files

## 5. Pereparing system for simulation

