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
   + [Hardware/Software requirements](#Hardware)
   + [Reference materials and further reading](#Reference_materials)
      + [Manuals and CheatCheets](#Manuals)
      + [Textbooks](#Textbooks)
      + [Other useful tutorials](#Other_tutorials)
2. [**System description and simulation strategy choice**](#System)
   + [Nucleosome core particle and its PDB structures](#NCP)
   + [Understanding your PDB structure](#Understanding_PDB)
   + [Dealing with flexible histone tails](#H_tails)
   + [Ionic conditions and choosing simulation box size](#Ions_box)
   + [Force field choice](#ForceField)
3. [**Installing software**](#installing_soft)
4. [**Obtaining force field files**](#Obtaining_FF_files)
5. [**Pereparing system for simulation**](#before_stimulation)
      



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

<a name="Prerequisites"/>

### Prerequisites
Basic knowledge of the following areas is required:
- Unix-like Operating System (e.g. Linux), recommended tutorial [here](http://swcarpentry.github.io/shell-novice/);
- Terminal / Command line interface and Bash scripting, recommended tutorial [here](http://swcarpentry.github.io/shell-novice/); 
- Basic Physical Chemistry and Biochemistry;
- PDB file structure, recommended tutorial [here](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction)

<a name="Hardware"/>

### Hardware/Software requirements for the tutorial
- Access to a Linux machine, ideally Ubuntu;
- Gromacs v. , [intallations instructions]();
- VMD v. , [intallations instructions]();
- Chimera v. , [intallations instructions]();
- Anaconda Python, [intallations instructions]();
...

<a name="Reference_materials"/>

### Reference materials and further reading

<a name="Manuals"/>

#### Manuals and CheatCheets
- Bash CheatSheet [[1]](https://gist.github.com/LeCoupa/122b12050f5fb267e75f) or [[2]](https://devhints.io/bash)
- [GROMACS manual](http://ftp.gromacs.org/pub/manual/manual-5.0.4.pdf)

<a name="Textbooks"/>

#### Textbooks 

- [Computer Simulation of Liquids](https://books.google.ru/books/about/Computer_Simulation_of_Liquids.html?id=O32VXB9e5P4C&redir_esc=y), Allen and Tildesley, 1989
- [Understanding Molecula Dynamics](https://www.sciencedirect.com/science/book/9780122673511), Frenkel and Smith, 2002
- [Molecular Driving Forces](https://books.google.ru/books/about/Molecular_Driving_Forces.html?id=hdeODhjp1bUC&redir_esc=y), Dill, 2003 
- [Intermolecular and Surface forces](https://www.sciencedirect.com/science/book/9780123751829), Israelachvili, 1985

<a name="Other_tutorials"/>

#### Other useful tutorials

- Quick and easy tutorial [Lysozyme in Water](http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/lysozyme/index.html)

<a name="System"/>

## 2. System description and simulation strategy choice

<a name="NCP"/>

### Nucleosome core particle and its PDB structures
Nucleosome is a basic unit of DNA package in eucaryotes. 
Nuclosome Core particle (NCP) consists of 1.67 left-handed super-helical turns of B-form DNA around an octamer of histone proteins. DNA wrap is about 167 bp. The octamer consists of 2 copies each of the core histones H2A, H2B, H3, and H4. 

Further information: [Nucleosome structural studies,](https://www.ncbi.nlm.nih.gov/pubmed/21176878), Davey, 2011

In this tutorial we are going to use 1KX5 PDB structure. You can try to find it on RSCB PDB by yourself or download from this [link](https://www.rcsb.org/structure/1kx5). We have chosen this structure because it has the best resolution - 1.94 A.  

We can 
[Check](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/resolution) what kind of problems you will get if you choose the structure with wrong resolution. 
структуры с каким разрешением стоит брать для МД - бокс - ссылку на другой файл в докс

Further information: [Protein crystallography for non-cristallographers](https://www.ncbi.nlm.nih.gov/pubmed/18034855), Wlodawer A., 2008

<a name="Understanding_PDB"/>

### Understanding your PDB structure 

Open your PDB structure in VMD using this command in command line
Make sure you are in the exact directory where 1kx5.pdb was downloaded. Then type:
> VMD 1kx5.pdb

Then you will see this:

<img src="../docs/understanding1.png" width="500">

*Fig.1. VMD presentation of the molecule*

For better representation I changed background from black to white which is not necessary for you. I also deleted axes. When you open your system in VMD the background is automaticaly black. 

For better and more understandable view open in VMD Main window:
> Graphics > Representations  

A new window called "Graphical Representation" will be opened. Then choose:
> Drawing method > NewCartoon 

You will see this: 

<a name="fig2"/>
<img src="../docs/understanding2.png" width="500">

*Fig.2. NewCartoon Drawing Method*


**Coloring by occupancy and B-factor**. 

At first, we need to understand what are occupancy and B-factor. 

Macromolecule crystall consists of many individual molecules. These molecules are packed into a symmetrical unit. In some crystalls there are slight conformational differences between the molecules. Scienmtist use the *occupancy* to estimate the amount of each conformation that is observed in the crystal. 

In the same window ("Graphical Representation"):
> Coloring Method > Occupancy

<img src="../docs/understanding3.png" width="500">

*Fig.3. Coloring by the Occupancy*

*B-factor* (often refered to as atomic displacement parameter) monitors the positional spread of each atom. It describes the displacement of the atomic positions from an average (mean) value.

For coloring by B-factor choose:
> Coloring Method > Beta

<img src="../docs/understanding4.png" width="500">

*Fig.4. Coloring by B-factor*

In this picture blue is an area with low B-factor. It means that these areas are more flexible than red ones. Consequently, the position of these sites probably cannot be trusted  

**DNA, Protein, Water and Ions**

Let's get back to the [Fig.2](#fig2) state with NewCartoon Drawing method. Now we're going to explore other molecules inside this system. We're going to start from water molecules. 

We are still working with Graphical Representation window.
First, you need to Create Rep. Tap on it and type "water" in Selected Atoms. Coloring method is VDW which is Van-der-Waals. You will see this:

<img src="../docs/water.png" width="500">

*Fig.5. Showing water molecules*

You will do the same to see the Ions except typing "ions" in Selected Atoms.

<img src="../docs/ions.png" width="500">

*Fig.6. Showing ions*

You can also choose DNA or protein the same way as described earlier. Now  you know how to understand your PDB system.

**The next step is to understand should we change water and ions condition**












наличие связанной воды - удалять или не удалять, как затем достраивать атомы водорода
разные ионы - что с ними делать


<a name="H_tails"/>

### Dealing with flexible histone tails

<a name="Ions_box"/>

### Ionic conditions and choosing simulation box size

<a name="ForceField"/>

### Force field choice

Ion parameters chosen as [(Yoo & Aksiementiev, JPC, 2012)](https://pubs.acs.org/doi/abs/10.1021/jz201501a)

<a name="installing_soft"/>

## 3. Installing software

<a name="Obtaining_FF_files"/>

## 4. Obtaining force field files

<a name="before_stimulation"/>

## 5. Pereparing system for simulation

