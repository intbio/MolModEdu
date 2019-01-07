# MolModEdu
Educational/Research Materials on Molecular Modeling: Tutorials, Code Examples and Templates, Resources, Course Materials

## Contents
1. [**Molecular Dynamics Simulations**](#MD) 
   + [Course and educational materials](#MDpresent)
   + [Tutorials](#MDtuturials)
   + [Code examples and templates](#MDcode)
        - [Gromacs](#MDcode_gromacs)
        - [NAMD](#MDcode_namd)
        - [MDanalysis](#MDcode_mdanalysis)
        
     
2. Miscelaneous
   + [Phenix](#phenix)



<a name="MD"/>

## Molecular Dynamics Simulations

<a name="MDpresent"/>

### MD: Courses and educational materials
- [Overview lecture by Alexey K. Shaytan](https://www.dropbox.com/s/y971h41by5wz0io/md_presentation.pptx?dl=0)

<a name="MDtutorials"/>

### MD: Tutotorials
- [A detailed Gromacs tutorial: simulating nucleosome dynamics. GROMACS/nucl](GROMACS/nucl)

<a name="MDcode"/>

### MD: Code examples and templates

<a name="MDcode_gromacs"/>

#### Gromacs
- The recommended way for setting up an in-house MD simulation is to use the extensive [gmx_template]( https://github.com/intbio/gmx_template) (currently available only locally), ideally in conjunction with the customly created ```moldyn``` conda environment on Newton linux cluster and bindings to Lomonosov-2 supercomputer.
- [A quick example of using Jupyter notebooks and MDAnalysis for preparing and running MD in GROMACS](MDanalysis/Nucleosome_dimer_MD_preparation.ipynb)
- [Molecular dynamics simulation of biomembranes at different temperature](https://github.com/intbio/MolModEdu/tree/master/GROMACS/biomembrane)

##### Gromacs Resources
- Gromacs force fileds https://github.com/intbio/gromacs_ff
- Typical Gromacs MDP files and protocols https://github.com/intbio/gmx_protocols
- Gromacs packed in conda package manager https://github.com/intbio/gromacs-conda 

<a name="MDcode_namd"/>

#### NAMD
- Example of a pipeline to simulate nucleosome dynamics and perform trajectory analysis from [JMB paper](https://www.ncbi.nlm.nih.gov/pubmed/26699921) - [here](MD/NAMD/nucl). 

<a name="MDcode_mdanalysis"/>

#### MDAnalysis
- [A quick example of using Jupyter notebooks and MDAnalysis for preparing and running MD in GROMACS](MDanalysis/Nucleosome_dimer_MD_preparation.ipynb)

<a name="phenix"/>

## PHENIX - Molecular modeling for X-ray crystallography
[https://www.phenix-online.org](https://www.phenix-online.org)
- [Minimizing molecular geometry using phenix](phenix/geo_minim.md)
