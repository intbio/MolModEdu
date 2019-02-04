# Tutorials on enhanced sampling techniques using PLUMED and Gromacs.

A number of  tutorials on enhanced sampling techniques in MD simulations using PLUMED and Gromacs.
All tutorials are provided as Jupyter notebooks.
## List of tutorials
1. [Understanding PLUMED and collective variables](plumed_intro.ipynb)
2. [Introduction to biasing simulations with PLUMED](plumed_bias.ipynb)
3. [Introduction to running Metadynamics](plumed_metad_intro.ipynb)
4. [Understanding Metadynamics parameters](plumed_metad_params.ipynb)
5. [Well tempered Metadynamics](plumed_WT_metad.ipynb)
6. [Advanced analysis with PLUMED (in Development)](plumed_analysis.ipynb)


## Prerequisites:
1. Understanding on how to use Gromacs (see other tutorials).
2. A computer running Linux or MacOS.
3. Basic familiarity with command line interface/Terminal is expected.
4. Install software (see below).

## Additional reading
- [A quick introduction to PLUMED](https://www.youtube.com/watch?v=PxJP16qNCYs)
- [Video lecture on MetaDynamics](https://www.youtube.com/watch?v=bZZggbV2r5E)
- [PLUMED manual](https://plumed.github.io/doc-v2.3/user-doc/html/_syntax.html)

## Getting you computer and software ready.
- You'll need a computer running Linux (preferably Ubuntu Linux) or MacOS.
- Install Miniconda Package Manager with Python 3:
   + Download archive from this link https://conda.io/miniconda.html
   + Open Terminal, navigate to the directory containing the downloaded file.
   + Make it executable ```chmod u+x DOWNLOADED_FILE.sh``` and run ```./DOWNLOADED_FILE.sh```
- Familiarize yourself with basic conda commands https://github.com/intbio/IT_notes/blob/master/conda.md
- Make a dedicated environment for our tutorial ```conda create --name MD```
- Activate the environment ```source activate MD```
- We will need to issue a number of commands to install the components needed for our tutorial:
   + `conda install -c conda-forge -c intbio gromacs=2018.4_plumed_2.5.0`
   + `conda install jupyter`
   + `conda install -c conda-forge mdanalysis wget nglview panedr`
   + `conda install -c intbio seq_tools`
- Test that jupyter notebook is working `jupyter notebook`

## Getting and running tutorial files.
- You will need to download and copy ipynb tutorial files from this directory to your computer and open them in jupyter notebook. Launch jupyter notebook with `jupyter notebook`.
- You can clear all previous output Cell->All Outputs->Clear and follow the steps in the tutorial file by running them (Press the run button).


