# Introductory tutorial on using PLUMED with Gromacs

This an introductory tutorial on using PLUMED with Gromacs.
Prerequisites:
1. Understanding of how to use Gromacs (see other tutorials).
2. A computer running Linux or MacOS.
3. Basic familiarity with command line interface/Terminal is expected.


Follow these preliminary steps to setup your computer.

## Step 1. Getting you computer and software ready.
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

## Step 2. Getting files.
- You will need to copy file from this directory to your computer, if you have not alread done so.
- One way to do it is to issue following commands
```
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/PLUMED/README.md
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/PLUMED/tutorial.ipynb
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/PLUMED/xvg_plot.py

```

## Step 3. Start the tutorial.
- Launch jupyter notebook `jupyter notebook`
- Open the [tutorial.ipynb](tutorial.ipynb) file
- You can clear all previous output Cell->All Outputs->Clear and follow the steps in that file by running them (Press the run button).
