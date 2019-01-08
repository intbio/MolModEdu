# Beginner's tutorial on molecular dynamics simulations in Gromacs

This a very simple tutorial with detailed instructions inteded to give the very basic yet practical skills in setting up MD simulations in Gromacs.
Prerequisites:
1. A computer running Linux or MacOS.
2. Basic familiarity with command line interface/Terminal is expected.


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
   + `conda install -c intbio gromacs=2018.3`
   + `conda install jupyter`
   + `conda install -c conda-forge mdanalysis wget nglview`
   + `conda install -c intbio seq_tools`
- Test that jupyter notebook is working `jupyter notebook`

## Step 2. Getting files.
- You will need to copy file from this directory to your computer, if you have not alread done so.
- One way to do it is to issue following commands
```
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/GROMACS/beginner/README.md
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/GROMACS/beginner/tutorial.ipynb
wget https://raw.githubusercontent.com/intbio/MolModEdu/master/GROMACS/beginner/xvg_plot.py

```

## Step 3. Start the tutorial.
- Launch jupyter notebook `jupyter notebook`
- Open the tutorial.ipynb file
- You can clear all prvious output Cell->All Outputs->Clear and follow the steps in that file by running them (Press the run button).