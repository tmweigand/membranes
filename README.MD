Collection of Membrane Anaylsis Tools 

LAMMPS Python Link: https://github.com/lammps/lammps/tree/develop/python


Installation Instructions for AmberTools (from https://ambermd.org/GetAmber.php#ambertools):

Most users should create a new conda "environment" to avoid conflicts with what you already have installed. To do this:

conda create --name AmberTools25 python=3.12 (Only done once.)

Next, install the code, which also only needs to be done once:

  conda activate AmberTools25
  conda config --add channels conda-forge
  conda config --set channel_priority strict
  conda install dacase::ambertools-dac=25


  To use the code:

  conda activate AmberTools25
  source $CONDA_PREFIX/amber.sh
        cd to some working folder, and run AmberTools programs
  conda deactivate


## Installation

### Step 1 — Create conda environment
```bash
conda env create -f environment.yml
conda activate membranes
```

### Step 2 — Build LAMMPS from source
```bash
bash src/scripts/install_lammps.sh
```

### Step 3 — Install the package
```bash
pip install -e .
```



Riley Examples
  https://zenodo.org/records/17274528


  variation from antechamber and what was used in paper because GAFF may be correct term but may be something else in mix. OPLS is what was used in postdoc. Using OPLS and a 4 point water molecule is one easy step up. Does Antechamber come with OPLS but lipargen does (https://traken.chem.yale.edu/ligpargen/). OPLS primarly for organics in liquid. 

  LigparGen and Antechamber can account account for reactons. Reduced to minimum components. Process unreacted TMC and MPD and then reacted Oligomer. Ensure that all data are the same except for reaction sites and 3 atoms awat from reactions sites (dihreals!). 