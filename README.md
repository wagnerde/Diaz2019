# Diaz2019
scRNA-seq data analysis for figures appearing in the following publication: In vitro characterization of the human segmentation clock. Margarete Diaz-Cuadros et. al. (2019)

## Installation
Please follow the following instructions to create a Python environment with required dependencies.  The first step is to use conda to create a python 3.6 environment:
```
conda create --name py36b python=3.6
```
Activate the environment and begin installing packages with conda and pip.  
```
source activate py36

conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
pip install scanpy
pip install fa2
```
# Install JupyterLab
conda install -c conda-forge jupyterlab

# Install leidenalg for Leiden clustering
conda config --add channels conda-forge
conda install leidenalg

# Install Scrublet
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .

# Install MulticoreTSNE
pip install MulticoreTSNE

# Install bbknn
conda install -c bioconda bbknn


```
