# Diaz2019
scRNA-seq data analysis for figures appearing in the following publication: In vitro characterization of the human segmentation clock. Margarete Diaz-Cuadros et. al. (2019)

## Installation
Please follow the following instructions to create a Python environment with required dependencies.  The first step is to use conda to create a python 3.6 environment:
```
conda create --name py36b python=3.6
```
Activate the environment:
```
source activate py36
```
Begin installing packages with conda and pip:
```  
conda install seaborn scikit-learn statsmodels numba pytables
```
```
conda install -c conda-forge python-igraph louvain jupyterlab leidenalg
```
```
conda install -c bioconda bbknn
```
```
pip install scanpy fa2 MulticoreTSNE
```
Install Scrublet:
```
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
```

