 # Diaz2019
scRNA-seq data analysis for figures appearing in the following publication: 
In vitro characterization of the human segmentation clock. Margarete Diaz-Cuadros et. al. (2019)

Please follow the following instructions to create a Python 3.6 environment with required dependencies:
```
# Python environment 'py36'
conda create --name py36 python=3.6

# Activate py36
source activate py36

# Install ScanPy & depenencies
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph louvain
pip install scanpy
pip install fa2

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
