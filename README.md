# Diaz2019
scRNA-seq data analysis for figures appearing in the following publication: In vitro characterization of the human segmentation clock. Margarete Diaz-Cuadros et. al. (2019)

## Installation
Start by cloning this repository using git:  
```
git clone https://github.com/wagnerde/Diaz2019.git
cd Diaz2019
```
scRNA-seq data files are too large to include within the Github repository, and must be downloaded separately.  Download and unzip the Diaz et. al. 2019 scRNA-seq data. The unzipped '\_rawData' directory should then reside in the same 'Diaz2019' directory:  
```
wget https://kleintools.hms.harvard.edu/paper_websites/diaz_2019/Diaz2019_inDropsCountsTables.zip --no-check-certificate
unzip Diaz2019_inDropsCountsTables.zip
```
Create a Python 3.6 conda environment to manage the required software packages:
```
conda create --name py36 python=3.6 -y
```
Activate the environment:
```
source activate py36
```
Begin installing packages with conda and pip:
```  
conda install -y -q seaborn scikit-learn statsmodels numba pytables
conda install -y -q -c conda-forge python-igraph louvain jupyterlab leidenalg
conda install -y -q -c bioconda bbknn
pip install scanpy fa2 MulticoreTSNE
```
Install Scrublet:
```
git clone https://github.com/AllonKleinLab/scrublet.git
cd scrublet
pip install -r requirements.txt
pip install --upgrade .
cd ..
```

## Run the Jupyter notebooks
Activate the py36 environment and run Jupyter Lab.  Within Jupyter lab, navigate to one of the dataset folders (mmE95, mmES, hsIPS), and open the .ipynb file.
```
source activate py36
jupyter lab
```


