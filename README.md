# Diaz2019
scRNA-seq data analysis for figures appearing in the following publication: In vitro characterization of the human segmentation clock. Margarete Diaz-Cuadros et. al. (2019)

## Installation
Start by cloning this repository using git:  
```
git clone https://github.com/wagnerde/Diaz2019.git
cd Diaz2019

```
Download and unzip the Diaz et. al. 2019 scRNA-seq data. The unzipped '\_rawData' directory should be placed into the same 'Diaz2019' directory:  
```
wget https://kleintools.hms.harvard.edu/paper_websites/diaz_2019/Diaz2019_inDropsCountsTables.zip --no-check-certificate
unzip Diaz2019_inDropsCountsTables.zip
```
Create a Python 3.6 conda environment to manage required packages:
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



