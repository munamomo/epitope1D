# epitope1D
Accurate taxonomy-aware b-cell linear epitope prediction.

epitope1D, an explainable machine learning method capable of accurately identifying linear B-cell epitopes, leveraging two new descriptors: 
- Graph-based signature representation of protein sequences;
- Organism Ontology information.

You can also use epitope1D in a friendly-user web interface, here: https://biosig.lab.uq.edu.au/epitope1d/

## Install

1) First clone this repository to you local machine.
2) We recommend the use of a virtual environment (anaconda, for instance) with the following packages:

channels:
  - conda-forge
  - bioconda
  - defaults
  - anaconda
dependencies:
  - python=3.6.13
  - joblib=1.2.0
  - pip
  - pandas
  - numpy
  - pip:
    - aaindex==0.0.2
    - scikit-learn==0.24.1

3) Download the .sav file of the trained ML model located at [Dropbox](https://www.dropbox.com/s/yn535ic814raypp/epitope1d_RF_model.sav). Since it's a large file, GitHub recommends Dropbox storage.
	

=> **To sum up**, the installation would go like this in a Linux command line (assuming you use conda package management system):
```
$ git clone https://github.com/munamomo/epitope1D.git
```
--> Open the newly created directory epitope1D, then proceed:
```
$ conda env create -f requirements.yml

$ conda activate epitope1d

$ cd ./epitope1d/models

$ wget https://www.dropbox.com/s/yn535ic814raypp/epitope1d_RF_model.sav
```
---------------------------------------------------------------------------------------------------------------------------------------------------------
## How to use it

epitope1D runs with python v.3.6. 
The following illustrates how to run it under the default folder.
```
python run_epitope1d.py -fasta example.fasta -organism Ribozyviria -window 25
```
Where:
- -fasta, Path to the input fasta file. Multiple peptides or protein sequences can be provided as input in fasta format file. Each sequence should be placed in a different line and must have a header, starting with '>'. Also, sequences should be of the same Taxonomy group and be at least 6 amino acid in length.
- -organism, Organism Taxonomy name. Only one per fasta file is accepted. See options bellow. Type: string.
- -window, Desired window size (amino acid length) to screen for linear epitopes, default = 25. Type: integer (from 6 to 25).

There are 20 possible organisms available, chose 1 per fasta file: 
```
Metamonada, Discoba, Sar, Viridiplantae, Opisthokonta, Terrabacteria_group, Proteobacteria, PVC_group,
Spirochaetes, FCB_group, Thermodesulfobacteria, Fusobacteria, Riboviria, 
Duplodnaviria, Monodnaviria, Varidnaviria, Ribozyviria, Anelloviridae, Naldaviricetes, Adnaviria.
```

---------------------------------------------------------------------------------------------------------------------------------------------------------
## Output

Every time you run epitope1D, it will create a random identifier number and folder with the corresponding output result file.
The output prediction is a .csv file with the following collumns:
```
- ID: sequence ordering of the input file considering the peptide spliting by window length.
- Fasta_header: informed fasta header from input file (information from white space are disregarded)
- Peptide: splited peptide, from original sequence, according to window length.
- Prediction: 0 or 1
- Score_Epitope: probability score (from 0 to 1. Threshold is 0.5)
- Classification: non-epitope (prediction 0) or epitope (prediction 1)
```
---------------------------------------------------------------------------------------------------------------------------------------------------------
## Citation
```
@article {da Silva2022,
	author = {da Silva, Bruna Moreira and Ascher, David B. and Pires, Douglas E. V.},
	title = {epitope1D: Accurate Taxonomy-Aware B-Cell Linear Epitope Prediction},
	year = {2022},
	journal = {bioRxiv}
 }
 ```
---------------------------------------------------------------------------------------------------------------------------------------------------------
