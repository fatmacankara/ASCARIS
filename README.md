# FeatureVectorGeneration

## Abstract 

Abstract comes here

## Development and DependenciesCancel changes 

- [Python 3.7.3](https://www.python.org/downloads/release/python-373/)
- [Pandas 1.1.4](https://pandas.pydata.org/pandas-docs/version/1.1.4/getting_started/install.html)
- [Numpy 1.19.5](https://numpy.org/devdocs/release/1.19.5-notes.html)
- [Ssbio 0.9.9.8.post1](https://pypi.org/project/ssbio/)
- [Freesasa 2.0.3.post7](https://pypi.org/project/freesasa/2.0.3.post7/)
- [Requests 2.22.0](https://pypi.org/project/requests/)
- [Biopython 1.78](https://biopython.org/docs/1.78/api/Bio.html)

## Descriptions of folders and files in the FeatureVectorGeneration repository

### input_files
Inside **input_files** folder, users can find files that are necessary to work the code. Below, explanations related to the files can be found.

- **domains.txt** : Includes InterPro domains simplified as in the following order *(tab separated)* --> 
  [uniprotID      domainID        domainStartPosition     domainEndPosition]
- **significant_domains** :  Selected domains from *domains.txt* file according to Fisher's Exact Test result. Fisher's Exact Test applied to all domains in the training test to assess their significance with respect to the the deleteriousness outcome. p_values is chosen as 0.01.
- **H_sapiens_interfacesHQ.txt** :  High confidence interfaces downloaded from [Interactome Insider](http://interactomeinsider.yulab.org/downloads.html) for *Homo sapiens*
- **index.json** : [Swiss-Model metadata](https://swissmodel.expasy.org/repository) for model information is downloaded to obtain Swiss-Models. Please unzip the file before using. 
- **UP000005640_9606_HUMAN** : [AlphaFold Human proteome predictions](http://ftp.ebi.ac.uk/pub/databases/alphafold/latest/) for structure predictions from Alphafold. Please unzip the file to 'alphafold_structures' folder before using. 
### datasets

Datasets that are used to create machine learning models are provided in **datasets** folder. 

- **training_uptodate_full.txt** : Full training dataset obtained from UniProt, PMD and ClinVar.
- **training_uptodate_full_2014selected.txt** : 2014 subset of the training set.
- **swiss_ready.txt** : Benchmark set obtained from SwissVar database.
- **varibench_ready.txt** : Benchmark set obtained from VariBench database.
- **psnp_ready.txt** : Benchmark set obtained from PredictSNP database.
- **test_MT_benchmark_datapoints_wo_training_datapoints** : Benchmark set obtained from MutationTaster data.


### usage

python main.py

### input arguments

- Enter Query DataPoint \n
  Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T))
  Option 2: Enter comma-separated file path
- Enter 1 to use PDB-ModBase-SwissModel structures, Enter 2 to use Alphafold structures



### sample_run

Contains results of a sample run from **sample_input.txt** input file. Example input file is shown below. Columns represent UniProt ID of the protein, wild type amino acid, position of the amino acid change and mutated amino acid, respectively. Input file must be given without a header.

```
P12694	C	264	W
P13637	T	613	M
P05067	I	716	V
P41180	E	604	K
P08123	G	646	C
P06731	I	80	V
P29474	D	298	E
Q16363	Y	498	H
P23560	V	66	M
Q00889	H	85	D
```
Files in **out_files** folder are created by running the script on **sample_input.txt** file.

- **feature_vector.txt** : Sample feature vector.
- **pdb_structures** : Contains downloaded structure files from PDB for input proteins when applicable.
- **swissmodel_structures** : Contains downloaded model files from SwissModel for input proteins when applicable.
- **modbase_structures** : Contains downloaded model files from ModBase for input proteins when applicable. Each file contains all models related to one protein.
- **modbase_structures_individual** : Contains downloaded model files from ModBase for input proteins when applicable. Each file contains individual models related to one protein.
- **alignment_files** : Contains alignment files of protein sequences with structure files. 
- **3D_alignment** : Contains alignment files of structure files. This step is performed in order to avoid missing residues in the PDB files.
- **freesasa_files** : Contains calculated FreeSASA values for each data point.

## Output File Legend


| Column name in the output file  | Description | 
| ------------- | ------------- |
| uniprotID_wildtype_aminoacidposition_mutatatedaminoacid | Datapoint Identifier |
| composition | Change in composition values upon aa change. Composition is defined as the atomic weight ratio of hetero (noncarbon) elements in end groups or rings to carbons in the side chain.  |
| polarity | Change in polarity values upon aa change. |
| volume | Change in volume values upon aa change. |
| granthamScore | Change in Grantham scores values upon aa change. The Grantham score predicts the distance between two amino acids. Higher Grantham scores are considered more deleterious. The distance scores published by Grantham range from 5 to 215.|
| domain | Domain IDs : Categorical  |
| domain_fisher | Domains that are not found to be significant in Fisher's Exact Test are labelled as DomainX for simplicity.|
| domaindistance3D | Euclidian distance between the domains from its closest residue and the mutation site. |
| disulfide, intMet, intramembrane, naturalVariant, dnaBinding, activeSite, nucleotideBinding, lipidation, site, transmembrane, crosslink, mutagenesis, strand, helix, turn, metalBindingi repeat, caBinding, topologicalDomain, bindingSite, region, signalPeptide, modifiedResidue, zincFinger, motif, coiledCoil, peptide, transitPeptide, glycosylation, propeptide | Minimum distance of annotations from the mutation on protein structure : Real values |
| disulfideBinary,  intMetBinary, intramembraneBinary, naturalVariantBinary, dnaBindingBinary, activeSiteBinary, nucleotideBindingBinary, lipidationBinary, siteBinary, transmembraneBinary, crosslinkBinary, mutagenesisBinary, strandBinary, helixBinary, turnBinary, metalBindingBinary, repeatBinary, caBindingBinary, topologicalDomainBinary, bindingSiteBinary, regionBinary, signalPeptideBinary, modifiedResidueBinary, zincFingerBinary, motifBinary, coiledCoilBinary, peptideBinary, transitPeptideBinary, glycosylationBinary, propeptideBinary | Binary labels for UniProt annotations: 0: not annotated, 1: annotation present, but mutation is not on annotation, 2: mutation is on the annotation :  Binary |
| sasa | SASA values :  Real values |
| threeState_trsh4_HQ| categories for SASA values: categorical (surface-core-interface) |


## Description of Output Vector

As can be seen in the figure below, dimensions 1-4 correspond to physicochemical property values, 5-6 correspond to domain-related information, 7-36 correspond to binary information of mutations with respect to their presence within annotation regions, 37-66 coorespond to Euclidian distance between mutation point and closest point of the annotation, 67-68 coorespond to information regarding mutation's position on the protein in terms of core, interface or surface.

<img width="1082" alt="Screen Shot 2021-06-24 at 12 32 04 AM" src="https://user-images.githubusercontent.com/26777185/123170836-a640c900-d483-11eb-90eb-473d826a2a75.png">


## Overall process

<img width="654" alt="Screen Shot 2021-07-06 at 6 07 30 PM" src="https://user-images.githubusercontent.com/26777185/124623980-22c0a800-de85-11eb-9b32-ff31a1f0a758.png">


