# FeatureVectorGeneration

## Abstract 

Genomic variations may cause deleterious effects on protein functionality and perturb biological processes. Elucidating the effects of variations is important for developing novel treatments for diseases of genetic origin. Computational approaches have been aiding the work in this field by modeling and analyzing the mutational landscape. However, new approaches are required for accurate/comprehensive representation and cutting-edge data-centric analysis of sequence variations. 

In this study, we propose a new method for featurizing single amino acid variations (SAVs) on proteins called ASCARIS (Annotation and StruCture-bAsed RepresentatIon of SAVs) to be used in data-driven modeling of variations for various purposes such as predicting the functions of protein variants or constructing multi-omics-based models. We evaluated variations’ function-related properties by utilizing a combination of sequence annotations from UniProt and 3-D structural information from PDB and AlphaFold-DB. For this, we extracted and analyzed the correspondence between the varied residue and 30 different sequence-based feature annotations (e.g., active/lipidation/glycosylation sites; calcium/metal/DNA binding, inter/transmembrane regions, etc.), together with structural features such as protein domains, the location of variation (e.g., core/interface/surface), and the change in physicochemical properties due to the variation. We also mapped mutated and annotated residues to the 3-D structures of corresponding proteins and calculated the spatial distances in-between since proximity (e.g., sharing an interface) may also effect functionality

We quantitatively investigated the relationships between each of these features and the consequences of variations, and finally constructed 68-dimensional feature vectors to represent SAVs in a large dataset composed of ~100,000 data points. To analyze potential applications of ASCARIS, we trained machine learning-based variant effect predictor models that utilise ASCARIS representations as input. We carried out both an ablation study and comparison against the state-of-the-art methods over well-known benchmark datasets. According to our results, our method displays competing and complementary performance against widely-used predictors. ASCARIS can be used either alone or in combination with other approaches, to universally represent SAVs from a functional perspective, for intensive data-driven analysis of genomic variations.


<img width="706" alt="Screen Shot 2022-06-08 at 7 14 23 PM" src="https://user-images.githubusercontent.com/26777185/172726336-8ccc2914-4253-4ba7-b534-3581526651e0.png">

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

- **swissmodel_structures.txt.zip** : Includes summary file for Swiss-Model structure summary. Swiss model summary (INDEX-metadata) files are downloaded separately for each organism from https://swissmodel.expasy.org/repository and all of them manually merged into a single file. Must be unzipped to input files folder prior to usage.
- **domains.txt** : Includes InterPro domains simplified as in the following order *(tab separated)* --> 
  [uniprotID      domainID        domainStartPosition     domainEndPosition]
- **significant_domains** :  Selected domains from *domains.txt* file according to Fisher's Exact Test result. Fisher's Exact Test applied to all domains in the training test to assess their significance with respect to the the deleteriousness outcome. p_values is chosen as 0.01.
- **H_sapiens_interfacesHQ.txt** :  High confidence interfaces downloaded from [Interactome Insider](http://interactomeinsider.yulab.org/downloads.html) for *Homo sapiens*
- **index.json** : [Swiss-Model metadata](https://swissmodel.expasy.org/repository) for model information is downloaded to obtain Swiss-Models. Zipped version can be found under zipped_/. 
- **alphafold_structures** : [AlphaFold Human proteome predictions](http://ftp.ebi.ac.uk/pub/databases/alphafold/latest/) for structure predictions from Alphafold. zipped_/UP000005640_9606_HUMAN.tar file is unzipped to alphafold_structures folder. 
- **alphafold_summary: Processed data for AlphaFold structures. Includes protein ideentifier, chain id, sequence, model count for each entry.


### datasets

Datasets that are used to create machine learning models are provided in **datasets** folder. 

- **training_uptodate_full.txt** : Full training dataset obtained from UniProt, PMD and ClinVar.
- **training_uptodate_full_2014selected.txt** : 2014 subset of the training set.
- **swiss_ready.txt** : Benchmark set obtained from SwissVar database.
- **varibench_ready.txt** : Benchmark set obtained from VariBench database.
- **psnp_ready.txt** : Benchmark set obtained from PredictSNP database.
- **test_MT_benchmark_datapoints_wo_training_datapoints** : Benchmark set obtained from MutationTaster data.


## Usage
```
python3 main.py -o 1 -i P13637-T-613-M<br/>
python3 main.py -o 2 -i 'P13637-T-613-M, Q9Y4W6-N-432-T, Q9Y4W6-N-432-T'<br/>
python3 main.py -o 2 -i sample_input.txt<br/>
```
### input arguments

-o :  input option. 1: Use PDB-ModBase-SwissModel structures, 2: Use AlphaFold Structures<br/>
-i :  input datapoint. Datapoint/input file name.<br/>
 *     >  Option 1: Comma-separated list of idenfiers (UniProt ID-wt residue-position-mutated residue (e.g. Q9Y4W6-N-432-T or Q9Y4W6-N-432-T, Q9Y4W6-N-432-T))  
 *     > Option 2: Enter tab-separated file path



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
- **pdb_structures** : Contains downloaded structure files from PDB for input proteins when applicable. If the user has a folder wherein PDB structures are stored, this folder might be used to decrease run time. In this case,please change the extension to .txt
- **swissmodel_structures** : Contains downloaded model files from SwissModel for input proteins when applicable.
- **modbase_structures** : Contains downloaded model files from ModBase for input proteins when applicable. Each file contains all models related to one protein.
- **modbase_structures_individual** : Contains downloaded model files from ModBase for input proteins when applicable. Each file contains individual models related to one protein.
- **alignment_files** : Contains alignment files of protein sequences with structure files. 
- **3D_alignment** : Contains alignment files of structure files. This step is performed in order to avoid missing residues in the PDB files.
- **freesasa_files** : Contains calculated FreeSASA values for each data point.




## Description of Output Vector

As can be seen in the figure below, dimensions 1-4 correspond to physicochemical property values, 5-6 correspond to domain-related information, 7-36 correspond to binary information of mutations with respect to their presence within annotation regions, 37-66 coorespond to Euclidian distance between mutation point and closest point of the annotation, 67-68 coorespond to information regarding mutation's position on the protein in terms of core, interface or surface.

<img width="1082" alt="Screen Shot 2021-06-24 at 12 32 04 AM" src="https://user-images.githubusercontent.com/26777185/123170836-a640c900-d483-11eb-90eb-473d826a2a75.png">


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

