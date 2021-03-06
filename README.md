# metaPathwayMap
Associating chemical class predictions with MetaCyc pathway models

## Current version
v1.0

No older versions

## Alternate strategy
The workflow below enables you to make your own pathway representations for a species present in MetaCyc/Plant Metabolic Network. If you do not need to make your own pathway representations and are satisfied with PlantCyc, BrachypodiumCyc and SolCyc, you can use the Sol Genomics Network tool (https://metapathwaymap.solgenomics.net). 

If you are using plant data, it is recommended that you try your data with both PlantCyc and your species-specific database, since PlantCyc has additional pathways and compounts that may not be present in your species of interest.

## Requirements
Unix<br/>
Python3<br/>
Following python modules: Numpy, Scipy, Pandas, Networkx<br/>
CANOPUS output (adducts.tsv file)

The get_similar_pathways_wrapper.py script should take ~3-5 minutes to run, while the metaPathwayMap script completes in <10 seconds on a standard computer. 

## How to get CANOPUS output
1. Get mgf format file from your LC-MS dataset. This can be obtained from any LC-MS data analysis software such as MS-DIAL and XCMS
2. Download and install the SIRIUS software (https://bio.informatik.uni-jena.de/software/sirius/)
3. Import the mgf file into SIRIUS
4. Click on Compute All. You may want to choose the following parameters: SIRIUS (Select elements: CHNOPS), CSI-FingerID (Bio Databases or MetaCyc or PlantCyc)
5. Export Summaries
6. Ensure that in the tab-delimited file output by SIRIUS (the adducts file), the Most Specific Class is in Column 6, and subsequent annotations are after that. 
7. Use this file as input for metaPathwayMap

## Instructions for using metaPathwayMap
1. Download the following files into the working directory:<br/>
  a. All scripts in the /scripts folder<br/>
  b. All files in the /OtherRequiredFiles folder. You will need to unzip the zipped file.<br/>
  c. Example Canopus output from /Canopus_predictions folder. Not required if you have your own predictions<br/>
  d. Pathway representation, either /PlantCyc or /SolCyc <br/>  
  ```bash
   mkdir _workingDirectory
   cp scripts/*.py OtherRequiredFiles/* Canopus_predictions/*.tsv PlantCyc/* _workingDirectory/
   cd _workingDirectory
   unzip ChEBI_126_classyfire_21_annotations.csv.zip   
 ```     
  

__NOTE: You need to do Step 2 only if you are working with a new PMN pathway model. For PlantCyc, BrachypodiumCyc and SolCyc, the necessary output files are already provided in the respective database folder. You need to do this only once for each new PMN pathway you work with.__

2. __(OPTIONAL)__ First process the Pathway database files to get them into the right format using the wrapper script. This script will create a bunch of intermediate files and pathway associations, some of which are needed by the next script. Running the script without any options will give you information about these options. <br/>
```bash
   python get_similar_pathways_wrapper.py -cmp compounds.dat -pwy pathways.dat -extra no -chemont ChemOnt_2_1.obo.mod.tax.organic -chebi ChEBI_126_classyfire_21_annotations.csv.mod
 ```      
3. You can view the .network file in Cytoscape if you want (Import from File...), to see pathway clusters<br/>

4. Run the metaPathwayMap script as follows. Running the script without any options will give you information about these options. <br/>
```bash
python metaPathwayMap.py -pwy pathways.dat.cid.all -canopus massbank_compound_canopus2.tab -jaccard 0.6
```

5. The script outputs status messages that tell you which file is the most relevant (*.top.format.abbr.tab)<br/>

## Explanation of script functions and output files
There are multiple output files produced through these scripts, including intermediate files showing how the data was processed. The following files are most consequential for the user

### get_similar_pathways_wrapper.py
Runs step1-4 one after the other. You DO NOT need to run individual scripts listed below. Skip directly to metaPathwayMap.py.

### get_similar_pathways_step1_getCompoundIDs.py
This script performs the following steps:
* Extracts Chebi IDs of compounds from the pathways.dat and compounds.dat files
* Produces an output with the extension ***cid***. If the option ***extra=yes*** is specified, then SMILES of compounds without ChEBI IDs are extracted into files with extensions ***nochebi.****

### get_similar_pathways_step2_getCompoundAnnot.py
This script performs the following steps:
* Reads ChemOnt annotations of individual ChEBI IDs from the ChemOnt files
* Associates these annotations with ChEBI IDs extracted from pathways files in step1
* Creates output files ****.all*** and ****.all.cpd***. The format of the ****.all*** file is as follows:

 | Column 1 | Column 2 | Column 3 | Column 4 |
 |----------|----------|----------|----------|
 | PathwayID | Pathway name | All compounds in that pathway in a Python list object | All unique ChemOntCategories in a Python list object |
 
 The compounds entries in column 3 are in the following format:
 
 ***CompoundName|CompoundChEBI-ID|ChemOnt|ChemOntName*** e.g. 'D-LACTATE|16004|CHEMONTID:0000129|Alcohols_and_polyols
 
 The format of the ****.all.cpd*** file is as follows:

 | Column 1 | Column 2 | Column 3 |
 |----------|----------|----------|
 | CompoundNameCompoundChEBI-ID | ChemOnt | All pathways the compound is present in, in a Python list object |
 

### get_similar_pathways_step3_getSimilarPathways.py
This script performs the following steps:
* Reads in the ****.all*** file produced by step2
* Computes Jaccard coefficient between all-vs-all pathway comparisons
* Prints a tab-delimited file  ****.pwy.dist*** with pairwise distances (1-Jaccard coefficient)
* Randomly resamples 1000 pathway pairs 10,000 times, and each time computes the mean Jaccard distance of the entire distribution. Based on the 10,000 mean values, obtains the 1st percentile of the distribution
* Using this threshold (boot) or another fixed threshold (fixed), filters the pwy.dist file, producing ****.pwy.dist.boot.fil*** or ****pwy.dist.fixed020.fil*** files
* Some additional ****.complex*** files are also produced, which list which all pathways are similar to a given pathway. The format of these files is:

 | Column 1 | Column 2 | Column 3 |
 |----------|----------|----------|
 | Pathway name | Number of similar pathways | List of similar pathways, based on threshold |
 
* Creates a ****.network*** file that can be used as input for Cytoscape
* Using the Python networkx module, identifies connected components
* The ****.components*** file has the following format

 | Column 1 | Column 2 | Column 3 |
 |----------|----------|----------|
 | Component-ID | Number of pathways | Pathway IDs in a Python list object |
 
### get_similar_pathways_step4_topClasses.py
This script performs the following steps:
* Reads in the ****.all*** and ****.components*** file produced by previous scripts
* Identifies enriched ChemOnt categories in each connected component, using the fisher.exact function of the Python module Scipy.stats
* Produces file ****.enriched*** and ****.enriched.complex***. The enriched file has the following format. Only p-vlaues < 0.01 are output. No multiple testing correction is performed as of now. 

 | Column 1 | Column 2 | Column 3 | Column 4 | Column 5 | Column 6 | Column 7 | Column 8 |
 |----------|----------|----------|----------|----------|----------|----------|----------|
 | Component-ID | ChemOnt ID | Chemont+PWY+ | Chemont+PWY- | ChemOnt-PWY+ | ChemOnt-PWY- | p-value (fisher.exact, alternative="greater") | PWYs with ChemOnt ID |


### metaPathwayMap.py
This script performs the following steps:
* Reads in output files produced by get_similar_pathways
* For CANOPUS predictions, only those that have at least three main class annotations are analyzed
* The script looks for every pathway compound where at least three of the main classes match. This limits the search substantially, potentially causing false positives and negatives
* Calculates Jaccard Coefficient (Tanimoto coefficient) between ALL (main+alternative) classes of the CANOPUS prediction and the pathway compound
* The pathways of the best matches are extracted
* All related pathways (network components) of each of the pathways are extracted
* All hits are outputted in various files.
* The output file ****.metaPathwayMap*** is the final file to look at. It is also the most condensed, and if the users desire additional information and/or a view of additional complexity, the ****.chemont.top.format*** can be looked at
* The components column refers to the component in the ****.boot.fil.components*** file created by get_similar_pathways
