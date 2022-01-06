# metaPathwayMap
Associating chemical class predictions with MetaCyc pathway models

### Requirements:
Unix<br/>
Python3<br/>
Following python modules: Numpy, Scipy, Pandas, Networkx<br/>

### Instructions for use:
1. Download the following files into the working directory:<br/>
  a. All scripts in the /scripts folder<br/>
  b. All files in the /OtherRequiredFiles folder. You will need to unzip the zipped file.<br/>
  c. Example Canopus output from /Canopus_predictions folder. Not required if you have your own predictions<br/>
  d. Pathway representation, either /PlantCyc or /SolCyc (only PlantCyc available for now)<br/>

2. First process the Pathway database files to get them into the right format using the wrapper script. This script will create a bunch of intermediate files and pathway associations, some of which are needed by the next script. Running the script without any options will give you information about these options. <br/>
```bash
   python get_similar_pathways_wrapper.py -cmp compounds.dat -pwy pathways.dat -extra no -chemont ChemOnt_2_1.obo.mod.tax.organic -chebi ChEBI_126_classyfire_21_annotations.csv.mod
 ```
      
3. You can view the .network file in cytoscape if you want, to see pathway clusters<br/>

4. Run the metaPathwayMap script as follows. Running the script without any options will give you information about these options. <br/>
```bash
python metaPathwayMap.py -pwy pathways.dat.cid.all -canopus massbank_compound_canopus2.tab -jaccard 0.7
```

5. The script outputs status messages that tell you which file is the most relevant (*.top.format.abbr.tab)<br/>
