import sys

#Import accessory scripts
import get_similar_pathways_step1_getCompoundIDs
import get_similar_pathways_step2_getCompoundAnnot
import get_similar_pathways_step3_getSimilarPathways
import get_similar_pathways_step4_topClasses

#Set up steps
step1=get_similar_pathways_step1_getCompoundIDs.step1()
step2=get_similar_pathways_step2_getCompoundAnnot.step2()
step3=get_similar_pathways_step3_getSimilarPathways.step3()
step4=get_similar_pathways_step4_topClasses.step4()

def help1():
    print ("-cmp:       compounds.dat")
    print ("-pwy:       pathways.dat")
    print ("-extra:     Do you want Inchikey/SMILES of compounds without ChEBI ID (yes/no)")
    print ("-chemont:   ChemOnt obo.mod.tax.organic")
    print ("-chebi:     ChEBI_126_classyfire_21_annotations.csv.mod")
    print ("-canopus:   CANOPUS annotations of the compounds of interest")

#Get user inputs
for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-cmp':
        f1=sys.argv[i+1]
    if sys.argv[i]=='-pwy':
        f2=sys.argv[i+1]
    if sys.argv[i]=='-extra':
        opt=sys.argv[i+1]
    if sys.argv[i]=='-chemont':
        chemont=sys.argv[i+1]
    if sys.argv[i]=='-chebi':
        chebi=sys.argv[i+1]
    if sys.argv[i]=='-canopus':
        canopus=sys.argv[i+1]    

#Run steps
if f1=="" or f2=="" or chemont=="" or chebi=="":
    help1()
else:
    if opt=="":
        opt='no'
    print ("#####")
    #Parse pathways    
    step1.getCompoundIDs(f1,f2,opt)
    print ("#####")
    
    #Get compound annotations
    f3=f2+".cid"
    step2.getCompoundAnnot(chemont,chebi,f3)
    print ("#####")

    #Get similar pathways
    f3=f2+".cid.all"
    step3.getSimilarPwys(f3)
    print ("#####")

    #Connect metabolomics compounds to pathways
    f4=f3+".pwy.dist.boot.fil.components"
    step4.mostFrequentClass(chemont,f3,f4)
    print ("####")
    


print ("All steps done!")
