import sys, os

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
    print ("-log:   Name for log file e.g. plantcyc.log")

#Get user inputs
f1=""; f2=""; f3=""; chemont=""; chebi=""; mainlog=""
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
    if sys.argv[i]=='-log':
        logn=sys.argv[i+1]
        mainlog=open(logn,'w')

#Run steps
if f1=="" or f2=="" or chemont=="" or chebi=="" or mainlog=="":
    help1()
    print ("ERROR: Please provide all parameters")
    sys.exit()
else:
    if opt=="":
        opt='no'

    #Get working directory
    pwd=os.getcwd()
    mainlog.write('WORKING DIRECTORY: {}\n'.format(pwd))
    mainlog.write('python {}\n'.format(' '.join(sys.argv)))
    print ("#####")
    
    #Parse pathways    
    step1.getCompoundIDs(f1,f2,opt,mainlog)    
    print ("#####"); mainlog.write("#####\n")
    
    #Get compound annotations
    f3=f2+".cid"
    step2.getCompoundAnnot(chemont,chebi,f3,mainlog)
    print ("#####"); mainlog.write("#####\n")

    #Get similar pathways
    f3=f2+".cid.all"
    step3.getSimilarPwys(f3,mainlog)
    print ("#####"); mainlog.write("#####\n")

    #Connect metabolomics compounds to pathways
    f4=f3+".pwy.dist.boot.fil.components"
    step4.mostFrequentClass(chemont,f3,f4,mainlog)
    print ("####"); mainlog.write("#####\n")

    print ("All steps done!")
    mainlog.write("All steps done!\n")
    mainlog.close()
    print ("Check log file {} for job log".format(logn))



