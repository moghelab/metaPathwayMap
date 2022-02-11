from __future__ import division
import sys, os, random
import numpy as np
from scipy.stats import fisher_exact
from collections import Counter

'''
From network components output of step3, get most frequent Chemont Class per component
Modify the *.components file to add most frequent class/classes in Index 2
'''

class step4:
    def __init__(self):
        pass

    def most_frequent(mylist):
        return max(set(mylist), key = mylist.count)

    def getEnriched(d1, d2, d3, out, dlog):
        allonts=list(d1.keys())
        allpwys=list(d2.keys())
        allcomp=list(d3.keys())

        #Count occurrence of each ont in each pwy
        myd={}; o=0
        for ont in allonts:
            o+=1
            if o%100==0:
                print ("Ontology categories tested: ", o)
                dlog.write ("Ontology categories tested: {}\n".format(o))
            for compx in allcomp:
                ont1=0; ont2=0; ont3=0; ont4=0
                for pwyx in allpwys:
                    ponts=d2[pwyx]; ontcount=ponts.count(ont) #Occurrence of ont in the pathway
                    noncount=len(ponts)-ontcount
                    if pwyx in d3[compx]:
                        ont1+=ontcount; ont2+=noncount
                    else:
                        ont3+=ontcount; ont4+=noncount

                #Perform fisher exact test for each pathway
                if ont1>5:
                    table=np.array([[ont1,ont2], [ont3,ont4]])
                    oddsr, p=fisher_exact(table, alternative='greater')
                    cpwy=d3[compx]
                    if p<0.01:
                        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'. \
                                   format(compx,ont,ont1,ont2,ont3,ont4,p,cpwy))
                        sont='{}|{}'.format(ont,p)
                        if compx not in myd:
                            myd[compx]=[sont]
                        else:
                            if sont not in myd[compx]:
                                myd[compx].append(sont)
        return out, myd, dlog

    def mostFrequentClass(self,f1,f2,f3,mylog):
        print ("Step 4: Getting enriched Ontologies for each network component...")
        mylog.write ("Step 4: Getting enriched Ontologies for each network component...\n")

        #First get all the Class annotations
        file1=open(f1,'r')
        line1=file1.readline()
        dict1={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                co=tab1[0]; lev1=tab1[1]
                if lev1=='Subclass':
                    dict1[co]=1
            line1=file1.readline()
        file1.close()

        #List of excluded Classes
        excl=['CHEMONTID:0000323|Organooxygen_compounds', 'CHEMONTID:0003940|Organic_oxides',
              'CHEMONTID:0000278|Organonitrogen_compounds', 'CHEMONTID:0000265|Carboxylic_acids_and_derivatives',
              'CHEMONTID:0004139|Azacyclic_compounds', 'CHEMONTID:0004140|Oxacyclic_compounds',
              'CHEMONTID:0002279|Benzene_and_substituted_derivatives', 'CHEMONTID:0004144|Heteroaromatic_compounds']

        #Second, get all Chemonts of each pathway
        file1=open(f2,'r') #pathway.dat.cid.all
        line1=file1.readline()
        dict2={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                pid=tab1[0]; pname=tab1[1]; str1='{}|{}'.format(pid,pname)
                subs=eval(tab1[2]); slist=[]
                for sub in subs:
                    sp=sub.split('|')
                    co='|'.join([sp[2],sp[3]])
                    if co in dict1:
                        slist.append(co)
                dict2[str1]=slist
            line1=file1.readline()
        file1.close()

        #Find representative name for each pathway Component in INP2
        file1=open(f3,'r') #components file
        out1=open(f3+".enriched",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        line1=file1.readline()
        dict3={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                comp=tab1[0]; ncomp=tab1[1]; plist=eval(tab1[2])
                dict3[comp]=plist
            line1=file1.readline()
        file1.close()

        #Calculate enrichment
        out1,cdict,mylog=step4.getEnriched(dict1,dict2,dict3,out1,mylog)
        out1.close()

        #Write Classes to another output
        file1=open(f3,'r') #components file
        out1=open(f3+".enriched.complex",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        line1=file1.readline()
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                comp=tab1[0]; ncomp=tab1[1]; plist=eval(tab1[2])
                if comp in cdict:
                    olist=cdict[comp]
                else:
                    olist=[]
                out1.write('{}\t{}\t{}\t{}\n'.format(comp,ncomp,olist,plist))
            line1=file1.readline()
        file1.close(); out1.close()
        print ("Most enriched class written for each component!")
        mylog.write ("Most enriched class written for each component!\n")        
        
        

if __name__ == '__main__':    
    print ("INP1: ChemOnt_2_1.obo.mod.tax.organic")
    print ("INP2: pathways.dat.cid.all")
    print ("INP3: pathways.dat.cid.all.pwy.dist.boot.fil.components")
    print ("INP4: name of the log file (e.g. name.log)")
    name1=sys.argv[1];  pfile=sys.argv[2]
    components=sys.argv[3]; mlog=sys.argv[4]
    mostFrequentClass(name1,pfile,components,mlog)
    
    print ("Done!")
        
