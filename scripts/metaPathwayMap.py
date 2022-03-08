#!/usr/bin/python
from __future__ import division
import sys, os
import numpy as np

def help1():
    print ("-pwy:       pathways.dat.cid.all (output of getSimilarPathways)")
    print ("-canopus:   CANOPUS annotations of the compounds of interest (adducts file)")
    print ("-jaccard:   Threshold of Jaccard to use for output predictions (DEFAULT: 0.7)")
    print ("            0.6 should be the lowest. Based on our observations, predictions above 0.7 are " \
                        "generally trustworthy.")
    print ("-log:       Log file name")
    
#Get user inputs
f1=""; canopus=""; logn=""; thresh=0.7
for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-pwy':
        f1=sys.argv[i+1]
        f2=f1+".pwy.dist.boot.fil.complex"
        f3=f1+".pwy.dist.boot.fil.components"
    if sys.argv[i]=='-canopus':
        canopus=sys.argv[i+1]
    if sys.argv[i]=='-jaccard':
        thresh=float(sys.argv[i+1])
    if sys.argv[i]=='-log':
        logn=sys.argv[i+1]

def writemsg(mlist,outfile):
    for msg in mlist:
        print (msg)
        outfile.write('{}\n'.format(msg))
    return outfile
    
if f1=="" or canopus=="" or logn=="":
    help1()
else:
    file1=open(f1, 'r') #PlantCyc compounds
    logfile=open(logn, 'w')
    logfile.write('python {}\n'.format(' '.join(sys.argv)))
    
    user=3 #structural neighbors to pull out. Default: 3 is sufficient
    line1=file1.readline()
    dict1={}; dict2={}
    while line1:
        if line1.startswith('#'):
            pass
        else:
            tab1=line1.strip().split('\t')
            pwy=tab1[0]; name=tab1[1]; cpds=eval(tab1[2])
            str1='{}|{}'.format(pwy,name)
            for cpd in cpds:
                sp=cpd.split('|')
                xname=sp[0]; id1=sp[1]
                chemont=sp[3]
                cname='|'.join([xname,id1])
                if cname not in dict1:
                    dict1[cname]=[chemont]
                else:
                    if chemont not in dict1[cname]:
                        dict1[cname].append(chemont)

                if cname not in dict2:
                    dict2[cname]=[str1]
                else:
                    if str1 not in dict2[cname]:
                        dict2[cname].append(str1)
        line1=file1.readline()
    file1.close()

    #Write to out
    out1=open(f1+".cpd",'w')
    out1.write('#python {}\n'.format(' '.join(sys.argv)))
    for cpd in dict1:
        pwy=dict2[cpd]; chemont=dict1[cpd]
        for ont in chemont:
            out1.write('{}\t{}\t{}\n'.format(cpd,ont,pwy))
    out1.close()

    #Read in the Complexity file
    file1=open(f2, 'r')
    line1=file1.readline()
    dict3={}
    while line1:
        if line1.startswith('#'):
            pass
        else:
            tab1=line1.strip().split('\t')
            pwy=tab1[0]; others=eval(tab1[2])
            dict3[pwy]=others
        line1=file1.readline()
    file1.close()

    #Get the component names for each pathway
    file1=open(f3, 'r') 
    line1=file1.readline()
    compdict={}
    while line1:
        if line1.startswith('#'):
            pass
        else:
            tab1=line1.strip().split('\t')
            comp=tab1[0]; pwys=eval(tab1[2])
            for pwy in pwys:
                if pwy not in compdict:
                    compdict[pwy]=comp
                else:
                    print ("Pathway repeat: ", pwy)
        line1=file1.readline()
    file1.close()

    #Open CANOPUS annotations
    file1=open(canopus, 'r')
    out1=open(canopus+".chemont",'w')
    out1.write('#python {}\n'.format(' '.join(sys.argv)))
    out1.write('#ExptID\tExpt_Formula\tAdduct\tMetaboliteMatch\tJaccard\tMetabolitesPathway\tNetwork_Component\t' \
               'RelatedPathways\tCANOPUS_mainPred\tChemOnt\n')

    out11=open(canopus+".chemont.abbr.tab",'w')
    out11.write('#python {}\n'.format(' '.join(sys.argv)))
    out11.write('#ExptID\tCANOPUS_MostSpecific\tMetaboliteMatch\tJaccard\tMetabolitesPathway\tNetwork_Component\n')

    out2=open(canopus+".chemont.top",'w')
    out2.write('#python {}\n'.format(' '.join(sys.argv)))
    out2.write('#ExptID\tExpt_Formula\tAdduct\tMetaboliteMatch\tJaccard\tNetwork_Component\t' \
               'RelatedPathways\tCANOPUS_mainPreds\tCANOPUS_allPreds\tChemOnt_dbStruc\n')

    out3=open(canopus+".chemont.top.format",'w')
    out3.write('#python {}\n'.format(' '.join(sys.argv)))
    out3.write('#ExptID\tExpt_Formula\tAdduct\tMetaboliteMatch\tJaccard\t' \
               'MetabolitesPathway\tNetwork_Component\tCANOPUS_mainPreds\tCANOPUS_allPreds\tChemOnt_dbStruc\n')
    out31=open(canopus+".chemont.top.format.abbr.tab",'w')
    out31.write('#python {}\n'.format(' '.join(sys.argv)))
    out31.write('#ExptID\tCANOPUS_MostSpecific\tMetaboliteMatch\tJaccard\tMetabolitesPathway\n')

    line1=file1.readline()
    m=0; n=0; z=0; idict={}; cdict={}; toppred={}; xdict={}
    removelist=['Chemical_entities', 'Chemical_entities"', '"Chemical_entities',
                'Organic_compounds', '"Organic_compounds', 'Organic_compounds"']    
    while line1:
        if line1.startswith('#') or line1.startswith('name'):
            pass
        else:
            m+=1
            line1=line1.replace('"','')
            tab1=line1.strip().split('\t')
            #print (tab1)
            xid=tab1[0]; form=tab1[1]; adduct=tab1[2]; mainclasses=tab1[4:9]; altclasses=tab1[9]
            #If same compound multiple times, rename id
            if xid not in xdict:
                xdict[xid]=1
                id1='{}-1'.format(xid)                
            else:
                xlen=xdict[xid]+1
                xdict[xid]=xlen
                id1='{}-{}'.format(xid, xlen)
                
            mostsp=tab1[4].replace(' ','_')
            if mostsp=='':
                mostsp='NA'
            #print (mainclasses)
            #print (altclasses)
            #sys.exit()

            #Process names of the altclasses
            annot1=[x.replace(' ','_') for x in altclasses.split('; ')]        
            for item in removelist:
                if item in annot1:
                    annot1.remove(item)
            
            #Make a separate list of just the main classes
            annot=[]
            for x in mainclasses:
                if x!='None' and x!='':
                    annot.append(x.replace(' ','_'))
            
            #At least three class predictions for the CANOPUS peak
            #because frequently there is 'None' in the class/subclass/level5 column
            if len(annot)>=3: 
                n+=1
                if m%500==0:
                    print ("Reading line: ", m, " cpd: ", n)

                cdict={}
                for cpd in dict1: #PlantCyc compounds
                    #Get chemont classification for each PlantCyc compound
                    #and see if the categories match
                    chemont=dict1[cpd]                                
                    y=0

                    #Narrow down region based on main classes
                    for ont in annot:
                        if ont in chemont:
                            y+=1
                    
                    
                    if abs(y-len(annot))<=3: #At least three good matches                        
                        #Get info about pathways
                        pwy=dict2[cpd]
                        idict[id1]=1

                        #Calculate Tanimoto based on all annotations
                        numerator = len(set(annot1) & set(chemont))
                        denominator = len(annot1)+len(chemont)-numerator
                        overlap=numerator/min([len(annot1), len(chemont)])
                        tcoeff = numerator/denominator                        
                        tanimoto = '{:2f}'.format(tcoeff)

                        #Get related pathways
                        for pwyg in pwy:
                            relpwys=[]
                            #pname=pwyg.split('|')[0]
                            pname=pwyg
                            if pname in dict3:
                                relpwys=dict3[pname]

                            #get component names
                            compt=compdict[pwyg]

                            #write to out                            
                            out1.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'. \
                                       format(id1, adduct, form, cpd, tanimoto, pwyg, compt, relpwys, annot, annot1, chemont))
                            out11.write('{}\t{}\t{}\t{:2f}\t{}\t{}\n'. \
                                       format(id1, mostsp, cpd, tcoeff, pwyg, compt))
                        z+=1

                        #Compile all predictions
                        if tcoeff not in cdict:
                            cdict[tcoeff]=[cpd]
                        else:
                            if cpd not in cdict[tcoeff]:
                                cdict[tcoeff].append(cpd)

                #Get top pathway associations for each user-provided compound
                #how many are retrieved depends on user selection on command line (set to default 3)
                if cdict!={}:      
                    #Get the top tanimoto coefficients
                    flist=[]
                    try:
                        flist0=sorted(list(cdict.keys()), reverse=True)[0:user]
                    except:
                        flist0=sorted(list(cdict.keys()), reverse=True)

                    #get tanimotos greater than user-specified threshold
                    for fx in flist0:
                        if fx>=thresh:
                            flist.append(fx)

                    #get pathway predictions for those tanimotos
                    for f in flist:
                        clist=cdict[f]
                        for cpdx in clist:
                            chemontx=dict1[cpdx]; pwyx=dict2[cpdx]

                            #Get related pwys                        
                            for pwyn in pwyx:
                                relpwys=[]
                                #pname=pwyn.split('|')[0]
                                pname=pwyn
                                if pname in dict3:
                                    relpwys=dict3[pname]

                                #get component names
                                compt=compdict[pwyn]
                                
                                out2.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'. \
                                           format(id1, adduct, form, cpdx, f, pwyn, compt, relpwys, annot, annot1, chemontx))
                                toppred[id1]=1

                                #Format for cleaner output
                                jrel=';'.join(relpwys); jannot=';'.join(annot)
                                jannot1=';'.join(annot1); jchemontx=';'.join(chemontx)
                                out3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'. \
                                           format(id1, adduct, form, cpdx, f, pwyn, compt, jrel, jannot, jannot1, jchemontx))
                                out31.write('{}\t{}\t{}\t{:2f}\t{}\t{}\n'. \
                                       format(id1, mostsp, cpdx, f, pwyn,compt))                                
        line1=file1.readline()
    file1.close(); out1.close(); out2.close(); out11.close(); out3.close(); out31.close()

    os.system(f"sort -k5nr -t$'\t' {canopus}.chemont.top >  {canopus}.chemont.top.sort")

    os.system(f'cp {canopus}.chemont.top.format.abbr.tab {canopus}.metaPathwayMap')

    print ("~~~~")
    print ("# of lines in CANOPUS file: ", m)
    print ("# of lines with >=3 CANOPUS levels: ", n)
    print ("# of compound match lines in OUT: ", z)
    print ("# of CANOPUS compounds with matching pathway annotations: ", len(idict.keys()))
    print ("# of CANOPUS compounds with high-conf pathway annotations: ", len(toppred.keys()))
    print ("~~~~")
    print (f"Easiest, most relevant output file to look at --> {canopus}.chemont.top.format.abbr.tab")
    print (f"Easiest, most relevant output file to look at (copy of above file) --> {canopus}.metaPathwayMap\n")
    print (f"Other files with excruciating details --> {canopus}.*")

    logfile.write ("~~~~\n")
    logfile.write (f"# of lines in CANOPUS file: {m}\n")
    logfile.write (f"# of lines with >=3 CANOPUS levels: {n}\n")
    logfile.write (f"# of compound match lines in OUT: {z}\n")
    logfile.write (f"# of CANOPUS compounds with matching pathway annotations: {len(idict.keys())}\n")
    logfile.write (f"# of CANOPUS compounds with high-conf pathway annotations (>threshold): {len(toppred.keys())}\n")
    logfile.write ("~~~~\n")
    logfile.write (f"Easiest, most relevant output file to look at --> {canopus}.chemont.top.format.abbr.tab\n")
    logfile.write (f"Easiest, most relevant output file to look at (copy of above file) --> {canopus}.metaPathwayMap\n")
    logfile.write (f"Other files with excruciating details --> {canopus}.*\n")

    

    print ("Done!")
    logfile.write ("Done!\n")
    logfile.close()

    
    
    
