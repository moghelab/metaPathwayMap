from __future__ import division
import sys, operator, os

class step1:
    def __init__(self):
        pass

    def getCompoundIDs(self,name1,name2,optional,mylog):
        print ("Step 1: Extracting compound info from Cyc files...")
        mylog.write ("Step 1: Extracting compound info from Cyc files...")

        #Process the compounds file
        file1=open(name1,'r',encoding='cp1252')
        line1=file1.readline()
        dict1={}; dict2={}; dict3={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                sp=line1.strip().split(' - ')
                if line1.startswith('UNIQUE-ID'):
                    cid=sp[1]            
                elif line1.startswith('DBLINKS'):
                    sp2=sp[1].split()                                
                    for i in range(0,len(sp2)):
                        item=sp2[i]
                        if ('CHEBI') in item:                    
                            chebi='CHEBI|{}'.format(sp2[i+1].split('"')[1])
                            if cid not in dict1:
                                dict1[cid]=[chebi]
                            else:
                                if chebi not in dict1[cid]:
                                    dict1[cid].append(chebi)
                                #print ("CID repeat:  ", cid, chebi, dict1[cid])
                                
                elif line1.startswith('SMILES'):
                    smiles=sp[1]
                    if cid not in dict2:
                        dict2[cid]=smiles
                    else:
                        print ("Repeat CID smiles: {}\n".format(cid))
                        mylog.write ("Repeat CID smiles: {}\n".format(cid))
                elif line1.startswith('INCHI-KEY'):
                    inchikey=sp[1]
                    if cid not in dict3:
                        dict3[cid]=inchikey
                    else:
                        print ("Repeat CID smiles: {}\n".format(cid))
                        mylog.write("Repeat CID smiles: {}\n".format(cid))
            line1=file1.readline()
        file1.close()

        #Process the pathway file
        file1=open(name2,'r',encoding='cp1252')
        out1=open(name2+".cid",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        if optional=='yes':
            out2=open(name2+".nochebi.smiles",'w') #Extra files
            out3=open(name2+".nochebi.inchikey",'w') #Extra files
            out4=open(name2+".chebi.smiles",'w')
            out5=open(name2+".chebi.smiles.2col.tsv",'w')
            out2.write('#CID\tPathway\tSMILES\n')
            out2.write('#CID\tPathway\tInchiKey\n')
        nochebi=0; yeschebi=0; donesmiles=[]; doneinchikey=[]
        line1=file1.readline()
        noids=['RXN','L2R','R2L','-PRIMARIES','DIRECTION','']
        while line1:
            if line1.startswith('#'):
                pass
            else:
                sp=line1.strip().split(' - ')
                #pwy="NA"; tps="NA"
                if line1.startswith('UNIQUE-ID'):
                    pwy=sp[1]; elist=[]; clist=[]
                elif line1.startswith('COMMON-NAME'):
                    tps=sp[1].replace(' ','_')
                elif line1.startswith('REACTION-LAYOUT'):
                    #Modify spaces and brackets to ^, since couldn't figure
                    #out how to use re.split function
                    spx=sp[1].replace(' ','^')
                    spx=spx.replace('(','^')
                    spx=spx.replace(')','^')            
                    sp2=spx.split('^')

                    #Get ChEBI ID
                    for item in sp2:
                        flag=0
                        if item in dict1:
                            chebilist=dict1[item]
                            for chebi in chebilist:
                                str1='{}|{}'.format(item,chebi)
                                if str1 not in clist:
                                    clist.append(str1)
                                    yeschebi+=1; flag=1
                                    
                            if flag==0 and optional=='yes':
                                if item in dict2:
                                    smiles=dict2[item]
                                    if smiles not in donesmiles:
                                        out4.write('{}\t{}\t{}\n'.format(item,pwy,smiles))
                                        out5.write('{}\t{}\n'.format(item,smiles))
                                        donesmiles.append(smiles)
                        else:
                            if optional=='yes':
                                if item.startswith('CPD'):                                                
                                    sp3=item.split(')')[0]                        
                                    if sp3 in dict2:
                                        smiles=dict2[sp3]
                                        if smiles not in donesmiles:
                                            nochebi+=1
                                            out2.write('{}\t{}\t{}\n'.format(sp3,pwy,smiles))
                                            donesmiles.append(smiles)
                                        else:
                                            pass

                                    if sp3 in dict3:
                                        inchikey=dict3[sp3].split('Key=')[1]
                                        if inchikey not in doneinchikey:                                    
                                            out3.write('{}\t{}\t{}\n'.format(sp3,pwy,inchikey))
                                            doneinchikey.append(inchikey)
                                        else:
                                            pass                        
                                else:
                                    pass
                elif line1.startswith('CITATIONS'):
                    sp2=sp[1].split(':')
                    for item in sp2:
                        if item.startswith('EV') and item not in elist:
                            elist.append(item)
                elif line1.startswith('//'):
                    out1.write('{}\t{}\t{}\t{}\n'.format(pwy,tps,clist,elist))
                else:
                    pass
            line1=file1.readline()
        file1.close(); out1.close()
        
        print ("Compounds with CHEBI ID in Cyc file: ", yeschebi)
        mylog.write("Compounds with CHEBI ID in Cyc file: {}\n".format(yeschebi))
        
        if optional=='yes':
            print ("Compounds with no CHEBI ID in Cyc file: ", nochebi)
            mylog.write("Compounds with no CHEBI ID in Cyc file: {}".format(nochebi))
            out2.close(); out3.close(); out4.close(); out5.close()
            os.system('cut -f 1,3 {}.nochebi.smiles > {}.nochebi.smiles.2col.tsv'. \
                      format(sys.argv[2],sys.argv[2]))
            os.system('cut -f 3 {}.nochebi.inchikey > {}.nochebi.inchikey.list'. \
                      format(sys.argv[2],sys.argv[2]))
        print ("Step 1 completed."); mylog.write("Step 1 completed.\n")

if __name__ == '__main__':    
    print ("INP1: compounds.dat")
    print ("INP2: pathways.dat")
    print ("INP3: Do you want Inchikey/SMILES of compounds without ChEBI ID (yes/no)")
    print ("INP4: full log file name (name.log)")
    name1=sys.argv[1]; name2=sys.argv[2],opt=sys.argv[3];mlog=sys.argv[4]
    step1.getCompoundIDs(name1,name2,opt,mlog)
    print ("Done!")
        
