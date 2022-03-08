from __future__ import division
import sys, operator, os

class step2:
    def __init__(self):
        pass

    def add2dict(D,v1,v2):
        if v1 not in D:
            D[v1]=[v2]
        else:
            if v2 not in D[v1]:
                D[v1].append(v2)
        return D

    def getCompoundAnnot(self,name1,name2,name3,mylog):
        print ("Step 2: Associating chemont categories to compounds...")
        mylog.write ("Step 2: Associating chemont categories to compounds...\n")
        file1=open(name1,'r')
        line1=file1.readline()        
        dict1={}; dict2={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                cat=tab1[0]; tp=tab1[1]                       
                dict1=step2.add2dict(dict1,tp,cat)
                if cat not in dict2:
                    dict2[cat]=tp
            line1=file1.readline()
        file1.close()
        
        file1=open(name2,'r') #Chebi - ChemOnt - Name
        line1=file1.readline()
        dict3={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                if len(tab1)>2:
                    chebi=tab1[0]; id1=tab1[1]; name=tab1[2].replace(' ','_')
                    str1='{}|{}'.format(id1,name)
                    dict3=step2.add2dict(dict3,chebi,str1)
            line1=file1.readline()
        file1.close()

        
        file1=open(name3,'r') #MetaCyc output of Step1
        level='all'
        out1=open(name3+".{}".format(level),'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        xdict={}; ydict={}
        line1=file1.readline()
        m=0; n=0; x=0; y=0
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                pwy=tab1[0]; cpds=eval(tab1[2])
                slist=[]; clist=[]
                if cpds==[]:
                    pass
                else:
                    for cmp in cpds:
                        ydict[cmp]=1
                        sp=cmp.split('|')
                        mid=sp[0]; source=sp[1]; cid=sp[2]

                        #Check if cid has a chemont category
                        if cid in dict3:
                            chemonts=dict3[cid]

                            #get the types of the chemont categories
                            for chemont in chemonts:
                                if chemont in dict2:
                                    ctype=dict2[chemont]
                                    flag=0
                                    if level=='all':
                                        flag=1
                                    else:                            
                                        if ctype==level:
                                            flag=1
                                            
                                    if flag==1:
                                        str1='{}|{}|{}'.format(mid,cid,chemont)
                                        xdict[cid]=1
                                        #str2='{}|{}'.format(cid,chemont)
                                        if str1 not in slist:
                                            slist.append(str1); y+=1
                                        if chemont not in clist:
                                            clist.append(chemont)
                    out1.write('{}\t{}\t{}\t{}\n'.format(pwy,tab1[1],slist,clist))
                    m+=1
            n+=1
            line1=file1.readline()
        file1.close(); out1.close()

        file1=open(name3+".all",'r') #Output of previous step
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
        out1=open(name3+".all.cpd",'w')
        out1.write('#python {}\n'.format(' '.join(sys.argv)))
        for cpd in dict1:
            pwy=dict2[cpd]; chemont=dict1[cpd]
            for ont in chemont:
                out1.write('{}\t{}\t{}\n'.format(cpd,ont,pwy))
        out1.close()
        
        msg1=("# of lines in step 1 output: {}\n".format(n))
        msg2=("# of lines in step 2 output: {}\n".format(m))
        msg3=("# of compounds in step 1 output: {}\n".format(len(ydict.keys())))
        msg4=("# of compounds in step 2 output with chemont IDs: {}\n".format(len(xdict.keys())))
        msg5=("# of chemont annotations total of all compounds: {}\n".format(y))
        msg6=("Step 2 completed.\n")
        print (msg1,msg2,msg3,msg4,msg5,msg6)
        mylog.write('{}{}{}{}{}{}'.format(msg1,msg2,msg3,msg4,msg5,msg6))

if __name__ == '__main__':    
    print ("INP1: ChemOnt obo.mod.tax")
    print ("INP2: Chebi -- ChemOnt (ChEBI_126_classyfire_21_annotations.csv.mod)")
    print ("INP3: metacyc output of step1")
    print ("INP4: full log file name (name.log)")
    name1=sys.argv[1]; name2=sys.argv[2],name3=sys.argv[3];mlog=sys.argv[4]
    step1.getCompoundIDs(name1,name2,opt,mlog)
    print ("Done!")
        
