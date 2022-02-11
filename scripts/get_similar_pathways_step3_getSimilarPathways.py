from __future__ import division
import sys, os, random
import numpy as np
import pandas as pd
from collections import Counter
import networkx as nx

class step3:
    def __init__(self):
        pass

    def makelog(msgx,outfile,val):
        print ('{}\n'.format(msgx))
        outfile.write('{}\n'.format(msgx))

    def bootstrap(mylist, func, level, resample, batch):
        vlist=[]
        for i in range(resample):
            xlist=random.sample(mylist, batch)
            myset = np.percentile(xlist, level)        
            vlist.append(myset)    
        print ("Mean bootstrap: ", np.mean(vlist))
        return func(vlist)

    def most_common(lst):
        return max(set(lst), key=lst.count)

    def getCatCounts(mylist2): #e.g. 'D-LACTATE|16004|CHEMONTID:0000129|Alcohols_and_polyols'
        myclist=[]
        for item in mylist2:
            sp=item.split('|')
            cat='{}|{}'.format(sp[0],sp[1])
            myclist.append(cat)
        maxitem=max(set(myclist), key=myclist.count)
        print (myclist)
        print (Counter(myclist))        
        print ("####")
        print (maxitem)
        sys.exit()

    def getSimilarPwys(self,name1,mylog):
        print ("Step 3: Estimating substrate similarities between pathways...")
        file1=open(name1,'r') #Output of step2
        line1=file1.readline()
        dict1={}
        while line1:
            if line1.startswith('#'):
                pass
            else:
                tab1=line1.strip().split('\t')
                pwyid=tab1[0]; name=tab1[1]; onts=eval(tab1[3])
                pwy='{}|{}'.format(pwyid,name)
                if pwy not in dict1:
                    dict1[pwy]=onts
                else:
                    print ("Repeat: ", pwy); mylog.write("Repeat: {}\n".format(pwy))
            line1=file1.readline()
        file1.close()

        out1=open(name1+".pwy.dist",'w')
        #out1.write('#python {}\n'.format(' '.join(sys.argv)))
        out1.write('#PWY1\tPWY2\t1-Jaccard\tChemonts_PWY1\tChemonts_PWY2\tUnion\tIntersect\n')
        out2=open(name1+".pwy.dist.details",'w')
        out2.write('#python {}\n'.format(' '.join(sys.argv)))
        out2.write('#PWY1\tPWY2\t1-Jaccard\tChemonts_PWY1\tChemonts_PWY2\tUnion\tIntersect\n')
        plist=list(dict1.keys())
        m=0; slist=[]
        for i in range(0,len(plist)):
            for j in range(0,len(plist)):
                list1=dict1[plist[i]]#; com1=  step3.getCatCounts(list1) 
                list2=dict1[plist[j]]#; com2=  step3.getCatCounts(list2)

                #Calculate Jaccard coefficient
                #For binary values, Tanimoto becomes Jaccard
                intersect = len(set(list1) & set(list2))
                union = len(list1) + len(list2) - intersect                
                #print (len(list1), len(list2), union, intersect)

                #Calculate cosine score
                #(produces same results as Jaccard, so not used)
                '''
                if len(list1)!=0 and len(list2)!=0:
                    cosine=(intersect/np.sqrt(len(list1)*len(list2)))
                else:
                    cosine=0
                '''
                if union>0:
                    jaccard=(intersect/union)                
                    dist = 1-jaccard
                    #dist=1-cosine
                    
                    #Write to out
                    out1.write('{}\t{}\t{:.4f}\t{}\t{}\t{}\t{}\n'. \
                               format(plist[i], plist[j], dist, len(list1), len(list2), \
                                      union, intersect))
                    out2.write('{}\t{}\t{:.4f}\t{}\t{}\n'. \
                               format(plist[i], plist[j], dist, list1, list2))
                    slist.append(dist)
                else:
                    pass

                
                m+=1
                if m%50000==0:
                    print ("Comparisons completed: ", m)
                    mylog.write ("Comparisons completed: {}\n".format(m))
        print ("Total comparisons (lines in OUT): ", m)        
        mylog.write ("Total comparisons (lines in OUT): {}\n".format(m))
        file1.close(); out1.close(); out2.close()

        #Filter first file by 10th percentile of random bootstrap threshold
        print ("Now performing bootstrapping to determine threshold...")
        mylog.write ("Now performing bootstrapping to determine threshold...\n")
        thresh=step3.bootstrap(slist, np.mean, 1, 10000, 1000)
        
        #thresh=0.42347 #5th percentile
        #thresh=0.258 #1st percentile
        print ("Filtering threshold = ", thresh)
        print ("Filtering output file...")
        mylog.write ("Filtering threshold = {}\n".format(thresh))
        mylog.write ("Filtering output file...\n")

        #Open and make necessary files
        file1=open(name1+".pwy.dist",'r')
        out1=open(name1+".pwy.dist.boot.fil",'w')
        out1.write('#threshold={}\n'.format(thresh))
        out11=open(name1+".pwy.dist.boot.fil.network",'w')
        out11.write('PWY1\tPWY2\tDist\n')

        #Fixed 0.20 threshold
        out2=open(name1+".pwy.dist.fixed020.fil",'w')
        out2.write('#threshold=0.20\n')
        out21=open(name1+".pwy.dist.fixed020.fil.network",'w')
        out21.write('PWY1\tPWY2\tDist\n')

        #Complex files
        out3=open(name1+".pwy.dist.boot.fil.complex",'w')
        out4=open(name1+".pwy.dist.fixed020.fil.complex",'w')

        #Create empty graph
        G=nx.Graph()
        out5=open(name1+".pwy.dist.boot.fil.components",'w')

        line1=file1.readline()
        m=0; n=0; n2=0; cdict={}; cdict2={}
        while line1:
            if line1.startswith('#'):
                out1.write(line1)                
            else:
                tab1=line1.strip().split('\t')
                p1=tab1[0]; p2=tab1[1]; val=tab1[2]
                if p1 not in cdict:
                    cdict[p1]=[p1]
                    cdict2[p1]=[p1]
                    
                if float(tab1[2])<=thresh: #Bootstrap threshold                   
                    out1.write(line1)                    
                    out11.write('{}\t{}\t{}\n'.format(p1,p2,val))
                    G.add_edge(p1,p2,length=val) #Add to network to get connected components
                    n+=1                    
                    if p2 not in cdict[p1]:
                        cdict[p1].append(p2)

                if float(tab1[2])<=0.20:                    
                    out2.write(line1)                    
                    out21.write('{}\t{}\t{}\n'.format(p1,p2,val))
                    n2+=1                    
                    if p2 not in cdict2[p1]:
                        cdict2[p1].append(p2)
                    #out2.write('{}\t{}\n'.format(p1,p2))
                    #out3.write('{}\t{}\t{}\n'.format(p1,p2,val))
            line1=file1.readline()
        file1.close(); out1.close() ; out2.close(); out21.close(); out11.close()
        print ("Filtered lines in BOOTSTRAP OUT: ", n); mylog.write("Filtered lines in BOOTSTRAP OUT: {}\n".format(n))
        print ("Filtered lines in FIXED020 OUT: ", n2); mylog.write("Filtered lines in FIXED020 OUT: {}\n".format(n))

        #Get complexity file for threshold based on random distribution
        for pwy in cdict:
            sims=cdict[pwy]
            out3.write('{}\t{}\t{}\n'.format(pwy,len(sims),sims))
        out3.close()
        os.system('sort -k2nr {}.pwy.dist.boot.fil.complex > ' \
                  '{}.pwy.dist.boot.fil.complex.sort'.format(name1,name1))

        #Get complexity file for threshold based on fixed distance of 0.20 (~2nd percentile of random)
        for pwy in cdict2:
            sims=cdict2[pwy]
            out4.write('{}\t{}\t{}\n'.format(pwy,len(sims),sims))
        out4.close()
        os.system('sort -k2nr {}.pwy.dist.fixed020.fil.complex > ' \
                  '{}.pwy.dist.fixed020.fil.complex.sort'.format(name1,name1))

        #Get connected components
        comp=nx.connected_components(G)
        gcount=0
        for component in comp:
            gcount+=1
            nline=list(component); nlen=len(nline)
            out5.write('Component-{}\t{}\t{}\n'.format(gcount,nlen,nline))
        out5.close()
            
        print ("Check outputs with the extensions .pwy.dist.*")
        print ("# of connected components: ", gcount)
        print ("Step 3 completed.")
        mylog.write ("Check outputs with the extensions .pwy.dist.*\n")
        mylog.write ("# of connected components: {}\n".format(gcount))
        mylog.write ("Step 3 completed.\n")

if __name__ == '__main__':    
    print ("module load python/3.9.6")
    print ("INP1: Output of step2 (*.all)")
    print ("INP2: name of the log file (e.g. name.log)")    
    name1=sys.argv[1]; mlog=sys.argv[2]
    step1.getSimilarPwys(name1,mlog)
    print ("Done!")
        
