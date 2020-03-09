# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 16:51:18 2019

@author: steph
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 20:49:12 2019

@author: steph
"""
import numpy as np
import pandas as pd
import re
import subprocess
import urllib
import urllib.request
import pyiptmnet.api as api
import os
import pickle

file = input("Filename:")
subs = pd.read_csv(file,sep="\t",header=None).iloc[:,:4].values
n = len(subs)
if not subs.shape[1] == 4:
    raise TypeError("Dimension of input file is incorrect (expected 4 columns, found {0})".format(subs.shape[1])) 
feature_matrix = np.zeros((n,18))



uniprot_IDs_named = pd.read_csv("id2acc.csv",index_col=0).iloc[:,:].values

gene_phospho_file = open("gene_phosphosites.txt").readlines()
gene_phospho = {}
for gene in gene_phospho_file:
    pattern = re.match("(\S+)\t\[([^\t]*)\]",gene)
    sites = pattern.group(2)
    if sites:
        sites = sites.split(", ")
        sites = [int(n) for n in sites]
    else:
        sites = []
    gene_phospho[pattern.group(1)] = sites

#define functions
def fasta(acc,filename=None):
    url = "https://www.uniprot.org/uniprot/{0}.fasta".format(acc)
    if filename:
        urllib.request.urlretrieve(url,filename)
    else:
        return urllib.request.urlopen(urllib.request.Request(url)).read()

def fastasub(acc,pos,aa1,aa2, filename):
    outputf = open(filename,"w")
    fastafile = str(fasta(acc))
    start = fastafile.index("\\n")
    fastafile = fastafile.replace("\\n","")
    if fastafile[pos+start-1]==aa1:
        listseq = [c for c in fastafile] 
        listseq[pos+start-1]=aa2
        listseq.insert(start,"\n")
        fastafile_sub = "".join(listseq)
        outputf.write(fastafile_sub)
        outputf.close()
        return fastafile_sub
    else:
        raise ValueError("input amino acid {0} does not match accession {1}".format(aa1,fastafile[pos+start-1]))
   
def wgt_avg(x,w):
    if x==[]:
        return 0
    elif sum(w)==0:
        return 0
    else:
        return np.average(x,weights=w)
         
def minus(a,b):
    if a-b < 1:
        return 1
    else:
        return a-b
    
def plus(a,b,length):
    if a+b > length:
        return length
    else:
        return a+b
    
def acc2pdb(acc):
    url = "http://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{0}".format(acc)
    try:
        pdbfile = urllib.request.urlopen(urllib.request.Request(url)).read()
    except urllib.error.HTTPError:
        return None
    pdb = re.findall("\"pdb_id\":\"([\da-z]{4})\"",str(pdbfile))
    unp_start = re.findall("\"unp_start\":(\d+),",str(pdbfile))
    unp_end = re.findall("\"unp_end\":(\d+),",str(pdbfile))
    chain_id = re.findall("\"chain_id\":\"([A-Z])\",",str(pdbfile))
    
    array = np.column_stack((pdb,unp_start,unp_end,chain_id)).astype(object) # convert to object to allow mixed dtypes
    int_cols = [1,2]
    for i in int_cols:
        for j in range(len(array)):
            array[j,i] = int(array[j,i])
    return array

def ring(pdb_matrix):   
    m = len(pdb_matrix)
    if m == 0:
        return 0
    length = max(pdb_matrix[:,2]) #the most C-terminal part of the protein we have data for
    pdb_site = [i for i in range(m) if pdb_matrix[i,1] <= minus(position,5) and pdb_matrix[i,2] >= plus(position,5,length)]    
    pdb_filtered = pdb_matrix[pdb_site,:]
    if len(pdb_site) == 0: #if can't get 5 either side, just get the one with the longest possible shortest flanking sequence
        pdb_site = [i for i in range(m) if pdb_matrix[i,1] <= position and pdb_matrix[i,2] >= position]
        pdb_filtered = pdb_matrix[pdb_site,:]
        if len(pdb_site) == 0:
            return 0
        else:
            flank = [min(position-pdb_filtered[i,1],pdb_filtered[i,2]-position) for i in range(len(pdb_site))]
            j = np.argmax(flank)
            pdb_best = pdb_filtered[j]
    else:
        pdb_best = pdb_filtered[0] 

    pdb_id = pdb_best[0]
    chain = pdb_best[3]
    cmd1 = ['wget', '-q', 'http://www.rcsb.org/pdb/files/{0}.pdb.gz'.format(pdb_id), "-O", "structure.pdb.gz"] #-q suppresses output
    cmd2 = ['gunzip', '-rf', 'structure.pdb.gz'.format(pdb_id)] #-rf suppresses overwrite message
    cmd3 = ['./bin/Ring', '-i', '../structure.pdb'.format(pdb_id), '-c', chain] #only look at chain-internal contacts. Necessary for heterooligomers

    out = open("ring_out.txt","w")
    subprocess.run(cmd1)
    subprocess.run(cmd2)
    os.chdir('dist')
    subprocess.run(cmd3, stdout=out)
    out.close()

    os.chdir("../")
    pdb_file = open("structure.pdb").read()
    x=re.search("DBREF\s+[\dA-Z]+\s+{0}\s+(\d+)\s".format(chain),pdb_file)
    if x:
        pdb_start = int(x.group(1))
    else:
        pdb_start = pdb_best[1]
    
    ring_output = open("ring_out.txt").readlines()
    try:
        header = ring_output.index("NodeId\tChain\tPosition\tResidue\tDssp\tDegree\tBfactor_CA\tx\ty\tz\tpdbFileName\tRapdf\tTap\n")
        ring_data = ring_output[header:]
        with open("ring_data.txt","w") as ring_datafile:
            for line in ring_data:
                ring_datafile.write(line)
        ring_data = pd.read_csv("ring_data.txt",sep="\t").iloc[:,:].values
        ring_data = pd.read_csv("ring_data.txt",sep="\t").iloc[:,:].values
        residues = list(ring_data[:,2])
        residues = [i+pdb_best[1]-pdb_start for i in residues]
        if position in residues:
            row = residues.index(position)
            dssp = ring_data[row,4]
        else:
            dssp="?"
        pres = [i for i in range(len(phosphosites)) if phosphosites[i] in residues]
        prows = [residues.index(phosphosites[i]) for i in pres]
        rapdf,degree=[0,0]
        if prows:
            try:
                rapdf = wgt_avg(ring_data[prows,11].astype(float),[weights[i] for i in pres]) #but actually, we will do a weighted mean.
                degree = wgt_avg(ring_data[prows,5].astype(float),[weights[i] for i in pres])
            except TypeError:
                print(accession,prows,ring_data[prows,11],weights)
        return [int(dssp=="H"),int(dssp=="S"),int(dssp==" "),rapdf,degree]
    except ValueError: #if RING fails to compute
        filter = [i for i in range(len(pdb_filtered)) if pdb_filtered[i,0]!=pdb_id]
        pdb_filtered = pdb_filtered[filter,:]
        return ring(pdb_filtered)
#EBC
   
ebc_table = pd.read_csv("PPIN_ebc.txt",sep="\t",header=None).iloc[:,:].values     
ebc_dict = {}
for row in ebc_table:
    ebc_dict[(row[0],row[1])]=row[2]
    ebc_dict[(row[1],row[0])]=row[2]
ebc_median = np.median(list(ebc_dict.values()))

#PolyPhen scores

cmd = ['curl', '-F', '_ggi_project=PPHWeb2', '-F', '_ggi_origin=query', '-F', 
       '_ggi_target_pipeline=1', '-F', 'MODELNAME=HumVar', '-F', 'UCSCDB=hg19',
       '-F', 'SNPFUNC=m', '-F', '_ggi_batch_file=@{0}'.format(file), 
       '-D', '-', 'http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi']

out = open("polyphen_output.txt","w")
subprocess.run(cmd, stdout=out)
out.close()

polyphen_output = open("polyphen_output.txt").read()
job_id = re.search("polyphenweb2=([^;]+);", polyphen_output).group(1)

cmd1 = ['sleep','30']
cmd2 = ['wget', 'http://genetics.bwh.harvard.edu/ggi/pph2/{0}/1/pph2-short.txt'.format(job_id), '-O', 'pph2-short.txt']

job_done = 0
while job_done==0:
    subprocess.run(cmd1)
    try:
        subprocess.run(cmd2)
        polyphen_scores = pd.read_csv("pph2-short.txt",sep="\t").iloc[:-5,:].values
        job_done = 1
    except pd.EmptyDataError:
        continue
        
for i in polyphen_scores:
    for j in [0,2,3]:
        i[j] = re.sub("\s","",i[j])


a = 0
for i in range(n):
    if all(subs[i]==polyphen_scores[a,:4]):
        try:
            feature_matrix[i,0]=polyphen_scores[a,10]
        except ValueError:
           feature_matrix[i,0]=None
        a+=1
    else:
        feature_matrix[i,0]=None


#NetPhorest
cmd1 = ['cat','gene_fasta.txt'] 
cmd2 = ['cat','gene_fasta_sub.txt']
cmd3 = ['./NetPhorest_human_2.1/netphorest']


# biophysical features

MWs = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147, "G":75,"H":155,
       "I":131,"L":131,"K":146, "M":149,"F":165,"P":115,"S":105,"T":119,"W":204,
       "Y":181,"V":117}
properties = {"A":"nonpolar","R":"basic","N":"polar","D":"acidic","C":"polar",
              "Q":"polar","E":"acidic","G":"nonpolar","H":"basic","I":"nonpolar",
              "L":"nonpolar","K":"basic","M":"nonpolar","F":"nonpolar","P":"nonpolar",
              "S":"polar","T":"polar","W":"nonpolar","Y":"polar","V":"nonpolar"}



#feature annotation

n = len(subs)
successes=[]
errors={}
for i in range(n):

    accession, position, input_aa, output_aa = subs[i]
    position=int(position)
    feature_matrix[i,1] = MWs[input_aa]
    feature_matrix[i,2] = MWs[output_aa]
    feature_matrix[i,3] = int(properties[input_aa]=="acidic")
    feature_matrix[i,4] = int(properties[input_aa]=="basic")
    feature_matrix[i,5] = int(properties[input_aa]=="polar")
    feature_matrix[i,6] = int(properties[output_aa]=="acidic")
    feature_matrix[i,7] = int(properties[output_aa]=="basic")
    feature_matrix[i,8] = int(properties[output_aa]=="polar")

#the NetPhorest stuff
    try:
        if not feature_matrix[i,0]:
            raise ValueError("No PolyPhen score could be found for this variant")
            
        j = list(uniprot_IDs_named[:,0]).index(accession)
        gene_name = uniprot_IDs_named[j,2]
        phosphosites = []
        if position in gene_phospho.get(gene_name,[]):
            feature_matrix[i,9] = 1
            feature_matrix[i,10] = -1 #so that #phosphoneighbours excludes direct phosphosite
        for k in range(position-5,position+6):
            if k in gene_phospho.get(gene_name,[]):
                feature_matrix[i,10]+=1
                phosphosites.append(k)
            
        if phosphosites:
            
            fasta(accession,"gene_fasta.txt")
            out = open("netpho_wt.csv","w")
            p = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
            subprocess.run(cmd3, stdin = p.stdout, stdout=out)
            
            fastasub(accession,position,input_aa,output_aa, "gene_fasta_sub.txt")
            out = open("netpho_sub.csv","w")
            p = subprocess.Popen(cmd2, stdout=subprocess.PIPE)
            subprocess.run(cmd3, stdin=p.stdout, stdout=out)
    
            netpho_wt = pd.read_csv("netpho_wt.csv",sep="\t").iloc[:,:].values
            netpho_sub = pd.read_csv("netpho_sub.csv",sep="\t").iloc[:,:].values
            ptmppi = api.get_ptm_dependent_ppi(accession, dict=True)
            ptmppi = [j for j in ptmppi if j['substrate']['uniprot_id'] ==accession]
            weights = []
            ppin_scores = []
            for k in phosphosites:
                score_wt = sum(netpho_wt[j,7] for j in range(netpho_wt.shape[0]) if netpho_wt[j,0]==k)
                score_sub = sum(netpho_sub[j,7] for j in range(netpho_sub.shape[0]) if netpho_sub[j,0]==k)
                weight = abs(score_sub - score_wt)
                weights.append(weight)
            
                ppi_mod = [j for j in ptmppi if re.search(str(k)+"$",j['site'])]
                interactants = set([j['interactant']['uniprot_id'] for j in ppi_mod])
            
                ppin_scores.append(0)
            
                for j in interactants:
                    ppi_mod_int = [k for k in ppi_mod if k['interactant']['uniprot_id']==j]
                    assoc = [k['association_type'] for k in ppi_mod_int]
                    if all(k=='association' for k in assoc):
                        dependence = 1
                    else:
                        dependence = sum(2*(k=="inhibited_association") + int(k=="increased_association" or k=="decreased_association") for k in assoc)/sum(k!="association" for k in assoc)
                
                    ebc = ebc_dict.get((accession,j),ebc_median) 
                    ppin_scores[-1] += ebc*dependence
            feature_matrix[i,17]=sum(weights)
            if ppin_scores and sum(weights)!=0:
                nps = wgt_avg(ppin_scores,weights)
                feature_matrix[i,11] = nps
    
        # RING
        try:
            pdb_matrix = acc2pdb(accession)
        except ValueError:
            feature_matrix[i,12:17] = np.nan
            print(accession)
            continue
        if type(pdb_matrix) == np.ndarray:
            ring_features = ring(pdb_matrix)
            feature_matrix[i,12:17] = ring_features
        else:
            feature_matrix[i,12:17] = np.nan
        successes.append(i)
    except Exception as e:
        errors[i]=e
        continue
    

error_file = open("results/error_log.txt","w")
for e in errors:
    i = subs[e]
    error_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(i[0],i[1],i[2],i[3],errors[e]))
error_file.close()


#import classifier and make predictions

clf = pickle.load(open("classifier","rb"))
scaler = pickle.load(open("clf_scaler","rb"))


X = feature_matrix[successes,:]

for i in range(len(X)):
    for j in range(X.shape[1]):
        if str(X[i,j])=="nan":
            X[i,j]=0.0
            
X[:,:3] = scaler.transform(X[:,:3])
X[:,10:12] = scaler.transform(X[:,10:12])
X[:,15:] = scaler.transform(X[:,15:])


pred = clf.predict(X)
pred = np.reshape(pred,(len(pred),1))
prob = clf.predict_proba(X)[:,1]

subs_final = subs[successes,:]
result = np.hstack((subs_final,pred,prob))
np.savetxt("results/scores.txt",result,delimiter="\t",fmt="%s")   