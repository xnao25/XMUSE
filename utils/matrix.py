import itertools
import numpy as np
import pandas as pd
from Bio import SeqIO

def generate_kmers(k,seq='ACGT'):
    return [''.join(list(x)) for x in itertools.product(seq, repeat=k)]

def ASC(sequences,k,labels,sliding_window=0,block=False,block_size=5,block_slide=False):
    '''
    Make position x k-mer matrix for Aligned Sequence Characterization LDA analysis. Sequences of the same tag will be grouped and the counts of different k-mers at each position will be summarized.
    '''
    newlabels=[]
    num_positions=len(sequences[0])-k-sliding_window+1
    possible_kmers=generate_kmers(k)
    combine_matrix=np.matrix([])
    for x in set(labels):
        output_matrix = np.zeros((num_positions, len(possible_kmers)))
        for y in range(num_positions):
            newlabels.append(x)
        for s in range(len(sequences)):
            if labels[s]!=x:
                continue
            seq=sequences[s]
            for i in range(num_positions):
                for j in range(sliding_window+1):
                    the_kmer=seq[i+j:i+j+k]
                    output_matrix[i,possible_kmers.index(the_kmer)]+=1
        if combine_matrix.shape==np.matrix([]).shape:
            combine_matrix=output_matrix
        else:
            combine_matrix=np.concatenate((combine_matrix,output_matrix),axis=0)
    return combine_matrix,newlabels,pd.DataFrame(combine_matrix,index=newlabels,columns=possible_kmers)

def UWSC(sequences,k,labels,sliding_window=0,subgroup_size=1,singlechar='ACGT'):
    '''
    Make sequence/subgroup - k-mer count matrix for Unaligned Whole Sequence Characterization LDA analysis. The k-mer counts for single sequences or subgroup of sequences will be summarized.
    '''
    possible_kmers=generate_kmers(k,singlechar)
    output_matrix=np.zeros((int(len(sequences)/subgroup_size),len(possible_kmers)))
    newlabels=[]
    dfseq=pd.DataFrame()
    dfseq['sequence']=sequences
    dfseq['label']=labels
    dfseq=dfseq.sort_values('label')
    dfseq.index=list(range(len(dfseq)))
    for s in range(len(sequences)):
        seq=dfseq.loc[s,'sequence']
        if s%subgroup_size==0:
            newlabels.append(dfseq.loc[s,'label'])
        for i in range(len(seq)-k-sliding_window+1):
            for j in range(sliding_window+1):
                the_kmer=seq[i+j:i+j+k]
                output_matrix[int(s/subgroup_size),possible_kmers.index(the_kmer)]+=1
    return output_matrix,newlabels,dfseq,pd.DataFrame(output_matrix,index=newlabels,columns=possible_kmers)

def AWSC(sequences,k,labels,sliding_window=0,subgroup_size=1,block=False,block_size=5,block_slide=False,defined_feature=[],startpos=0):
    '''
    Make sequence/subgroup - k-mer count matrix for Aligned Whole Sequence Characterization LDA analysis. Instead of using k-mer counts, here we used positional k-mer counts, e.g. AG_at_0 for every possible k-mer at every possible position.
    '''
    possible_kmers=generate_kmers(k)
    newlabels=[]
    scount=0
    if block==False:
        possible_positional_kmers=[x+'_at_'+str(y+startpos) for x in possible_kmers for y in range(len(sequences[0])-k-sliding_window+1)]
        output_matrix=np.zeros((int(len(sequences)/subgroup_size),len(possible_kmers)*(len(sequences[0])-k-sliding_window+1)))
        #print(output_matrix.shape)
        for s in range(len(sequences)):
            seq=sequences[s]
            if scount%subgroup_size==0:
                newlabels.append(labels[scount])
            scount+=1
            for i in range(len(seq)-k-sliding_window+1):
                for j in range(sliding_window+1):
                    the_positional_kmer=seq[i+j:i+j+k]+'_at_'+str(i+startpos)
                    output_matrix[int(s/subgroup_size),possible_positional_kmers.index(the_positional_kmer)]+=1
    else:
        possible_positional_kmers=[x+'_at_'+str(y) for x in possible_kmers for y in range(int(len(sequences[0])/block_size))]
        output_matrix=np.zeros((int(len(sequences)/subgroup_size),len(possible_positional_kmers)))
        for s in range(len(sequences)):
            seq = sequences[s]
            if scount%subgroup_size==0:
                newlabels.append(labels[scount])
            for i in range(len(seq) - k + 1):
                the_positional_kmer = seq[i:i+k] + '_at_' + str(int(i)/block_size)
                output_matrix[int(s/subgroup_size), possible_positional_kmers.index(the_positional_kmer)] += 1
    return output_matrix,newlabels,possible_positional_kmers,pd.DataFrame(output_matrix,index=newlabels,columns=possible_positional_kmers)

def distinctive1(dftopic):#K_L divergence
    dfdis=pd.DataFrame()
    for kmer in dftopic.index:
        for topic in dftopic.columns:
            other_topic=[x for x in dftopic.columns if x!=topic]
            kmer_topic_distinct=np.min([dftopic.loc[kmer,topic]*np.log(dftopic.loc[kmer,topic]/dftopic.loc[kmer,x])+dftopic.loc[kmer,x]-dftopic.loc[kmer,topic] for x in other_topic])
            dfdis.loc[kmer,topic]=kmer_topic_distinct
    return dfdis

def poisson_KL(x,y):
    return x*np.log(x/y)+y-x

def distinctive_JS(dftopic):#Jensen-Shannon
    dfdis=pd.DataFrame()
    for kmer in dftopic.index:
        for topic in dftopic.columns:
            other_topic=[x for x in dftopic.columns if x!=topic]
            JS_divergence=[]
            for x in other_topic:
                M=(dftopic.loc[kmer,topic]+dftopic.loc[kmer,x])/2
                D_PM=poisson_KL(dftopic.loc[kmer,topic],M)
                D_QM=poisson_KL(dftopic.loc[kmer,x],M)
                JS_divergence.append(0.5*D_PM+0.5*D_QM)
            dfdis.loc[kmer,topic]=np.max(JS_divergence)
    return dfdis

def drivingKmer(dftopic_,divergence='KL'):
    dftop=pd.DataFrame()
    for x in dftopic_.columns:
        df0=dftopic_.sort_values(x,ascending=False) if divergence=='KL' else dftopic_.sort_values(x,ascending=True)
        dftop[x]=[y if np.isnan(df0.loc[y,x])==False else 'NA' for y in df0.index ]
    return dftop

def drivingKmer_bydistinct(dfdistinc_,dftopic_,multiple_topic=False,divergence='KL'):
    dfdistinc_0=dfdistinc_.copy()
    dfdistinc_1 = dfdistinc_.copy()
    for i in dfdistinc_.index:
        for c in dfdistinc_.columns:
            if multiple_topic:
                if dftopic_.loc[i, c] < np.mean(dftopic_.loc[i]):
                    dfdistinc_0.loc[i, c] = np.nan
                if dftopic_.loc[i, c] > np.mean(dftopic_.loc[i]):
                    dfdistinc_1.loc[i, c] = np.nan
            else:
                if dftopic_.loc[i,c]<max(dftopic_.loc[i]):
                    dfdistinc_0.loc[i,c]=np.nan
                if dftopic_.loc[i,c]!=min(dftopic_.loc[i]):
                    dfdistinc_1.loc[i, c] = np.nan
    return drivingKmer(dfdistinc_0,divergence),drivingKmer(dfdistinc_1,divergence)

def getfa_seqs(fapath,getlist=False): #get sequences from fasta files and output a dictionary with keys as titles and values as sequences
    outseq={}
    outseqlist=[]
    with open(fapath,'r') as handle:
        for record in SeqIO.parse(handle,'fasta'):
            if getlist:
                outseqlist.append(str(record.seq).upper())
            else:
                outseq[record.id]=str(record.seq).upper()
    if getlist:
        return outseqlist
    return outseq
