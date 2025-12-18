import random
import os,sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import LatentDirichletAllocation
from sklearn.preprocessing import StandardScaler

from utils import matrix,plot

class muse:
    def __init__(self, mode, k, sequences, labels, sliding_window=0, subgroup_size=1, topic_size=6,rs=0,decay=.7,matrix=None,startpos=0,singlechar='ACGT'):
        '''
        Initialize the muse class.
        mode: choose from ASC, AWSC, UWSC for different nucleotide analysis problems
        k: size of the k-mer
        sequences: list of sequences
        labels: list of sequence tags, must be the same length as sequences
        sliding_window: counting k-mers using a sliding_window, which will increase the number of k-mers in each sequence/position
        subgroup_size: for AWSC and UWSC modes only to use groups of sequences as single samples to fit LDA. This may get better performance of LDA since it increases the k-mer counts in each sample.
        topic_size: number of topics to be made
        rs: random state number, default 0
        decay: decay used for LDA fitting, default 0.7
        matrix: input k-mer count matrix for UWSC only. Instead of generating input matrix by this module, you can provide your own k-mer count matrix for specific purposes, e.g. counting k-mers in reading frames (which does not allow moving along the sequences by 1 but 3 nuc each time)
        startpos: the position of the starting point of the input sequences for ASC and AWSC only. By using this parameter, the first nucleotide in the input sequences will be refered as the startpos-th nucleotide, e.g. -2. 
        '''
        self.mode=mode
        self.original_sequences=sequences
        self.original_labels=labels
        self.labels=labels
        self.k=k
        self.sliding_window=sliding_window
        self.subgroup_size=subgroup_size
        self.topic_num=topic_size
        self.random_state=rs
        self.decay=decay
        self.test_matrix=[]
        self.test_label=[]
        self.matrix=matrix
        self.startpos=startpos
        self.singlechar=singlechar
        np.random.seed(self.random_state)
        os.environ['PYTHONHASHSEED'] = str(self.random_state)
        random.seed(self.random_state)
    
    def matrix_generate(self):
        if self.mode=='ASC':
            self.matrix,self.labels,self.detailtable=matrix.ASC(self.original_sequences, self.k, self.original_labels, self.sliding_window)
        elif self.mode=='UWSC':
            self.matrix,self.labels,self.dfseq,self.detailtable=matrix.UWSC(self.original_sequences, self.k, self.original_labels, self.sliding_window, self.subgroup_size,singlechar=self.singlechar)
        elif self.mode=='AWSC':
            self.matrix,self.labels,self.positionalkmers,self.detailtable=matrix.AWSC(self.original_sequences, self.k, self.original_labels, self.sliding_window, self.subgroup_size,startpos=self.startpos)
        self.sparcity=(self.matrix > 0).sum()/self.matrix.size
    
    def LDA(self,recursive=False):
        """
        Make input k-mer count matrix for model fitting step.
        """
        os.environ['PYTHONHASHSEED'] = str(self.random_state)
        random.seed(self.random_state)
        np.random.seed(self.random_state)
        self.lda=LatentDirichletAllocation(n_components=self.topic_num,learning_decay=self.decay,random_state=self.random_state).fit(self.matrix)
        self.transdata=self.lda.transform(self.matrix)
        self.score=self.lda.score(self.matrix)
        self.topics=self.lda.components_/self.lda.components_.sum(axis=1)[:, np.newaxis]
        kmers = matrix.generate_kmers(self.k,self.singlechar)
        self.sampletopic=pd.DataFrame(self.transdata,columns=['topic_' + str(x) for x in list(range(1, self.topic_num + 1))])
        self.sampletopic['label']=self.labels
        #calculate topic prior
        self.prior=np.sum(self.lda.components_,1)/np.sum(self.lda.components_)
        if self.mode!='AWSC':
            self.topicdf = pd.DataFrame(self.topics, columns=kmers,
                                index=['topic_' + str(x) for x in list(range(1, self.topic_num + 1))]).transpose()
        else:
            self.topicdf = pd.DataFrame(self.topics,columns=self.positionalkmers,
                                        index=['topic_' + str(x) for x in list(range(1, self.topic_num + 1))]).transpose()
        if recursive:
            pass
        #calculate expectation for training samples
        expectation_val = np.matmul(self.matrix / self.matrix.sum(1)[:, None],
                                    self.topicdf.to_numpy() / self.topicdf.to_numpy().sum(1)[:, None])
        self.train_expectation = pd.DataFrame(expectation_val, columns=['topic_' + str(i + 1) for i in range(self.topic_num)])

        expectation_val_train = np.matmul(self.matrix / self.matrix.sum(1)[:, None],
                                          self.topicdf.to_numpy() / self.topicdf.to_numpy().sum(1)[:, None])
        self.expectation_train = pd.DataFrame(expectation_val_train,
                                              columns=['topic_' + str(i + 1) for i in range(self.topic_num)])
        self.likelihood_train = pd.DataFrame(np.matmul(self.matrix / self.matrix.sum(1)[:, None],
                                          self.topicdf.to_numpy()), columns=['topic_' + str(i + 1) for i in range(self.topic_num)])
        return  
    
    def transform_newdata(self,input_sequences,input_labels,input_matrix=[]):
        '''
        Transform new sequences to topic distribution using fitted model.
        '''
        if input_matrix!=[]:
            self.test_original_matrix = input_matrix
            self.test_matrix = self.lda.transform(input_matrix)
            self.test_label = input_labels
            self.sampletopic_test = pd.DataFrame(self.test_matrix,
                                                 columns=['topic_' + str(i + 1) for i in range(self.topic_num)])
            self.sampletopic_test['label'] = self.test_label
            return
        if self.mode=='ASC':
            newmatrix,newlabels,self.testdetail=matrix.ASC(input_sequences, self.k, input_labels, self.sliding_window,block=self.wblock,block_size=self.wblocksize,block_slide=self.wblockslide)
        elif self.mode=='UWSC':
            newmatrix,newlabels,self.dfseq_test,self.testdetail=matrix.UWSC(input_sequences, self.k, input_labels, self.sliding_window, 1)
        elif self.mode=='AWSC':
            newmatrix,newlabels,self.newpositionalkmers,self.testdetail=matrix.AWSC(input_sequences, self.k, input_labels , self.sliding_window, 1)
        self.test_original_matrix=newmatrix
        self.test_matrix=self.lda.transform(newmatrix)
        self.test_label=newlabels
        self.sampletopic_test=pd.DataFrame(self.test_matrix,columns=['topic_'+str(i+1) for i in range(self.topic_num)])
        self.sampletopic_test['label']=self.test_label
        return
    
    def drivingKmer(self,multitopic=False,lowtp_remove=True,divergence='KL'):
        dftopic_ = self.topicdf
        if divergence=='regular':
            #compare w/ background and log
            original_frac=self.matrix.sum(0)/self.matrix.sum()
            new_frac=(dftopic_.transpose()/original_frac).transpose()
        if lowtp_remove:
            rm_topic=[x for x in dftopic_.columns if round(dftopic_[x].min(),6)==round(dftopic_[x].max(),6)]
            dftopic_=dftopic_.drop(labels=rm_topic,axis=1)
        #print(dftopic_)
        #if KL==False:
        dftop=pd.DataFrame()
        for x in dftopic_.columns:
            df0=dftopic_.sort_values(x,ascending=False)
            dftop[x]=df0.index
        self.dKmer_raw=dftop
        if divergence=='regular':
            self.dKmer=pd.DataFrame()
            self.dKmer_neg=pd.DataFrame()
            for x in new_frac.columns:
                df0=new_frac.sort_values(x,ascending=False)
                self.dKmer[x]=df0.index
                df0=new_frac.sort_values(x,ascending=True)
                self.dKmer_neg[x]=df0.index
        else:
            dfdistinct=matrix.distinctive1(dftopic_) if divergence=='KL' else matrix.distinctive_JS(dftopic_)
            self.kmerdist=dfdistinct
            dfdk_distinct=matrix.drivingKmer_bydistinct(dfdistinct,dftopic_,multiple_topic=multitopic,divergence=divergence)
            self.dKmer,self.dKmer_neg=dfdk_distinct
        return
    
    def probSeqTopic(self,random_by='self'):
        if random_by not in ['self','uniform']:
            print("The parameter 'random_by' can be either 'uniform' or 'self'. Abborting...")
            return
        self.likelihood=pd.DataFrame(index=list(range(self.test_original_matrix.shape[0])),columns=['topic_'+str(i+1) for i in list(range(self.topic_num))]).fillna(0)
        self.joint = pd.DataFrame(index=list(range(self.test_original_matrix.shape[0])),columns=['topic_' + str(i + 1) for i in list(range(self.topic_num))]).fillna(0)

        expectation_val = np.matmul(self.test_original_matrix / self.test_original_matrix.sum(1)[:, None],
                                    self.topicdf.to_numpy() / self.topicdf.to_numpy().sum(1)[:, None])
        self.expectation=pd.DataFrame(expectation_val,columns=['topic_'+str(i+1) for i in range(self.topic_num)])

        if random_by=='uniform':
            self.random_generate=np.sum(self.test_original_matrix,1)*np.log10(1/(self.topics.shape[1]))#(1/self.k)**np.sum(self.test_original_matrix,1)
        else:
            observed_freq=self.matrix.sum(0)/self.matrix.sum()
            self.random_generate=np.sum(self.test_original_matrix*np.array([np.log10(x) if x!=0 else 0 for x in observed_freq ]),1)
        for i in range(self.test_original_matrix.shape[0]):
            for x in self.joint.columns:
                valuearray=self.test_original_matrix[i,]*np.log10(np.array(self.topicdf[x]))
                self.likelihood.loc[i,x]=np.sum(valuearray)
                self.joint.loc[i,x]=self.likelihood.loc[i,x]+np.log10(self.prior[list(self.joint.columns).index(x)])
        self.normalized_likelihood=pd.DataFrame(self.likelihood.to_numpy()-np.log10((10**self.likelihood.to_numpy()).sum(1)[:,None]),columns=['topic_' + str(i + 1) for i in range(self.topic_num)])
        self.logratio_likelihood = pd.DataFrame(self.likelihood.to_numpy() - self.random_generate[:, None],columns=['topic_' + str(i + 1) for i in range(self.topic_num)])
        scaler_likelihood = StandardScaler()
        scaler_likelihood.fit(self.likelihood.to_numpy().transpose())
        self.scaled_likelihood=pd.DataFrame(scaler_likelihood.transform(self.likelihood.to_numpy().transpose()).transpose(),columns=['topic_'+str(i+1) for i in range(self.topic_num)])

        return
    
    def structurePlot(self,nrow=1,sortby=False,fontsize=17,plot_test=False):
        if plot_test:
            inmatrix=self.test_matrix
            inlabel=self.test_label
        else:
            inmatrix=self.transdata
            inlabel=self.labels
        f,ax=plot.horizontal_structureplot1(inmatrix,inlabel,nrow=nrow,sortby=sortby,fontsize=fontsize)
        return f,ax
    
    def drivingKmerPlot(self,knum=5,blank=0.05,fontsize=15):
        f,ax=plot.plotdriving_kmer(self,knum=knum,blank=blank,font_size=fontsize)
        return f,ax

'''
functions for cmd pipeline
'''

def data_prepare(fasta_path,label_path):
    #can be one fasta + a list of tags or multiple fasta + multiple tags one for each file
    if ',' in fasta_path:
        fastas=fasta_path.split(',')
        sequence_list=[]
        label_list=[]
        labels=label_path.split(',')
        fcount=0
        for fa in fastas:
            if fa=='':
                continue
            sequence_sub=matrix.getfa_seqs(fa,getlist=True)
            sequence_list+=sequence_sub
            labels+=[labels[fcount] for i in range(len(sequence_sub))]
            fcount+=1
    else:
        sequence_list=matrix.getfa_seqs(fasta_path,getlist=True)
        with open(label_path,'r') as handle:
            label_list=[x for x in handle.read().split('\n') if x!='']
    return sequence_list,label_list

def run_model(mode,k,sequences,labels,topic_num,sliding_window=0,subgroup_size=1,startpos=0):
    model=muse(mode=mode,k=k,sequences=sequences,labels=labels,topic_size=topic_num,sliding_window=sliding_window,startpos=startpos,subgroup_size=subgroup_size)
    model.matrix_generate()
    model.LDA()
    model.drivingKmer()
    return model

def load_model(model_path):
    with open(model_path,'rb') as handle:
        model=pickle.load(handle)
    return model

def run_xmuse(fasta_path,label_path,mode,k,topic_num,sliding_window,subgroup_size,startpos,modelpath='',outpath=''):
    '''
    1. prepare data from input
    2. fit model or transform by provided model
    3. save result and the model
    '''
    #1. prepare data from input
    sequences,labels=data_prepare(fasta_path=fasta_path,label_path=label_path)
    #2. fit model or transform new sequence
    if modelpath!='':
        #load model
        print('Loading model...')
        model=load_model(model_path=modelpath)
        mode=model.mode
        print('Transforming data...')
        model.transform_newdata(sequences,labels)
    else:
        print('Building model...')
        model=run_model(mode=mode,k=k,sequences=sequences,labels=labels,topic_num=topic_num,sliding_window=sliding_window,subgroup_size=subgroup_size,startpos=startpos)
    #3. save result and the model
    if outpath!='':
        if modelpath!='':
            with open(os.path.join(outpath, mode+'_model_updated.pickle'),'wb') as handle:
                pickle.dump(model,handle)
            #save the testing figure
            print('Generating structure plot...')
            f,ax=model.structurePlot(plot_test=True)
            plt.savefig(os.path.join(outpath, mode+'_test_structureplot.pdf'),dpi=300,bbox_inches='tight')
            print('Writing tables to disk...')
            dfdist=pd.DataFrame(model.test_matrix,columns=['topic_'+str(i+1) for i in range(model.topic_num)])
            dfdist['label']=model.test_label
            dfdist.to_csv(os.path.join(outpath, mode+'_test_sequence_topic_dist.csv'))
        else:
            with open(os.path.join(outpath, mode+'_model.pickle'),'wb') as handle:
                pickle.dump(model,handle)
            #save the fitting figure
            print('Generating structure plot...')
            f,ax=model.structurePlot()
            plt.savefig(os.path.join(outpath, mode+'_fit_structureplot.pdf'),dpi=300,bbox_inches='tight')
            f,ax=model.drivingKmerPlot()
            plt.savefig(os.path.join(outpath, mode+'_fit_drivingkmerplot.pdf'),dpi=300,bbox_inches='tight')
            print("Writing tables to disk...")
            model.topicdf.to_csv(os.path.join(outpath, mode+'_fit_topic_kmer_dist.csv'))
            dfdist=pd.DataFrame(model.transdata,columns=['topic_'+str(i+1) for i in range(model.topic_num)])
            dfdist['label']=model.labels
            dfdist.to_csv(os.path.join(outpath, mode+'_fit_sequence_topic_dist.csv'))
            model.dKmer.to_csv(os.path.join(outpath, mode+'_fit_driving_kmers.csv'))
    print('Finished!')
    return
