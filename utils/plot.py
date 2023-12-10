import math
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as patches

def horizontal_structureplot1(LDAmatrix1,LDAtag,nrow=1,startpalette=0,fontsize=20,rotate=False,sortby='False'):
    plt.style.use('default')
    plt.rcParams.update({'font.size': fontsize,'font.family':'Arial'})
    topics=['topic_'+str(i+1) for i in range(LDAmatrix1.shape[1])]
    if sortby and sortby in topics:
        dfx=pd.DataFrame(LDAmatrix1,columns=topics)
        dfx['tag']=LDAtag
        dfx=dfx.sort_values(sortby,ascending=False)
        LDAmatrix1=dfx[topics].to_numpy()
        LDAtag=list(dfx['tag'])
    
    LDAmatrix=np.asarray([x[::-1] for x in LDAmatrix1])
    #LDAmatrix=LDAmatrix1
    tag_types=list(sorted(list(set(LDAtag))))
    category_names = list(str(x) for x in range(1,1+LDAmatrix.shape[1]))[::-1]
    category_colors = plt.get_cmap('Paired')(
        np.linspace(0., (LDAmatrix.shape[1]+startpalette)/12-0.05, LDAmatrix.shape[1]+startpalette))[startpalette:][::-1]
    #fig, axs = plt.subplots(int(np.ceil(len(tag_types)/col)),col,figsize=(6*1.618,10))
    if nrow==1:
        fig, axs = plt.subplots(len(tag_types), 1, figsize=(8, 4 * len(tag_types)))
    else:
        if rotate:
            fig, axs = plt.subplots(nrow,math.ceil(len(tag_types)/nrow), figsize=(8* len(tag_types)/nrow, (8/rotate) * len(tag_types)/nrow))
        else:
            fig, axs = plt.subplots(nrow,math.ceil(len(tag_types)/nrow), figsize=(8* len(tag_types)/nrow, 6 * len(tag_types)/nrow))
    #plt.rcParams.update({'font.size': fontsize,'font.family':'Arial'})
    
    for i in range(len(tag_types)):
        if nrow==1:
            ax = axs[i] if len(tag_types)>1 else axs
        else:
            ax=axs[i%nrow][int(i/nrow)]
        #ax=axs[i]
        tag_ = tag_types[i]
        subidx = [j for j in range(len(LDAtag)) if LDAtag[j]==tag_]
        sub_array = LDAmatrix[subidx,]
        ax.set_xlabel("Sample ID")
        ax.set_title(tag_)
        #ax.text(-96,0.6,tag_,color='white')
        ax.set_ylabel('Topic Membership')
        ax.set_ylim(0, 1)
        ax.set_xlim(-0.5,sub_array.shape[0]-0.5)
        accumulate=[]
        for j in range(LDAmatrix.shape[1]):
            if j==0:
                ax.bar([x for x in range(sub_array.shape[0])], sub_array[:,j], 1, label=category_names[j], color=category_colors[j])
                accumulate=sub_array[:,j]
            else:
                ax.bar([x for x in range(sub_array.shape[0])], sub_array[:,j], 1, label=category_names[j], color=category_colors[j], bottom=accumulate)
                accumulate=accumulate+sub_array[:,j]
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        #ax.set_xticks([])
    if nrow==1:
        axs[0].legend(ncol=1, bbox_to_anchor=(-0.15, 1),loc='best',prop={'size':fontsize/1.2})
    else:
        axs[0][0].legend(ncol=1, bbox_to_anchor=(-0.15, 1),loc='best',prop={'size':fontsize/1.2})
    #plt.rcParams.update({'font.size': fontsize,'font.family':'Arial'})
    fig.tight_layout()
    
    return fig,axs

def plotdriving_kmer(model,knum=5,blank=0.03,font_size=13):
    plt.rcParams.update({'font.size': font_size,'font.family':'Arial'})
    fig, axs = plt.subplots(figsize=((0.25+(knum-1)*blank*knum)*6,((model.topic_num-1)*0.1+blank+0.1)*6))
    
    for i in range(model.topic_num):
        rect = patches.Rectangle((0, 0+i*0.1), 0.1, 0.1, linewidth=1, edgecolor='none', facecolor=sns.color_palette('Paired')[i])
        axs.add_patch(rect)
        axs.text(0.15,0+i*0.1+blank,str(i+1)+': ')
        for j in range(knum):
            axs.text(0.25+j*blank*knum,0+i*0.1+blank,model.dKmer['topic_'+str(i+1)][j]+'')

    plt.axis('off')
    plt.xlim(0,0.25+j*blank*(knum))
    plt.ylim(0,i*0.1+blank+0.1)
    plt.rcParams.update({'font.size': 13+0.5*knum,'font.family':'Arial'})
    return fig,axs
