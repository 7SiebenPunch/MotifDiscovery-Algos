import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from sax import SAX_subsequences
from datetime import datetime
import stumpy
import matrixprofile as mp
from grammar import GrammarInduction
from mkalgo.mk import mk, mk_eab



def stomp(ts, path=None, windowsize=100, topk=3):
    """
    Extract and visualize Motif in time series using STMOP algorithm
    ts: DataFrame, two columns, one for time and one for variables
    windowsize:Time Window
    topk:The topk motifs to extract and visualize
    """
    #starttime = datetime.now()
    ts1 = ts[ts.columns[1]].values
    profile = mp.algorithms.stomp(ts1, windowsize)
    profile = mp.discover.motifs(profile, k=topk)


#output motifs
    for i in range(topk):   
        m_idx = profile['motifs'][i]['motifs'][0]
        motifs = ts.iloc[m_idx :m_idx + windowsize]
        motifs.to_csv(f'{path}{ts.columns[1]}motif{i+1}_Stomp.csv',index = False)
        m_idxneib = profile['motifs'][i]['motifs'][1]
        

        _, ax = plt.subplots(figsize=(14, 6)) 
        ax.set_title(f'{ts.columns[1]}-' + 'Motif %i' % (i + 1) + ' STOMP')       
        ax.plot(ts[ts.columns[1]].values, alpha=0.7, c='k')
       # ax.set_ylabel(ts.columns[1], fontsize='20')
        ax.plot(ts[ts.columns[1]].iloc[profile['motifs'][i]['motifs'][0]:profile['motifs'][i]['motifs'][0] + windowsize], c='red', alpha=0.9)
        ax.plot(ts[ts.columns[1]].iloc[profile['motifs'][i]['motifs'][1]:profile['motifs'][i]['motifs'][1] + windowsize], c='red', alpha=0.9)
        ax.set_xlabel('Time')
       # plt.savefig(f'{path}{ts.columns[1]}motif{i+1}_Stomp.jpg')
   
    #endtime = datetime.now()
    #print(endtime - starttime)
    return 


    
    

def scrimp_2plus(ts, path=None, windowsize=100, topk=3):
    #starttime = datetime.now()
    """
    Extract and visualize Motif in time series using Scrimp++ algorithm
    ts: DataFrame, two columns, one for time and one for variables
    windowsize:Time Window
    topk:The topk motifs to extract and visualize
    """
    ts1 = ts[ts.columns[1]].values
    profile = mp.algorithms.scrimp.scrimp_plus_plus(ts1, windowsize)
    profile = mp.discover.motifs(profile, k=topk)
    # mp.visualize(profile)
    
    # output motifs
    for i in range(topk):   
        m_idx = profile['motifs'][i]['motifs'][0]
        motifs = ts.iloc[m_idx :m_idx + windowsize]
        motifs.to_csv(f'{path}{ts.columns[1]}motif{i+1}_Scrimp.csv',index = False)
  
        

        _, ax = plt.subplots(figsize=(14, 6)) 
        ax.set_title(f'{ts.columns[1]}-'+'Motif %i' % (i + 1) + ' SCRIMP++')       
        ax.plot(ts[ts.columns[1]].values, alpha=0.7, c='k')
       #  ax.set_ylabel(ts.columns[1], fontsize='20')
        ax.plot(ts[ts.columns[1]].iloc[profile['motifs'][i]['motifs'][0]:profile['motifs'][i]['motifs'][0] + windowsize], c='red', alpha=0.9)
        ax.plot(ts[ts.columns[1]].iloc[profile['motifs'][i]['motifs'][1]:profile['motifs'][i]['motifs'][1] + windowsize], c='red', alpha=0.9)
        ax.set_xlabel('Time')        
    # endtime = datetime.now()
    # print(endtime - starttime)
    # plt.savefig(f'{path}{ts.columns[1]}motif{i+1}_Scrimp.jpg')
   
    return


def stump(ts, path=None, windowsize=100, topk=3):
    """
    Extract and visualize Motif in time series using Stump(SCAMP) algorithm
    ts: DataFrame, two columns, one for time and one for variables
    windowsize:Time Window
    topk:The topk motifs to extract and visualize
    """
    #starttime = datetime.now()

    m = windowsize
    mp = stumpy.stump(ts[ts.columns[1]], m)
    motif_idx = np.argsort(mp[:, 0])[0]
    nearest_neighbor_idx = mp[motif_idx, 1]

    for i in range(topk):
        fig, ax = plt.subplots(figsize=(14, 6), sharex=True, gridspec_kw={'hspace': 0})
        plt.suptitle( f'{ts.columns[1]}-' + 'Motif %i' % (i + 1) + ' STUMP')
      #  plt.suptitle('Motif %i' % (i + 1), fontsize='30')
        motif_idx = np.argsort(mp[:, i])[0]
        nearest_neighbor_idx = mp[motif_idx, 1]
        
        motif = ts.iloc[motif_idx:motif_idx + m]
        motif.to_csv(f'{path}{ts.columns[1]}motif{i+1}_Stump.csv',index = False)
        
        ax.plot(ts[ts.columns[1]].values, alpha=0.7, c='k')
       #  ax.set_ylabel(ts.columns[1], fontsize='20')
        ax.plot(ts[ts.columns[1]].iloc[motif_idx:motif_idx + m], c='red', alpha=0.9)
        ax.plot(ts[ts.columns[1]].iloc[nearest_neighbor_idx:nearest_neighbor_idx + m], c='red', alpha=0.9)
        ax.set_xlabel('Time')
#       plt.savefig(f'{path}{ts.columns[1]}motif{i+1}_Stump.jpg')
        plt.show()

    #endtime = datetime.now()
    #print(endtime - starttime)

    return



def gpu_stump(ts, path=None, windowsize=100, topk=3):
    #starttime = datetime.now()
    """
    Extract and visualize Motif in time series using gpu_Stump algorithm
    ts: DataFrame, two columns, one for time and one for variables
    windowsize:Time Window
    topk:The topk motifs to extract and visualize    
    """
    m = windowsize
    mp = stumpy.gpu_stump(ts[ts.columns[1]], m)
    motif_idx = np.argsort(mp[:, 0])[0]
    nearest_neighbor_idx = mp[motif_idx, 1]
    for i in range(topk):
        fig, ax = plt.subplots(figsize=(14, 6), sharex=True, gridspec_kw={'hspace': 0})
        plt.suptitle( f'{ts.columns[1]}-' + 'Motif %i' % (i + 1) + ' GPU_STUMP')
        motif_idx = np.argsort(mp[:, i])[0]
        nearest_neighbor_idx = mp[motif_idx, 1]
        
        motif = ts.iloc[motif_idx:motif_idx + m]
        motif.to_csv(f'{path}{ts.columns[1]}motif{i+1}_GStump.csv',index = False)
        
        ax.plot(ts[ts.columns[1]].values, alpha=0.7, c='k')
        ax.plot(ts[ts.columns[1]].iloc[motif_idx:motif_idx + m],  c='red', alpha=0.9)
        ax.plot(ts[ts.columns[1]].iloc[nearest_neighbor_idx:nearest_neighbor_idx + m], c='red', alpha=0.9)
        ax.set_xlabel('Time')
        
#        plt.savefig(f'{path}{ts.columns[1]}motif{i+1}_GStump.jpg')
   
        plt.show()

    #endtime = datetime.now()
    #print(endtime - starttime)

    return


def grammarintroduction(df, path=None, w=6, n=100, k=10, topk=3):
    
    """
    Extracting and visualizing Motifs in time series using Grammarintroduction
    ts：DataFrame, two columns, one for time and one for variables
    w:Subsequence character length
    n:subsequences length
    k:stride of k
    topk:The topk motifs to extract and visualize
    """
    
    # starttime = datetime.now()
    gi = GrammarInduction(df, w=w, n=n, k=k)
    gi.process()
    gi.show_rules(i_col=0, ordered=True)

    i_col = 0
    i_rule = np.arange(topk)
    for i in i_rule:
        # Encode the reduced representation in suffix trees
        if not gi.trees:
            gi.trees = gi.generate_suffix_trees(gi.reduced)
        # Retrieve correct variables for the search
        tree = gi.trees[i_col]
        motif = gi.ranked_rules[i_col].iloc[i, 1]
        mapping = gi.mapping[gi.df.columns[i_col + 1]]

        # Search for the motif using the tree
        _, starting_positions = tree.find_motifs(motif)

        """
        We must correct the starting positions of the motif to 
        take into account the size of a symbol (+ 1 because of 
        the space separation between two successive symbols)
        Note: we have to do this because the tree encodes each 
        character as a unique symbol. But the way the tree is
        implemented could cause some other issues.
        E.g.: 'ABC ACB' is encoded as 'A', 'B', 'C', 'A','C''B'
        in the tree rather than 'ABC', 'ACB'. This may create
        erroneous matches when searching for patterns, for 
        instance: 'BCA' could be detected in the previous 
        string even though it should not be.
        """
        starting_positions = [pos // (len(motif.split(' ')[0]) + 1)
                              for pos in starting_positions]

        # Reverse the numerosity reduction
        motif_indices = gi.find_true_motif_positions(mapping,
                                                     motif,
                                                     starting_positions)

        # Find the time values associated with these indices
        if not gi.sax_time_vectors:
            gi.sax_time_vectors = [df['t'].tolist()
                                   for df in gi.sax.df_SAX]

        time_indices = gi.find_time_values(gi.sax_time_vectors,
                                           motif_indices)

        """
        Plots the motifs.
        """

        cut_dataframes = []
        for indices in time_indices:
            cut_dataframes.append(
                df.loc[(df['t'] >= indices[0]) &
                       (df['t'] < indices[1])]
            )
        motifdf = DataFrame(cut_dataframes[1])
        motifdf.to_csv(f'{path}{df.columns[1]}motif{i+1}_GI.csv',index = False)
        # Plot
        _, ax = plt.subplots(figsize=(14, 6))

        ax.plot(df.iloc[:, 0], df.iloc[:, i_col + 1], c='k', alpha=0.7)

        for d in cut_dataframes:
            ax.plot(d.iloc[:, 0], d.iloc[:, i_col + 1], c='red', alpha=0.9)

        ax.set_title(f'{d.columns[i_col + 1]}-' + 'Motif %i' % (i + 1) + ' GRAMMARINTRODUCTION')
        ax.set_xlabel('Time', fontsize=15)
#        plt.savefig(f'{path}{ts.columns[1]}motif{i+1}_GI.jpg')
   
    # endtime = datetime.now()
    # print(endtime - starttime)
    return


def mkal(ts, path=None, windowsize=100, metric='euclidean'):

    """
    Search for the most significant motif in the time series using MK
    ts: DataFrame, two columns, one for time and one for variables
    windowsize: subsequence length
    metric:Distance calculation method: 'dtw' or 'euclidean'
    """
    #starttime = datetime.now()

    # Use the index matrix to speed up the calculation process
    sax = SAX_subsequences(ts, w=8, n=100, k=2, alphabet_type='letters')

    data = ts[ts.columns[1]].values.tolist()
    obj = mk(l=windowsize, metric=metric,r=200)
    motif_a, motif_b = obj.search(data)



    # Plot
    cut_dataframes = []
    motif = []
    cut_dataframes.append(
        ts.loc[(ts.index >= motif_a['begin']) &
               (ts.index < motif_a['end'])]
    )   
    
    cut_dataframes.append(
        ts.loc[(ts.index >= motif_b['begin']) &
               (ts.index < motif_b['end'])]
    )
    
    _, ax = plt.subplots(figsize=(14, 6))
    ax.plot(ts.iloc[:, 0], ts.iloc[:, 1], c='k', alpha=0.7)
    for d in cut_dataframes:
        plt.plot(d.iloc[:, 0], d.iloc[:, 1], c='red', alpha=0.9)
    ax.set_title(f'{d.columns[1]}-' + 'Motif 1' + ' MK')
    ax.set_xlabel('Time', fontsize=15)
    motifdf = DataFrame(cut_dataframes[1])
    motifdf.to_csv(f'{path}{d.columns[1]}motif_MK.csv', index = False)
    #plt.savefig('MK.jpg')
    #endtime = datetime.now()
    #print(endtime - starttime)

    return 
