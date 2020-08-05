
"""
Created on Mon Jul 20 10:33:31 2020

@author: farihatanjin
"""

import pandas as pd
import scanpy as sc
import numpy as np 


def cell_grouping(condition):

    adata = sc.read_csv('scRecover+scImpute_'+  condition +'_condition.csv').transpose()
    

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    sc.tl.pca(adata, svd_solver='arpack')
    
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    sc.tl.umap(adata)
                                    
    
    sc.tl.leiden(adata)
    sc.pl.umap(adata,color='leiden')
    
    raw = pd.DataFrame(data=adata.X, columns=adata.var_names)
    
    avg_gata=np.average(raw['Gata2'].to_numpy())
    avg_sox=np.average(raw['Sox2'].to_numpy())
    avg_zic=np.average(raw['Zic3'].to_numpy())
    
    
    labels=[]
    
    
    for i in range(0,len(raw)):
        if raw['Gata2'][i] > avg_gata:
            labels.append('2c')
        elif raw['Sox2'][i] > avg_sox:
            labels.append('naive')
        elif raw['Zic3'][i] > avg_zic:
            labels.append('primed')
        else:
            labels.append('unknown')
    
    raw.set_index(adata.obs_names)
    
    adata.obs.leiden=labels
    
    
    sc.pl.umap(adata,color='leiden')
    
    adata_2c = adata[adata.obs.leiden=='2c']
    adata_naive=adata[adata.obs.leiden=='naive']
    adata_primed = adata[adata.obs.leiden=='primed']
    
    
    raw_2c = pd.DataFrame(adata_2c.X, columns=adata_2c.var_names)
    raw_naive = pd.DataFrame(adata_naive.X, columns=adata_naive.var_names)
    raw_primed = pd.DataFrame(adata_primed.X, columns=adata_primed.var_names)
    
    raw_2c = raw_2c.transpose()
    raw_naive = raw_naive.transpose()
    raw_primed = raw_primed.transpose()
    
    
    
    raw_2c.to_csv(condition + '_2c.csv')
    raw_naive.to_csv(condition + '_naive.csv')
    raw_primed.to_csv(condition + '_primed.csv')
    
    
    
conditions = ['2iLif', 'CiTC', 'LB']
for x in conditions: cell_grouping(x)


