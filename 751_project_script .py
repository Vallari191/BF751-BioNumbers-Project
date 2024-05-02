#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('pip', 'install celltypist')


# In[2]:


# use for cell annotation
import celltypist
from celltypist import models


# In[ ]:





# In[12]:


import scanpy as sc
import scvi


# In[31]:


import anndata as ad


# In[13]:


import numpy as np


# In[23]:


data1 = sc.read_text("/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071/GSM4453576_P1_exp.txt").T
data2 = sc.read_text("/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071/GSM4453577_P2_exp.txt").T
data3 = sc.read_text("/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071/GSM4453578_P3_exp.txt").T


# In[27]:


data3.obs


# In[28]:


data2


# In[29]:


data3


# In[32]:


adatas = [data1, data2, data3]
adatas = ad.concat(adatas, join="outer")


# In[33]:


adatas


# In[34]:


#transform the data format
from scipy.sparse import csr_matrix
dir = '/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071/'
adatas.X = csr_matrix(adatas.X)
adatas.X
#store the data
adatas.write_h5ad(dir+'combined.h5ad')


# In[51]:


#read the data
dir = '/Users/irisxu/Desktop/BU/Spring/BF751 molecular bio/project/GSE148071/'
adata = sc.read_h5ad(dir+'combined.h5ad')
adata


# In[52]:


adata.obs


# In[53]:


#get rid of cells with fewer than 200 genes
sc.pp.filter_cells(adata, min_genes=200) 
#get rid of cells with more than 5000 genes
sc.pp.filter_cells(adata, max_genes=5000)
#get rid of genes that are found in fewer than 3 cells
sc.pp.filter_genes(adata, min_cells=3) 


# In[54]:


adata


# In[55]:


#filter mtDNA > 30%
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) 
adata = adata[adata.obs.pct_counts_mt<30,:]


# In[56]:


adata


# In[57]:


#removed cells that have too many expression genes
adata = adata[adata.obs.n_genes_by_counts<2500,:] #The number of genes with at least 1 count in a cell.


# In[58]:


adata


# In[59]:


adata.obs.head()


# In[60]:


#1.the number of genes expressed in the count matrix; 2.the total counts per cell; 3.the percentage of counts in mitochondrial genes
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)


# ## Normalization and log transformation

# In[61]:


#just to see the difference before and after normalization
adata.X.sum(axis = 1)


# In[62]:


#normalization + log transformation
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
#copy the raw data
adata.raw = adata


# In[63]:


adata.X.sum(axis = 1)


# In[64]:


import matplotlib.pyplot as plt


# In[65]:


#identify and visualize high variable gene
sc.pp.highly_variable_genes(adata, n_top_genes = 2000) #original paper select 600
sc.pl.highly_variable_genes(adata)
plt.savefig(dir+"03-highly_variable_genes.png")


# In[66]:


adata = adata[:, adata.var.highly_variable]
adata


# ## PCA and clustering

# In[67]:


#remove covariables
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
#scale each gene to unit variance
sc.pp.scale(adata, max_value=10)
#PCA reduction
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50) #Elbow method
plt.savefig(dir+"04-pca_variance.png")


# In[68]:


#Compute a neighborhood graph
sc.pp.neighbors(adata, n_pcs=20) 
#t-SNE; set n_pcs = 10
sc.tl.tsne(adata,n_pcs=20) 

#t-SNE plotting
sc.pl.tsne(adata, save='10pc_tsne_plot.png')


# In[73]:


#clustering using Leiden graph-clustering method
sc.tl.leiden(adata, resolution=0.4 )
adata.obs['leiden']


# In[74]:


#plot the clusters
sc.pl.tsne(adata, color = 'leiden')


# ## Differential analysis & find marker gene

# In[75]:


#using wilcoxon method to do the analysis
sc.tl.rank_genes_groups(adata, groupby = 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.settings.verbosity = 2  # reduce the verbosity
#sc.get.rank_genes_groups_df(adata, group="0")


# In[ ]:





# In[ ]:




