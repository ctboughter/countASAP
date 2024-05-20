# Need to run this in Ipython only.
import pandas
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib import rc
import glob
import gc
import umap
import sklearn.cluster as cluster
from sklearn.decomposition import PCA
import MEMSEAL
import anndata
from muon import MuData
from muon import prot as pt
    

################################
# This bit is for that figure formatting. Change font and font size if desired
font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 20}
COLOR = 'black'
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = COLOR
rcParams['xtick.color'] = COLOR
rcParams['ytick.color'] = COLOR

rc('font', **font)

# Lastly this custom colormap is for 
import matplotlib as mpl
upper = mpl.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
################################

# THESE ARE WHAT YOU CHANGE TO READ IN DIFFERENT FILES
day=1; fileNum = 0
# Choose if you want to drop highly correlated vectors/genes
drop_corr = False
# Choose if you want to display PCA or UMAP dimension reduction
dimRed = 'umap'
# Choose where you want to send outputs
finDir='./'

# LOAD IN YOUR ATACSEQ processed data (processed using the included scripts)
atac_mat = pandas.read_csv('rawDat/mat_'+str(day)+'0'+str(fileNum+1)+'.csv')
atac_names = pandas.read_csv('rawDat/gene_meta_'+str(day)+'0'+str(fileNum+1)+'.csv')

# We need to process our "atac_names" to pull out just headers for columns:
gene_names = atac_names['gene_info']
copyNA = gene_names[gene_names.isna()].index

# This is probably not good practice but it works!!
# Move over chromosomal names to replace NAN where there aren't precise gene names...
gene_names.loc[copyNA] = atac_names['regions'].loc[copyNA]
# Save these just in case....
chr_regions = atac_names['regions']

atac_mat.index = gene_names

# LOAD IN ASAPseq data, processed using CountASAP and HTOdemux
asap_mat = pandas.read_csv('rawDat/d'+str(day)+'0'+str(fileNum+1)+'_asap_count.csv',index_col=0)
hto_assign = pandas.read_csv('rawDat/asap'+str(day)+'0'+str(fileNum+1)+'_assign.csv')
# We need to read in the HTO normalized counts to get back out the 
hto_umis = pandas.read_csv('rawDat/asap'+str(day)+'0'+str(fileNum+1)+'_normCount.csv')
htoF = hto_assign['x']
htoF.index = hto_umis.columns[1:]

# Concatenate our HTOs into our ASAPseq matrix, just to be sure...
test_asap = pandas.concat([htoF,asap_mat],axis=1)
if len(htoF) == len(test_asap):
    print('hto and asap indices match')
# SO NOTE, in the original analysis, the ADT and SCT (CITE and GEX) were both arranged
# such that the gene names were in the indices and the barcodes were in the columns
# WE ACTUALLY DO NOT WANT OUR HTOS IN OUR FINAL ASAP MAT....
# And then all of this is for the data normalization through muon
ann_raw = anndata.AnnData(asap_mat)
mdata = MuData({"A": ann_raw})
pt.pp.clr(mdata['A'])
muon_norm = mdata['A'].to_df()
fin_asap = np.transpose(muon_norm)

# ALRIGHT NOW EVERYTHING SHOULD BE LOADED UP... GO INTO PROCESSING..
# Double check that the barcodes are in the exact same order (they should be)

# So unfortunately we need to play some dirty
# tricks here to match up the columns (annoying)
col_df1 = pandas.DataFrame(fin_asap.columns)
col_df2 = pandas.DataFrame(atac_mat.columns)
holdit=[]
for col1 in col_df1.values:
    if len(col_df2[col_df2[0] == col1[0]]) != 0:
        holdit = holdit + [col1[0]]

#full_df = pandas.concat([fin_asap[holdit],atac_mat[holdit]],axis=0)
full_df = atac_mat[holdit]
if drop_corr:
    parsed_mat = MEMSEAL.drop_corr(full_df)
else:
    parsed_mat = full_df

# This *is* a parallel step, so probably should request multiple
# compute nodes for this on LOCUS...
reducer = umap.UMAP(n_components=3, n_neighbors = 100)
if dimRed == 'umap':
    init_red = reducer.fit_transform(np.transpose(parsed_mat.values))
elif dimRed == 'pca':
    # Technically this is UMAP AND PCA
    pca = PCA(n_components=10)#, svd_solver='full')
    pc_vects=pca.fit_transform(np.transpose(parsed_mat.values))
    init_red = reducer.fit_transform(pc_vects)

# Could try other clustering algorithms, but we use kmeans
#clusts = cluster.DBSCAN(min_samples=100,eps=0.5).fit_predict(init_red)
clusts = cluster.KMeans(n_clusters=4).fit_predict(init_red)
# Have this as the standard for now
#clusts = cluster.DBSCAN(min_samples=25).fit_predict(init_umap)
clust_df = pandas.DataFrame(clusts)
clust_df.index = parsed_mat.columns

# Plot all of the things that we are seeing...
# Make sure that we have good purity
clust_input=init_red#.values

fig,ax = pl.subplots(1,1,figsize = (12,8))
#x= pl.scatter(clust_input[:,0],clust_input[:,1],edgecolors='black',c=np.sum(parsed_mat.loc[['S100a9']].values,axis=0),cmap=cm.viridis,vmin=0,vmax=3)
#pl.scatter(clust_input[:,0],clust_input[:,1],edgecolors='black',c=parsed_mat.loc[['Yam1']].iloc[0],cmap=cm.viridis)
x = pl.scatter(clust_input[:,0],clust_input[:,1],edgecolors='black',c=clust_df.values,cmap=cm.viridis)
#pl.colorbar(x)

pl.xlabel('UMAP1')
pl.ylabel('UMAP2')
pl.savefig('figure4.png',format='png',dpi=600)
pl.close()

####################################################
# This section for plotting bar graphs:
cite_cd4 = []; cd4 = []; cite_cd8=[];cd8 = []
cite_cd11 = []; s100a9=[]; cite_H2= []; cite_cd38=[]
cite_IA = []; cite_cd3 = []; cd3 = []

std_cite_cd4 = []; std_cd4 = []; std_cite_cd8=[]; std_cd8 = []
std_cite_cd11 = []; std_s100a9=[]; std_cite_H2= []; std_cite_cd38=[]
std_cite_IA = []; std_cite_cd3 = []; std_cd3 = []

pax5 = [];hsf2=[]
for i in clust_df.drop_duplicates().sort_values(0).values:
    subF = parsed_mat[clust_df[clust_df[0]==i[0]].index]
    subF_asap = ASAPF[clust_df[clust_df[0]==i[0]].index]

    cite_cd4 = cite_cd4 + [np.average(subF_asap.loc['CD4'].values)]
    cd4 = cd4 + [np.average(subF.loc['Cd4'].values)]
    cite_cd8 = cite_cd8 + [np.average(subF_asap.loc['CD8a'].values)]
    cd8 = cd8 + [np.average(subF.loc['Cd8a'].values)]
    cite_cd3 = cite_cd3 + [np.average(subF_asap.loc['CD3'].values)]
    cd3 = cd3 + [np.average(subF.loc['Cd3e'].values)]
    cite_IA = cite_IA + [np.average(subF_asap.loc['I.A_I.E'].values)]
    cite_H2 = cite_H2 + [np.average(subF.loc['H2-DMb2'].values)]
    cite_cd38 = cite_cd38 + [np.average(subF_asap.loc['CD38'].values)]
    cite_cd11 = cite_cd11 + [np.average(subF_asap.loc['CD11b'].values)]
    s100a9 = s100a9 + [np.average(subF.loc['S100a9'].values)]
    pax5=pax5+[np.average(subF.loc['Pax5'].values)]
    hsf2=hsf2+[np.average(subF.loc['Hsf2'].values)]

pl.figure(figsize=(12,8))
pl.bar(np.arange(4)+0.0,cd4/max(cd4),width=1/10,color='firebrick')
#pl.bar(np.arange(5)+0.6,cd8/max(cd8),width=1/10,color='red')
#pl.bar(np.arange(5)+0.0,cd3/max(cd3),width=1/10,color='red',edgecolor='black')
pl.bar(np.arange(4)+0.1,cite_cd4/max(cite_cd4),width=1/10,color='pink',hatch='/',edgecolor='black')
#pl.bar(np.arange(5)+0.7,cite_cd8/max(cite_cd8),width=1/10,color='red',hatch='o')
#pl.bar(np.arange(5)+0.1,cite_cd3/max(cite_cd3),width=1/10,color='red',hatch='O',edgecolor='black')

pl.bar(np.arange(4)+0.2,pax5/max(pax5),width=1/10,color='royalblue',edgecolor='black')
pl.bar(np.arange(4)+0.3,cite_cd38/max(cite_cd38),width=1/10,color='cyan',edgecolor='black',hatch='/')

pl.bar(np.arange(4)+0.4,s100a9/max(s100a9),width=1/10,color='purple',edgecolor='black')
pl.bar(np.arange(4)+0.5,cite_cd11/max(cite_cd11),width=1/10,color='violet',edgecolor='black',hatch='/')

pl.legend(['CD4','ASAP CD4','Pax5','ASAP CD38','S100A9','ASAP CD11'],ncol=3,bbox_to_anchor=(0.37,1.01))
pl.ylabel('Average Expression (Normalized)')
pl.savefig('asap_expression_atacOnly.pdf')
pl.close()