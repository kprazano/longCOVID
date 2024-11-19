#Preparing the AnnData files from csv count files

import pandas as pd
import anndata as ad

#Read the data and transpose so that columns are cells and rows are genes
sample = pd.read_csv("<sample>.csv", header=None, skiprows=7) #skipping comments in the first 7 rows
sample = sample.T
sample = pd.DataFrame(sample)

#Create and check the AnnData object
adata = ad.AnnData(sample)
print(adata)

#Save
adata.write("sample.h5ad")
