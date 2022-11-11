import scanpy as sc
import scvi
import seaborn as sns

results_file = './output/rmDoublrawDataset5000.h5ad'

adata = sc.read("./dataSaveOriginal/rawDataset.h5ad")
adata.var_names_make_unique()
adata.X = adata.X.astype('float64')

print(adata)

scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()

solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train()

df = solo.predict()
df['prediction'] = solo.predict(soft = False)

df.index = df.index.map(lambda x: x[:-2])

print(df)

print(df.groupby('prediction').count())

df['dif'] = df.doublet - df.singlet
print(df)

df.to_csv('doubletData.csv', index=False)
print("Saved to CSV")

sns.displot(df[df.prediction == 'doublet'], x = 'dif')

doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
print(doublets)

adata.obs['doublet'] = adata.obs.index.isin(doublets.index)

print(adata.obs)

adata.write(results_file)

adata = adata[~adata.obs.doublet]

print(adata)
