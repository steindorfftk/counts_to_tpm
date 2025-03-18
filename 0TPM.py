import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

tab_file = 'SRR5926242.tabular'

df = pd.read_csv(tab_file,sep='\t')
dfb = pd.read_csv('lenghts.txt',sep=',')

df = df.dropna(subset=["GeneSymbol"])

# Get unique values
unique_df = set(df['GeneSymbol'])
unique_dfb = set(dfb.iloc[:, 0])  # First column of dfb

# Filter df
df = df[df['GeneSymbol'].isin(unique_dfb)]

df = df.merge(dfb, left_on='GeneSymbol', right_on='gene', how='left')

df.drop(columns=['gene'], inplace=True)

df["RPK"] = (df["counts"] / df['len'])

rpk_sum = df['RPK'].sum()
scaling = rpk_sum / 1000000

df["TPM"] = (df["RPK"] / scaling)

total_counts = df['counts'].sum()

df["RPM"] = (df["counts"] / total_counts) * 10**6

df = df.sort_values(by='TPM', ascending=False)

df.to_csv(f'{tab_file[:-8]}_TPM.tabular')
