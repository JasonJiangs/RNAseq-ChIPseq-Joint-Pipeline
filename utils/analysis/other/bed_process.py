import pandas as pd

# read bed file
bed = pd.read_csv('result/promoters_regions_UD10kb.bed', sep='\t', header=None)

# replace all negative values with 1 in second column
bed[1] = bed[1].apply(lambda x: 1 if x < 1 else x)

# drop last column
bed = bed.drop(columns=[3])

# save
bed.to_csv('result/promoters_regions_UD10kb_Mod.bed', sep='\t', header=False, index=False)
