import pandas as pd

##########
#import all data
##########
ela_archaea=pd.read_csv("Strict_US_trained_Features_Archaea.tsv",sep="\t")
ela_preterm=pd.read_csv("Strict_US_trained_Features.tsv",sep="\t")
print(len(ela_archaea["feature"]))
print(len(ela_preterm["feature"]))
overlap=ela_archaea[ela_archaea["feature"].isin(ela_preterm["feature"])]
#overlap.to_csv("overlap_taxa_ElasticNet.tsv",sep="\t",index=True)
mas_preterm=pd.read_csv("preterm_sig.tsv",sep=",")
mas_archaea=pd.read_csv("archaea_sig.tsv",sep=",")

###########
#Rename features
###########
ela_archaea["feature_rename"] = ela_archaea["feature"].str.replace(";", ".", regex=False)
print(ela_archaea["feature_rename"].head)
print(ela_preterm["feature"].head)
print(mas_archaea.iloc[:,0].head)
print(mas_preterm.iloc[:,0].head)


###########
#Venn diagram
###########
from venn import venn
import matplotlib.pyplot as plt
sets={
    "Archaea_EN":set(ela_archaea["feature_rename"]),
    "Preterm_EN":set(ela_preterm["feature"]),
    "Archaea_MaAsLin2":set(mas_archaea.iloc[:,0]),
    "Preterm_MaAsLin2":set(mas_preterm.iloc[:,0])
}
venn(sets)
intersection = sets["Archaea_MaAsLin2"] & sets["Preterm_EN"]
print(intersection)
plt.show()
