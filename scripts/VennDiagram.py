# %% imports 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import pandas as pd
import numpy as np
import os

os.chdir("/home/david/Documents/Python")

# %%fig1
data = pd.read_csv("regulation.csv")
plt.style.use("ggplot")
fig, axes = plt.subplots(1, 4)


i = 0

maximum = data["GeneCount"].max()

for pest in data["Part"].unique():
    pestdata = data[data["Part"] == pest]
    for part in pestdata["Pest"].unique():
        partdata = pestdata[pestdata["Pest"] == part]
        x = partdata["Names"]
        y = partdata["GeneCount"]
        axes[i].bar(x,y, color=["r", "b"])
        axes[i].set_xticks(partdata["Names"])
        axes[i].set_ylim(0, maximum)
        if i != 0:
            axes[i].get_yaxis().set_visible(False)
        axes[i].set_xlabel(part + "\n" + pest)
        axes[0].set_ylabel("Gene Count")

        i += 1
fig.suptitle("Gene Regulation")
fig.savefig("fig.png", bbox_inches="tight", dpi=250)

#%% fig2
plt.style.use("ggplot")

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles


deet = pd.read_csv("Genes/deet.csv", header=None)
deetcont = pd.read_csv("Genes/deetcont.csv", header=None)
perm = pd.read_csv("Genes/perm.csv", header=None)
permcont = pd.read_csv("Genes/permcont.csv", header=None)
genes = pd.read_csv("Genes/listofgenes.csv", header=None)
def compare(s1, s2):
    c = len(list(set(s1.iloc[:, 0]) & set(s2.iloc[:, 0])))
    return(c)

diff = compare(deet, genes)
len(deet)
fig, axes = plt.subplots(2)
v1 = venn2_unweighted(subsets=(len(deet), len(deetcont), len(genes)-len(deetcont)-len(deet)), set_labels=("Up Regulated", "Down Regulated"), ax=axes[0])
v2 = venn2_unweighted(subsets=(len(perm), len(permcont), len(genes)-len(permcont)-len(perm)), set_labels=("Up Regulated", "Down Regulated"), ax=axes[1])
axes[0].set_title("Deet")
axes[1].set_title("Perm")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/DeetnPerm.png", bbox_inxhes="tight", dpi=250)
plt.show()

#%% fig3

genes = pd.read_csv("Genes/listofgenes.csv", header=None)
body = pd.read_csv("Genes/Body.csv")
leg = pd.read_csv("Genes/Leg.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(body), len(leg),len(genes)-len(leg) - len(body)), set_labels=("Body", "Leg"), ax=axes)

axes.set_title("Leg vs Body")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/BodyvsLeg.png", bbox_inxhes="tight", dpi=250)


plt.show()

#%% fig4

genes = pd.read_csv("Genes/listofgenes.csv", header=None)
body = pd.read_csv("Genes/Body.csv")
upDeetBody = pd.read_csv("Genes/Genelists/upDeetBody.csv")
downDeetBody = pd.read_csv("Genes/Genelists/downDeetBody.csv")


fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upDeetBody), len(downDeetBody), len(body)- len(upDeetBody)-len(downDeetBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Body Genes")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/DeetBody.png")
plt.show()


#%% DeetLeg
leg = pd.read_csv("Genes/Leg.csv")
upDeetLeg = pd.read_csv("Genes/Genelists/upDeetLeg.csv")
downDeetLeg = pd.read_csv("Genes/Genelists/downDeetLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upDeetLeg), len(downDeetLeg), len(leg)-len(upDeetLeg)-len(downDeetLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Leg Genes")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/DeetLeg.png")


#%% PermLeg
leg = pd.read_csv("Genes/Leg.csv")
upPermLeg = pd.read_csv("Genes/Genelists/upPermLeg.csv")
downPermLeg = pd.read_csv("Genes/Genelists/downPermLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPermLeg), len(downPermLeg), len(leg)-len(downPermLeg)-len(upPermLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Leg Genes")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/PermLeg.png")


#%% PermBody

body = pd.read_csv("Genes/Body.csv")
upPermBody = pd.read_csv("Genes/Genelists/upPermBody.csv")
downPermBody = pd.read_csv("Genes/Genelists/downPermLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPermBody), len(downPermBody), len(body)-len(upPermBody)-len(downPermBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Body Genes")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/PermBody.png")


#%% Perm Expresiiion

genes = pd.read_csv("Genes/listofgenes.csv")
upPerm = pd.read_csv("Genes/Genelists/upPerm.csv")
downPerm = pd.read_csv("Genes/Genelists/downPerm.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPerm), len(downPerm), len(genes)-len(upPerm)-len(downPerm)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression")
fig.tight_layout()
fig.savefig("VennDiagrams/RevisedDiagrams/Perm.png")

