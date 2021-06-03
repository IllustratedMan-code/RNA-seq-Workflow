import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
import pandas as pd
import numpy as np
import os
plt.style.use("ggplot")
os.chdir("/home/david/Documents/BenoitLab/RNA-seq")

data = pd.read_csv("edgeR-Genes/venn-regulation.csv")
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
        axes[i].bar(x, y, color=["r", "b"])
        axes[i].set_xticks(partdata["Names"])
        axes[i].set_ylim(0, maximum)
        if i != 0:
            axes[i].get_yaxis().set_visible(False)
        axes[i].set_xlabel(part + "\n" + pest)
        axes[0].set_ylabel("Gene Count")

        i += 1
fig.suptitle("Gene Regulation")
fig.savefig("figures/EdgeR/gene-regulation-bar-plot.png", bbox_inches="tight", dpi=250)
plt.show()

plt.style.use("ggplot")

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles


deet = pd.read_csv("edgeR-Genes/genelists/upDeet.csv", header=None)
deetcont = pd.read_csv("edgeR-Genes/genelists/downDeet.csv", header=None)
perm = pd.read_csv("edgeR-Genes/genelists/upPerm.csv", header=None)
permcont = pd.read_csv("edgeR-Genes/genelists/downPerm.csv", header=None)
genes = pd.read_csv("edgeR-Genes/genelists/listofgenes.csv", header=None)
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
fig.savefig("figures/EdgeR/DeetAndPermRegulation.png", bbox_inches="tight", dpi=250)
plt.show()

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
body = pd.read_csv("GeneReference/Body.csv")
leg = pd.read_csv("GeneReference/Leg.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(body), len(leg),len(genes)-len(leg) - len(body)), set_labels=("Body", "Leg"), ax=axes)
v.get_patch_by_id("10").set_color("blue")
v.get_patch_by_id("11").set_color("blue")
v.get_patch_by_id("01").set_color("grey")
v.get_patch_by_id("01").set_alpha(0.5)
axes.set_title("Leg vs Body")
fig.tight_layout()
fig.savefig("figures/figure2/Deseq-BodyvsLeg.png", bbox_inches="tight", dpi=250)
fig.savefig("figures/figure2/Deseq-BodyvsLeg.pdf", bbox_inches="tight", dpi=250)

plt.show()

genes = pd.read_csv("edgeR-Genes/genelists/listofgenes.csv", header=None)
body = pd.read_csv("edgeR-Genes/bodyExpr.csv")
upDeetBody = pd.read_csv("edgeR-Genes/genelists/upDeetBody.csv")
downDeetBody = pd.read_csv("edgeR-Genes/genelists/downDeetBody.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upDeetBody), len(downDeetBody), len(body)- len(upDeetBody)-len(downDeetBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Body Genes")
fig.tight_layout()
fig.savefig("figures/figure3/DeetBody.png")
plt.show()

leg = pd.read_csv("edgeR-Genes/legExpr.csv")
upDeetLeg = pd.read_csv("edgeR-Genes/genelists/upDeetLeg.csv")
downDeetLeg = pd.read_csv("edgeR-Genes/genelists/downDeetLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upDeetLeg), len(downDeetLeg), len(leg)-len(upDeetLeg)-len(downDeetLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Leg Genes")
fig.tight_layout()
fig.savefig("figures/figure3/DeetLeg.png")

leg = pd.read_csv("edgeR-Genes/legExpr.csv")
upPermLeg = pd.read_csv("edgeR-Genes/genelists/upPermLeg.csv")
downPermLeg = pd.read_csv("edgeR-Genes/genelists/downPermLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPermLeg), len(downPermLeg), len(leg)-len(downPermLeg)-len(upPermLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Leg Genes")
fig.tight_layout()
fig.savefig("figures/EdgeR/PermLeg.png")

body = pd.read_csv("edgeR-Genes/bodyExpr.csv")
upPermBody = pd.read_csv("edgeR-Genes/genelists/upPermBody.csv")
downPermBody = pd.read_csv("edgeR-Genes/genelists/downPermLeg.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPermBody), len(downPermBody), len(body)-len(upPermBody)-len(downPermBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Body Genes")
fig.tight_layout()
fig.savefig("figures/EdgeR/PermBody.png")

genes = pd.read_csv("GeneReference/listofgenes.csv")
upPerm = pd.read_csv("edgeR-Genes/genelists/upPerm.csv")
downPerm = pd.read_csv("edgeR-Genes/genelists/downPerm.csv")

fig, axes = plt.subplots()
v = venn2_unweighted(subsets=(len(upPerm), len(downPerm), len(genes)-len(upPerm)-len(downPerm)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)

axes.set_title("Perm Expression")
fig.tight_layout()
fig.savefig("figures/EdgeR/Perm.png")

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
body = pd.read_csv("edgeR-Genes/bodyExpr.csv")
leg = pd.read_csv("edgeR-Genes/legExpr.csv")
fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(body), len(leg),len(genes)-len(leg) - len(body)), set_labels=("Body", "Leg"), ax=axes)
v.get_patch_by_id("10").set_color("blue")
v.get_patch_by_id("11").set_color("blue")
v.get_patch_by_id("01").set_color("grey")
v.get_patch_by_id("01").set_alpha(0.5)
axes.set_title("Leg vs Body")
fig.tight_layout()
fig.savefig("figures/figure2/EdgeR-BodyvsLeg.png", bbox_inches="tight", dpi=250)
fig.savefig("figures/figure2/EdgeR-BodyvsLeg.pdf", bbox_inches="tight", dpi=250)

plt.show()

plt.style.use("ggplot")

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles


deet = pd.read_csv("Deseq-Genes/upDeet.csv")
deetcont = pd.read_csv("Deseq-Genes/downDeet.csv")
perm = pd.read_csv("Deseq-Genes/upPerm.csv")
permcont = pd.read_csv("Deseq-Genes/downPerm.csv")
genes = pd.read_csv("GeneReference/listofgenes.csv")


def compare(s1, s2):
    c = len(list(set(s1.iloc[:, 0]) & set(s2.iloc[:, 0])))
    return c


diff = compare(deet, genes)
len(deet)
fig, axes = plt.subplots(2)
v1 = venn2_unweighted(
    subsets=(len(deet), len(deetcont), len(genes) - len(deetcont) - len(deet)),
    set_labels=("Up Regulated", "Down Regulated"),
    ax=axes[0],
)
v2 = venn2_unweighted(
    subsets=(len(perm), len(permcont), len(genes) - len(permcont) - len(perm)),
    set_labels=("Up Regulated", "Down Regulated"),
    ax=axes[1],
)
axes[0].set_title("Deet")
axes[1].set_title("Perm")
fig.tight_layout()
fig.savefig("figures/Deseq2/DeetAndPermRegulation.png", bbox_inches="tight", dpi=250)
fig.savefig("figures/Deseq2/DeetAndPermRegulation.pdf", bbox_inches="tight", dpi=250)
plt.show()

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
body = pd.read_csv("Deseq-Genes/Body.csv")
upDeetBody = pd.read_csv("Deseq-Genes/upDeetBody.csv")
downDeetBody = pd.read_csv("Deseq-Genes/downDeetBody.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upDeetBody), len(downDeetBody), len(body)- len(upDeetBody)-len(downDeetBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Body Genes")
fig.tight_layout()
fig.savefig("figures/figure3/DeetBodyDeseq.png")
fig.savefig("figures/figure3/DeetBodyDeseq.pdf")
plt.show()

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
body = pd.read_csv("Deseq-Genes/Body.csv")
upPermBody = pd.read_csv("Deseq-Genes/upPermBody.csv")
downPermBody = pd.read_csv("Deseq-Genes/downPermBody.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upPermBody), len(downPermBody), len(body)- len(upPermBody)-len(downPermBody)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Body Genes")
fig.tight_layout()
fig.savefig("figures/figure3/PermBodyDeseq.png")
fig.savefig("figures/figure3/PermBodyDeseq.pdf")
plt.show()

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
Leg = pd.read_csv("Deseq-Genes/Leg.csv")
upDeetLeg = pd.read_csv("Deseq-Genes/upDeetLeg.csv")
downDeetLeg = pd.read_csv("Deseq-Genes/downDeetLeg.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upDeetLeg), len(downDeetLeg), len(Leg)- len(upDeetLeg)-len(downDeetLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Deet Expression in Leg Genes")
fig.tight_layout()
fig.savefig("figures/figure3/DeetLegDeseq.png")
fig.savefig("figures/figure3/DeetLegDeseq.pdf")
plt.show()

genes = pd.read_csv("GeneReference/listofgenes.csv", header=None)
Leg = pd.read_csv("Deseq-Genes/Leg.csv")
upPermLeg = pd.read_csv("Deseq-Genes/upPermLeg.csv")
downPermLeg = pd.read_csv("Deseq-Genes/downPermLeg.csv")

fig, axes = plt.subplots()

v = venn2_unweighted(subsets=(len(upPermLeg), len(downPermLeg), len(Leg)- len(upPermLeg)-len(downPermLeg)), set_labels=("Up Regulated", "Down Regulated"), ax=axes)
axes.set_title("Perm Expression in Leg Genes")
fig.tight_layout()
fig.savefig("figures/figure3/PermLegDeseq.png")
fig.savefig("figures/figure3/PermLegDeseq.pdf")
plt.show()
