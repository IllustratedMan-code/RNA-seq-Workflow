import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

plt.rcParams["font.family"] = "Bahnschrift"
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    "", ["dodgerblue", "lightgrey"]
)

dotplotdata = pd.read_csv("../GeneReference/dotplotdata.csv")  # infers header existance

pesticide = dotplotdata.Category
pathway = dotplotdata.Pathway
enrichmentratio = dotplotdata.Enrichment
pval = dotplotdata["P-Value"]

fig, ax = plt.subplots(1)
p = ax.scatter(
    pathway, pesticide, c=pval, s=enrichmentratio * 3, cmap=cmap, vmin=0, vmax=0.05
)
plt.legend(loc="lower left", markerscale=2.0, scatterpoints=1, fontsize=10)
plt.xticks(rotation="vertical")

# Legend:
handles, labels = p.legend_elements(prop="sizes", alpha=0.6, num=4, func=lambda x : x/3)
legend2 = ax.legend(
    handles,
    labels,
    loc="upper right",
    title="Enrichment Ratio",
    bbox_to_anchor=(1.4, 1.1),
)

# Colorbar:
cbar = fig.colorbar(p, shrink=0.5)
cbar.ax.set_ylabel("p Value", rotation=270, labelpad=15)

# Title, labels, and visual modifications:
plt.title("Metabolic Pathway Enrichment")
plt.xlabel("Metabolic Pathway")
plt.ylabel("Pesticide")
ax.grid(axis="x", linestyle="dotted")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.show()

plt.savefig('dotplot.pdf', dpi=300, bbox_inches='tight')
plt.savefig('dotplot.jpg', dpi=300, bbox_inches='tight')
