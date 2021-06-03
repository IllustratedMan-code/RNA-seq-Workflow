import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib
import pandas as pd

plt.rcParams["font.family"] = "Bahnschrift"
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    "", ["dodgerblue", "lightgrey"]
)
# pandas section
dotplotdata = pd.read_csv(
    "../GX2C9-Z99E-7HeneReference/dotplotdata.csv"
)  # infers header existance
#
# these are all pandas dataframes as well
pesticide = dotplotdata.Category
pathway = dotplotdata.Pathway
enrichmentratio = dotplotdata.Enrichment
pval = dotplotdata["P-Value"]
# have to index it like this because of the dash in P-Value
# matplotlib section
# The reason you couldn't multiply the enrichment ratio here before is because you were using python lists
# instead of a numpy array or pandas dataframe. Python lists will just replicate the elements [1 , 2] * 2 -> [1, 2, 1, 2]
fig, ax = plt.subplots(1)
p = ax.scatter(
    pathway, pesticide, c=pval, s=enrichmentratio * 3, cmap=cmap, vmin=0, vmax=0.05
)
plt.legend(loc="lower left", markerscale=2.0, scatterpoints=1, fontsize=10)
plt.xticks(rotation="vertical")
# Legend:
handles, labels = p.legend_elements(
    prop="sizes",
    alpha=0.6,
    num=4,
    func=lambda x: x / 3
    # use the lambda function to reverse any scaling done to the data. Erichmentratio is still multiplied by 3 within the p scatter object
    # Lamdas are just short function definitions with a slightly different syntax
    # you could use a regular function here instead
)
legend2 = ax.legend(
    handles,  # handles are not needed
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

plt.savefig("dotplot.pdf", dpi=300, bbox_inches="tight")
# plt.savefig('dotplot.jpg', dpi=300, bbox_inches='tight')
plt.show()

# Comments:
# The figure becomes janky when creating 2 subplots, so I left it as 1 large plot that I split manually
# using image/pdf editor at the end. So, make sure you also add the Y-axis "group" labels (subplot1 & 2,
# such as "Upregulated" and "Downregulated") when editing the figure!
