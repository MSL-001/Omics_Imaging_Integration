import matplotlib.pyplot as plt
import math
from PIL import Image
import numpy as np
from matplotlib.pyplot import savefig

color = "Reds"

met_ff_height = np.asarray(Image.open("C:/Users/matti/Downloads/met_ff_height_heatmap_"+color+".png"))
met_vol_height = np.asarray(Image.open("C:/Users/matti/Downloads/met_vol_height_heatmap_"+color+".png"))
prot_ff_height = np.asarray(Image.open("C:/Users/matti/Downloads/prot_ff_height_heatmap_"+color+".png"))
prot_vol_height = np.asarray(Image.open("C:/Users/matti/Downloads/prot_vol_height_heatmap_"+color+".png"))

met_ff = np.asarray(Image.open("C:/Users/matti/Downloads/met_ff_heatmap_"+color+".png"))
met_vol = np.asarray(Image.open("C:/Users/matti/Downloads/met_vol_heatmap_"+color+".png"))
prot_ff = np.asarray(Image.open("C:/Users/matti/Downloads/prot_ff_heatmap_"+color+".png"))
prot_vol = np.asarray(Image.open("C:/Users/matti/Downloads/prot_vol_heatmap_"+color+".png"))

fig, axs = plt.subplots(2, 4,figsize=(12, 3), dpi=300)
plt.subplots_adjust(wspace=0, hspace=0)
axs[0, 0].imshow(prot_vol)
axs[0, 1].imshow(prot_ff)
axs[0, 2].imshow(prot_vol_height)
axs[0, 3].imshow(prot_ff_height)
axs[1, 0].imshow(met_vol)
axs[1, 1].imshow(met_ff)
axs[1, 2].imshow(met_vol_height)
axs[1, 3].imshow(met_ff_height)

for ax in axs.flat:
    for ax in axs.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)


for ax in axs.flat:
    ax.label_outer()
    ax.xaxis.set_label_position('top')

row_labels = ["Proteomics", "Metabolomics"]

column_labels = ["VOL", "FF", "VOL", "FF"]

for i, label in enumerate(row_labels):
    axs[i, 0].set_ylabel(label, fontsize=10, rotation=90)

for i, label in enumerate(column_labels):
    axs[0, i].set_xlabel(label, fontsize=10)

# after creating fig, axs

fig.canvas.draw()  # ensures positions are updated

# Get positions of the top row axes
pos1 = axs[0, 0].get_position()
pos2 = axs[0, 1].get_position()
pos3 = axs[0, 2].get_position()
pos4 = axs[0, 3].get_position()

# Centers of column groups
center_12 = (pos1.x0 + pos2.x1) / 2
center_34 = (pos3.x0 + pos4.x1) / 2

# Y position slightly above the top row
y = pos1.y1 + 0.07

# Add group labels
fig.text(center_12, y, "", ha='center', va='bottom', fontsize=16)
fig.text(center_34, y, "Height Adjusted", ha='center', va='bottom', fontsize=16)



plt.savefig(
    "C:/Users/matti/Downloads/combined_heatmaps_"+color+".png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0.05
)

plt.show()
