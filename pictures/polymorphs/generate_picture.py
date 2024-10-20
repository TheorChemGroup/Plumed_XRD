import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
import string
from matplotlib import cm, colors, patches

def readFile(name: str, func=float):
    x, y = [], []
    with open(name, 'r') as file:
        for line in file:
            if '#' not in line:
                line = [func(i) for i in line.split()]
                x.append(line[0])
                y.append(line[1])
    return x, y
    
    
def countPeaks(name: str, value: float):
    number_of_big_peaks = 0
    with open(name, 'r') as file:
        for line in file:
            if '#' not in line and float(line.split()[5]) > 10:
                number_of_big_peaks += 1
    return number_of_big_peaks


def drawEnclosingRect(x0, y0, x1, y1, x_offset=0.01, y_offset=0.02, color="black", linewidth=3):
    rect = patches.Rectangle((x0-x_offset, y0-y_offset), x1-x0 + x_offset*2, y1-y0 + y_offset*2, edgecolor=color, linewidth=linewidth, fill=False)
    fig.add_artist(rect)
    return rect


diff_color="navy"
border_color="dimgray"
norm = colors.Normalize(0, 100)
border_width = 2
peak_value = 10
cmap = cm.get_cmap("Blues")


theta_mono, int_mono = readFile("mono.dat")
theta_ortho, int_ortho = readFile("ortho.dat")

im_mono = mpimg.imread("mono.png")
im_ortho = mpimg.imread("ortho.png")

fig, axs = plt.subplot_mosaic("AB;CD", figsize=(6.4,4.8))
plt.tight_layout()

for axis, letter in zip(axs.values(), string.ascii_lowercase):
    axis.annotate(f"{letter})", xycoords="axes fraction", xy=(0.02, 0.9), fontsize=12)

axs["A"].plot(theta_mono, int_mono, color=diff_color)
axs["A"].set_xlabel("2θ, degrees", labelpad=1)
axs["A"].set_ylabel("Intensity", labelpad=-4)

axs["B"].imshow(im_mono)
axs["B"].set_xticks([])                                    
axs["B"].set_yticks([])  
axs["B"].set_xlabel("view along the {100} plane\nSpace group: P21/c")

axs["C"].plot(theta_ortho, int_ortho, color=diff_color)
axs["C"].set_xlabel("2θ, degrees", labelpad=1)
axs["C"].set_ylabel("Intensity", labelpad=-4)

axs["D"].imshow(im_ortho)
axs["D"].set_xticks([])                                    
axs["D"].set_yticks([])  
axs["D"].set_xlabel("view along the {010} plane\nSpace group: Pbca")

boxes = {key: value.get_position() for key, value in axs.items()}

boxes["C"].y0 -= 0.11
boxes["C"].y1 -= 0.11

boxes["B"].x0 -= 0.03
boxes["B"].x1 -= 0.03
boxes["D"].x0 -= 0.03
boxes["D"].x1 -= 0.03

boxes["B"].y0 = boxes["A"].y0
boxes["B"].y1 = boxes["A"].y1
boxes["D"].y0 = boxes["C"].y0
boxes["D"].y1 = boxes["C"].y1


rectA = drawEnclosingRect(boxes["A"].x0-0.07, boxes["B"].y0-0.08, boxes["D"].x1+0.01, boxes["A"].y1, color=border_color)
rectC = drawEnclosingRect(boxes["C"].x0-0.07, boxes["D"].y0-0.08, boxes["D"].x1+0.01, boxes["C"].y1, color=border_color)

for axis in axs:
    axs[axis].set_position(boxes[axis])

axs["A"].annotate("Monoclinic polymorph (Form I)", xy=(0.5, 1.156), xycoords="figure fraction", fontsize=16, ha="center", va="center", weight="bold")
axs["C"].annotate("Orthorhombic polymorph (Form II)", xy=(0.5, 1.167/2-0.022), xycoords="figure fraction", fontsize=16, ha="center", va="center", weight="bold")

fig.savefig("diff.svg", bbox_inches='tight', dpi=3000)
