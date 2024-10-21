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
cmap = cm.get_cmap("Blues")
norm = colors.Normalize(0, 100)
border_width = 2
peak_value = 10

time, energy = [], []
with open("global_diff", 'r') as file:
    for i, line in enumerate(file):
        line = line.split()
        time.append((i+1)*2*25/1000)
        energy.append(float(line[0]))
        

patterns = []
with open('diffractions_push', 'r') as file:
    count = -1
    for line in file:
        if 'FRAME' in line:
            count += 1
            patterns.append([])
            continue
        patterns[count].append(float(line))	



fig, axs = plt.subplot_mosaic("AB;CC;DD;EF", figsize=(9.6, 12.8))
plt.tight_layout()

for axis, letter in zip(axs.values(), string.ascii_lowercase):
    axis.annotate(f"{letter})", xycoords="axes fraction", xy=(0.02, 0.9), fontsize=12)

theta_end, int_end = readFile("pushed.dat")
theta_ref, int_ref = readFile("mono.dat")

ref_big_peaks = countPeaks("ref_peaks.txt", peak_value)
end_big_peaks = countPeaks("end_peaks.txt", peak_value)

im_push = mpimg.imread('pushed.png')
im_ref = mpimg.imread('mono.png')

axs["A"].imshow(im_ref)
axs["A"].set_xticks([])                                    
axs["A"].set_yticks([])
axs["A"].set_xlabel("view along the\n{100} plane")

axs["B"].plot(theta_ref, int_ref, color=diff_color)
axs["B"].set_xlabel("2θ, degrees", labelpad=1)
axs["B"].set_ylabel("Intensity", labelpad=-4)
axs["B"].annotate(f"{ref_big_peaks} peaks \nwith intensity > {peak_value}", xy=(0.6, 0.78), xycoords="axes fraction")

axs["C"].imshow(patterns, cmap=cmap, interpolation='spline36', aspect='auto', extent=[5, 50, 6.6, 0])
axs["C"].set_xlabel("2θ, degrees", labelpad=1)
axs["C"].set_ylabel("Simulation time, ps", labelpad=10)
axs["C"].annotate("", xy=(-0.04, -0.02), xytext=(-0.04, 1.02), xycoords='axes fraction', arrowprops=dict(arrowstyle="->"))

axs["D"].plot(time, energy, color=diff_color)
axs["D"].set_xlabel("Simulation time, ps")
axs["D"].set_ylabel("Difference", labelpad=1)

axs["E"].imshow(im_push)
axs["E"].set_xticks([])                                    
axs["E"].set_yticks([])
axs["E"].set_xlabel("view along the\n{100} plane")

axs["F"].plot(theta_end, int_end, color=diff_color)
axs["F"].set_xlabel("2θ, degrees", labelpad=1)
axs["F"].set_ylabel("Intensity", labelpad=-4)
axs["F"].annotate(f"{end_big_peaks} peaks \nwith intensity > {peak_value}", xy=(0.6, 0.78), xycoords="axes fraction")




boxes = {key: value.get_position() for key, value in axs.items()}
shift = -0.00865
boxes["B"].y0 -= shift
boxes["B"].y1 += shift
axs["B"].set_position(boxes["B"])
boxes["F"].y0 -= shift
boxes["F"].y1 += shift
axs["F"].set_position(boxes["F"])
boxes["C"].x0 += 0.05
boxes["C"].y0 -= 0.01
boxes["C"].y1 -= 0.01
axs["C"].set_position(boxes["C"])
boxes["D"].x0 += 0.05
boxes["D"].x1 -= 0.1
boxes["D"].y0 -= 0.01
boxes["D"].y1 -= 0.01
axs["D"].set_position(boxes["D"])


shift = 0.04
boxes["E"].y0 -= shift
boxes["E"].y1 -= shift
boxes["F"].y0 -= shift
boxes["F"].y1 -= shift
axs["E"].set_position(boxes["E"])
axs["F"].set_position(boxes["F"])


rectA = drawEnclosingRect(boxes["A"].x0, boxes["A"].y0-0.015, boxes["B"].x1, boxes["A"].y1-0.0115, color=border_color, linewidth=border_width)
rectF = drawEnclosingRect(boxes["E"].x0, boxes["E"].y0-0.015, boxes["F"].x1, boxes["E"].y1-0.0115, color=border_color, linewidth=border_width)
rectAbox = rectA.get_bbox()
rectFbox = rectF.get_bbox()
axs["A"].annotate("Initial structure", xy=(0.483, 1.00 + shift), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")
# axs["D"].annotate("Final structure", xy=(0.483, (1.02 + 0.07 + 0.08 + 0.04)/2 + 0.015), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")
axs["F"].annotate("Final structure", xy=(0.483, (1.00+shift)/4), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")

style = "Simple, tail_width=3, head_width=8, head_length=12"
kw = dict(arrowstyle=style, color=border_color)
arrow1 = patches.FancyArrowPatch((rectAbox.x1-0.05, rectAbox.y0+0.002), (rectAbox.x1-0.05, rectFbox.y1), **kw)
fig.add_artist(arrow1)

axs["C"].annotate("pushing process", xy=(1.22, -0.18), xycoords="axes fraction", ha="center", va="center", rotation=90, weight="bold", fontsize=14)

style = "Simple, tail_width=2, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color=border_color)
arrow2 = patches.FancyArrowPatch((boxes["D"].x0+0.05, boxes["D"].y0+0.01), (boxes["D"].x1-0.03, boxes["D"].y1-0.105), connectionstyle="arc3,rad=-.1", **kw)
fig.add_artist(arrow2)
axs["D"].annotate("pushing process", xy=(0.5, 0.3), xycoords="axes fraction", ha="center", va="center", rotation=10, weight="bold", fontsize=12)


fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=[axs["C"]], label="Intensity")

fig.savefig("push.svg", bbox_inches='tight', dpi=3000)
