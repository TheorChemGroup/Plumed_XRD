import string
import queue
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
import matplotlib
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
border_width = 2
peak_value = 10

time, cv = readFile("../../calculations/pull/COLVAR")
cv = list(map(np.log, cv))

theta_found, int_found = readFile("pulled.dat")
theta_ref, int_ref = readFile("ortho.dat")
theta_start, int_start = readFile("start.dat")

ref_big_peaks = countPeaks("ref_peaks.txt", peak_value)
cur_big_peaks = countPeaks("cur_peaks.txt", peak_value)
start_big_peaks = countPeaks("start_peaks.txt", peak_value)        

im_found = mpimg.imread('pulled.png')
im_ref = mpimg.imread('ortho.png')
im_start = mpimg.imread('start.png')

fig, axs = plt.subplot_mosaic("AB;CC;DE;FG", figsize=(9.6, 9.6))
plt.tight_layout()


for axis, letter in zip(axs.values(), string.ascii_lowercase):
    axis.annotate(f"{letter})", xycoords="axes fraction", xy=(0.02, 0.9), fontsize=12)

axs["A"].imshow(im_start, aspect='auto')
axs["A"].set_xticks([])                                    
axs["A"].set_yticks([])
axs["A"].set_xlabel("view along the\n{010} plane")

axs["B"].plot(theta_start, int_start, color=diff_color)
axs["B"].set_xlabel("2θ, degrees", labelpad=1)
axs["B"].set_ylabel("Intensity", labelpad=-4)
axs["B"].annotate(f"{start_big_peaks} peaks \nwith intensity > {peak_value}", xy=(0.5, 0.78), xycoords="axes fraction")

axs["C"].plot(time, cv, color=diff_color)
axs["C"].set_xlabel("Simulation time, ps")
axs["C"].set_ylabel("ln(CV)", labelpad=-1)

axs["D"].imshow(im_found, aspect='auto')
axs["D"].set_xticks([])                                    
axs["D"].set_yticks([])
axs["D"].set_xlabel("view along the\n{010} plane")

axs["E"].plot(theta_found, int_found, color=diff_color)
axs["E"].set_xlabel("2θ, degrees", labelpad=1)
axs["E"].set_ylabel("Intensity", labelpad=-4)
axs["E"].annotate(f"{cur_big_peaks} peaks \nwith intensity > {peak_value}", xy=(0.5, 0.78), xycoords="axes fraction")

axs["F"].imshow(im_ref, aspect='auto')
axs["F"].set_xticks([])                                       
axs["F"].set_yticks([])
axs["F"].set_xlabel("view along the\n{010} plane")

axs["G"].plot(theta_ref, int_ref, color=diff_color)
axs["G"].set_xlabel("2θ, degrees", labelpad=1)
axs["G"].set_ylabel("Intensity", labelpad=-4)
axs["G"].annotate(f"{ref_big_peaks} peaks \nwith intensity > {peak_value}", xy=(0.5, 0.78), xycoords="axes fraction")

boxes = {key: value.get_position() for key, value in axs.items()}

boxes["C"].x0 += 0.048
boxes["C"].y0 -= 0.03
boxes["C"].y1 -= 0.03
axs["C"].set_position(boxes["C"])

shift = 0.03
for letter in list(axs.keys())[3:]:
    boxes[letter].y0 -= shift
    boxes[letter].y1 -= shift
    axs[letter].set_position(boxes[letter])

rectA = drawEnclosingRect(boxes["A"].x0, boxes["A"].y0-0.03, boxes["B"].x1, boxes["A"].y1-0.01, linewidth=border_width, color=border_color)
rectAbox = rectA.get_bbox()


style = "Simple, tail_width=3, head_width=8, head_length=12"
kw = dict(arrowstyle=style, color=border_color)
arrow1 = patches.FancyArrowPatch((boxes["A"].x0+0.09, rectAbox.y0+0.002), (boxes["A"].x0+0.09, rectAbox.y0-0.035), **kw)
fig.add_artist(arrow1)

style = "Simple, tail_width=2, head_width=4, head_length=8"
kw = dict(arrowstyle=style, color=border_color)
arrow2 = patches.FancyArrowPatch((boxes["A"].x0+0.12, rectAbox.y0-0.036), (boxes["C"].x1-0.05, boxes["C"].y0+0.07), connectionstyle="arc3,rad=.1", **kw)
fig.add_artist(arrow2)
axs["C"].annotate("pulling process", xy=(0.5, 0.55), xycoords="axes fraction", ha="center", va="center", rotation=-10, weight="bold", fontsize=12)

shift = 0.07
for letter in list(axs.keys())[3:]:
    boxes[letter].y0 -= shift
    boxes[letter].y1 -= shift
    axs[letter].set_position(boxes[letter])

style = "Simple, tail_width=3, head_width=8, head_length=12"
kw = dict(arrowstyle=style, color=border_color)
arrow3 = patches.FancyArrowPatch((boxes["C"].x1-0.04, boxes["C"].y0+0.04), (boxes["C"].x1-0.04, boxes["E"].y1+0.009), **kw)
fig.add_artist(arrow3)


rectD = drawEnclosingRect(boxes["D"].x0, boxes["D"].y0-0.03, boxes["E"].x1, boxes["D"].y1-0.01, linewidth=border_width, color=border_color)

shift = 0.08
for letter in list(axs.keys())[5:]:
    boxes[letter].y0 -= shift
    boxes[letter].y1 -= shift
    axs[letter].set_position(boxes[letter])

rectF = drawEnclosingRect(boxes["F"].x0, boxes["F"].y0-0.03, boxes["G"].x1, boxes["F"].y1-0.01, linewidth=border_width, color=border_color)


axs["A"].annotate("Initial structure", xy=(0.483, 1.02 + 0.07 + 0.08 + 0.04), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")
axs["D"].annotate("Final structure", xy=(0.483, (1.02 + 0.07 + 0.08 + 0.04)/2 + 0.015), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")
axs["F"].annotate("Reference structure", xy=(0.483, (1.02 + 0.07 + 0.08 + 0.04)/4 - 0.01), xycoords="figure fraction", fontsize=18, ha="center", va="center", weight="bold")



fig.savefig("pull.svg", bbox_inches='tight', dpi=3000)
