#!/usr/bin/python
#

from pylab import *
import glob
import numpy as np
import sys
from optparse import OptionParser

parser = OptionParser(usage="draw boxplot")
parser.add_option("-1", "--label1", help="name for legend1", dest="label1", default="RPV3C")
parser.add_option("-2", "--label2", help="name for legend2", dest="label2", default="RPV3E")


(opts, args) = parser.parse_args()

g1 = glob.glob(args[0])
g2 = glob.glob(args[1])


# function for setting the colors of the box plots pairs
spacing = 7

def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')

# Some fake data to plot

fig = figure()
ax = axes()
hold(True)

ymin = 1e4
ymax = 0

ndata = len(g1)
xticks = []
xticks_labels = []
for idx in range(0, ndata):
  f1 = g1[idx]
  f2 = g2[idx]

# first boxplot pair
  data1 = np.loadtxt(f1)
  data1 = data1[:,1:].flatten()

  data2 = np.loadtxt(f2)
  data2 = data2[:,1:].flatten()

  print np.mean(data1), np.std(data1)
  print np.mean(data2), np.std(data2)

  ymin = min(np.amin(np.hstack((data1,data2))) - 100, ymin)
  ymax = max(np.amax(np.hstack((data1,data2))) + 100, ymax)

  A= [ data1, data2 ]

  bp = boxplot(A, positions = [1.0/3*spacing+idx*spacing, 2.0/3*spacing+idx*spacing], widths = 0.6/3*spacing)
  setBoxColors(bp)
  xticks.append(1.5/3*spacing + idx * spacing)
  xticks_labels.append("R_%02d" % (idx + 1))

# second boxplot pair
#bp = boxplot(B, positions = [4, 5], widths = 0.6)
#setBoxColors(bp)

# thrid boxplot pair
#bp = boxplot(C, positions = [7, 8], widths = 0.6)
#setBoxColors(bp)

# set axes limits and labels
xlim(0,max(xticks)+1.0/3*spacing)
ylim(ymin-100,ymax+100)
#ax.set_xticklabels(['A', 'B', 'C'])
ax.set_xticks(xticks)
ax.set_xticklabels(xticks_labels)

# draw temporary red and blue lines and use them to create a legend
hB, = plot([1,1],'b-')
hR, = plot([1,1],'r-')
legend((hB, hR),(opts.label1, opts.label2))
hB.set_visible(False)
hR.set_visible(False)

savefig(args[2])
show()

