#! /usr/bin/python3

import fileinput
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def list_flt(vals):
    a = []
    for val in vals:
        a.append(float(val))
    return a

def list_avg(vals):
    sum = 0.0
    for val in vals:
        sum += float(val)
    return sum/len(vals)


plt.interactive(False)
i=0
prevprod = 1000000.0
makeitblack = False
madeitblack = False
for line in fileinput.input(sys.argv[1]):
    cols = line.split(',');
    if (cols[0] == 'LAM'):
        continue
    majr = float(cols[4])
    minr = float(cols[7])
    a = float(cols[10]) * 180 / math.pi
    z = float(cols[11])
    n = int(cols[12])
    xvals = list_flt(cols[13:13+n])
    yvals = list_flt(cols[13+n:13+2*n])
    xavg = list_avg(xvals)
    yavg = list_avg(yvals)

    # make an ellipse object
    factor = 2 / math.sqrt(n)
    #factor *= 5.991 # 95%
    factor *= 4.605 # 90%
    e = Ellipse(xy=[xavg,yavg], width=factor*majr, height=factor*minr, angle=a)
    thisprod = majr * minr
    if thisprod>prevprod:
        if not madeitblack: 
            makeitblack = True
    if (makeitblack):
        e.set_alpha(1.0) # 0=white, 1=black
    else:
        e.set_alpha(.25) # 0=white, 1=black
    e.set_facecolor([1,1,1]) # white insides
    prevprod = thisprod

    # set up a plot; I'm not sure what 'figure' and 'subplot' are
    plt.figure(1).clear()
    plt.subplot(111).add_artist(e)
    plt.subplot(111).set_xlim(-15,15)
    plt.subplot(111).set_ylim(-15,15)

    # plot the stuff
    plt.plot(xvals, yvals, 'kx') # black x's
    plt.plot(0, 0, 'r+')
    #plt.show()

    nholds = 1
    if makeitblack:
        nholds = 10
    for holdi in (range(nholds)):
        stri = str(i)
        i += 1
        while len(stri) < 4:
            stri = '0' + stri
        plt.savefig('ell'+stri+'.png');
        print('ell'+stri+'.png');
    if makeitblack:
        makeitblack = False
        madeitblack = True
