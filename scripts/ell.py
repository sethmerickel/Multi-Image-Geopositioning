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
freeze = False
froze  = False
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
    factor = 4.605;
    #factor = 5.991 # 95%
    e_big = Ellipse(xy=[xavg,yavg], width=factor*majr, height=factor*minr, angle=a)
    e_big.set_alpha(0.25)        # grey line
    e_big.set_facecolor([1,1,1]) # white insides

    factor /= math.sqrt(n)
    e_sma = Ellipse(xy=[xavg,yavg], width=factor*majr, height=factor*minr, angle=a)
    e_sma.set_alpha(1.0) # 0=white, 1=black
    e_sma.set_facecolor([1,1,1]) # white insides

    # set up a plot; I'm not sure what 'figure' and 'subplot' are
    plt.figure(1).clear()
    plt.subplot(111).add_artist(e_big)
    plt.subplot(111).add_artist(e_sma)
    plt.subplot(111).set_xlim(-15,15)
    plt.subplot(111).set_ylim(-15,15)

    # plot the stuff
    plt.plot(xvals, yvals, 'kx') # black x's
    plt.plot(0, 0, 'r+')
    plt.plot(xavg, yavg, 'b+')
    #plt.show()

    thisprod = majr * minr
    if thisprod>prevprod and not froze:
        freeze = True
    prevprod = thisprod

    nholds = 1
    if freeze:
        nholds = 10
    for holdi in (range(nholds)):
        stri = str(i)
        i += 1
        while len(stri) < 4:
            stri = '0' + stri
        plt.savefig('ell'+stri+'.png');
        print('ell'+stri+'.png');
    if freeze:
        freeze = False
        froze  = True
