#! /usr/bin/python3

import fileinput
import sys
import math
import matplotlib.pyplot as plt


svxs  = []; svys  = []; svzs = []
mvxs  = []; mvys  = []; mvzs = []
for line in fileinput.input(sys.argv[1]):
    cols = line.split(',')
    if (cols[0] == 'Nbig'):
        i=0
        continue
    if  (cols[0] == 'Niter'):
        i=1
        continue
    nstr = cols[0+i]
    svx = float(cols[18+i]); svy = float(cols[19+i]); svz = float(cols[20+i])
    mvx = float(cols[27+i]); mvy = float(cols[28+i]); mvz = float(cols[29+i])
    svxs.append(svx); svys.append(svy); svzs.append(svz);
    mvxs.append(mvx); mvys.append(mvy); mvzs.append(mvz);



pltno = 0
plt.interactive(False)

lim=0.15
pltno += 1
plt.figure(pltno)
plt.title('Horizontal Variances, N=400,200,100')
plt.xlabel('MIG Variance (m^2)')
plt.ylabel('Projected Hourglass Variance (m^2)')
plt.axes().set_xlim([0.0, lim])
plt.axes().set_ylim([0.0, lim])
plt.axes().set_aspect('equal')
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.plot(mvxs, svxs, 'k>', fillstyle='none', label='X')
plt.plot(mvys, svys, 'k^', fillstyle='none', label='Y')
plt.plot([0.0,lim], [0.0,lim], 'k-');
plt.legend(numpoints=1, loc='upper left')
plt.savefig('../papers/fig_self_var_xy.png')

lim=0.7
pltno += 1
plt.figure(pltno)
plt.title('Vertical Variances, N=400,200,100')
plt.xlabel('MIG Variance (m^2)')
plt.ylabel('Projected Hourglass Variance (m^2)')
plt.axes().set_xlim([0.0, lim])
plt.axes().set_ylim([0.0, lim])
plt.axes().set_aspect('equal')
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.plot(mvzs, svzs, 'kD', fillstyle='none', label='Z', markersize=5)
plt.plot([0.0,lim], [0,lim], 'k-');
plt.legend(numpoints=1, loc='upper left')
plt.savefig('../papers/fig_self_var_z.png')



plt.show()
