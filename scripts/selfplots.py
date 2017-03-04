#! /usr/bin/python3

import fileinput
import sys
import math
import matplotlib.pyplot as plt

svxs   = []; svys   = [];  svzs = []
mvxs   = []; mvys   = [];  mvzs = []
svxs4  = []; svys4  = []; svzs4 = []
mvxs4  = []; mvys4  = []; mvzs4 = []
svxs2  = []; svys2  = []; svzs2 = []
mvxs2  = []; mvys2  = []; mvzs2 = []
svxs1  = []; svys1  = []; svzs1 = []
mvxs1  = []; mvys1  = []; mvzs1 = []
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

    mvx = float(cols[33+i]); mvy = float(cols[34+i]); mvz = float(cols[35+i])
    if (nstr == '400'):
        svxs4.append(svx); svys4.append(svy); svzs4.append(svz);
        mvxs4.append(mvx); mvys4.append(mvy); mvzs4.append(mvz);
    elif (nstr == '200'):
        svxs2.append(svx); svys2.append(svy); svzs2.append(svz);
        mvxs2.append(mvx); mvys2.append(mvy); mvzs2.append(mvz);
    elif (nstr == '100'):
        svxs1.append(svx); svys1.append(svy); svzs1.append(svz);
        mvxs1.append(mvx); mvys1.append(mvy); mvzs1.append(mvz);



pltno = 0
plt.interactive(False)

lim=0.15
pltno += 1
plt.figure(pltno)
plt.title('Horizontal Variances, N=400,200,100')
plt.xlabel('Self-Projected MIG Variance (m^2)')
plt.ylabel('Self-Projected Hourglass Variance (m^2)')
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
plt.xlabel('Self-Projected MIG Variance (m^2)')
plt.ylabel('Self-Projected Hourglass Variance (m^2)')
plt.axes().set_xlim([0.0, lim])
plt.axes().set_ylim([0.0, lim])
plt.axes().set_aspect('equal')
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.plot(mvzs, svzs, 'kD', fillstyle='none', label='Z', markersize=5)
plt.plot([0.0,lim], [0,lim], 'k-');
plt.legend(numpoints=1, loc='upper left')
plt.savefig('../papers/fig_self_var_z.png')

lim=0.15
pltno += 1
plt.figure(pltno)
plt.title('Horizontal Variances, N=100,200,400')
plt.xlabel('Self-Projected MIG Variance (m^2)')
plt.ylabel('Self-Projected Hourglass Variance (m^2)')
plt.axes().set_xlim([0.0, lim])
plt.axes().set_ylim([0.0, lim])
plt.axes().set_aspect('equal')
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.plot(mvxs1, svxs1, 'kx', fillstyle='none', label='100')
plt.plot(mvxs2, svxs2, 'k.', fillstyle='none', label='200')
plt.plot(mvxs4, svxs4, 'k+', fillstyle='none', label='400')
#plt.plot(mvys, svys, 'k^', fillstyle='none', label='Y')
plt.plot([0.0,lim], [0.0,lim], 'k-');
plt.legend(numpoints=1, loc='upper left')
plt.plot([0.02,0.02], [0,1], 'k-',
         [0.04,0.04], [0,1], 'k-',
         [0.08,0.08], [0,1], 'k-')
plt.savefig('../papers/fig_self_var_mig_x.png')

lim=0.7
pltno += 1
plt.figure(pltno)
plt.title('Vertical Variances, N=100,200,400')
plt.xlabel('Self-Projected MIG Variance (m^2)')
plt.ylabel('Self-Projected Hourglass Variance (m^2)')
plt.axes().set_xlim([0.0, lim])
plt.axes().set_ylim([0.0, lim])
plt.axes().set_aspect('equal')
#plt.gca().set_xscale('log')
#plt.gca().set_yscale('log')
plt.plot(mvzs1, svzs1, 'kx', fillstyle='none', label='100')
plt.plot(mvzs2, svzs2, 'k.', fillstyle='none', label='200')
plt.plot(mvzs4, svzs4, 'k+', fillstyle='none', label='400')
#plt.plot(mvzs, svzs, 'kD', fillstyle='none', label='Z', markersize=5)
plt.plot([0.0,lim], [0,lim], 'k-');
plt.legend(numpoints=1, loc='upper left')
plt.plot([0.07,0.07], [0,1], 'k-',
         [0.15,0.15], [0,1], 'k-',
         [0.3,0.3], [0,1], 'k-')
plt.savefig('../papers/fig_self_var_mig_z.png')



plt.show()
