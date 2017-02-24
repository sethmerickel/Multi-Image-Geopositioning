#! /usr/bin/python3

import fileinput
import sys
import math
import matplotlib.pyplot as plt

def dh(xyz):
    x = xyz[0]
    y = xyz[1]
    dx = float(x) -  455000
    dy = float(y) - 3984000
    dh = math.sqrt(dx*dx + dy*dy)
    return dh


threshn = 100
allns = []; allrefs = []   
smamxs = []; smamys = []; smamzs = []
smapxs = []; smapys = []; smapzs = []
bigmxs = []; bigmys = []; bigmzs = []
bigpxs = []; bigpys = []; bigpzs = []

ns     = []
avgxs  = []; avgys  = []; avgzs = []
pce90s = []; ple90s = []; pce90s_scl = []; ple90s_scl = []
mce90s = []; mle90s = []; mce90s_scl = []; mle90s_scl = []
for line in fileinput.input(sys.argv[1]):
    cols = line.split(',')
    n = int(cols[1])
    if (cols[0] == 'ALL'):
        allns.append(n)
        allrefs.append(float(cols[2]))
        mx=float(cols[3]); my=float(cols[4]); mz=float(cols[5])
        px=float(cols[6]); py=float(cols[7]); pz=float(cols[8])
        if (n > threshn):
            bigmxs.append(mx); bigmys.append(my); bigmzs.append(mz)
            bigpxs.append(px); bigpys.append(py); bigpzs.append(pz)
        else:
            smamxs.append(mx); smamys.append(my); smamzs.append(mz)
            smapxs.append(px); smapys.append(py); smapzs.append(pz)
        continue

    # else cols[0] == 'AGG'
    ns.append(int(cols[1]))
    avgxs.append(float(cols[2]))
    avgys.append(float(cols[3]))
    avgzs.append(float(cols[4]))
    pce90s.append(float(cols[5]))
    ple90s.append(float(cols[6]))
    sn = math.sqrt(n)
    pce90s_scl.append(float(cols[5]) * sn)
    ple90s_scl.append(float(cols[6]) * sn)
    mce90 = float(cols[7])
    mle90 = float(cols[8])
    fpc = math.sqrt((1000-n)/(1000-1))
    if (fpc > 0):
        mce90 /= fpc
        mle90 /= fpc
    else:
        mce90 = float(cols[5])
        mle90 = float(cols[6])
    mce90s.append(mce90)
    mle90s.append(mle90)
    mce90s_scl.append(mce90 * sn)
    mle90s_scl.append(mle90 * sn)




pltno = 0
plt.interactive(False)

pltno += 1
plt.figure(pltno)
plt.gcf().set_size_inches(8,4)
plt.xlabel('Number of images')
plt.ylabel('meters')
plt.axes().set_ylim([-.3, .3])
plt.plot(ns, avgxs, 'k>', fillstyle='none', label='Avg X err')
plt.plot(ns, avgys, 'k^', fillstyle='none', label='Avg Y err')
plt.plot(ns, avgzs, 'kD', fillstyle='none', label='Avg Z err', markersize=5)
plt.legend(numpoints=1)
plt.plot([0,1000], [mx,mx], 'k-',
         [0,1000], [my,my], 'k-',
         [0,1000], [mz,mz], 'k-')
plt.savefig('../papers/fig_avgxyz.png', bbox_inches='tight')

pltno += 1
plt.figure(pltno)
plt.axes().set_aspect(80)
plt.axes().tick_params(left='off')
plt.plot(allns, allrefs, 'k,', label='Each pixel is one refvar')
#plt.legend(numpoints=1)
plt.savefig('../papers/fig_refvar.png', bbox_inches='tight')

pltno += 1
plt.figure(pltno)
plt.gcf().set_size_inches(8,4)
plt.ylabel('meters')
plt.gca().set_yscale('log')
plt.plot(ns, pce90s, 'k+', label='Pred CE90')
plt.plot(ns, ple90s, 'kx', label='Pred LE90')
plt.plot(ns, mce90s, 'k>', label='Meas LE90*FPC', fillstyle='none')
plt.plot(ns, mle90s, 'kD', label='Meas LE90*FPC', fillstyle='none')
plt.legend(numpoints=1)
plt.savefig('../papers/fig_pred_meas.png', bbox_inches='tight')

pltno += 1
plt.figure(pltno)
plt.gcf().set_size_inches(8,4)
plt.axes().set_yticklabels([])
plt.axes().tick_params(left='off')
plt.axes().set_ylim([5, 15])
plt.plot(ns, pce90s_scl, 'k+', label='Pred CE90*sqrtN')
plt.plot(ns, ple90s_scl, 'kx', label='Pred LE90*sqrtN')
plt.plot(ns, mce90s_scl, 'k>', label='Meas CE90*sqrtN*FPC', fillstyle='none')
plt.plot(ns, mle90s_scl, 'kD', label='Meas CE90*sqrtN*FPC', fillstyle='none')
plt.legend(numpoints=1, loc='upper center')
plt.savefig('../papers/fig_pred_meas_sqrtn.png', bbox_inches='tight')

pltno += 1; 
plt.figure(pltno)
plt.gcf().set_dpi(160)
plt.gcf().set_size_inches(8,2.5)
lim=15
p=plt.subplot(131); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(smamxs, smapxs, 'ko', [-20,20], [-20,20], 'k-', markersize=2)
plt.xticks([-15,-5,0,5,15])

p=plt.subplot(132); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(smamys, smapys, 'ko', [-20,20], [-20,20], 'k-', markersize=2)
plt.xticks([-15,-5,0,5,15])

lim=20
p=plt.subplot(133); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(smamzs, smapzs, 'ko', [-20,20], [-20,20], 'k-', markersize=2)
plt.xticks([-20,-10,0,10,20])
plt.savefig('../papers/fig_mig_ply_xyz_sma.png', bbox_inches='tight')

pltno += 1; 
plt.figure(pltno); 
plt.gcf().set_dpi(160)
plt.gcf().set_size_inches(8,2.5)
lim=1.0
p=plt.subplot(131); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(bigmxs, bigpxs, 'k,', [-20,20], [-20,20], 'k-')
#plt.xticks([-15,-5,0,5,15])
p=plt.subplot(132); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(bigmys, bigpys, 'k,', [-20,20], [-20,20], 'k-')
#plt.xticks([-15,-5,0,5,15])

lim=2.0
p=plt.subplot(133); p.set_xlim(-lim,lim); p.set_ylim(-lim,lim)
p.plot(bigmzs, bigpzs, 'k,', [-20,20], [-20,20], 'k-')
plt.gcf().set_size_inches(8,2.5)
plt.xticks([-2,-1,0,1,2])
plt.savefig('../papers/fig_mig_ply_xyz_big.png', bbox_inches='tight')


#plt.show()

