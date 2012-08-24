# auswertung der neuen vergleichsrechnung stepstaycont

# nice pdf plots in TU style
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

import matplotlib.pyplot as plt
import numpy as np

# read files
t1arraycont = np.genfromtxt("contdens,T1,T1err,beta,betaerr.txt").transpose()
t1arraystay = np.genfromtxt("staydens,T1,T1err,beta,betaerr.txt").transpose()
t1arraystep = np.genfromtxt("stepdens,T1,T1err,beta,betaerr.txt").transpose()
print t1arraystay

# figure for plot
fig= plt.figure(frameon=False)
ax = fig.add_subplot(111)

# plot
ax.plot(t1arraycont[0],t1arraycont[1],"g-")
#ax.errorbar(t1arraycont[0],1/t1arraycont[1],t1arraycont[2]/t1arraycont[1]**2,fmt="go",label="continuous")
ax.errorbar(t1arraycont[0],t1arraycont[1],t1arraycont[2],fmt="go",label="continuous")


ax.plot(t1arraystay[0],t1arraystay[1],"r-")
#ax.errorbar(t1arraystay[0],1/t1arraystay[1],t1arraystay[2]/t1arraystay[1]**2,fmt="ro",label="stay")
ax.errorbar(t1arraystay[0],t1arraystay[1],t1arraystay[2],fmt="ro",label="stay")

ax.plot(t1arraystep[0],t1arraystep[1],"b-")
#ax.errorbar(t1arraystep[0],1/t1arraystep[1],t1arraystep[2]/t1arraystep[1]**2,fmt="bo",label="step")
ax.errorbar(t1arraystep[0],t1arraystep[1],t1arraystep[2],fmt="bo",label="step")

#plot power law
plx = np.array([2e23,8e25])
ply = 1e25*plx**-1
ax.plot(plx,ply,"k-",label='$m = -1$')

# plot options
ax.legend(loc=3,prop=dict(size=17))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel("$T\ped{1n}$/s")
ax.set_xlabel("$n\\ped F$/m$^{-3}$")
ax.text(1e25,5e0,"$T\ped{1n}\\propto n\\ped F^{-1}$",color="k")
# distances
ax.text(1.1e23,4e0,"$d/b = 16$",color="k",size='smaller')
ax.text(7e23,4e-1,"$d/b = 7.5$",color="k",size='smaller')
ax.text(7e24,4e-2,"$d/b = 3.5$",color="k",size='smaller')
ax.text(4e25,3e-3,"$d/b = 1.6$",color="k",size='smaller')

fig.savefig("vergleichT1zeiten.pdf")
import sys
sys.exit()

# bar and label
barrange = np.array([5e22,2e24])
barval = 1.5e-3*(barrange/4e22)**1
ax.plot(barrange,barval,"r",label="Potenzgesetz")

# save plot

fig.savefig("./t1auswertung6vardens.pdf")

