# simple test of relaxsim.py module with classes

import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *



size = .5e-8
C = 1e-54
D = 1e-16
b = 1e-9
evo_time = 3e+1
walks = 500

fig = plt.figure()
ax = fig.add_subplot(111)


mean_dist = 1.85e-9
density = 1/(mean_dist**3*np.pi*4/3.)

ceninst = RelaxCenters(size,density,C,D,b,name ="compare")

ax = fig.ax = fig.gca(projection='3d')

ax.scatter(*(ceninst.center_positions.transpose()), c='b')

cen = ceninst.center_positions.flatten()
dx = ceninst.b/2.           # coarse grid for plotting
dt = dx**2./(6*ceninst.D)
#print "dt",dt,"    dx",dx
xr = np.arange(0,+ceninst.size+dx,dx)       # regular grid for sampling
yr = np.arange(0,+ceninst.size+dx,dx)
xlen = len(xr)
#print "number of grid points for plotting centers:",xlen**2
# make grid (xr and yr now of 2d with squared as many elements, giving x,y coordinates respectively)
xr,yr = np.meshgrid(xr,yr)
xr = xr.flatten()      # flatten for loop
yr = yr.flatten()
rr = np.empty(xr.shape,dtype=xr.dtype)
# begin rate scanning
#progress = Progress(len(xr))
array_for_rate = np.array([ceninst.b,ceninst.size]) # saving time making the array before loop
middle_point = 0.5*ceninst.size*0                   # same here
for point in range(len(xr)):
    rr[point] = accum_dist_C(cen,array_for_rate,np.array([xr[point],yr[point],middle_point]))
    if rr[point] == -1: #accum_dist_C returns -1, if inside radius: set to value on border
        rr[point] = ceninst.b**-6.
    #progress.increment()
    #progress.print_status_line()
rr = 1.-(dt*ceninst.C*rr)        # apply relaxation rate
rr = rr.clip(0,1)             # clamp the array (not really necessary...)
rrmin = min(rr)
# reshape to 2d for plotting
rr = rr.reshape((xlen,xlen))
yr = yr.reshape((xlen,xlen))
xr = xr.reshape((xlen,xlen))
# prepare plot

levels = np.linspace(rrmin,1,num=7) # levels for the contour plot

cen_plot = ax.contourf(xr,yr,rr, levels=levels, alpha=0.6,offset=middle_point)

#print "dt",dt,"    dx",dx
xr = np.arange(0,+ceninst.size+dx,dx)       # regular grid for sampling
yr = np.arange(0,+ceninst.size+dx,dx)
xlen = len(xr)
#print "number of grid points for plotting centers:",xlen**2
# make grid (xr and yr now of 2d with squared as many elements, giving x,y coordinates respectively)
xr,yr = np.meshgrid(xr,yr)
xr = xr.flatten()      # flatten for loop
yr = yr.flatten()
rr = np.empty(xr.shape,dtype=xr.dtype)
# begin rate scanning
#progress = Progress(len(xr))
array_for_rate = np.array([ceninst.b,ceninst.size]) # saving time making the array before loop
middle_point = 0.7*ceninst.size                   # same here
for point in range(len(xr)):
    rr[point] = accum_dist_C(cen,array_for_rate,np.array([xr[point],yr[point],middle_point]))
    if rr[point] == -1: #accum_dist_C returns -1, if inside radius: set to value on border
        rr[point] = ceninst.b**-6.
    #progress.increment()
    #progress.print_status_line()
rr = 1.-(dt*ceninst.C*rr)        # apply relaxation rate
rr = rr.clip(0,1)             # clamp the array (not really necessary...)
rrmin = min(rr)
# reshape to 2d for plotting
rr = rr.reshape((xlen,xlen))
yr = yr.reshape((xlen,xlen))
xr = xr.reshape((xlen,xlen))
# prepare plot

levels = np.linspace(rrmin,1,num=7) # levels for the contour plot

cen_plot = ax.contourf(xr,yr,rr, levels=levels, alpha=0.6,offset=middle_point)

ax.set_xlabel('$x\SI{}{/\meter}$')
ax.set_xlim(0, size)
ax.set_ylabel('$y\SI{}{/\meter}$')
ax.set_ylim(0, size)
ax.set_zlabel('$z\SI{}{/\meter}$')
ax.set_zlim(0, size)
ax.legend()

fig.savefig("centers.pdf",format="pdf")