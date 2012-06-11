# varying the B field

import matplotlib
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


def Cconst(bfield,tau):
    """return C from magnetic field"""
    const = 6.549e-44       # prefactor before J(omega) of no dimension
    gamma = 2.516e8         # gyromagnetic ratio of flourine
    omega = gamma*bfield
    return const*tau/(1+(omega*tau)**2)

size = 15e-8
C = 1e-54
D = 1e-16
b = 1e-9
evo_time = 1e+2
mean_dist = 1e-8
density = 1/(mean_dist**3*np.pi*4/3.)
walks = 500

fig = plt.figure()
ax = fig.add_subplot(111)

T1dens = []

fig = plt.figure()
ax = fig.add_subplot(111)

tau=1e-7

for bfield in np.arange(1,5,0.2):
    C = Cconst(bfield,tau)
    
    ceninst = RelaxCenters(size,density,C,D,b,name ="varB1")
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks)
    
    fu = "$\\SI{"+"{:.2}".format(bfield)+"}{T}$"
    
    result = RelaxResult(experiment=experiment,logscale=True)
    result.plot_data(axes=ax,label=fu,logx=True)
    filename = result.write_hdf()
    
    T1dens.append([bfield, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1dens =  np.array(T1dens).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.errorbar(T1dens[0],1/T1dens[1],yerr=T1dens[2]/T1dens[1]**2,label="$T_1^{-1}$/s$^{-1}$")
ax2.errorbar(T1dens[0],T1dens[3],yerr=T1dens[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("C / m$^{6}$s$^{-1}$")

ax2.legend()

print "saving figures as pdf"
fig.savefig('magnetizations.pdf',format='pdf')
fig2.savefig('t1beta.pdf',format='pdf')