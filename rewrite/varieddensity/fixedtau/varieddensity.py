# varying the density after finding tau from ../findtau/findtaubg.py
# hoping to get similar curve to klemptphd fig 6.6

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})
import matplotlib.pyplot as plt

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *

size = 15e-8*10
C = 1e-54
D = 2.026e-16
b = 1.914e-9
evo_time = 1e+2
walks = 50

bfield=1.22
tau=53e-8*250

print "tau:",tau

def Cconst(bfield,tau):
    """return C from magnetic field"""
    const = 6.549e-44       # prefactor before J(omega) of no dimension
    gamma = 2.516e8         # gyromagnetic ratio of flourine
    omega = gamma*bfield
    return const*tau/(1+(omega*tau)**2)
C = Cconst(bfield,tau)  # C = 5.79e-55
bg_rate = 2.5e-2        # bgrate from klemptphd p 40 (T1=40 s (F) unbestrahltes LiF)

T1dens = []
fig = plt.figure()
ax = fig.add_subplot(111)

for density in relaxsim._log_range(3e22,2e24,5)/100:
#for mean_dist in [4e-9]:
    mean_dist = (density*np.pi*4/3.)**(-1/3.)
    print "density",density
    print "mean_dist",mean_dist
    print "mean_dist/b",mean_dist/b

    ceninst = RelaxCenters(size,density,C,D,b,name ="1varden"+"{:07.1}".format(density))
    #fig,axes = ceninst.plot(fname="")
    #fig,axes = ceninst.plot(fname="centerplot")

    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks,background=bg_rate)
    
    fu = "$\\SI{"+"{:.2}".format(density)+"}{m^{-3}}$"
    print fu
    
    result = RelaxResult(experiment=experiment,logscale=True)
    result.plot_data(activewalkers=True,axes=ax,label=fu,logx=True)
    filename = result.write_hdf(result.name)
        
    T1dens.append([result.density, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1dens =  np.array(T1dens).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.errorbar(T1dens[0],1/T1dens[1],yerr=T1dens[2]/T1dens[1]**2,label="$T_1^{-1}$/s$^{-1}$")
#ax2.errorbar(T1dens[0],T1dens[3],yerr=T1dens[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("density / m$^{-3}$")

ax2.legend()

fig.savefig('./magnetizationsbg.pdf',format='pdf')
fig2.savefig('./t1betabg.pdf',format='pdf')