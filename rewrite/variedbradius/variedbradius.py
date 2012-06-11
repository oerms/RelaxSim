# finding the correlation time tau from klemptphd p 43f.
#  at 1.22 T B0 field, 4.3e23 /m^3 centers density, we should get a T1 of 1e2 s (fig 6.6, p46)

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

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

size = 10e-8
# constants from stork2010
D = 3.02e-17
a = 2.843e-10
bfield=1.22
tau=1e-4
b=1.5*a
C = Cconst(bfield,tau)
print "C  :",C,"m^6/s"

#b = 1.914e-9
evo_time = 1e+3
#mean_dist = 1e-8
#density = 1/(mean_dist**3*np.pi*4/3.)
density = 4.3e23        # mean_dist approx 8.22e-9m
walks = 150

fig = plt.figure()
ax = fig.add_subplot(111)

T1tau = []

for density in relaxsim._log_range(.4e23,20e23,num=5):
    print "density:", density
    
    ceninst = RelaxCenters(size,density,C,D,b,name ="variedb"+"{:07.2}den".format(b)+"{:07.2}".format(density))
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks)
    
    lab = "$\\SI{"+"{:07.2}".format(density)+"}{s}$"
    
    result = RelaxResult(experiment=experiment,logscale=True)
    result.plot_data(activewalkers=True,axes=ax,label=lab,logx=True)
    filename = result.write_hdf()
    
    T1tau.append([density, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1tau =  np.array(T1tau).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.errorbar(T1tau[0],1/T1tau[1],yerr=T1tau[2]/T1tau[1]**2,label="$T_1^{-1}$/s$^{-1}$")
#ax2.errorbar(T1tau[0],T1tau[1],yerr=T1tau[2],label="$T_1$/s")
#ax2.errorbar(T1tau[0],T1tau[3],yerr=T1tau[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("density / m$^{-3}$")

ax2.legend(loc=0)

print "saving figures as pdf"
fig.savefig('magnb{:07.2}.pdf'.format(b),format='pdf')
fig2.savefig('t1b{:07.2}.pdf'.format(b),format='pdf')

# we find from t1beta.pdf: tau=1.2e-6 s so that T1 = 9e1 s