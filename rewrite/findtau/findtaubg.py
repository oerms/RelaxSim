# finding the correlation time tau from klemptphd p 43f.
#  at 1.22 T B0 field, 4.3e23 /m^3 centers density, we should get a T1 of 1e1 s for F at RT (fig 6.6, p46)

# RESULT w/ bg rate: tau = 2e-7 s yields good agreement

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

size = 20e-8
# constants from exel table
D = 2.026e-16
b = 1.914e-9
evo_time = 1e+2
#mean_dist = 1e-8
#density = 1/(mean_dist**3*np.pi*4/3.)
bfield=1.22
density = 17e23        # mean_dist approx 8.22e-9m
bg_rate = 2.5e-2        # bgrate from klemptphd p 40 (T1=40 s (F) unbestrahltes LiF)
walks = 200  

fig = plt.figure()
ax = fig.add_subplot(111)

T1tau = []

fig = plt.figure()
ax = fig.add_subplot(111)

for tau in relaxsim._log_range(1e-8,1e-7,num=8):
    print "tau:",tau,"s"
    C = Cconst(bfield,tau)
    print "C  :",C,"m^6/s"
    ceninst = RelaxCenters(size,density,C,D,b,name ="vartaubg"+"{:.2}".format(tau))
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks,background=bg_rate)
    
    lab = "$\\SI{"+"{:.2}".format(tau)+"}{s}$"
    
    result = RelaxResult(experiment=experiment,logscale=True)
    result.plot_data(activewalkers=True,axes=ax,label=lab,logx=True)
    filename = result.write_hdf()
    
    T1tau.append([tau, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1tau =  np.array(T1tau).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.errorbar(T1tau[0],1/T1tau[1],yerr=T1tau[2]/T1tau[1]**2,label="$T_1^{-1}$/s$^{-1}$")
#ax2.errorbar(T1tau[0],T1tau[1],yerr=T1tau[2],label="$T_1$/s")
ax2.errorbar(T1tau[0],T1tau[3],yerr=T1tau[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("$\\tau$ / s")

ax2.legend(loc=0)

print "saving figures as pdf"
fig.savefig('magnetizationsbg.pdf',format='pdf')
fig2.savefig('t1betabg.pdf',format='pdf')

# we find from t1beta.pdf: tau=53e-8 s so that T1 = 1e1 s