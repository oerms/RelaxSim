# simple test of relaxsim.py module with classes

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

size = 1e-8
#C = 5.79e-55            # after fixing tau at B0 = 1.22 T at room temperature
tau=.1e-6
bg_rate = 1/40.
D = 2.026e-16
b = 1.914e-9
evo_time = 1e+2         # max with 4 gb ram is 6e+3 sec evol time!!
walks = 10

density = 5e23
mean_dist = (density*np.pi*4/3.)**(-1/3.)
print "density",density
print "mean_dist",mean_dist
print "mean_dist/b",mean_dist/b

T1dens = []
fig = plt.figure()  # for magnetizations
ax = fig.add_subplot(111)

fig2 = plt.figure()  # for rates
ax2 = fig2.add_subplot(111)
#ax2.errorbar(T1dens[0],T1dens[3],yerr=T1dens[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("$B_0$ / T")
ax2.set_ylabel("$T_1^{-1}$/s$^{-1}$")

#for rad_ratio in np.array([0.1,0.6,1.2]):
rad_ratio = 0.3
T1densran = []
T1densdet = []

#for bfield in np.array( [0.6, 1.22, 3.8] ):
for bfield in np.array( [1.22] ):
    C = Cconst(bfield,tau)
    C=1e-57
    C=0
    print "halo_rad/size:",rad_ratio
    print "bfield       :",bfield
    print "C            :",C

    ceninst = RelaxCenters(size,density,C,D,b,distribution='hom',halo_rad=rad_ratio*size,name ="bgtest1")
    #ceninst.plot(fname="centers{}.pdf".format(rad_ratio))
    #sys.exit()
    label = "$\\SI{"+"{:.2}".format(bfield)+"}{T}$ "

    ranexperiment = RelaxExperiment(ceninst,evo_time,method="det")
    ranexperiment.run_experiment(walks=walks,background=bg_rate)
    ranresult = RelaxResult(experiment=ranexperiment,logscale=True)
    ranresult.plot_data(activewalkers=False,axes=ax,label=label+"det",logx=True)
    filename = ranresult.write_hdf(ranresult.name)
    T1densran.append([bfield, ranresult.fitT1, ranresult.fitT1err, ranresult.fitbeta, ranresult.fitbetaerr])

T1densran =  np.array(T1densran).transpose()
ax2.errorbar(T1densran[0],1/T1densran[1],yerr=T1densran[2]/T1densran[1]**2,label="$\\sigma/l ="+"{:.2}".format(rad_ratio)+"$ ran")

ax2.legend()

fig.savefig('bgmagdetfast.pdf',format='pdf')
fig2.savefig('bgdett1.pdf',format='pdf')