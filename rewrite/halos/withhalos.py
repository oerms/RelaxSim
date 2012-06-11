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

size = 20e-8
#C = 5.79e-55            # after fixing tau at B0 = 1.22 T at room temperature
tau=1.2e-6
D = 2.026e-16
b = 1.914e-9
evo_time = 1e+3
walks = 200

density = 1e22
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

for rad_ratio in [0.1,0.3,.6,0.9,1.2]:
    T1dens = []
    for bfield in np.linspace(0.1,10.,num=20):
        print "halo_rad/size:",rad_ratio
        print "bfield       :",bfield
        C = Cconst(bfield,tau)

        ceninst = RelaxCenters(size,density,C,D,b,distribution='clu',halo_rad=rad_ratio*size,name ="varhaloratio"+"{:.2}".format(rad_ratio)+"bfield{:.2}".format(bfield))
        #ceninst.plot(fname="centers.pdf")
        #sys.exit()
        experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
        experiment.run_experiment(walks=walks)
        
        label = "$\\SI{"+"{:.2}".format(density)+"}{m^{-3}}$"
        
        result = RelaxResult(experiment=experiment,logscale=True)
        result.plot_data(axes=ax,label="",logx=True)
        filename = result.write_hdf(result.name)
        print "T1   :",result.fitT1
        print "beta :",result.fitbeta
        T1dens.append([bfield, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
    T1dens =  np.array(T1dens).transpose()
    ax2.errorbar(T1dens[0],1/T1dens[1],yerr=T1dens[2]/T1dens[1]**2,label="$\\sigma/l ="+"{:.2}".format(rad_ratio)+"$")

ax2.legend()

fig.savefig('magnetizations.pdf',format='pdf')
fig2.savefig('t1beta.pdf',format='pdf')