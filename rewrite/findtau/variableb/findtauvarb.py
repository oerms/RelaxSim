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

size = 1e-8
# constants from exel table
D = 2.026e-16
b = 1.914e-9
evo_time = 1e+3

fig = plt.figure()
ax = fig.add_subplot(111)

T1tau = []

D = Dconst()    # 
print "D    :",D
a = 2.847e-10   # lattice constant
density = 4e26
mean_dist = avdist(density)

print "den:",density
print "avd:",mean_dist
evo_time = 1e2
walks = 500

#plotting b constants
bfield = 1.22

for tau in relaxsim.log_range(5e-10,5e-8,num=5):
    print "tau:",tau,"s"
    mue = avmup(bfield,tau)
    b = bradius(mue)
    C = Cconst(bfield,tau)
    print "C    :",C,"m^6/s"
    print "b    :",b,"m"
    print "b/a  :",b/2.843e-10
    print "C/b^6:",C/b**6,"1/s"
    ceninst = RelaxCenters(size,density,C,D,b,name ="2vartau"+"{:.2}".format(tau))
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks)
    
    lab = "$\\SI{"+"{: 3.0f}".format(tau/1e-8)+"e-8}{s}$"
    
    result = RelaxResult(experiment=experiment,logscale=True)
    result.plot_data(activewalkers=True,axes=ax,label=lab,logx=True)
    filename = result.write_hdf()
    
    T1tau.append([tau, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1tau =  np.array(T1tau).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
#ax2.errorbar(T1tau[0],1/T1tau[1],yerr=T1tau[2]/T1tau[1]**2,label="$T_1^{-1}$/s$^{-1}$")
ax2.errorbar(T1tau[0],T1tau[1],yerr=T1tau[2],label="$T_1$/s")
ax2.errorbar(T1tau[0],T1tau[3],yerr=T1tau[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("$\\tau$ / s")

leg = ax2.legend(loc=0)
#shift = max([t.get_window_extent().width for t in leg.get_texts()])
#for t in leg.get_texts():
    #t.set_ha('right') # ha is alias for horizontalalignment
    #t.set_position((shift,0))

print "saving figures as pdf"
fig.savefig('magvarbsmalltau.pdf',format='pdf')
fig2.savefig('t1varbsmalltau.pdf',format='pdf')

# we find from t1beta.pdf: tau=1.2e-6 s so that T1 = 9e1 s