# simple test of relaxsim.py module with classes

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *


size = 6e-8
C = 1e-54
D = 1e-16
b = 1e-9
evo_time = 2e+2
walks = 500

fig = plt.figure()
ax = fig.add_subplot(111)

T1dens = []

for mean_dist in [1e-9,5e-9,11e-9,16e-9]:
#mean_dist = 5e-9
    density = 1/(mean_dist**3*np.pi*4/3.)

    ceninst = RelaxCenters(size,density,C,D,b,name ="compare2")
    #fig,axes = ceninst.plot(fname="centerplot")

    detexp = RelaxExperiment(ceninst,evo_time,method="deterministic")
    detexp.run_experiment(walks=walks)
    detresult = RelaxResult(experiment=detexp,logscale=True)
    detresult.plot_data(axes=ax,label="determ",logx=True)
    filename = detresult.write_hdf()

    rwexp = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    rwexp.run_experiment(walks=walks)
    rwresult = RelaxResult(experiment=rwexp,logscale=True)
    rwresult.plot_data(axes=ax,label="random",logx=True)
    filename = rwresult.write_hdf()
    
    T1dens.append([result.density, result.fitT1,result.fitT1err, result.fitbeta, result.fitbetaerr])


fig.savefig("comparemagn.pdf",format="pdf")
T1dens =  np.array(T1dens).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.errorbar(T1dens[0],1/T1dens[1],yerr=T1dens[2]/T1dens[1]**2,label="$T_1^{-1}$/s$^{-1}$")
ax2.errorbar(T1dens[0],T1dens[3],yerr=T1dens[4],label="$\\beta$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("density / m$^{-3}$")
ax2.legend()
fig2.savefig('comparet1beta.pdf',format='pdf')