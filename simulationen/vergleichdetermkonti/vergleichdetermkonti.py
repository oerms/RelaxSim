# simple test of relaxsim.py module with classes

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

import relaxsimv1
from relaxsimv1 import *
reload(relaxsimv1)
from relaxsimv1 import *

# RENAME ME PLEASE!!
name = "compare_"

size = 3e-8
evo_time = 1e+2
walks = 200
tau=1e-4
density = 1e24

fig = plt.figure()
ax = fig.add_subplot(111)

T1densdet = []
T1densrw = []

for bfield in np.array([1e-3,1e-2,1e-1,1e0,1e1]):
    # constants
    #bfield = 1.2
    print 'field is now:', bfield, " T"
    
    #tau = relaxsimv1.gettau(density,bfield)
    C,D,b = relaxsimv1.getCDb(bfield,tau)

    ceninst = RelaxCenters(size,density,C,D,b,name=name+"{:.2}".format(density))
    print ceninst
    #fig,axes = ceninst.plot(fname="centerplot")

    detexp = RelaxExperiment(ceninst,evo_time,method="deterministic")
    detexp.run_experiment()
    detresult = RelaxResult(experiment=detexp,logscale=True)
    detresult.plot_data(axes=ax,label="determ "+"$B_0 = \SI{"+"{:.2e}".format(bfield)+"}{T}$",logx=True)
    filename = detresult.write_hdf()

    rwexp = RelaxExperiment(ceninst,evo_time,method="randomwalks",walk_type='continuous')
    rwexp.run_experiment(walks=walks)
    rwresult = RelaxResult(experiment=rwexp,logscale=True)
    rwresult.plot_data(axes=ax,ploterror=True,label="random "+"$B_0 = \SI{"+"{:.2e}".format(bfield)+"}{T}$",logx=True)
    filename = rwresult.write_hdf()

    T1densdet.append([detresult.density, detresult.fitT1, detresult.fitT1err, detresult.fitbeta, detresult.fitbetaerr])
    T1densrw.append([ rwresult.density,  rwresult.fitT1,  rwresult.fitT1err,  rwresult.fitbeta,  rwresult.fitbetaerr])


ax.legend(loc=3,prop=dict(size=12))
fig.savefig(name+"magn.pdf",format="pdf")

T1densdet =  np.array(T1densdet).transpose()
T1densrw =  np.array(T1densrw).transpose()
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(T1densdet[0],1/T1densdet[1],'rs',label="deterministisch")
ax2.plot(T1densrw[0],1/T1densrw[1],'ko',label="Random Walks")
#ax2.errorbar(T1dens[0],T1dens[3],yerr=T1dens[4],label="$\\beta$")
ax2.set_title("\\Large Walker bleiben stehen.")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel("$B_0$ / T")
ax2.set_ylabel("$T_1^{-1}$ / s$^{-1}$")
ax2.set_xlim(8e-2,2e1)
ax2.legend(prop=dict(size=15))
fig2.savefig(name+'t1.pdf',format='pdf')

np.savetxt("det,T1,T1err,beta,betaerr",T1densdet.transpose())
np.savetxt("det,T1,T1err,beta,betaerr",T1densrw.transpose())

