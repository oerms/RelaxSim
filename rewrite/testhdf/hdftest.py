# simple test of relaxsim.py module with classes

import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *


size = 3e-8
C = 1e-54
D = 1e-16
b = 1e-9
mean_dist = 6e-9
density = 1/(mean_dist**3*np.pi*4/3.)
evo_time = 2e+1
walks = 30

ceninst = RelaxCenters(size,density,C,D,b,name ="test")
fig,axes = ceninst.plot(fname="")
#fig,axes = ceninst.plot(fname="centerplot")

experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
experiment.run_experiment(walks=walks)

result = RelaxResult(experiment=experiment,logscale=True)
fig2,axes2 = result.plot_data(label="",logx=True)
filename = result.write_hdf()

#fig.show()
fig2.show()
#fig2.savefig(filename[-4]+".pdf", format='pdf', facecolor='none', edgecolor='none')
