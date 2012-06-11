# simple test of relaxsim.py module with classes

import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *


size = 5e-8
C = 1e-54
D = 1e-16
b = 1e-9
evo_time = 1e+1
walks = 1

T1dens = []

fig = plt.figure()
ax = fig.add_subplot(111)

#for mean_dist in [2e-9,3e-9,4e-9,5e-9,6e-9,7e-9,8e-9,9e-9,10e-9]:
mean_dist = 10e-9
density = 1/(mean_dist**3*np.pi*4/3.)

ceninst = RelaxCenters(size,density,C,D,b,name ="varden")
fig,axes = ceninst.plot(fname="")
#fig,axes = ceninst.plot(fname="centerplot")

experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
def prof():
    experiment.run_experiment(walks=walks,plotaxes=axes)

import cProfile
cProfile.run('prof()','profout.cpr')

result = RelaxResult(experiment=experiment,logscale=True)
fig2, ax2 = result.plot_data(label="",logx=True)
filename = result.write_hdf()

fig.show()
fig2.show()