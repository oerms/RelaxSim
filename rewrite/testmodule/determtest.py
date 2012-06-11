# simple test of relaxsim.py module with classes
import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *

size = 1e-8
C = 1e-54
D = 1e-16
b = 1e-9
mean_dist = 3e-9
density = 1/(mean_dist**3*np.pi*4/3.)
evo_time = 1e2

print density


ceninst = RelaxCenters(size,density,C,D,b,name ="test")

experiment = RelaxExperiment(ceninst,evo_time,method="deterministic")

experiment.run_experiment()
result = RelaxResult(experiment=experiment,logscale=True)

experiment.run_experiment(withM=True,redo=True)
result2 = RelaxResult(experiment=experiment,logscale=True)

fig2,ax = result.plot_data(label="with")
fig,ax = result2.plot_data(axes=ax,label="without")
fig2.show()