# simple test of relaxsim.py module with classes
import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *


size = 1e-8
C = 1e-54
D = 1e-16
b = 1e-9
mean_dist = 2e-9
density = 1/(mean_dist**3*np.pi*4/3.)
evo_time = 1e2
walks = 200

def doexpran():
    ceninst = RelaxCenters(size,density,C,D,b,name ="ran")
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks)

    result = RelaxResult(experiment=experiment)
    fig2,adf = result.plot_data(label="testlabel")
    fig2.show()
    
import cProfile
cProfile.run("doexpran()","statoutran")
sys.exit()
def doexpdet():
    ceninst = RelaxCenters(size,density,C,D,b,name ="det")

    experiment = RelaxExperiment(ceninst,evo_time,method="deterministic")
    experiment.run_experiment(walks=walks)

    result = RelaxResult(experiment=experiment)
    fig2 = result.plot_data(label="asdf1234!!")

cProfile.run("doexpdet()","statoutdet")
