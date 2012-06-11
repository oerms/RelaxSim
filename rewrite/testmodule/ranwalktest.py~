# simple test of relaxsim.py module with classes
import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *

mean_dist = 2e-9
density = 1/(mean_dist**3*np.pi*4/3.)

ceninst = RelaxCenters(2e-8,density,1e-54,1e-16,1e-9,name ="test")
#               size, density, C, D, b, name=""):

fig = ceninst.plot(fname="centerplot")

import cProfile
experiment = RelaxExperiment(ceninst,1e-1,method="randomwalks")
#                      centers, evo time
experiment.run_experiment(walks=50)

result = RelaxResult(experiment=experiment)
fig2 = result.plot_data(label="asdf1234!!")
#fig.show()
fig2.show()
