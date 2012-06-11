# simple test of relaxsim.py module with classes

import sys
sys.path.append('..')

import relaxsim
from relaxsim import *
reload(relaxsim)
from relaxsim import *

mean_dist = 5e-9
density = 1/(mean_dist**3*np.pi*4/3.)

ceninst = RelaxCenters(5e-8,density,1e-54,1e-16,1e-9,name ="test")
#               size, density, C, D, b, name=""):
print(ceninst)
fig,ax = ceninst.plot(fname="centerplot")
fig.show()