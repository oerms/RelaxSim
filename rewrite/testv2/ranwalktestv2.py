# simple test of relaxsim.py module with classes

# nice pdf plots in TU style
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
rc('text', usetex=True)
rc('text.latex',preamble='\usepackage[charter]{mathdesign}\n\usepackage{phystex_base}')
rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
rc('mathtext', **{'it':'Charter', 'fontset':'custom'})

# import and reload module
import relaxsimv2
from relaxsimv2 import *
reload(relaxsimv2)
from relaxsimv2 import *

name = "ranwalktestv2_a_"     # name of simulation (ALWAYS rename please!)
bfield  = 0.1
tau     = 1e-4      # stork2010 @ 0.33 T: tau ~ 1e-4, @ 1.22T: tau ~ 2e-6 (with tau ~ B^-3)
density = 1e24
walks   = 10
evo_time = 1e+2

## ## ## ## ## ## ## ## ## ## ## 
import cProfile
def profilefunc():
    size = 5e-8
    C,D,b = getCDb(bfield,tau)
    b /= 2

    fig = plt.figure()
    ax = fig.add_subplot(111)

    mean_dist = (density*np.pi*4/3.)**(-1/3.)
    print "density    ",density
    print "mean_dist  ",mean_dist
    print "radius b   ",b
    print "mean_dist/b",mean_dist/b

    ceninst = RelaxCenters(size,density,C,D,b,name=name+"{:.2}".format(density))
    print ceninst
    experiment = RelaxExperiment(ceninst,evo_time,method="randomwalks")
    experiment.run_experiment(walks=walks)
    result = RelaxResult(experiment=experiment,logscale=True)

    plotlabel = "$\\SI{"+"{:.2}".format(density)+"}{m^{-3}}$"
    plotlabel = ""
    result.plot_data(activewalkers=True,ploterror=False,axes=ax,label=plotlabel,logx=True,color="b")

    filename = result.write_hdf()

    fig.savefig('./'+name+'magn.pdf',format='pdf')
    
cProfile.run("profilefunc()","statoutran")