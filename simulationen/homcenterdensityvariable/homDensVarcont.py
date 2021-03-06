#import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import relaxsim
from relaxsim import  *
reload(relaxsim)
from relaxsim import  *

name = "cont_homogeneous"#+"{}".format(i)

denssteps = 8
densvec = log_range(1e23,1e26,denssteps)

numberwalkers = 200
bfield = 1e-2  			# 0,125 < B < 6
tau = 1e-4				#
evo_time = 1e2
size = 5e-8

T1dens = []

for dens in range(len(densvec)):
    density = densvec[dens]         # NumberofCenters / size^3 , 10^22 < density < 10^26 

    C,D,b = getCDb(bfield,tau)
    
    thisname = name+','+'dens'+"{:.4}".format(density)+',magn.pdf'

    ceninst = RelaxCenters(size,C,D,b,density, tau=tau, bfield=bfield,distribution="hom",halo_rad=None,name=thisname,\
    centerpositions=None,trackvec = None)

    experiment = RelaxExperiment(ceninst, evo_time, dxb_ratio=1/6., method="randomwalks",fake=False,walk_type="continuous")
    experiment.run_experiment(walks=numberwalkers,plotaxes=None,background = 0,threshhold=1e-4)

    result = RelaxResult(experiment=experiment,logscale=True)
    filename = result.write_hdf()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plotlabel = "Centers: "+"{}".format(ceninst.numberofcentersreal)+\
    "\n"+"density: "+"{:.2}".format(float(ceninst.density))+\
    "\n"+"Walkers: "+"{}".format(numberwalkers)
    
    result.plot_data(activewalkers=True,ploterror=True,axes=ax,label=plotlabel,logx=True)
    fig.savefig('./'+thisname)

    result.plot_data3d(label=name+','+'dens:'+"{:.4}".format(density),showplot="no",plottitle='',legendsize=17) 

    T1dens.append([result.density, result.fitT1, result.fitT1err, result.fitbeta, result.fitbetaerr])
    
T1dens = np.array(T1dens)

np.savetxt("contdens,T1,T1err,beta,betaerr.txt",T1dens)

