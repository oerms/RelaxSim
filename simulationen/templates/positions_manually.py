#import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import relaxsim
from relaxsim import  *
reload(relaxsim)
from relaxsim import  *

bfieldsteps = 2
blogvec = log_range(1e-3,1e-1,bfieldsteps)
evologvec = log_range(1e0,1e1,bfieldsteps)

for bf in range(len(blogvec)):
    name = "Positioned"#+"{}".format(i)
    bfield = blogvec[bf]   			# 0,125 < B < 6
    tau = 1e-4				#
    numberwalkers = 3
    evo_time = evologvec[bf]
    size = 2e-8			#  

    #set start  of relative location
    i = 0.1
    #set end of relative location
    end = 0.4
    #number of steps
    steps = 2
    #calculating steplength
    if steps > 1:
        step = (end - i) / (steps-1)
    else:
        end = i
        step = i
        
    C,D,b = getCDb(bfield,tau)
    
    while i <= end:
        vec = [[i,i,i],[1-i,i,i],[i,1-i,i],[i,i,1-i],[1-i,1-i,i],[1-i,i,1-i],[i,1-i,1-i],[1-i,1-i,1-i]]
        centerpositions = np.array(vec)

        ceninst = RelaxCenters(size,C,D,b,density=None,tau=tau,bfield=bfield, distribution="clu",halo_rad=None,name=name,var="{}".format(i),\
        centerpositions=centerpositions,trackvec = None)

        experiment = RelaxExperiment(ceninst, evo_time, dxb_ratio=1/6., method="randomwalks",fake=False,walk_type="continuous")
        experiment.run_experiment(walks=numberwalkers,plotaxes=None,background = 0,threshhold=1e-4)

        result = RelaxResult(experiment=experiment,logscale=True)
        filename = result.write_hdf()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        plotlabel = "Centers: "+"{}".format(ceninst.numberofcentersreal)+\
        "\n"+"density: "+"{:.2}".format(float(ceninst.density))+\
        "\n"+"Walkers: "+"{}".format(numberwalkers)#+" rad$^{-3}$"
        
        result.plot_data(activewalkers=True,ploterror=False,axes=ax,label=plotlabel,logx=True)
        fig.savefig('./'+name+','+'bfield:'+"{:.3}".format(bfield)+','+'i:'+"{:.3}".format(i)+',magn.pdf',format='pdf')

        result.plot_data3d(label=name+','+'bfield:'+"{:.3}".format(bfield)+','+'i:'+"{:.3}".format(i),showplot="no",plottitle='',legendsize=17)
    
        
        i = i + step

