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
    name = "Clustered"#+"{}".format(i)
    bfield = blogvec[bf]   			# 0,125 < B < 6
    tau = 1e-4				#
    numberwalkers = 5
    evo_time = evologvec[bf]
    size = 5e-8			#  
    density = 1e24          # NumberofCenters / size^3 , 10^22 < density < 10^26 

    #set minimum halo_radius
    halo_radMin = 1e-9
    #set maximum halo_radius
    halo_radMax = 5e-9
    #number of steps (number of radii)
    radiussteps = 3
    #calculating steplangth
    if radiussteps > 1:
        endradius = halo_radMax
        step = (halo_radMax - halo_radMin) / (radiussteps-1)
    else:
        endradius = halo_radMin
        step = (halo_radMax - halo_radMin)
        
    #initialize start radius
    rad = halo_radMin

    C,D,b = getCDb(bfield,tau)

    while rad <= endradius:

        ceninst = RelaxCenters(size,C,D,b,density, tau=tau, bfield=bfield,distribution="clu",halo_rad=rad,name=name,var="{}".format(rad),\
        centerpositions=None,trackvec = [[0.1,0.3],[0.7,0.1]])

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
        fig.savefig('./'+name+','+'bfield:'+"{}".format(bfield)+','+'radius:'+"{}".format(rad)+',magn.pdf',format='pdf')

        result.plot_data3d(label=name+','+'bfield:'+"{}".format(bfield)+','+'radius:'+"{}".format(rad),showplot="no",plottitle='',legendsize=17) 
        
        
        rad = rad + step