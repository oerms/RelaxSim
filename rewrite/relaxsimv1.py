"""relaxsim.py is a python module for simulating nuclear relaxation by paramagnetic centers

Workflow:
create a RelaxCenters instance
give it to a RelaxExperiment instance
run the experiment
create a RelaxResult instance and read the data from the experiment
plot and/or write results to HDF file
"""

__author__="Simon Quittek"
__date__="2012/03/18"
__version__="0.1"

import time as t
import sys
import os
import errno
import math
import numpy as np
import numpy.random as rnd; rnd.seed()
import scipy.odr as odr # orthogonal data regression
import matplotlib
#matplotlib.use('Agg')
#from matplotlib import rc
#rc('text', usetex=True)
#rc('text.latex',preamble='\usepackage[charter]{mathdesign}')
#rc('font',**{'family':'serif','serif':['Charter'], 'size':16})
#rc('mathtext', **{'it':'Charter', 'fontset':'custom'})
import matplotlib.pyplot as plt
from matplotlib import cm
# inline C code
from instant import inline_with_numpy
# a progress bar from the internets
from progress import *
# save/load objects to file
import pickle


## ## ## ## HELPER FUNCTIONS BEGIN ## ## ## ##
def gaussn(x,sigma):
    """gauss function normalized to maximum = one"""
    return np.exp(-(x/float(sigma))**2/2.)
    
def lorentzn(x,sigma):
    """lorentz curve normalized to maximum = one"""
    return sigma**2/(float(sigma)**2+x**2)
    
def overlapg(x,sigma,sigma2=None):
    """overlap of two gaussians with sigma and sigma2/twice sigma"""
    if sigma2 == None:
        sigmanew = np.sqrt(2)*sigma
    else:
        sigmanew = np.sqrt(sigma**2+sigma2**2)
    return gaussn(x,sigmanew)

def overlapl(x,sigma,sigma2=None):
    """overlap of two lorentzians with sigma and sigma2/twice sigma"""
    if sigma2 == None:
        sigmanew = 2*sigma
    else:
        sigmanew = sigma+sigma2
    return gaussn(x,sigmanew)

def dens(mean_dist):
    """calculate density from average distance"""
    return 1/(mean_dist**3*np.pi*4/3.)

def avdist(density):
    """calculate average distance from density"""
    return 1/(density**3*np.pi*4/3.)**(1/3.)
    
def coth(x):
    """hyperbolic cotangent"""
    return 1/np.tanh(x)
    
def csch(x):
    """hyperbolic cosecant"""
    return 1/np.sinh(x)
    
def BS(x,S):
    """The Brillouin function"""
    return (2*S+1)/(2*S) * coth((2*S+1)/(2*S)*x) - 1/(2*S) * coth(x/(2*S))

def BSdash(x,S):
    """first derivative of the Brillouin function"""
    return -((2*S+1)/(2*S))**2 * csch((2*S+1)/(2*S)*x)**2 + 1/(4*S*S) * csch(x/(2*S))**2

def avmup(bfield,tau,temperature=300,nucleus="F",T2n=''):
    """calculate average mu of impurity according to rorschach"""
    mue = 9.285e-24         # electronic magnetic moment in J/T
    kB  = 1.381e-23         # Boltzmann constant in J/K
    xforBS  = mue * bfield/(kB*temperature)
    if nucleus == "F":
        gamma   = 2.516e8       # gyromagnetic ratio of flourine
        S       = 1/2.          # flourine spin quantum number
        munuc   = 1.327e-26     # magnetic moment of fluorine
        T2nuc   = 53.5e-6       # T2 time of fluorine in LiF (stork2010)
    elif nucleus == "Li":
        print "TODO"
        # TODO
    else:
        raise ValueError("No other element than fluorine implemented.")
    if T2n != '':
        T2nuc = T2n
    parenth = BS(xforBS,S)**2 + BSdash(xforBS,S)*2/np.pi * np.arctan(2*np.pi*tau/T2nuc)
    return mue*np.sqrt(parenth)

def bradius(avmup,nucleus="F",lattconst=2.847e-10):
    """bradius from average mup"""
    if nucleus == "F":
        munuc   = 1.327e-26     # magnetic moment of fluorine
    elif nucleus == "Li":
        print "TODO"
        # TODO
    else:
        raise ValueError("No other element than fluorine implemented.")
    return pow(3*avmup/munuc,1/4.)*lattconst

def Cconst(bfield,tau,nucleus="F"):
    """relaxation constant C from magnetic field and tau"""
    if nucleus == "F":
        const = 6.549e-44       # prefactor before J(omega) of no dimension
        gamma = 2.516e8         # gyromagnetic ratio of flourine
    elif nucleus == "Li":
        print "TODO"
        # TODO
    else:
        raise ValueError("No other element than fluorine implemented.")
    omega = gamma*bfield
    return const*tau/(1+(omega*tau)**2)

def Dconst(a=2.847e-10,nucleus="F",T2=''):
    """diffusion constant D"""
    if nucleus == "F":
        T2nuc = 53.5e-6
    elif nucleus == "Li":
        print "TODO"
        # TODO
    else:
        raise ValueError("No other element thatn fluorine implemented.")
    if T2 != '':
        T2nuc = T2
    return a**2/(50*T2nuc)

def getCDb(bfield,tau,nucleus="F",**kwargs):
    """return (C,D,b) from input
    **kwargs are given to avmup"""
    return (Cconst(bfield,tau,nucleus=nucleus),Dconst(),bradius(avmup(bfield,tau,nucleus=nucleus,**kwargs)))

def log_range(start, stop, num=50):
    """return numpy array with equidistant values on a logarithmic scale
    
    input:
        start   starting value
        stop    stopping value (included)
        num     number of values (optional default=50)"""
    if (start<=0 or stop<=0 or num<1):
        raise ValueError("start, stop must be positive and num must be >=1")
    return np.logspace(np.log10(start),np.log10(stop), num=num)

def _upd_avg_var(mk,vk,k,newx,samcov=True):
    """Update the sample-average mk and variance vk with new value newx.
    The sample consists of k   values before updating,
    the sample consists of k+1 values after  updating.

    Input:
    mk, vk, k       sample average, variance, number of samples
    newx            new sample datum

    samcov == True  Update using the corrected sample-variance.
                    Otherwise using the standard deviation.
                    
    Return:
    mk, vk, k       same as input, but now k = k+1"""
    #check for valid lengths
    
    #if ( bool(type(mk) != int) != bool(type(mk) != float) != bool(type(mk) != complex)): # this is an xor!
        #if ( len(mk) != len(vk) ):
            #raise RWError(3,"when updating, mk/vk lengths do not match!")
        #if ( len(mk) != len(newx) ):
            #raise RWError(3,"when updating, mk/newx lengths do not match!")
            
    # convert to np.arrays for addition and exponentiation, k to float
    mk   = np.array(mk)
    vk   = np.array(vk)
    newx = np.array(newx)
    k    = float(k)
    # update variance
    if samcov == True:
        vk = (k-1)*vk/k + (newx - mk)**2/(k+1)
        #print 'max difference',max((newx-mk)**2)
    else :
        vk = (k*vk + k*(newx - mk)**2/(k+1))/(k+1)
    # update average
    mk = (k*mk+newx)/(k+1)
    # return new mean, variance, and sample size for next update
    return mk,vk,k+1


"""Compile the C code for random walks and bind to variable. Do same with python version.

Return the accumulated distances^-6 from position to all centers. Return -1 if inside b_rad

Input:
centers flattened array of positions
consts  array of b_rad (when inside return -1), size (for periodic boundaries)
pos     current position

works with real positions as well as with coordinates j,l,m"""
accum_code = """
double accum_dist_C (int lencen, double* centers, int constlen, double* consts, int dim, double* pos ){
    double rel = 0;
    double dx; double dy; double dz;
    double distsq;
    for (int i=0; i<lencen/3; i++) {
        //accumulate rate by adding 1/r^6, keeping periodic boundary conditions in mind
        //have to use local header with cmath include and fabs!
        dx = (fabs(pos[0]-centers[3*i  ])) < fabs(consts[1]-fabs(pos[0]-centers[3*i  ])) ?\
             (fabs(pos[0]-centers[3*i  ])) : fabs(consts[1]-fabs(pos[0]-centers[3*i  ]));
        dy = (fabs(pos[1]-centers[3*i+1])) < fabs(consts[1]-fabs(pos[1]-centers[3*i+1])) ?\
             (fabs(pos[1]-centers[3*i+1])) : fabs(consts[1]-fabs(pos[1]-centers[3*i+1]));
        dz = (fabs(pos[2]-centers[3*i+2])) < fabs(consts[1]-fabs(pos[2]-centers[3*i+2])) ?\
             (fabs(pos[2]-centers[3*i+2])) : fabs(consts[1]-fabs(pos[2]-centers[3*i+2]));
        distsq = dx*dx + dy*dy + dz*dz;
        if (distsq < consts[0]*consts[0]) {
            return -1;   // turn stepping off
            }
        rel += pow(distsq,-3);
        }
    return rel;
}
"""
# this compiles the code into a built-in function when importing the module, need "local_headers.h" file
fileh = open("local_headers.h","w")
fileh.write("#include <cmath>")
fileh.close()
accum_dist_C = inline_with_numpy(accum_code, arrays = [['lencen', 'centers'],['constlen','consts'],['dim','pos']],local_headers=["local_headers.h"])
import os
os.remove("local_headers.h")

def accum_dist_P(centers, const, pos,exit = True):
    """accumulate relaxation (sum r**-6)
    rel_cen         flattened array of x,y,z positions of relax centers
    brad            array (!) of brad (only one element)
    pos             array of x,y,z positions of walker
    exit = True     if True, returns -1 when inside brad
                    if False, returns sum r**-6

    The C version of the same function (exit=True fixed!) is called accum_dist_C and accepts the same arguments."""
    rel = 0
    for i in range(len(rel_cen)/3):
        distsq = \
    (min( abs(pos[0]-centers[3*i  ]), const[1]-abs(pos[0]-centers[3*i  ])))**2\
   +(min( abs(pos[1]-centers[3*i+1]), const[1]-abs(pos[1]-centers[3*i+1])))**2\
   +(min( abs(pos[2]-centers[3*i+2]), const[1]-abs(pos[2]-centers[3*i+2])))**2
        if distsq < brad[0]**2 and exit :
            return -1
        rel += distsq**-3
    return rel

"""
fold_back_C folds back any coordinates that are outside [0,size] in place

call:   fold_back_C(np.array(coordinates),np.array([size]))

"""
fold_code="""
void fold_back_C (int lencoo, double* coords, int constlen, double* consts){
    
    for (int i=0; i<lencoo; i++) {
        if (coords[i] > consts[0]) { coords[i] -= consts[0]*floor(coords[i]/consts[0]); }
        if (coords[i] <        0 ) { coords[i] += consts[0]*ceil(-coords[i]/consts[0]); }
        }
    return;
}
"""
fold_back_C = inline_with_numpy(fold_code, arrays = [['lencoo', 'coords'],['constlen','consts']])


fold_code="""
}
"""
def fold_back_P(coords, consts):
    """
    fold back any coordinates that are outside [0,size] conforming to periodic boundaries

    call:   fold_back_C(np.array(coordinates),np.array([size]))
    
    The C version of the same function is called fold_back_C and accepts the same arguments."""
    for i in range(lencoo):
        if (coords[i] > consts[0]):
            coords[i] -= consts[0]*floor(coords[i]/consts[0])
        if (coords[i] <        0 ):
            coords[i] += consts[0]*ceil(-coords[i]/consts[0])
    return
    
## ## ## ## HELPER FUNCTIONS END ## ## ## ##


## ## ## ## CLASS DEFINITIONS BEGIN ## ## ## ##

class RelaxError(Exception):
    """A simple Error class for the relaxation calculations"""
    def __init__(self,err_code,err):
        self.err = err
        self.err_code = err_code
        
    def __str__(self):
        return self.err


class RelaxCenters():
    """A class for the "sample": the positions of the paramagnetic centers, and other constants."""
    def __init__(self, size, density, C, D, b, distribution="homogeneous",halo_rad=5e-9,name=""):
        """constructor

        input:
            density     centers per m^2         float
            size        side length in meter    float
            C           m^6/s                   float
            D           m^2/s                   float
            b           m                       float
        optional input:
            distribution    "homogenous", "clustered", etc.
                
        output:
            self.centers    array of real centers' positions"""
        t_inicen = t.time()
        self.name       = name+"cen"
        self.density    = density
        self.mean_dist = (self.density*math.pi*4/3.)**(-1/3.)
        self.size       = size
        self.C          = C
        self.D          = D
        self.b          = b
        self.distribution= distribution
        # numbers of centers in center box
        self.center_number = int(density*size**3)
        print "Number of centers:",self.center_number
                
        if self.distribution == "homogeneous" or self.distribution == 'hom':
            self.name = self.name + "hom" + "cen"
            self.center_positions = self.size*rnd.rand(self.center_number,3)
        elif self.distribution == 'clustered' or self.distribution == 'clu':
            x = rnd.normal(self.size/2.,halo_rad,self.center_number)
            y = rnd.normal(self.size/2.,halo_rad,self.center_number)
            fold_back_C(x,np.array([self.size]))
            fold_back_C(y,np.array([self.size]))
            print size,max(x)
            print size,min(x)
            z = rnd.rand(self.center_number)*self.size
            self.center_positions = np.array([x,y,z]).transpose()
        else:
            print "Distribution not supported yet, generating homogeneous distribution."
            
        #print "generated center distribution. number of centers:",num_cen
        print "initiated centers in {0:.0f}s.".format(t_inicen-t.time())
    
    def __str__(self,show='all'):
        """return string for printing"""
        return "RelaxCenters instance of name "+str(self.name)+\
            "\ndensity:      "+str(self.density)+\
            "\nsize:         "+str(self.size)+\
            "\nC:            "+str(self.C)+\
            "\nD:            "+str(self.D)+\
            "\nb:            "+str(self.b)+\
            "\ndistribution: "+str(self.distribution)+\
            "\nn/o centers:  "+str(self.center_number)+"\n"
        
    def dump(self,filename="",path="./",overwrite=False):
        """dump the centers to a file
        optional input:
            filename    filename of dump, defaults to self.name
            path        dump directory relative to ".", defaults to "./"
            overwrite   bool if existing file should be overwritten, defaults to False
        """
        if filename=="":
            filename = self.name
        path = os.path.join(os.getcwd(),path)
        filename = os.path.join(path,filename,".relcen")
        import errno
        # test if directory can be created
        try :
            os.makedirs(path)
            print "creating dump directory",directory
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise
            print "existing dump directory for Centers instance",directory
        if ( not os.path.exists(filename) ) or overwrite == True:
            filehandler = open(path,"wb")
            pickle.dump(self,filehandler)
            filehandler.close()
        else:
            raise RelaxError(3,"RelaxCenters.dump() failed because file existing and not to be overwritten.")
        
    def plot(self, figure=None, axes=None, fname="",integrated=False,**kwargs):
        """plot a 2D projection of relaxation rate with cross section through z=1/2*size
        input:
            figure      figure to plot in
            axes        axes to plot in
            fname       supply string if you want to save figure
            integrated  if True, no axes attributes like title, labels, will be set
            **kwargs    sent to axes.contourf(**kwargs)
            """
        t_pltcen_b = t.time()
        cen = self.center_positions.flatten()
        dx = self.b/2.           # coarse grid for plotting
        dt = dx**2./(6*self.D)
        #print "dt",dt,"    dx",dx
        xr = np.arange(0,+self.size+dx,dx)       # regular grid for sampling
        yr = np.arange(0,+self.size+dx,dx)
        xlen = len(xr)
        #print "number of grid points for plotting centers:",xlen**2
        # make grid (xr and yr now of 2d with squared as many elements, giving x,y coordinates respectively)
        xr,yr = np.meshgrid(xr,yr)
        xr = xr.flatten()      # flatten for loop
        yr = yr.flatten()
        rr = np.empty(xr.shape,dtype=xr.dtype)
        # begin rate scanning
        array_for_rate = np.array([self.b,self.size]) # saving time making the array before loop
        middle_point = 0.5*self.size                   # same here
        for point in range(len(xr)):
            rr[point] = accum_dist_C(cen,array_for_rate,np.array([xr[point],yr[point],middle_point]))
            if rr[point] == -1: #accum_dist_C returns -1, if inside radius: set to value on border
                rr[point] = self.b**-6.
        rr = 1.-(dt*self.C*rr)        # apply relaxation rate
        rr = rr.clip(0,1)             # clamp the array (not really necessary...)
        rrmin = min(rr)
        # reshape to 2d for plotting
        rr = rr.reshape((xlen,xlen))
        yr = yr.reshape((xlen,xlen))
        xr = xr.reshape((xlen,xlen))
        # prepare plot
        if figure == None and axes == None:
            figure = plt.figure()
        if axes == None and figure != None:
            axes = figure.add_subplot(111)
        levels = np.linspace(rrmin,1,num=7) # levels for the contour plot
        cen_plot = axes.contourf(xr,yr,rr, levels=levels, rstride=1, cstride=1,**kwargs)
        if not integrated:
            bar = plt.colorbar(cen_plot)      # colorbar at right side of plot
            import matplotlib.ticker as tick
            bar.formatter = tick.FormatStrFormatter('%6f')
            bar.update_ticks()
            axes.set_title("Relaxation factor $f_R$ in $z=0$ plane ($m^{n+1}=f_R m^n$)\n\n")
            axes.set_xlabel("$x$")
            axes.set_ylabel("$y$")
            axes.set_aspect('equal', adjustable='box') # make x/y axes of same aspect and adjust whole box to fit the plit (no white spaces)
        if fname != "":
            plt.savefig(fname[:-3]+".pdf", format='pdf')
        #print "plotted centers in {0:.0f}s.".format(t_pltcen_b-t.time())
        if figure != None:
            return figure,axes
        else:
            return True


def load_centers(filename,path="./"):
    """load a RelaxCenters instance from a pickle
    
    input:
        filename    filename of pickle
        path        directory of file relative to ".",defaults to "./" 
    
    return:
        a Centers instance if successfull
        otherwise: raise error"""
    path = os.path.join(os.getcwd(),path,filename)
    if os.path.isfile(path):
        filehandler = open(path,"rb")
    else:
        raise RelaxError(4,"load_centers() did not receive a valid filename.")
    centers = pickle.load(filehandler)
    filehandler.close()
    if isinstance(centers,RelaxCenters):
        return centers
    else:
        raise RelaxError(2,"Loaded pickle is not a RelaxCenters instance.")


class RelaxResult():
    """A class for a relaxation experiment result."""
    def __init__(self, timesamples=1000, experiment=None, name="",**kwargs):
        """constructor
        optional input:
            timesamples     for interpolated data
            experiment      directly call read_experiment with a RelaxExperiment instance
            name            name of the result (for saving to file)
            **kwargs        send to RelaxResult.read_experiment() call, when called with experiment"""
        self.timesamples   = timesamples
        self.dataread   = False
        self.name       = name
        if experiment != None:
            self.read_experiment(experiment,**kwargs)
    
    def read_experiment(self, experiment, overwrite=False, logscale=False):
        """write data to internal variables, interpolating to less points, fitting to exponential decay
        
        input:
            experiment      an experiment instance
                       
        optional input:
            overwrite       bool if existing data in RelaxResult instance should be overwritten
            logscale        interpolate the experiment data logarithmic with t_ev/steps many points, linear if False
        """
        if not isinstance(experiment,RelaxExperiment):
            raise RelaxError(5,"experiment in RelaxResult.read_experiment() is not an RelaxExperiment instance.")
        
        if experiment.finished != True:
            raise RelaxError(6,"experiment in RelaxResult.read_experiment() is not finished.")
        if self.dataread==True and overwrite==False:
            raise RelaxError(6,"Some data has already been read into RelaxResult instance. Use optional argument to overwrite.")
        
        self.experiment = experiment
        self.evolution_time = experiment.evolution_time
        self.samplestep = self.evolution_time/self.timesamples
        self.method     = experiment.method
        if self.name == "":
            self.name   = experiment.name[:-3] + "res"
        self.size       = experiment.size
        self.b          = experiment.b
        self.C          = experiment.C
        self.D          = experiment.D
        self.density            = experiment.density
        self.center_positions   = experiment.center_positions
        
        if self.method == "randomwalks":
            self.walkers_number = experiment.walkers_number
        
        self.logscale = logscale
            
        import scipy
        import scipy.interpolate
        # set up new timearray
        if logscale == False:
            self.timearray = np.array([i*self.samplestep for i in range(self.timesamples)])
        else:
            self.timearray = log_range(experiment.timearray[1],self.evolution_time-self.samplestep,self.timesamples)
        # interpolate old timearray with an interpolating function
        magninterpfunc = scipy.interpolate.interp1d(experiment.timearray,experiment.magn)
        # interpolate magnetization
        self.magnarray  = magninterpfunc(self.timearray)
        # interpolate other data for random walks
        if self.method == "randomwalks" or self.method=="ran":
            magnerrinterpfunc   = scipy.interpolate.interp1d(experiment.timearray,experiment.magnerr)
            self.magnerrarray   = magnerrinterpfunc(self.timearray)
            actwalinterpfunc    = scipy.interpolate.interp1d(experiment.timearray,experiment.activewalkers)
            self.activewalkers  = actwalinterpfunc(self.timearray)
        elif self.method == "deterministic" or self.method=="det":
            self.magnerrarray   = None
            self.activewalkers  = None
        else:
            raise RelaxError(12,"Method of experiment not known: "+self.method)
        self.dataread = True
        print "Data has successfully been read into RelaxResult. Now fitting data to exponential decay."
        
        # fit data to exponential decay
        def strexpdecay(par,x):
            return np.exp(-(x/par[0])**par[1])
        strexpmodel = odr.Model(strexpdecay)
        # with error or without?
        odrdata = odr.RealData(x=self.timearray,y=self.magnarray)
        par0 = [ 1., 1.]
        odrinst = odr.ODR(odrdata, strexpmodel, beta0=par0,maxit=1000)
        odrout = odrinst.run()
        
        parfit = odrout.beta
        parerr = odrout.sd_beta
        #print odrout.stopreason
        self.fitT1 = parfit[0]
        self.fitbeta = parfit[1]
        self.fitT1err = parerr[0]
        self.fitbetaerr = parerr[1]
        self.fitresidualvar = odrout.res_var
        
        self.fitarray = strexpdecay(parfit,self.timearray)    
        
        print "Fit yields:   T1 =",self.fitT1,"+/-",self.fitT1err
        print "Fit yields: beta =",self.fitbeta,"+/-",self.fitbetaerr
        
    def update_data(self,experiment): 
        """update existing data with another data set from new experiment
        ## DO NOT USE THIS, NOT DONE ##"""
        #TODO
        #shall be made working for deterministic (more center distributions with same density,...)
        #as well as for random walks (more walks, more center distributions with same density,...)
        #failsafe: check for same values of C,D,b,density,etc?
        if not isinstance(experiment,RelaxExperiment):
            raise RelaxError(5,"experiment in RelaxResult.read_experiment() is not an RelaxExperiment instance.")
        print "Data has successfully been updated into RelaxResult."
                
    def write_hdf(self,filename=""):
        """write the results to a hdffile
        optional input:
        filename    filename of dump, defaults to self.name, increases counter if file exists
        """
        if filename=="":
            filename = self.name[:-3]
        
        filename,fname = os.path.split(filename)
        pname = filename

        # prepare directory
        directories=[]
        currentpath = os.getcwd()
        while filename!="":
            filename,nextdir=os.path.split(filename)
            directories.append(nextdir)
        directories.reverse()
        for directory in directories:
            try:
                os.makedirs(directory)
            except OSError as err:
                if err.errno != errno.EEXIST:
                    raise
            os.chdir(directory)

        # prepare filename
        if len(fname.split("."))>1:
            fext = "."+fname.split(".")[-1]
            fbody = fname[:-len(fext)]
            if not (fext == ".h5" or fext == ".hdf" or fext == ".hdf5"):
                fext=".hdf"
                fbody=fname
        else:
            fbody = fname
            fext=".hdf"

        # check for old names, rename
        count=1
        if os.path.exists(fbody+fext):
            while count<100:
                if not os.path.exists(fbody+"({:02})".format(count)+fext):
                    fname=fbody+"({:02})".format(count)+fext
                    break
                count += 1
            else:
                raise RelaxError(15,"Saving HDF files: too many files with same name!")
        else:
            fname = fbody+fext

        os.chdir(currentpath)
        filename = os.path.join(pname,fname)
        print "Saving Result in",filename
        
        import inspect
        import h5py
        h5fh = h5py.File(filename)
        h5strtype = h5py.special_dtype(vlen=str)        # need to define a string type
        
        # saving module and script
        moddata = h5fh.create_dataset('module source',(1,),dtype=h5strtype)
        modpath = inspect.getabsfile(self.__class__)
        modfh = open(modpath)
        moddata[0] = modfh.read()
        moddata.attrs["module path"] = modpath
        modfh.close()

        scrdata = h5fh.create_dataset('script source',(1,),dtype=h5strtype)
        scrpath = os.path.abspath(sys.modules['__main__'].__file__)
        scrfh = open(scrpath)
        scrdata[0] = scrfh.read()
        scrdata.attrs["script path"] = scrpath
        scrfh.close()
        
        # saving input data
        inpdata = h5fh.create_dataset('input data',(1,),dtype=h5strtype)
        inpdata[0] = "See attributes! All values in basic SI units."
        inpdata.attrs["size"]       = self.size
        inpdata.attrs["density"]    = self.density
        inpdata.attrs["distribution"] = self.experiment.centers.distribution
        inpdata.attrs["mean centers distance"] = (self.density*math.pi*4/3.)**(-1/3.)
        inpdata.attrs["number of centers"] = self.density*self.size**3
        inpdata.attrs["C"]          = self.C
        inpdata.attrs["D"]          = self.D
        inpdata.attrs["b"]          = self.b
        inpdata.attrs["method"]     = self.method
        inpdata.attrs["evolution time"] = self.evolution_time
        inpdata.attrs["finishing time of simulation"] = time.strftime("%Y/%m/%d %H:%M:%S")
        if self.method == "randomwalks":
            inpdata.attrs["number of walkers"] = self.walkers_number
        
        # saving output
        outdata = h5fh.create_dataset('output data',(len(self.timearray),4),dtype='d')
        outdata[:,0] = self.timearray
        outdata[:,1] = self.magnarray
        if self.method == "randomwalks":
            outdata[:,2] = self.magnerrarray
            outdata[:,3] = self.activewalkers
        else:
            outdata[:,2:] = -1
        outdata.attrs["column data"] = "0: timearray, 1: magnarray, 2: magnerrarray, 3: activewalkers"
        outdata.attrs["data in log_range"] = str(self.logscale)
        
        
        if self.method == "randomwalks":
            outdata.attrs["step length"] = self.experiment.dx
            outdata.attrs["dt"] = self.experiment.dt
            outdata.attrs["number of steps"] = self.experiment.steps_number
            outdata.attrs["estimated MSD"] = self.experiment.MSD
        else:
            outdata.attrs["dx"] = self.experiment.dx
            outdata.attrs["dt"] = self.experiment.dt
            outdata.attrs["lattice points"] = self.experiment.lattice_points
            outdata.attrs["time steps"] = self.experiment.time_steps
        
        outdata.attrs["T1 time"] = self.fitT1
        outdata.attrs["T1 error"] = self.fitT1err
        outdata.attrs["stretching beta"] = self.fitbeta
        outdata.attrs["stretching beta error"] = self.fitbetaerr
        outdata.attrs["residual variance"] = self.fitresidualvar
        
        # saving plots
        # TODO!!
        
        h5fh.close()
        return filename
        
    def plot_data(self,figure=None,axes=None,activewalkers=False,ploterror=True,plotfit=True,plotmagn=True,logx=False,logy=False,label="",plottitle='',legendsize=20,**kwargs):
        """Plot the data using matplotlib.
        optional input:
            figure          a matplotlib figure
            axes            plot axis
            activewalkers   plot active walkers
            ploterror       plot the error
            plotfit         plot exponential fit
            plotmagn        plot data (if False, you MUST specify a color, if you don't want all lines to be of different color!)
            logx, logy      logarithmic x, y scale 
            label           label for plot
            plottitle       title of plot
            legendsize      size of legend, default: 20
            **kwargs          optional arguments passed to axes.plot
        return:
            (figure,axes)   as tuple
        """
        import matplotlib.pyplot as plt
        if figure==None and axes==None:
            figure = plt.figure()
        if axes == None:
            axes = figure.add_subplot(111)
        # make color of error/walkers same color as magnetization
        color = None
        if "color" in kwargs:
            color = kwargs.pop("color")
        # plot magnetization
        if plotmagn:
            if color != None :
                if plotfit:             # point markers are too small for legend
                    plmag    = axes.plot(self.timearray,self.magnarray,',',color=color,**kwargs)
                else:
                    plmag    = axes.plot(self.timearray,self.magnarray,',',color=color,label=label,**kwargs)
            else :
                if plotfit:
                    plmag    = axes.plot(self.timearray,self.magnarray,',',**kwargs)
                else:
                    plmag    = axes.plot(self.timearray,self.magnarray,',',label=label,**kwargs)
            color = plmag[0].get_color()
        # plot error if applicable
        if self.magnerrarray!=None and ploterror:
            plmagperr = axes.plot(self.timearray,self.magnarray+self.magnerrarray,'--',color=color,**kwargs)
            plmagmerr = axes.plot(self.timearray,self.magnarray-self.magnerrarray,'--',color=color,**kwargs)
        # plot activewalkers if applicable
        if activewalkers and self.activewalkers!=None:
            plactwalk = axes.plot(self.timearray,self.activewalkers,'-',color=color,**kwargs)
        if plotfit:
            plfit = axes.plot(self.timearray,self.fitarray,'-',label=label, color=color)
        # adjust plot options
        if logx == True or self.logscale == True:
            axes.set_xscale('log')
        if logy == True:
            axes.set_yscale('log')
        axes.set_xlabel("$t/$s")
        axes.set_ylabel("$m(t)/m_0$")
        axes.set_ylim(-0.1,1.1)
        axes.set_xlim(0,self.timearray[-1])
        axes.set_title(plottitle)
        axes.legend(loc=0,prop=dict(size=legendsize))
        return figure,axes
            
    def __str__(self):
        """return string representation of RelaxResult instance."""
        return self.name+" has self.centers information:\n"+self.centers.__str__()


class RelaxExperiment():
    """A class for the actual "Experiment" on a center distribution.
    Instances are meant to be transient and not saved/duplicated."""
    def __init__(self, centers, evolution_time, dxb_ratio=1/3., method="randomwalks",**kwargs):
        """constructor
        
        input:
            centers         a instance of the Centers class
            evolution_time  in second
            
        optional input:
            dxb_ratio       ratio of step size/mesh size over quenching radius b
            method          "randomwalks"/"ran" or "deterministic"/"det"
            **kwargs        passed on to specific initializer
            
            deterministic:
                dt          set a timestep manually
            randomwalks:
                none
        """
        if not isinstance(centers,RelaxCenters):
            raise RelaxError(1,"center input in RelaxExperiment instance is not of class RelaxCenters!")
        else:
            self.centers        = centers
        #copy data from centers
        self.name       = centers.name[:-3] + "exp"
        self.size       = self.centers.size
        self.C          = self.centers.C
        self.D          = self.centers.D
        self.b          = self.centers.b
        self.dx         = self.b*dxb_ratio
        self.density            = self.centers.density
        self.center_positions   = centers.center_positions
        
        self.evolution_time = evolution_time
        self.finished       = False
        
        # check values
        if self.size < 10*self.dx or self.b < self.dx:
            raise RelaxError(15,"RelaxExperiment.init_deterministic(): dx too large")
        
        self.method = method
        if self.method=="randomwalks" or self.method == "ran":
            self._init_randomwalks(**kwargs)
        elif method=="deterministic" or self.method == "det":
            self._init_deterministic(**kwargs)
        else:
            raise RelaxError(12,"Method of experiment not known. Choose 'deterministic' or 'randomwalks'.")        

    def run_experiment(self,redo=False,**kwargs):
        """run the experiment with data provided by the constructor.
        Keyword arguments are passed to experimenter.
        
        optional input:
            redo        boolean of finished experiment should be redone
            **kwargs    passed on to specific initializer:
            
            deterministic:
                background      apply background relaxation rate,
                                 0 for no relaxation,
                                 ~7e-2 /s as usual rate (4mm sample klempt2003 paper)
                
            random walks:
                walks       number of walks to be done
                plotaxes    to plot the actual walk, give an axis
                background  apply background relaxation, 0 for no relaxation,
                            usual rate: ~7e-2 /s (4mm sample klempt2003 paper)
        """
        if self.finished == True:
            if redo != True:
                raise RelaxError(13,"Experiment already finished. Use optional argument to redo experiment.")
            else:
                self.__init__(self.centers,self.evolution_time,method=self.method)
        if self.method=="randomwalks" or self.method == "ran":
            self._run_randomwalks(**kwargs)
        elif self.method=="deterministic" or self.method == "det":
            self._run_deterministic(**kwargs)
        else:
            raise RelaxError(12,"Method of experiment not known.")
        self.finished = True
        
    def _init_deterministic(self, dt=0):
        """Test values and creating matrices from it.

        input:
            dt:         evolution timestep in seconds
        """
        t_detini_b = t.time()

        if dt == 0:
            self.dt = 1./(12.*self.D/(self.dx*self.dx) + self.C*(self.b/2)**(-6))
            print "dt: 1/(",12.*self.D/(self.dx*self.dx),"+",self.C*(self.b/2)**(-6),") s"
            print "dt=",self.dt,"s"
        else :
            self.dt = dt
            
        if self.evolution_time < self.dt:
            raise RelaxError(16,"RelaxExperiment.init_deterministic(): dt too long")

        self.lattice_points = int(round(self.size/self.dx))
        self.time_steps     = int(round(self.evolution_time/self.dt))

        print "Initializing deterministic simulation with"
        print "system size:           ", self.size, "(",self.lattice_points,"pts in one dim,",7*self.lattice_points**3,"in matrix)"
        print "evolution time:        ", self.evolution_time,"s (",self.time_steps,"steps)"
        print "   dt =",self.dt,"   dx =",self.dx
        
        # helper functions To Index and From Index
        def ti(j,l,m):
            return j%self.lattice_points +\
                    (l%self.lattice_points)*self.lattice_points +\
                    (m%self.lattice_points)*self.lattice_points**2
        def fo(i):
            return i%self.lattice_points, int(i/self.lattice_points), int(i/self.lattice_points**2)

        # the magnetization, discretized on a lattice, initialized to one (m0)
        self.lattice = np.ones(self.lattice_points**3)
        # vector timearray holding the time steps
        self.timearray = np.arange(0,self.evolution_time+self.dt/2.,self.dt)
        
        if len(self.timearray) != self.time_steps+1:
            raise RelaxError(18,"Failure creating timearray[i]\nlen(self.t)     = "+str(len(self.t))+"\nself.time_steps = "+str(self.time_steps))

        # convert quenching radius and centers array from percentages to indices of the lattice
        self.b_l       = self.b / self.dx
        self.centers_l = (self.center_positions/self.dx).flatten()

        print "Creating D and C vectors."
        # one-dim  diffusion vector D[lat_pts**3] 
        self.D_vector = np.ones(self.lattice_points**3)
        # one-dim relaxation vector C[lat_pts**3] 
        self.C_vector = np.zeros(self.lattice_points**3)
        # one-dim magnetization weight vector M[lat_pts**3] (which points to take into magnetization calculation)
        self.M_vector = np.ones(self.lattice_points**3)
        self.Meqzero = 0   # counter of zero entries of M
        
        # fill D,C,B vectors
        array_for_rate = np.array([self.b_l,self.lattice_points],dtype="d") # saving time making the array before loop
        progress = Progress(self.lattice_points)
        for j in range(self.lattice_points):
            for l in range(self.lattice_points):
                for m in range(self.lattice_points):
                    ind_jlm = ti(j,l,m) # current index, no need to calculate every time
                    self.C_vector[ind_jlm] = accum_dist_C(self.centers_l,array_for_rate,np.array([j,l,m],dtype="d"))
                    if self.C_vector[ind_jlm] == -1 :
                        # if inside any quenching radius, set C_vec to value on border, diffusion to zero
                        self.C_vector[ind_jlm] = self.b_l**-6
                        self.D_vector[ind_jlm] = 0
                        self.M_vector[ind_jlm] = 0
                        self.Meqzero += 1
            progress.increment()
            progress.print_status_line()
        # amend vectors by constants
        self.D_vector *= self.D*self.dt*self.dx**-2
        self.C_vector *= self.C*self.dt*self.dx**-6

        # RHS matrix B, according to p 15, eq (IX)
        print "Creating matrix B."
        import scipy.sparse as spsp
        from scipy.sparse import linalg
        self.B_matrix = spsp.lil_matrix((self.lattice_points**3, self.lattice_points**3))
        progress = Progress(self.lattice_points)
        for j in range(self.lattice_points):
            for l in range(self.lattice_points):
                for m in range(self.lattice_points):  # prefactors are already in D_vec/C_vec
                    ind_jlm = ti(j,l,m)               # these two lines should speed up
                    D_vec_entry = self.D_vector[ind_jlm] #  the process significantly
                    self.B_matrix[ind_jlm,ind_jlm]   = 1 -6*D_vec_entry - self.C_vector[ind_jlm]
                    self.B_matrix[ind_jlm,ti(j,l-1,m)] = D_vec_entry
                    self.B_matrix[ind_jlm,ti(j-1,l,m)] = D_vec_entry
                    self.B_matrix[ind_jlm,ti(j+1,l,m)] = D_vec_entry
                    self.B_matrix[ind_jlm,ti(j,l+1,m)] = D_vec_entry
                    self.B_matrix[ind_jlm,ti(j,l,m-1)] = D_vec_entry
                    self.B_matrix[ind_jlm,ti(j,l,m+1)] = D_vec_entry

            progress.increment()
            progress.print_status_line()

        print "converting matrix B to CSR format"
        self.B_matrix = self.B_matrix.tocsr()

        print "deterministic simulation initialized in {0:.0f}s.".format(t.time() - t_detini_b)
        
    def _run_deterministic(self, background = 0, withM=False):
        """ Evolve the system. Write magnetization per step in array.
        optional argument:
            background      apply background relaxation rate,
                             0 for no relaxation,
                             ~7e-2 /s as usual rate (4mm sample klempt2003 paper)
        """
        t_detrun_b = t.time()
        print "Evolving system with", self.time_steps, "timesteps"

        def calc_magn(lattice,M_vector):
            # this works just fine, magnetization near centers decays as if on border
            # (->faster than anything else) and does not contribute after several steps,
            # and due to small quenched space compared to bulk volume
            if withM:
                return np.dot(lattice,M_vector)/(len(lattice)-self.Meqzero)
            else:
                return float(lattice.sum())/len(lattice)

        # prepare array of magnetization
        self.magn      = np.zeros(self.time_steps+1)
        self.magn[0]   = calc_magn(self.lattice,self.M_vector)
        
        progress = Progress(self.time_steps/1000)
        for step in range(self.time_steps):
            self.lattice = self.B_matrix*self.lattice     # calculate time step
            self.magn[step+1] = calc_magn(self.lattice,self.M_vector)   # calculate magnetization
            if step%1000 == 0:
                progress.increment()
                progress.print_status_line()

        #apply background relaxation
        if background != 0:
            perstep = exp(-self.dt*background)
            self.magn *= perstep**np.arange(self.time_steps+1)

        print "ran deterministic experiment in {0:.0f}s.".format(t.time() - t_detrun_b)

    def _init_randomwalks(self):
        """Initialize random walks"""
        self.dt   = self.dx**2./(6*self.D)
        if self.evolution_time < self.dt:
            raise RelaxError(16,"RelaxExperiment.init_randomwalks(): dt > t_evo")
        if self.dt*self.C/self.b**6>1:
            print "WARNING: dt and dx were reset to meet stability at radius b!"
            self.dt = self.b**6/self.C*0.9
            self.dx = np.sqrt(6*self.D*self.dt)
        self.steps_number = int(self.evolution_time/self.dt)
        self.MSD  = math.sqrt(6*self.D*self.evolution_time) 
        
        # plotting some information
        print "initializing walks"
        print "number of steps:    {}".format(self.steps_number)
        print "length of steps:    {:.2} m".format(self.dx)
        print "time step:          {:.2} s".format(self.dt)
        print "total time:         {:.2} s".format(self.evolution_time)
        print "estimated MSD:      {:.2} m".format(self.MSD)
        print "simulation box size:{:.2} m".format(self.size)
        
        self.magn           = np.zeros( self.steps_number+1   ,dtype='d')
        self.magnvar        = np.zeros( self.steps_number+1   ,dtype='d')
        self.magnerr        = np.zeros( self.steps_number+1   ,dtype='d')
        self.activewalkers  = np.zeros( self.steps_number+1   ,dtype='d')
        self.activewalkersvar = np.zeros( self.steps_number+1 ,dtype='d')
        self.timearray      = np.arange(0, self.dt*(self.steps_number+1/2.), self.dt)
    
    
    def _run_randomwalks(self,walks=10,plotaxes=None,background = 0):
        """Run random walks.
        
        optional arguments
            walks           number of walks to be done
            plotaxes        to plot the actual walk, give an axis
            background      apply background relaxation rate,
                             0 for no relaxation,
                             ~7e-2 /s as usual rate (4mm sample klempt2003 paper)
        """

        self.walkers_number = walks
        print "number of walkers:  ", walks
        
        def random_spherical_coord():
            """generate random position on unit 3D sphere around origin"""
            while True :
                #step = rnd.rand(3)*2-1
                step = rnd.uniform(-1,1,3)
                norm = reduce(lambda x,y:x+y**2,step)
                if norm <= 1 :
                    break
            return step/norm
        
        magntmp = np.zeros( self.steps_number+1   ,dtype='d')   # temp magnetization
        postmp  = np.zeros((self.steps_number+1)*3,dtype='d')   # temp position
        actwalktmp =  np.zeros(self.steps_number+1,dtype='d')   # temp active walker         
                
        def relax_rw(position, magnetization):
            """Uniform random walk with relaxating magnetization.
            Return last step, where the walker comes to rest
            
            Input:
                position        position vector             len = 3*(steps+1)
                magnetization   magnetization               len = steps+1
            
            magnetization   contains the magnetization
            position        contains the walker's positions"""
            
            array_for_rate = np.array([self.b,self.size],dtype="d")
            centers = self.center_positions.flatten()
            # start the steps
            for step in range(self.steps_number):
                # accumulate distance**-6
                dist6 = accum_dist_C(centers, array_for_rate, position[step:step+3])
                if dist6 == -1 :
                    dist6 = 1.-self.dt*self.C*self.b**-6            # calculate relaxation rate
                    position[3*step+3::3] = position[3*step  ]      # do not step anymore
                    position[3*step+4::3] = position[3*step+1]
                    position[3*step+5::3] = position[3*step+2]
                    magnetization[step+1:] = magnetization[step]* dist6**np.arange(self.steps_number-step)     # relaxate by simple multiplication        
                    break                                           # break out of stepping loop
                else:
                    dist6 = 1.-self.dt*self.C*dist6                 # calculate relaxation rate
                    magnetization[step+1] = magnetization[step]*dist6      # apply relaxation rate (clamp included)
                    step3 = step*3
                    step_coords = random_spherical_coord()*self.dx
                    position[step3+3:step3+6] = position[step3:step3+3]+step_coords # make the step
                    #correct for periodic boundaries
                    if position[step3+3]>self.size:
                        position[step3+3] -= self.size
                    if position[step3+3]<0:
                        position[step3+3] += self.size
                    if position[step3+4]>self.size:
                        position[step3+4] -= self.size
                    if position[step3+4]<0:
                        position[step3+4] += self.size
                    if position[step3+5]>self.size:
                        position[step3+5] -= self.size
                    if position[step3+5]<0:
                        position[step3+5] += self.size
            return step       # returns last step
        
        #start the walks
        progress = Progress(walks)
        for walk in range(walks):
            # reset first three values to empty temp arrays
            magntmp[0]    = 1
            postmp[:3]    = rnd.rand(3)*self.size # random starting position within size
            actwalktmp[:] = 1
            #print "doing walk",walk+1
            last_step = relax_rw(postmp, magntmp)
            actwalktmp[last_step+1:] = 0
            if plotaxes != None :
                # plot the walk (only x,y dimensions)
                rwplot = plotaxes.plot(postmp[0::3],postmp[1::3],'k-')
            
            if walk == 0:  # if first walk, assign directly
                self.magn[:]        = magntmp[:]
                self.activewalkers  = actwalktmp[:]
            else : # from second on, update values
                self.magn, self.magnvar, devnull =\
                        _upd_avg_var( self.magn, self.magnvar, walk, magntmp )
                self.activewalkers, devnull, devnull =\
                        _upd_avg_var( self.activewalkers, self.activewalkersvar, walk, actwalktmp )
            
            progress.increment()
            progress.print_status_line()
        
        # apply background relaxation
        if background != 0:
            perstep = np.exp(-self.dt*background)
            self.magn *= perstep**np.arange(self.steps_number+1)
            
        self.magnerr   = np.sqrt( self.magnvar/walks )
        #print 'walks_done:',walkers_number
        #print "max of magnvar",max(self.magnvar)
        
        return 0






