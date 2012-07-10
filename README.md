RelaxSim
========

Nuclear Magnetization Relaxation Simulation by paramagnetic impurities in LiF and other Fluorine Systems

planned additions
-----------------

* continuous quenching
* omit steps if far from centers
* make walk class?!
* do not remember every step => save RAM (in fake simulation as well!)
* clean up variables (is C really a local var of Experiment?)
* time simulations (save start/stop time in hdf)
* implement " fake simulation" for fast execution while debugging
* deterministic method for magnetization of single spins (big rate equation) (-> BPP)

* * *

changelog
---------

v4:
additions:

* optional dict for `RelaxResult.write_hdf()`
* switch `do_step` in `RelaxExperiment._run_randomwalks()` to switch stepping on/off after touching a quenching area

v3:
First commit.
For changes see file `relaxsim.py`

