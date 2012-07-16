RelaxSim
========

Nuclear Magnetization Relaxation Simulation by paramagnetic impurities in LiF and other Fluorine Systems

planned additions
-----------------
(sorted by urgency)

* do not remember every step => save RAM (in fake simulation as well!)
* continuous quenching
* omit steps if far from centers (can be part of continuous quenching)
* clean up variables (is `C` really a local var of Experiment?)
* make walk class?!
* time simulations (save start/stop time in hdf)
* deterministic method for magnetization of single spins (big rate equation) (=> BPP)

* * *

changelog
---------

v5 (not tagged yet):
additions:
* implement "fake simulation" for fast execution while debugging:
`fake=True` in `RelaxExperiment()`

v4:
additions:

* optional dict for `RelaxResult.write_hdf()`
* switch `do_step` in `RelaxExperiment._run_randomwalks()` to switch stepping on/off after touching a quenching area

v3:
First commit.
For changes see file `relaxsim.py`

