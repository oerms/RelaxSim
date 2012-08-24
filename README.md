RelaxSim
========

Nuclear Magnetization Relaxation Simulation by paramagnetic impurities in LiF and other Fluorine Systems.

The following changelog only lists changes in the module.
Additionally there are files for simulating certain situations.
The resulting `.hdf` and `.pdf` files are not part of the git repository.

planned additions and changes
-----------------------------
(kinda sorted by urgency)

* deterministic method for magnetization of single spins (big rate equation) (=> BPP)
* do not remember every step => save RAM (in fake simulation as well!)
* clean up variables (e.g. is `C` really a local var of Experiment?)
* clean up error codes in `RelaxError` class
* correct `free_steps` if there are too many free steps for the walker
* generally rework free steps at `stay` method
* save center number in hdf
* move `plot_datad3` (a center plot) out of `RelaxResult` into `RelaxCenters`
* test if continuous diffusion in deterministic method is magnetization conserving
* import results from HDF file for easy analysis
* make somehow clear, that: number of steps in HDF file != number of steps in simulation

* * *

changelog
---------

v6 (nothing tagged yet):  
additions:
* continuous quenching in deterministic method
* free steps now in every walk method, show percentage at end of walk (might be too verbose)
* setting centerpositions manually now possible, additional 3D center plot option

fixes/changes:
* free steps did add not actual steps of length ~1e-10 but of lengths 1
* added `<cmath>` when compiling `fold_back_C` and `find_nearest_C`
* `density` is now optional argument to `RelaxCenters()` (not neccessary when setting centers manually)
* fixed a bug in `stay` method: `array_for_rate` had wrong dimension, which led to segfaults

v5:  
additions:
* continuous quenching
* omit steps if far from centers (not only in continuous quenching)
* implement "fake simulation" for fast execution while debugging: `fake=True` in `RelaxExperiment()`
* add function `quenched_diffusion()` for continous diffusion quenching with radial dependence
* add function `quenched_step()` for continous quenched step: use an ellipse of step lengths
* new C function for searching two nearest centers

fixes/changes:
* `overlapl()` used gauss function instead of lorentzian => fix also plot `continbrad_thesis_overlmsdradii_neu.pdf`
* when fitting T1 and beta could be negative => modified `strexpdecay()` in `RelaxResult.read_experiment()`
* make parameter `b` optional input in `pull_center()`, default to `size/2` if not given
* `RelaxCenters()` now takes `tau` and `bfield` as optional parameters and `RelaxResult` saves it
* `RelaxExperiment._init_randomwalks()` now takes the parameter about what type of walk. Should be passed on via `**kwargs` from `RelaxExperiment()`.

v4:  
additions:

* optional dict for `RelaxResult.write_hdf()`
* switch `do_step` in `RelaxExperiment._run_randomwalks()` to switch stepping on/off after touching a quenching area

v3:  
First commit.
For changes see file `relaxsim.py`

