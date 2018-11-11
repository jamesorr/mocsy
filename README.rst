mocsy
=====

Routines to model ocean carbonate system thermodynamics

Synopsis: mocsy is a Fortran 95 package designed to compute all
carbonate system variables from total dissolved inorganic carbon (DIC)
and total alkalinity, particularly from models. It updates previous
OCMIP code, avoids 3 common model approximations, and offers the
best-practice constants as well as more recent options. It agrees with
CO2SYS within 0.005%.

The mocsy package is described by Orr and Epitalon (2015) and has been
compared to other packages that compute marine carbonate chemistry by
(Orr et al., 2015).  More recently, new routinnes were added to
propagate uncertainties and compute sensitivities of derived variables
to input variables (Orr et al., 2018)

**Documentation and Examples**

* Documentation: http://ocmip5.ipsl.jussieu.fr/mocsy/
* Example scripts: see *examples* directory
* Jupyter notebooks (interactive examples): see  *notebooks* directory

**REFERENCES**

Orr, J. C. and Epitalon, J.-M. (2015) Improved routines to model the
ocean carbonate system: mocsy 2.0, Geosci. Model Dev., 8, 485-499,
https://doi.org/10.5194/gmd-8-485-2015 .

Orr, J. C., J.-P. Gattuso, and J.-M. Epitalon (2015) Comparison of ten
packages that compute ocean carbonate chemistry, Biogeosciences, 12,
1483–1510, https://doi.org/10.5194/bg-12-1483-2015 .

Orr, J.C., J.-M. Epitalon, A. G. Dickson, and J.-P. Gattuso (2018) Routine
uncertainty propagation for the marine carbon dioxide system, in prep. for
Mar. Chem., in press, https://doi.org/10.1016/j.marchem.2018.10.006 .

