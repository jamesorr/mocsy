import numpy as np
from mocsy import *

temp    = np.array([-1, 12, 25])  
salt    = temp*0 + 35.

phi0_cfc11 =    gasx.solgas('cfc11', temp, salt)
phi0_cfc12 =    gasx.solgas('cfc12', temp, salt)
phi0_sf6   =    gasx.solgas('sf6',   temp, salt)
phi0_co2   =    gasx.solgas('co2',   temp, salt)
phi0_n2o   =    gasx.solgas('n2o',   temp, salt)

print "phi0_cfc11 =", phi0_cfc11
print "phi0_cfc12 =", phi0_cfc12
print "phi0_sf6   =", phi0_sf6
print "phi0_co2   =", phi0_co2
print "phi0_n2o   =", phi0_n2o


