import numpy as np
from mocsy import *

temp    = np.array([-1, 12, 25])  
salt    = temp*0 + 35.

phi0_cfc11 =    gasx.phizero('cfc11', temp, salt)
phi0_cfc12 =    gasx.phizero('cfc12', temp, salt)
phi0_sf6   =    gasx.phizero('sf6',   temp, salt)
phi0_co2   =    gasx.phizero('co2',   temp, salt)
phi0_n2o   =    gasx.phizero('n2o',   temp, salt)

print "phi0_cfc11 =", phi0_cfc11
print "phi0_cfc12 =", phi0_cfc12
print "phi0_sf6   =", phi0_sf6
print "phi0_co2   =", phi0_co2
print "phi0_n2o   =", phi0_n2o


