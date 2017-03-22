
# coding: utf-8

# change this as appropriate
mocsy_dir = "/home/username/Documents/mocsy"

## Preliminaries
import os, sys
import numpy as np
sys.path.append (mocsy_dir)
import mocsy



# Brief program to test subroutine vars()
# ---------------------------------------

## input arrays
temp = np.array ((18.0,))
sal = np.array ((35.0,))
alk = np.array ((2300.e-6,))
DIC = np.array ((2000.e-6,))
sil = np.array ((60.e-6,))
phos = np.array ((2.e-6,))
patm = np.array ((1.0,))
depth = np.array ((0.0,))
lat = np.array ((0.0,))

# computed variables at 1 input point
ph = np.zeros ((1,)).astype('f')
pco2 = np.zeros ((1,)).astype('f')
fco2 = np.zeros ((1,)).astype('f')
co2 = np.zeros ((1,)).astype('f')
hco3 = np.zeros ((1,)).astype('f')
co3 = np.zeros ((1,)).astype('f')
omegaa = np.zeros ((1,)).astype('f')
omegac = np.zeros ((1,)).astype('f')
betad = np.zeros ((1,)).astype('f')
rhosw = np.zeros ((1,)).astype('f')
p = np.zeros ((1,)).astype('f')
tempis = np.zeros ((1,)).astype('f')

# Call subroutine 'vars()' with TEOS-10 imput temperature and salinity.
ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis = \
    mocsy.mvars.vars(temp,sal, alk, DIC,
                    sil,phos,patm,depth,lat,optcon="mol/kg",optt='Tcsv',optp='db',
                    optb='u74', optk1k2='l', optkf='dg', opts='Sabs')

# Results
print "ph",   ph[0]
print "pco2", pco2[0]
print "fco2", fco2[0]
print "co2",  co2[0]
print "hco3", hco3[0]
print "co3",  co3[0]
print "omegaA", omegaa[0]
print "omegaC", omegac[0]

