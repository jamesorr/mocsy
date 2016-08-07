import numpy as np
from mocsy import *

invar = ['alk','dic','pho','sil','tem','sal','k0 ','k1 ','k2 ','kb ','kw ','ka ','kc ']

optCON  = 'mol/kg' 
optT    = 'Tinsitu'
optP    = 'm'      
optB    = 'l10'
#optB    = 'u74'
optK1K2 = 'l'
optKf   = 'dg'

temp   = 18.0  
sal    = 35.0  
alk    = 2300.e-6
dic    = 2000.e-6
phos   = 2.0e-6   
#phos   = 0.0e-6   
sil   = 60.0e-6  
#sil    = 0.0e-6  
depth  = 0.
Patm   = 1.0        
lat    = 0.

[ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis] =    mvars.vars(
temp, sal, alk, dic, sil, phos, Patm, depth, lat,
optCON, optT, optP, optb=optB, optk1k2=optK1K2, optkf=optKf    )

# [H+] concentration
H = 10.0**(-ph[0])

print "derivnum: Absolute derivatives" 
print "             dh_dx         dpco2_dx       dfco2_dx         dco2_dx      dhco3_dx       dco3_dx       dOmegaA_dx     dOmegaC_dx"

for i in range(13):
    [dh_dx, dpco2_dx, dfco2_dx, dco2_dx, dhco3_dx, dco3_dx, dOmegaA_dx, dOmegaC_dx] = \
    mderivnum.derivnum (
    temp, sal, alk, dic, sil, phos, Patm, depth, lat, invar[i],
    optCON, optT, optP, optb=optB, optk1k2=optK1K2, optkf=optKf    )
    #
    print "%3s  :  %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e" %  \
    (invar[i], dh_dx[0], dpco2_dx[0], dfco2_dx[0], dco2_dx[0], dhco3_dx[0], dco3_dx[0], 
    dOmegaA_dx[0], dOmegaC_dx[0] )

