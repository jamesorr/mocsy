import numpy as np
from mocsy import *

invar = ['Alk ','DIC ','Phos','Sil ','T   ','S   ']

optCON  = 'mol/kg' 
optT    = 'Tinsitu'
optP    = 'm'      
optB    = 'l10'
#optB    = 'u74'
optK1K2 = 'l'
optKf   = 'dg'

temp   = 18.0  
sal    = 35.0  
alk    = 2300.*1.e-6
dic    = 2000.*1.e-6
sil    = 60.0e-6  
#sil    = 0.0e-6
phos   = 2.0e-6   
#phos   = 0.0e-6
depth  = 0.
Patm   = 1.0        
lat    = 0.

[ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis] =    mvars.vars(
temp, sal, alk, dic, sil, phos, Patm, depth, lat,
optCON, optT, optP, optb=optB, optk1k2=optK1K2, optkf=optKf    )

[ph_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv] = mderivauto.derivauto(
temp, sal, alk, dic, sil, phos, Patm, depth, lat, 
optCON, optT, optP, optb=optB, optk1k2=optK1K2, optkf=optKf    )

# [H+] concentration
H = 10.0**(-ph[0])
# derivative of [H+] deduced from that of pH
H_deriv = - H * ph_deriv * np.log(10.0)

invar = ['alk', 'dic', 'pho', 'sil', 'tem', 'sal']

# output results
#print "      dH/dx         dpCO2/dx       dfCO2/dx       d[CO2*]/dx      d[HCO3-]/dx    d[CO3--]/dx     dOmegaA/dx     dOmegaC/dx"
print "derivauto: Absolute derivatives" 
print "             dh_dx         dpco2_dx       dfco2_dx         dco2_dx      dhco3_dx       dco3_dx       dOmegaA_dx     dOmegaC_dx"
for i in range(6):
       print "%3s  :  %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e" %  \
       (invar[i], H_deriv[i,0], pco2_deriv[i,0], fco2_deriv[i,0], co2_deriv[i,0], hco3_deriv[i,0], co3_deriv[i,0], omegaa_deriv[i,0], omegac_deriv[i,0])
       # print invar[i], H_deriv[i,0], pco2_deriv[i,0], fco2_deriv[i,0], co2_deriv[i,0], hco3_deriv[i,0], co3_deriv[i,0], omegaa_deriv[i,0], omegac_deriv[i,0]

# output results: another way 
#print np.hstack ((H_deriv, pco2_deriv, fco2_deriv, co2_deriv, hco3_deriv, co3_deriv, omegaa_deriv, omegac_deriv))
