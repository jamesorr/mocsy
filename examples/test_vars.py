
# coding: utf-8

# change this as appropriate
mocsy_dnad_dir = "/home/jean-marie/Documents/DNAD/mocsy"

## Preliminaries
import os, sys
import numpy as np
sys.path.append (mocsy_dnad_dir)
import mocsy



# Program that computes automatic derivative of (pCO2 [uatm] as a function of DIC [umol/kg])
# ------------------------------------------------------------------------------------------

## input arrays
#temp = np.array ((20.0,))
#sal = np.array ((35.0,))
#alk = np.array ((2300 * 1.e-6,))
#DIC = np.array ((2000 * 1.e-6,))
#sil = np.array ((0.0,))
#phos = np.array ((0.0,))
#patm = np.array ((1.0,))
#depth = np.array ((0.0,))
#lat = np.array ((0.0,))

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

# derivatives w/ respect to 6 input variables
ph_deriv    = np.zeros((6,)).astype('f')
pco2_deriv  = np.zeros((6,)).astype('f')
fco2_deriv  = np.zeros((6,)).astype('f')
co2_deriv   = np.zeros((6,)).astype('f')
hco3_deriv  = np.zeros((6,)).astype('f')
co3_deriv   = np.zeros((6,)).astype('f') 
omegaa_deriv  = np.zeros((6,)).astype('f') 
omegac_deriv  = np.zeros((6,)).astype('f')

ph_deriv    = pco2_deriv.reshape ((6,1), order='F')
pco2_deriv  = pco2_deriv.reshape ((6,1), order='F')
fco2_deriv  = pco2_deriv.reshape ((6,1), order='F')
co2_deriv   = pco2_deriv.reshape ((6,1), order='F')
hco3_deriv  = pco2_deriv.reshape ((6,1), order='F')
co3_deriv   = pco2_deriv.reshape ((6,1), order='F') 
omegaa_deriv  = pco2_deriv.reshape ((6,1), order='F') 
omegac_deriv  = pco2_deriv.reshape ((6,1), order='F')



## Call subroutine 'vars()' with one optionnal output parameter (pco2_deriv)
#mocsy.mvars.vars(ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis,temp,sal,
                    #alk, DIC,
                    #sil,phos,patm,depth,lat,n=1,optcon="mol/kg",optt='Tinsitu',optp='db',
                    #pco2_deriv=pco2_deriv )
## Call subroutine 'vars()' with all optionnal output parameters (xxx_deriv, yyy_deriv, ...)
#mocsy.mvars.vars(ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis,temp,sal,
                    #alk, DIC,
                    #sil,phos,patm,depth,lat,n=1,optcon="mol/kg",optt='Tinsitu',optp='db',
                    #ph_deriv=ph_deriv, pco2_deriv=pco2_deriv, fco2_deriv=fco2_deriv, co2_deriv=co2_deriv,
                    #hco3_deriv=hco3_deriv, co3_deriv=co3_deriv, omegaa_deriv=omegaa_deriv, omegac_deriv=omegac_deriv)

# Call subroutine 'vars()' with one optionnal output parameter (pco2_deriv)
mocsy.mvars.vars(ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis,temp=20,sal=35,
                    alk=2300e-6, dic=2000*1e-6,
                    sil=0,phos=0,patm=1,depth=0,lat=0,n=1,optcon="mol/kg",optt='Tinsitu',optp='db',
                    pco2_deriv=pco2_deriv )
# Result OK
print "pco2_deriv[1,0]", pco2_deriv[1,0]

# Call subroutine 'vars()' with all optionnal output parameters (xxx_deriv, yyy_deriv, ...)
mocsy.mvars.vars(ph,pco2,fco2,co2,hco3,co3,omegaa,omegac,betad,rhosw,p,tempis,temp=20,sal=35,
                    alk=2300e-6, dic=2000*1e-6,
                    sil=0,phos=0,patm=1,depth=0,lat=0,n=1,optcon="mol/kg",optt='Tinsitu',optp='db',
                    ph_deriv=ph_deriv, pco2_deriv=pco2_deriv, fco2_deriv=fco2_deriv, co2_deriv=co2_deriv,
                    hco3_deriv=hco3_deriv, co3_deriv=co3_deriv, omegaa_deriv=omegaa_deriv, omegac_deriv=omegac_deriv)

# Result for pco2_deriv is wrong
print "pco2_deriv[1,0]", pco2_deriv[1,0]

