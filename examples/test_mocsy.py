# -*- coding: utf-8 -*-
import numpy as np
import mocsy

# Define input data (typical values at depth from 0 to 5000 meters)
temp = np.repeat(2.0, 6).astype('float32')
depth = np.arange (0, 6000, 1000).astype('float32')
sal = np.repeat(35.0, 6).astype('float32')
alk = np.repeat(2295.*1.e-6, 6).astype('float32')
dic = np.repeat(2154.*1.e-6, 6).astype('float32')
sil = phos = np.repeat(0.0, 6).astype('float32')
Patm = np.repeat(1.0, 6).astype('float32')
optK1K2 = 'l'

# Create output arrays
# --------------------

# computed variables at 6 input points
lat    = np.zeros((6,)).astype('float32')
ph     = np.zeros((6,)).astype('float32')
pco2   = np.zeros((6,)).astype('float32')
fco2   = np.zeros((6,)).astype('float32')
co2    = np.zeros((6,)).astype('float32')
hco3   = np.zeros((6,)).astype('float32')
co3    = np.zeros((6,)).astype('float32')
OmegaA = np.zeros((6,)).astype('float32')
OmegaC = np.zeros((6,)).astype('float32')
BetaD  = np.zeros((6,)).astype('float32')
rhoSW  = np.zeros((6,)).astype('float32')
p      = np.zeros((6,)).astype('float32')
tempis = np.zeros((6,)).astype('float32')
# values of derivatives w/ respect to 6 input variables and at 6 input points
ph_deriv     = np.zeros((6*6,)).astype('float32')
pco2_deriv   = np.zeros((6*6,)).astype('float32')
OmegaA_deriv = np.zeros((6*6,)).astype('float32')
ph_deriv     = ph_deriv.reshape ((6,6), order='F')
pco2_deriv   = pco2_deriv.reshape ((6,6), order='F')
OmegaA_deriv = OmegaA_deriv.reshape ((6,6), order='F')

# Run mocsy
mocsy.mvars.vars (ph, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, rhoSW, p, tempis, temp, sal, alk, dic, sil, phos, Patm, depth, lat, 6,
    optcon='mol/kg', optt='Tinsitu', optp='db', optb='l10', optk1k2=optK1K2, optkf='dg',
    ph_deriv=ph_deriv, pco2_deriv=pco2_deriv, omegaa_deriv=OmegaA_deriv )

# print mocsy results
# -------------------

print "pH     pCO2   fCO2     CO2*       HCO3-       CO32-      OmegaA OmegaC  R    Density Press  Temperature"
for i in range (0, 6):
    print ph[i], pco2[i], fco2[i], co2[i], hco3[i], co3[i], OmegaA[i], OmegaC[i], BetaD[i], rhoSW[i], p[i], tempis[i]

# Compute buffer factors from Egleston
# pco2_deriv[2,:] are derivatives of pCO2 w/ respect to DIC
# pco2_deriv[1,:] are        ...          w/ respect to Alk
gamma_DIC = pco2 / pco2_deriv[1,:]
gamma_Alk = pco2 / pco2_deriv[0,:]

beta_DIC  = -1. / (np.log(10.) * ph_deriv[1,:])
beta_Alk  = -1. / (np.log(10.) * ph_deriv[0,:])

# Here, we use Omega of Aragonite (use of Calcite would have been equaly valid)
omega_DIC = OmegaA / OmegaA_deriv[1,:]
omega_Alk = OmegaA / OmegaA_deriv[0,:]

print ""
print "gamma_DIC     gamma_Alk     beta_DIC     beta_Alk    omega_DIC    omega_Alk"
for i in range (0, 6):
    print gamma_DIC[i], gamma_Alk[i], beta_DIC[i], beta_Alk[i], omega_DIC[i], omega_Alk[i]

# Print derivatives of pH with respect to phosphate, silicate, temperature and salinity
print ""
print "dpH/dPhos  dpH/dSil  dpH/dT   dpH/dS"
print "pH/µMol     pH/µMol  pH/°C   pH/psu"
for i in range (0, 6):
    print ph_deriv[2,i], ph_deriv[3,i], ph_deriv[4,i], ph_deriv[5,i]
