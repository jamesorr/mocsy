# !/usr/bin/env python                                                    
#                                                                         
"""@package docstring
 \file mocsy_GLODAP.py
  \BRIEF Python routine computing 3-D global gridded carb system vars with GLODAP & WOA2009 input
"""
#                                                                         
# <*>=                                                                    
# <Module Documentation>=                                                 

"""
    From observed  DIC and Alkalinity, Temperature and Salinity, and Silica and Phosphate:
    compute the following variables:
        - pH (Total scale),
        - ocean partial pressure of CO2 (uatm),
        - CO2 fugacity (uatm),
        - CO2*, HCO3- and CO3-- concentrations (mol/m^3),
        - Omega-C and Omega-A, saturation state of calcite and aragonite,
        - BetaD, homogeneous buffer factor (a.k.a. the Revelle Factor)
    Save them in a NetCDF file.

    This python script calls a library (shared object file built with f2py) that exploits 
    fortran subroutines (mocsy.f90) that compute ocean carbonate chemistry variables.
"""

#=============================================================================
#=============================================================================

# <Import library modules>=                                               

import numpy as np
import cdms2 as cdms, MV2 as MV
#import cdutil
import sys, os
import types

#print "Imported libraries"

sys.path.append ('.') 
import mocsy

#=============================================================================
#=============================================================================

# <Constant definitions>=                                                 

# Define aliases name for all dimension types: time, level, longitude, latitude
# cdms.axis.time_aliases.extend (('t', 'time', 'ttime', 'taxis', 'tt'))
cdms.axis.time_aliases.extend (('t', 'ttime', 'taxis', 'tt', 'TIME_COUNTER'))
cdms.axis.level_aliases.extend (('z', 'depth', 'zdepth', 'zaxis', 'zz', 'zlevel', 'depth2','DEPTHT'))
cdms.axis.longitude_aliases.extend (('x', 'xlon', 'xaxis', 'xx','X'))
cdms.axis.latitude_aliases.extend (('y', 'ylat', 'yaxis', 'yy','Y'))
#print "Defined aliases for axes"

#=============================================================================
#=============================================================================

# <Sub-routines of this module>=                                          
# None

#=============================================================================
#=============================================================================

# <Main program>=                                                         

#=============================================================================
#=============================================================================

#print "Before usage documentation"
if __name__ == "__main__":

    # <Read command line parameters>=                                         

    usage = """
        Usage: python mocsy_GLODAP.py [preindustrial_option]

        preindus-option:   Option for computing from preindustrial DIC or total DIC
                            0 = total (1994), i.e., modern era, this is the default
                            1 = preindustrial (1765), after removing anthropogenic DIC from total DIC

        From observed DIC & Alk (GLODAP), T, S, PO4, and SiO2 (WOA2009)
        compute following 3D variables:
            - pH,
            - CO2*, HCO3- and CO3-- concentrations,
            - OmegaA & OmegaC (saturation states for aragonite and calcite),
            - Revelle factor (BetaD or homogeneous buffer factor)
        Save them in appropriate netCDF files
        preindustrial_option = 0 
    """ + __doc__

    nb_arguments = len (sys.argv)
    if nb_arguments in (1, 2):
        if nb_arguments == 2:
            if sys.argv[1] in ("0", "1"):
                preindustrial_option = int (sys.argv[1])
            else:
                print usage
                sys.exit (1)
    else:
        print usage
        sys.exit (1)


    #-----------------------------------------------------------------------------

#   Data directories and filenames: whereabouts and names 
#   Local directory:
#   OCMIP5_DAT='/prodigfs/OCMIP5/DATA'
#   For remote access:
    OCMIP5_DAT='http://dods.ipsl.jussieu.fr/cgi-bin/nph-dods/ocmip/phase5/DATA'

    GLODAPsubdir = 'gridded/glodap'
    WOA2009subdir = 'gridded/WOA2009'

    DIC_dsubdir = GLODAPsubdir
    DIC_dfile = 'TCO2.nc'
    DIC_dvar = 'TCO2'

    Alk_dsubdir = GLODAPsubdir
    Alk_dfile = 'TALK.nc'
    Alk_dvar = 'TALK'

#   In situ temperature
    Temp_dsubdir = WOA2009subdir
    Temp_dfile = "temperature_annual_1deg.nc"
    Temp_dvar = 't_an'

    Salt_dsubdir = WOA2009subdir
    Salt_dfile = "salinity_annual_1deg.nc"
    Salt_dvar = 's_an' 

    PO4_dsubdir = WOA2009subdir
    PO4_dfile = "phosphate_annual_1deg.nc"
    PO4_dvar = 'p_an'

    Si_dsubdir = WOA2009subdir
    Si_dfile = "silicate_annual_1deg.nc"
    Si_dvar = 'i_an'

    DIC_datfile  = OCMIP5_DAT  + "/" + DIC_dsubdir  + "/" + DIC_dfile
    Alk_datfile  = OCMIP5_DAT  + "/" + Alk_dsubdir  + "/" + Alk_dfile
    Temp_datfile = OCMIP5_DAT  + "/" + Temp_dsubdir + "/" + Temp_dfile
    Salt_datfile = OCMIP5_DAT  + "/" + Salt_dsubdir + "/" + Salt_dfile
    PO4_datfile  = OCMIP5_DAT  + "/" + PO4_dsubdir  + "/" + PO4_dfile
    Si_datfile   = OCMIP5_DAT  + "/" + Si_dsubdir   + "/" + Si_dfile

#   Get Depth and Latitude information from GRID file
#   These are "reshaped" later so that their shapes conform with other input variables
    j = cdms.open (DIC_datfile)
    dicdat = j(DIC_dvar)
    depth1D = dicdat.getLevel().getValue()
    lat2D   = dicdat.getLatitude().getValue()

    f = cdms.open (DIC_datfile)
    dicd = f(DIC_dvar)
    dic_axes=dicd.getAxisList()        #Stores axes for later use (in all output variables)
    dic_grid = dicd.getGrid()

    if preindustrial_option:
        # Compute preindustrial DIC
        aDIC_dsubdir = 'gridded/glodap'
        aDIC_dfile = 'ANTCO2.nc'
        aDIC_dvar = 'ANTCO2'
        aDIC_datfile  = OCMIP5_DAT  + "/" + aDIC_dsubdir  + "/" + aDIC_dfile
        h = cdms.open (aDIC_datfile)
        anth_dic = h(aDIC_dvar)
        dicd -= anth_dic
        del(anth_dic)
        h.close()

    dicd = dicd.filled(1e+20)
    dicd = cdms.MV2.masked_equal(dicd, dicd[0,0,0])
    f.close()
    
    f = cdms.open (Alk_datfile)
    alkd = f(Alk_dvar)
    alkd = alkd.filled(1e+20)
    alkd = cdms.MV2.masked_equal(alkd, alkd[0,0,0])
    f.close()

    f = cdms.open (Temp_datfile)
    td = f(Temp_dvar)
    td = td.filled(1e+20)
    td = cdms.MV2.masked_equal(td, td[0,0,0])
    f.close()

    f = cdms.open (Salt_datfile)
    sd = f(Salt_dvar)
    sd = sd.filled(1e+20)
    sd = cdms.MV2.masked_equal(sd, sd[0,0,0])
    f.close()

    f = cdms.open (PO4_datfile)
    po4d = f(PO4_dvar)
    po4d = po4d.filled(1e+20)
    po4d = cdms.MV2.masked_equal(po4d, po4d[0,0,0])
    f.close()

    f = cdms.open (Si_datfile)
    sio2d = f(Si_dvar)
    sio2d = sio2d.filled(1e+20)
    sio2d = cdms.MV2.masked_equal(sio2d, sio2d[0,0,0])
    f.close()

    print "Completed reading input files"

    #-----------------------------------------------------------------------------
    depth = depth1D
    lat = lat2D
    #
    #-----------------------------------------------------------------------------
    # Read grid dimensions from file variable (which has a dummy time dimension)
    level_nb , lat_nb, lon_nb = dicd.shape
    #
    #-----------------------------------------------------------------------------

    # Adjust shape and size of 2 variables (from GRID file): depth, lat
    depth = np.resize (depth, [lon_nb, lat_nb, level_nb])
    lat   = np.resize (lat,   [lon_nb, level_nb, lat_nb])
    depth = np.transpose (depth)
    lat   = np.transpose (lat,axes=(1,2,0))
    # explanation for tupple (1,2,0)   
    # 1: put current 2nd dim (level) in first position
    # 2: put current 3rd dim (lat) in  secodn position
    # 0: put current 1st dim (lon) in   third position
    # Make a flat array of each of the same 2 variables:

    depth = np.ravel (depth)
    lat   = np.ravel (lat)


    #-----------------------------------------------------------------------------

    # Convert input data to flat arrays
    # Mask DIC where any of the input variables is masked
    common_mask = np.zeros(dicd.shape)
    for v in [dicd, alkd, sd, td, sio2d, po4d]:
        v = cdms.MV2.reshape(v, dicd.shape)
    	m=v.mask
        #print "var shape = ", v.shape
    	if (m is not None):
        		common_mask=np.logical_or(common_mask,m)
    dicd = cdms.MV2.masked_where(common_mask,dicd)
    #
    # Make a flat array of 3D variables: t, s, alk, sio2 & po4
    td = np.ravel (td.getValue())
    sd = np.ravel (sd.getValue())
    alkd = np.ravel (alkd.getValue())
    po4d = np.ravel (po4d.getValue())
    sio2d = np.ravel (sio2d.getValue())

    # Convert units
    # Alk: from ueq/kg to mol/kg
    alkd = alkd * 1.0e-6
    # DIC: from umol/kg to mol/kg
    dicd = dicd * 1.0e-6
    # Silica and Phosphate: convert from umol/L to mol/kg
    pfd = mocsy.mdepth2press(depth,lat)    #Compute pressure (db) from depth (m) and latitude
    rhois = mocsy.mrhoinsitu(sd, td, pfd)  #Compute in situ density (kg/m3)
    sio2d = sio2d[:] / (rhois * 1.0e+3)
    po4d  = po4d[:] / (rhois * 1.0e+3)

    # Atmospheric pressure (1 atm everywhere)
    patmd = td[:]*0 + 1

    # Flatten array, missing values filled with 0
    # dic_array = np.ravel (dicd.filled(0.))
    dicd = dicd.filled(0.)
    #
    #-----------------------------------------------------------------------------
    # Compute other carbonate system variables (pH, CO2*, HCO3-, CO3--, Omegas, R (BetaD), (& in-situ rho, pressure, & T)
    print "Compute carbonate chemistry"
    pH, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, DENis, p, Tis = (
    mocsy.mvars (td, sd, alkd, dicd.ravel(), sio2d, po4d, patmd,  depth, lat,
               optcon='mol/kg', optt='Tinsitu', optp='m', optb="u74", optk1k2='l', optkf="dg", optgas='Pinsitu')
    )
    del (td)
    del (sd)
    # del (alkd)
    # del (dic_array)
    del (sio2d)
    del (po4d)
    #del (p)
    del (depth)
    del (lat)
    #
    #-----------------------------------------------------------------------------
    # Make them variables of correct shape
    pH   = np.reshape (pH, dicd.shape)
    pco2  = np.reshape (pco2, dicd.shape)
    fco2  = np.reshape (fco2, dicd.shape)
    co2  = np.reshape (co2, dicd.shape)
    hco3 = np.reshape (hco3, dicd.shape)
    co3  = np.reshape (co3, dicd.shape)
    OmegaA = np.reshape (OmegaA, dicd.shape)
    OmegaC = np.reshape (OmegaC, dicd.shape)
    BetaD  = np.reshape (BetaD, dicd.shape)
    Alk =  np.reshape (alkd, dicd.shape)
    DIC =  np.reshape (dicd.ravel(), dicd.shape)
    #
    DENis   = np.reshape (DENis, dicd.shape)
    p       = np.reshape (p, dicd.shape)
    Tis  = np.reshape (Tis, dicd.shape)

    print 'Convert CO2, HCO3, CO3, Alk, & DIC from mol/kg to umol/kg (standard data units)'
    co2 = co2 * 1e+6
    co3 = co3 * 1e+6
    hco3 = hco3 * 1e+6
    Alk = Alk * 1e+6
    DIC = DIC * 1e+6
    #
    #-----------------------------------------------------------------------------
    # Caution: fill_value parameter is mandatory !!!
    pH   = cdms.createVariable (pH, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    pco2  = cdms.createVariable (pco2, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    fco2  = cdms.createVariable (fco2, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    co2  = cdms.createVariable (co2, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    hco3 = cdms.createVariable (hco3, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    co3  = cdms.createVariable (co3, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    OmegaA = cdms.createVariable (OmegaA, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    OmegaC = cdms.createVariable (OmegaC, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    BetaD  = cdms.createVariable (BetaD, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    Alk  = cdms.createVariable (Alk, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    DIC  = cdms.createVariable (DIC, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    DENis  = cdms.createVariable (DENis, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    p      = cdms.createVariable (p, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    Tis  = cdms.createVariable (Tis, mask=common_mask, fill_value=1e20, axes=dic_axes, grid=dic_grid)
    #
    #-----------------------------------------------------------------------------
    # Save them under correct name
    pH.id = "pH"
    pH.long_name = "pH"
    pH.units = "pH units"
    del(pH.name)
    #    
    pco2.id = "pCO2"
    pco2.long_name = "partial pressure of CO2"
    pco2.units = "uatm"
    del(pco2.name)
 	#
    fco2.id = "fCO2"
    fco2.long_name = "fugacity of CO2"
    fco2.units = "uatm"
    del(fco2.name)
    #
    co2.id = "CO2"
    co2.long_name = "CO2* concentration"
    co2.units = "umol kg^-1"
    del(co2.name)
    #
    hco3.id = "HCO3"
    hco3.long_name = "Bicarbonate ion concentration"
    hco3.units = "umol kg^-1"
    del(hco3.name)
    #
    co3.id = "CO3"
    co3.long_name = "Carbonate ion concentration"
    co3.units = "umol kg^-1"
    del(co3.name)
    #
    OmegaA.id = "OmegaA"
    OmegaA.long_name = "Saturation state of Aragonite"
    OmegaA.units = "None"
    del(OmegaA.name)
    #
    OmegaC.id = "OmegaC"
    OmegaC.long_name = "Saturation state of Calcite"
    OmegaC.units = "None"
    del(OmegaC.name)
    #
    BetaD.id = "BetaD"
    BetaD.long_name = "Revelle factor (Homogeneous buffer factor)"
    BetaD.units = "None"
    del(BetaD.name)
    #
    DENis.id = "DENis"
    DENis.long_name = "in situ Density"
    DENis.units = "kg m^-3"
    del(DENis.name)
    #
    p.id = "p"
    p.long_name = "Pressure"
    p.units = "db"
    del(p.name)
    #
    Tis.id = "Tis"
    Tis.long_name = "in situ Temperature"
    Tis.units = "C"
    del(Tis.name)
    #
    DIC.id = "DIC"
    DIC.long_name = "Dissolved Inorganic Carbon"
    DIC.units = "umol kg^-1"
    del(DIC.name)
    #
    Alk.id = "Alk"
    Alk.long_name = "Total Alkalinity"
    Alk.units = "ueq kg^-1"
    del(Alk.name)
    #
    # dic.id = "DIC"
    # dic.long_name = "Dissolved Inorganic Carbon"
    # dic.units = "mol/m^3"
    #
    #-----------------------------------------------------------------------------
    # Save results in output file
    dirnameo = '/prodigfs/OCMIP5/DATA/gridded/glodap/derived2'
    # dirnameo = '.'
    #version=os.system("date +'v20%y%m%d'")
    version=os.popen("date +'v20%y%m%d'").read()[0:9]
    if preindustrial_option:
        filenameo = "GLODAP_pH_Carb_1765_3Da_" + str(version) + ".nc"
    else:
        filenameo = "GLODAP_pH_Carb_1994_3Da_" + str(version) + ".nc"
    #
    dirfilenameo = dirnameo + "/" + filenameo
    dirfileo = cdms.open (dirfilenameo, "w")

    for variable in (pH, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, DENis, Alk, DIC):
    	dirfileo.write (variable)
    	print '  *',variable.id, ':', dirfilenameo

    dirfileo.close()
    j.close()
