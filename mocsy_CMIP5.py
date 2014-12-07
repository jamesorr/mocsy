# Adapted for numpy/ma/cdms2 by convertcdms.py
# !/usr/bin/env python                                                    
#                                                                         
#                                                                         
# <*>=                                                                    
# <Module Documentation>=                                                 

"""
    From modelled DIC and Alkalinity, Temperature and Salinity, and Silica and Phosphate:
    compute the following variables:
        - pH (Total scale),
        - ocean partial pressure of CO2 (uatm),
        - CO2 fugacity (uatm),
        - CO2*, HCO3- and CO3-- concentrations (mol/m^3),
        - Omega-C and Omega-A, saturation state of calcite and aragonite,
        - BetaD, homogeneous buffer factor (a.k.a. the Revelle Factor)
    Save them in a NetCDF file.

    This python script calls a Fortran subroutine (mocsy), with recent developments by J.C. Orr.
    It was derived from original code used by O. Aumont for carb. chem. in OPA-PISCES model,
    part of which was based on code from E. Maier-Reimer (HAMMOC3); other parts were taken from OCMIP2.
    That code was first extended and corrected by J.-M. Epitalon & J. Orr, and also 
    made consistent with J-P Gattuso's Seacarb R software (in 2004). This effort provided bug
    fixes to seacarb (sent to J.P. Gattuso who implemented them in the
    subsequent seacarb version).  Originally, this python-fortran code was
    used for the OCMIP-2 & GLODAP analysis described in (Orr et al., 2005,
    Nature). It has subsequently undergone improvements and functionalities added
    by J. Orr.  Extensive tests 2012 and 2013 show it produces
    results that are essentially identical to those from CO2sys and seacarb 
    (within round-off error).
"""

#=============================================================================
#=============================================================================

# <Import library modules>=                                               

import numpy as np
import cdms2 as cdms, MV2 as MV
#import cdutil
import sys, os
import types

sys.path.append ('.')    
import mocsy

#=============================================================================
#=============================================================================

# <Constant definitions>=                                                 

# Define aliases name for all dimension types: time, level, longitude, latitude
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
        Usage: mocsy_OCMIP5.py <DIC-file> <DIC-var> <Alk-file> <Alk-var> <Temp-file> <Temp-var> \\
                           <Salt-file> <Salt-var> <PO4-file> <PO4-var> <Si-file> <Si-var> 
                           
        Ouput file is created & stored following OCMIP5 (CMIP5) archive nomenclature 
        (under /prodigfsOCMIP5/derived_CMIP5). Its name is derived from DIC filename:
         * dir as for input dir but "dissic" replaed by output var names
         * filename handled likewise with "dissic" replaced by output var names"
    """ + __doc__

    nb_arguments = len (sys.argv)
    if nb_arguments == 16:
        DIC_file  = sys.argv[1]
        DIC_var   = sys.argv[2]
        Alk_file  = sys.argv[3]
        Alk_var   = sys.argv[4]
        Temp_file = sys.argv[5]
        Temp_var  = sys.argv[6]
        Salt_file = sys.argv[7]
        Salt_var  = sys.argv[8]
        PO4_file  = sys.argv[9]
        PO4_var   = sys.argv[10]
        Si_file   = sys.argv[11]
        Si_var    = sys.argv[12]
        stimebeg   = sys.argv[13]
        stimeendref = sys.argv[14]
        stimeendmin = sys.argv[15]
        print "Number arguments: ", nb_arguments, " OK!"
    else:
        print usage
        print "Number arguments: ", nb_arguments, " != 16"
        sys.exit (1)

#   Convert input arguments from string to integer
    timebeg=int(stimebeg)
    timeendref=int(stimeendref)
    timeendmin=int(stimeendmin)

#   Separate filename (only) and dirname (only) from DIC_file (which contains both)
    g = cdms.open (DIC_file)
    dicfile = os.path.basename(DIC_file)
    dicdir  = os.path.dirname(DIC_file)
    outdir  = dicdir

#   Extract model name from filename
    filnameparts=dicfile.split('_')
    model=filnameparts[2]
    expmt=filnameparts[3]
    dirparts=dicdir.split('/')
    institute=dirparts[4]
#   g.close()

#   Get Depth and Latitude information from GRID file
#   These are "reshaped" later so that their shapes conform with other input variables
#   dicall = g(DIC_var)        #whole array
#   depth1D = dicall.getLevel().getValue()
#   lat2D   = dicall.getLatitude().getValue()
    dic0 = g(DIC_var,time=slice(0,1))
    depth1D = dic0.getLevel().getValue()
    lat2D   = dic0.getLatitude().getValue()


#   Now start long sequence of 2 imbedded for loops: 1) increments timestep and 2) increments output variable
#   ---------------------------------------------------------------------------------------------------------
#   For loop to treat results 1 year at a time 
#   (N.B. cannot treat all years at once: results in Segmentation fault)
 
#   years=xrange(1860,2100)
#   yearbeg=years[0]
#   yearend=years[-1]
#   for year in years:
#	print 'Year:', year

#   Get time axis & its length, starting value, and ending value
    tax = g.getAxis('time')
    taxlen = len(tax)
    taxbeg = tax[0]
    taxend = tax[taxlen-1]

#   Get corresponding time indices
#   itims = xrange(0,1)
#   itims = xrange(taxlen)
    itims = xrange(timeendmin+1-timebeg)
#   itimbeg = itims[0]
#   itimend = itims[-1]
#   g.close()

    for itim in itims:
	print 'Index, timestep:', itim, tax[itim]
    	f = cdms.open (DIC_file)
    	dic = f(DIC_var,time=slice(itim,itim+1))
        dic_axes=dic.getAxisList()        #Stores axes for later use (in all output variables)
        dic_grid = dic.getGrid()
	dic=dic.filled(1e+20)
        dic = cdms.MV2.masked_equal(dic,dic[0,0,0,0])
    	f.close()
    	#
    	f = cdms.open (Alk_file)
    	alk = f(Alk_var,time=slice(itim,itim+1))
	alk=alk.filled(1e+20)
        alk = cdms.MV2.masked_equal(alk,alk[0,0,0,0])
    	f.close()
    	#
    	f = cdms.open (Temp_file)
    	t = f(Temp_var,time=slice(itim,itim+1))
        t = t - 273.15    #Convert from K to C
	t=t.filled(1e+20)
        t = cdms.MV2.masked_equal(t,t[0,0,0,0])
    	f.close()
    	#
    	f = cdms.open (Salt_file)
    	s = f(Salt_var,time=slice(itim,itim+1))
	s=s.filled(1e+20)
        s = cdms.MV2.masked_equal(s,s[0,0,0,0])
    	f.close()
    	#
        # c HadCM3LC model does not have PO4 file; COSMOS HISTA2U run has no available PO4 file
        po4 = t * 0.
	sio2 = t * 0.
        # Correct for error in MIROC output (Alk is too low by a factor of two!) - still the case on 11 June 2012
	if  (institute == 'MIROC'):
            alk = 2.*alk
	#if  (model == 'HadCM3LC'):
	#	po4 = t * 0.
	#elif  (model == 'COSMOS') & (expmt == 'HISTA2U'):
	#	po4 = t * 0.
	#else:
	#        f = cdms.open (PO4_file)
    	#        po4 = f(PO4_var,time=slice(itim,itim+1))
	#        po4=po4.filled(1e+20)
        #        po4 = cdms.MV2.masked_equal(po4,po4[0,0,0,0])
        #        f.close()
    	#
	# 5 models do not have SiO2 (only IPSL-CM4 and CCSM3 have this varible)
	#if (model == 'BCM-C') | (model == 'CSM1-4') | (model == 'COSMOS') | (model == 'HadCM3LC') | (model == 'UVIC2-8'):
	#	sio2 = t * 0.
	#else:
   	#	f = cdms.open (Si_file)
	#	sio2 = f(Si_var,time=(tax[itim],tax[itim]))
	#	sio2=sio2.filled(1e+20)
        #        sio2 = cdms.MV2.masked_equal(sio2,sio2[0,0,0,0])
    	#	f.close()
	#
    	#-----------------------------------------------------------------------------
    	depth = depth1D
        lat = lat2D
	#
    	#-----------------------------------------------------------------------------
    	# Read grid dimensions from file variable (which has a dummy time dimension)
    	time_nb, level_nb , lat_nb, lon_nb = dic.shape
	#
    	#-----------------------------------------------------------------------------
    	# Convert input data to flat arrays
    	# Mask DIC where any of the input variables is masked
    	# common_mask = reduce (np.logical_or, (dic.mask, alk.mask, s.mask, t.mask, sio2.mask, po4.mask ))
    	# dic = cdms.MV2.array (dic, mask=common_mask)
    	common_mask = np.zeros(dic.shape)
    	for v in [dic, alk, s, t, sio2, po4]:
        	m=v.mask
        	if (m is not None):
            		common_mask=np.logical_or(common_mask,m)
    	dic = cdms.MV2.masked_where(common_mask,dic)
	#
     	# Adjust shape and size of 2 variables (from GRID file): depth, lat
    	depth = np.resize (depth, [lon_nb, lat_nb, level_nb, time_nb])
    	lat   = np.resize (lat,   [lon_nb, lat_nb, level_nb, time_nb])
    	depth = np.transpose (depth)
    	lat   = np.transpose (lat)
    	# Make a flat array of each of the same 2 variables:
    	depth = np.ravel (depth)
    	lat   = np.ravel (lat)
    	# Make a flat array of 3D variables: t, s, alk, sio2 & po4
    	t = np.ravel (t.getValue())
    	s = np.ravel (s.getValue())
    	alk = np.ravel (alk.getValue())
    	po4 = np.ravel (po4.getValue())
    	sio2 = np.ravel (sio2.getValue())
    	# Flatten array, missing values filled with 0
    	# dic_array = np.ravel (dic.filled(0.))
	dic = dic.filled(0.)
        patm = t[:]*0 + 1
	#
    	#-----------------------------------------------------------------------------
        # Computed other carbonate system variables (pH, CO2*, HCO3-, CO3--, Omegas, R (BetaD), (& in-situ rho, pressure, & T)
        pH, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, DENis, p, Tis = (
        mocsy.mvars (t, s, alk, dic.ravel(), sio2, po4, patm, depth, lat,
                     optcon='mol/m3', optt='Tpot', optp='m', optb="u74", optk1k2='l', optkf="dg", optgas="Pinsitu")
        )
    	del (t)
    	del (s)
    	# del (alk)
    	# del (dic_array)
    	del (sio2)
    	del (po4)
    	#del (p)
        del (depth)
        del (lat)
    	del (patm)
	#
    	#-----------------------------------------------------------------------------
    	# Make the variables have correct shape
    	pH   = np.reshape (pH, dic.shape)
   	pco2  = np.reshape (pco2, dic.shape)
    	fco2  = np.reshape (fco2, dic.shape)
    	co2  = np.reshape (co2, dic.shape)
    	hco3 = np.reshape (hco3, dic.shape)
    	co3  = np.reshape (co3, dic.shape)
    	OmegaA = np.reshape (OmegaA, dic.shape)
    	OmegaC = np.reshape (OmegaC, dic.shape)
    	BetaD  = np.reshape (BetaD, dic.shape)
    	Alk =  np.reshape (alk, dic.shape)
    	# DIC =  np.reshape (dic_array, dic.shape)
    	DIC =  np.reshape (dic.ravel(), dic.shape)
        #
    	DENis   = np.reshape (DENis, dic.shape)
    	p       = np.reshape (p, dic.shape)
    	Tis  = np.reshape (Tis, dic.shape)
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
    	fco2.long_name = "CO2 fugacity"
    	fco2.units = "uatm"
    	del(fco2.name)
    	#
    	co2.id = "CO2"
    	co2.long_name = "CO2* concentration"
    	co2.units = "mol/m^3"
    	del(co2.name)
    	#
    	hco3.id = "HCO3"
    	hco3.long_name = "Bicarbonate ion concentration"
    	hco3.units = "mol/m^3"
    	del(hco3.name)
    	#
    	co3.id = "CO3"
    	co3.long_name = "Carbonate ion concentration"
    	co3.units = "mol/m^3"
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
    	BetaD.long_name = "Homogeneous buffer factor (Revelle factor)"
    	BetaD.units = "None"
    	del(BetaD.name)
    	#
    	DENis.id = "DENis"
    	DENis.long_name = "in-situ Density"
    	DENis.units = "kg/m^3"
    	del(DENis.name)
    	#
    	p.id = "p"
    	p.long_name = "Pressure"
    	p.units = "db"
    	del(p.name)
    	#
    	Tis.id = "Tis"
    	Tis.long_name = "in-situ Temperature"
    	Tis.units = "C"
    	del(Tis.name)
    	#
    	DIC.id = "DIC"
    	DIC.long_name = "Dissolved Inorganic Carbon"
    	DIC.units = "mol/m^3"
    	del(DIC.name)
    	#
    	Alk.id = "Alk"
    	Alk.long_name = "Total Alkalinity"
    	Alk.units = "eq/m^3"
    	del(Alk.name)
    	#
    	# dic.id = "DIC"
    	# dic.long_name = "Dissolved Inorganic Carbon"
    	# dic.units = "mol/m^3"
	#
    	#-----------------------------------------------------------------------------
    	# Save results in output file
	#
    	# Loop to write each variable into separate files w/ dir and filename consistent with input files
    	# for variable in (pH, co3):
    	for variable in (pH, pco2, fco2, co2, hco3, co3, OmegaA, OmegaC, BetaD, DENis):
        	varfile = dicfile.replace('dissic',variable.id)
        	varfile = varfile.replace(stimeendref,stimeendmin)
        	vardir  = outdir.replace('dissic',variable.id)
        	os.system("mkdir -p " + vardir)
        	output_file_name = vardir + "/" + varfile
		#if itim==itimbeg:
		if itim==0:
			# dict = {}                             #Create empty dictionary
			# dict[variable.id]=variable.id+'file'  #Add key:val to dictionary (for each variable)
        	        # output_file = cdms.open (dict[variable.id], "w")
        	        output_file = cdms.open (output_file_name, "w")
        	        #output_file = cdms.open (output_file_name, "a")
		else:
        	        output_file = cdms.open (output_file_name, "a")
		#
        	output_file.write (variable)
		print '  *',variable.id, ':',  output_file_name
        	output_file.close()

        #g.close()
