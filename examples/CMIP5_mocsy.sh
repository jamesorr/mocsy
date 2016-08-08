#!/bin/bash

#==================================================
#
# Usage: CMIP5_mocsy.sh Institute Model Experiment Frequency Ensemble Version 
#
# Ex: ./CMIP5_mocsy.sh IPSL IPSL-CM5A-LR historical yr r1i1p1 latest 
#
#==================================================
# Original by James Orr, LSCE/CEA-IPSL, 24 Feb 2012
#==================================================


  Institute=$1
  Model=$2
  Experiment=$3
  Frequency=$4
  Ensemble=$5
  Version=$6
# Variable=$7

# Seed variable (upon which filenames will be searched and for which files for other vars will be built)
  Variable=dissic

  shopt -s expand_aliases
# alias cdo="/home/brocksce/cdo/bin/cdo"
# WORKDIR=/data/jomce/work ; export WORKDIR

#  files=( `grep -E ".*/*CMIP5/derived_CMIP5/${Institute}/${Model}/${Experiment}/${Frequency}/.*/.*/${Ensemble}/${Version}/${Variable}/.*\.nc" /home/brocksce/lists/prodigfs.list | sort -t'/' -k12.2n` )
files=(`find /prodigfs/OCMIP5/derived_CMIP5/${Institute}/${Model}/${Experiment}/${Frequency} -name "*.nc" | grep ${Ensemble} | grep ${Version} | grep ${Variable} | sort -t'/' -k12.2n`)

  echo "files: " $files

# Variable names: (default names of subdirectories and variables in input files)
# ------------------------------------------------------------------------------
  alkv="talk"
  dicv="dissic"
  po4v="po4"
  tempv="thetao"
  saltv="so"
  sio2v="si"

  echo "Variable names: " $alkv, $dicv, $tempv, $saltv, $po4v, $sio2v

  #Define function to seek out files with slightly different temporal subset (end year), if necessary
  # - use in next for loop
  function searchfile ()
  {
       dirfile=${1} ; var=${2}
       echo "dirfile: " $dirfile  
       echo "var: " $var
       #Get reference start & end years (dissic file)
       let nfile=$nfile+1
       echo "nfile: " $nfile
       temporalsubsetref=`basename ${dirfiledic} | awk -F'[_|.]' '{print $6}'`
       timebegref=`echo ${temporalsubsetref} | awk -F'-' '{print $1}'`
       timeendref=`echo ${temporalsubsetref} | awk -F'-' '{print $2}'`
       if [ $nfile -eq 1 ] ; then
          timeend=$timeendref
          timeendmin=$timeendref
          timebeg=$timebegref
          echo "timebeg, timeendref, timeendmin: " $timebeg $timeendref $timeendmin
       fi
       #Test to see if file exists
       #ls $dirfile ; err=$?
       #if [ $err -eq 0 ] ; then
       if [ -e $dirfile ] ; then
         #Exit status OK: proceed with default file (it exists)
         dirfilenew=$dirfile
       else
         #Exit status BAD: search for file, with perhaps different year range
         # echo "Error code:" $err
         #Get name of Ensemble directory by going up directory tree 2 steps
         dirvarname=`dirname $dirfile`
         dirvsname=`dirname $dirvarname`
         dirensname=`dirname $dirvsname`
         #Find any netCDF file under the ensemble directory for same variable
         #
         dirfiles=`find $dirensname -wholename "*$var\/$var_*.nc"`
         numfiles=`echo $dirfilenew | wc -w`
         echo "Files: " $numfiles
         dirfilenew=`echo $dirfiles | awk '{print $1}'`
         echo "dirfilenew: " $dirfilenew
         temporalsubset=`basename ${dirfilenew} | awk -F'[_|.]' '{print $6}'`
         timebeg=`echo ${temporalsubset} | awk -F'-' '{print $1}'`
         timeend=`echo ${temporalsubset} | awk -F'-' '{print $2}'`
         if [ $timeend -lt $timeendmin ] ; then
            timeendmin=$timeend
            echo "timeendmin: " $timeendmin
         fi
       fi
  }

  for file in ${files[*]} ; do
      echo "file:"$n ${file}
      dirfiledic=$file
      dirdic=`dirname ${file}`
      diralk=`echo ${dirdic}  | sed -e "s/${dicv}/${alkv}/g"`
      dirtemp=`echo ${dirdic} | sed -e "s/${dicv}/${tempv}/g"  -e "s/ocnBgchem/ocean/g"`
      dirsalt=`echo ${dirdic} | sed -e "s/${dicv}/${saltv}/g"  -e "s/ocnBgchem/ocean/g"`
      dirpo4=`echo ${dirdic}  | sed -e "s/${dicv}/${po4v}/g"`
      dirsio2=`echo ${dirdic} | sed -e "s/${dicv}/${sio2v}/g"`

      filedic=`basename ${file}`
      filealk=`echo ${filedic}  | sed -e "s/${dicv}/${alkv}/g"`
      filetemp=`echo ${filedic} | sed -e "s/${dicv}/${tempv}/g"`
      filesalt=`echo ${filedic} | sed -e "s/${dicv}/${saltv}/g"`
      filepo4=`echo ${filedic}  | sed -e "s/${dicv}/${po4v}/g"`
      filesio2=`echo ${filedic} | sed -e "s/${dicv}/${sio2v}/g"`

      #dirfiledic=$dirdic/$filedic
      dirfilealk=$diralk/$filealk
      dirfiletemp=$dirtemp/$filetemp
      dirfilesalt=$dirsalt/$filesalt
      dirfilepo4=$dirpo4/$filepo4
      dirfilesio2=$dirsio2/$filesio2
      # For now, assume no models carry PO4 and SiO2 
      # Perhaps the right way to go for 2 reasons:
      # 1) better comparison since most models don't have PO4; also many don't have SiO2
      # 2) models such as IPSL, MOHC, and MPI-M carry only "practical alkalinity", not total alkalinity
      # In regards to (2), it may still be useful to make calcs with PO4 & SiO2 for
      #    - improved data-model comparison?
      #    - to study future trends in acidification due to changes in these nutrients
      dirfilepo4=None
      dirfilesio2=None
      po4v=None
      sio2v=None

      #Compute carbonate chemistry
      if [[ -e ${dirfiledic} &&  -e ${dirfilealk} &&  -e ${dirfiletemp} &&  -e ${dirfilesalt} ]] ; then
        echo "All 4 standard files present: proceed to compute carbonate chemistry"
      else
        echo "Missing files: (at least one of following 4 files)"
        #ls -l ${dirdic}/${filedic} ${diralk}/${filealk} ${dirtemp}/${filetemp} ${dirsalt}/${filesalt}
        ls -l ${dirfiledic}    #This file exists by definition, but the 3 others may not
        ls -l ${dirfilealk}
        ls -l ${dirfiletemp}
        ls -l ${dirfilesalt}
      fi

      #Get times of start & end, and replace default files if they don't exist 
      nfile=0                                               #Counter for searchfile function
      searchfile $dirfiledic $dicv ; dirfiledic=$dirfilenew
      searchfile $dirfilealk $alkv ; dirfilealk=$dirfilenew
      searchfile $dirfiletemp $tempv ; dirfiletemp=$dirfilenew
      searchfile $dirfilesalt $saltv ; dirfilesalt=$dirfilenew

      echo "timebeg, timeendref, timeendmin: " $timebeg $timeendref $timeendmin
      \rm cmd_dirfile
      echo cdat ./mocsy_CMIP5.py ${dirfiledic} ${dicv} ${dirfilealk} ${alkv} ${dirfiletemp} ${tempv} ${dirfilesalt} ${saltv} ${dirfilepo4} ${po4v} ${dirfilesio2}  ${sio2v} ${timebeg} ${timeendref} ${timeendmin}
      echo cdat ./mocsy_CMIP5.py ${dirfiledic} ${dicv} ${dirfilealk} ${alkv} ${dirfiletemp} ${tempv} ${dirfilesalt} ${saltv} ${dirfilepo4} ${po4v} ${dirfilesio2}  ${sio2v} ${timebeg} ${timeendref} ${timeendmin} >> cmd_dirfile
      cdat ./mocsy_CMIP5.py ${dirfiledic} ${dicv} ${dirfilealk} ${alkv} ${dirfiletemp} ${tempv} ${dirfilesalt} ${saltv} ${dirfilepo4} ${po4v} ${dirfilesio2}  ${sio2v} ${timebeg} ${timeendref} ${timeendmin} 
  done


exit 0
