'''
Graham Kerr
June 2022

O I in Flares

NAME: IRIS_SG_initial.py

PURPOSE: To perform some intial analysis on IRIS SG data, reading it into memory, 
         normalising exposure time, and performing quicklooks. 

INPUTS: None (probably best to open this up and run it line by line in the terminal/ipython)

OUTPUTS: Some .pkl files (these can get large, so maybe isolate the flare times before saving)

NOTES:

'''

###################
### SOME SET UP 
###################

## Import Modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# %matplotlib tk
import sys
import os
import re
from scipy import io
import pickle
import copy
from iris_lmsalpy import extract_irisL2data
import iris_lmsalpy
import pandas as pd

'''
Read in the data 
> - IRIS can operate in multiple modes: raster scanning (each raster has multiple slit positions) or sit-and-stare (each raster has one slit position).
> - If in sit-and-stare mode there is one raster file, which contains the time-series of the sit-and-stare observation.
> - If there are multiple slit positions then each raster file contains one set of slit positions, the next file contains the repeat scan through each slit position, and so-on.
> - In the latter case for a very long observation we might want to only grab some portion of the rasters (ie. when the flare is), which is why I have the option to grab the first and final rasters to study below. The default is all of them
'''
##
## The SG filenames
## 

dir1 = '/Users/gskerr1/Documents/Research/Melissa_OI_IRIS/2014_09_10_1130/'
file_search = r'iris_l2_20140910_112825_3860259453_raster_.*\.fits'

#os.chdir(dir1)
## Search for all raster files that start with "file_search"
rasterfiles = [f for f in os.listdir(dir1) if re.match(file_search, f)]
rasterfiles.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

# First raster file to study
rfi = 0
# rfi = 73 

# Final raster file to study
rff = len(rasterfiles)
# rff = 289 #17UT

rfinds = np.arange(rfi,rff)
nraster = (rfinds.shape)[0]

##
## > - Lets check what data we have in one of the raster files
##
## 
## Print the fits files info
##
tmp = extract_irisL2data.info_fits(dir1+rasterfiles[rfi])


## > - Extract the linelist info, and select the windows we want to study.
## > - For now I've selected Mg II NUV, Si IV and the O I window (that also has Fe XXI).
## > - The cell below grabs the index that corresponds to the lines of interest
##
## Show the line info and winids, and grab the ones we want
##
linelist = ['Mg II k 2796','O I 1356', 'Si IV 1403']
lines=extract_irisL2data.show_lines(dir1+rasterfiles[rfi])
lineid = np.zeros(len(linelist),dtype=int)
for i in range (len(linelist)):
    lineid[i] = (np.where(lines == linelist[i]))[0]

# > - Using those indices we create the raster object
# > - Note that for now we are only doing this for the first raster in our set so we can get some generic info... of course, if this is a sit-and-stare obs that is all we have.

########
## Lets have a quicklook at the object
########
%matplotlib tk
iris_sg_init.quick_look()
## run this to go back to plotting within the notebook
%matplotlib inline

# > - The data are saved as memory map objects, but if we want to modify them (e.g. normalise by exposure time etc.,) then we need to have them in the computer memory properly.
iris_sg_init.flush()

'''
Grab the slit positions, and headers
> - Below I have typed out one block for each window we are interested in here, but to grab the others in the observation too just copy, paste, and re-name.
> - There is probably a more convenient way to collate this info, but for now i'm doing it this way.
> - This will collate all the raster headers.. or in the case of sit-and-stare, just grab the single values. 
> - For sit-and-stare each variable will be [E, 1]... where E = number of exposures, ie. time or position on the Sun.
> - For multiple slit positions this will be [N, R]... where N = number of slit positons, and R = raster number (ie. time)
> - ... these definitons will get easier to understand the more you play with the data
'''

ind_mgii = (np.where(lines == 'Mg II k 2796'))[0][0]+1
ind_oi = (np.where(lines == 'O I 1356'))[0][0]+1
ind_siiv1403 = (np.where(lines == 'Si IV 1403'))[0][0]+1

##
## Grab the headers and slit positions
##

if nraster == 1:

    ## Extension indices
    ind_mgii = (np.where(lines == 'Mg II k 2796'))[0][0]+1
    ind_oi = (np.where(lines == 'O I 1356'))[0][0]+1
    ind_siiv1403 = (np.where(lines == 'Si IV 1403'))[0][0]+1


    # Create header to hold bulk of info from extension 0
    hdr_primary = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[0], extension = 0)

    # Hold the slit positions
    rast_dims = iris_sg_init.raster['Mg II k 2796']['_raster__dim_data']
    ntimes = rast_dims[1]
    slitxpos = np.zeros([rast_dims[1]], dtype=np.float64)
    slitypos = np.zeros([rast_dims[0], ntimes], dtype=np.float64)
    
    
    minx_coord = np.zeros(ntimes, dtype=np.float64)
    maxx_coord = np.zeros(ntimes, dtype=np.float64)
    miny_coord = np.zeros(ntimes, dtype=np.float64)
    maxy_coord = np.zeros(ntimes, dtype=np.float64)


    ## Create some dictionaries to hold the data. 
    ## There will be some repeated info, but I like having a header associated with each data cube

    hdr_mgiik = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(nraster,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }
    hdr_oi = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(nraster,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }
    hdr_siiv1403 = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(ntimes,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }
    
    hdr_mgiik_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[0]], extension = ind_mgii)
    hdr_oi_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[0]], extension = ind_oi)
    hdr_siiv1403_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[0]], extension = ind_siiv1403)
    hdr_raster_vals = iris_lmsalpy.extract_irisL2data.only_data(dir1+rasterfiles[rfinds[0]], extension = -2)
    iris_sg = extract_irisL2data.load(dir1+rasterfiles[rfinds[0]], window_info = lines[lineid])

    hdr_mgiik['NAXIS'][0] =  hdr_mgiik_tmp['NAXIS']
    hdr_mgiik['NAXIS1'][0] =  hdr_mgiik_tmp['NAXIS1']
    hdr_mgiik['NAXIS2'][0] =  hdr_mgiik_tmp['NAXIS2']
    hdr_mgiik['NAXIS3'][0] =  hdr_mgiik_tmp['NAXIS3']
    hdr_mgiik['CDELT1'][0] =  hdr_mgiik_tmp['CDELT1']
    hdr_mgiik['CDELT2'][0] =  hdr_mgiik_tmp['CDELT2']
    hdr_mgiik['CDELT3'][0] =  hdr_mgiik_tmp['CDELT3']
    hdr_mgiik['CRPIX1'][0] =  hdr_mgiik_tmp['CRPIX1']
    hdr_mgiik['CRPIX2'][0] =  hdr_mgiik_tmp['CRPIX2']
    hdr_mgiik['CRPIX2'][0] =  hdr_mgiik_tmp['CRPIX3']
    hdr_mgiik['CRVAL1'][0] =  hdr_mgiik_tmp['CRVAL1']
    hdr_mgiik['CRVAL2'][0] =  hdr_mgiik_tmp['CRVAL2']
    hdr_mgiik['CRVAL3'][0] =  hdr_mgiik_tmp['CRVAL3']
    hdr_mgiik['CUNIT1'][0] =  hdr_mgiik_tmp['CUNIT1']
    hdr_mgiik['CUNIT2'][0] =  hdr_mgiik_tmp['CUNIT2']
    hdr_mgiik['CUNIT3'][0] =  hdr_mgiik_tmp['CUNIT3']
    hdr_mgiik['exptimen'][:,0] = hdr_raster_vals[:,4]
    hdr_mgiik['exptimef'][:,0] = hdr_raster_vals[:,3]
    hdr_mgiik['xcen'][:,0] = hdr_raster_vals[:,13]
    hdr_mgiik['ycen'][:,0] = hdr_raster_vals[:,14]
   # hdr_mgiik['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']

    hdr_oi['NAXIS'][0] =  hdr_oi_tmp['NAXIS']
    hdr_oi['NAXIS1'][0] =  hdr_oi_tmp['NAXIS1']
    hdr_oi['NAXIS2'][0] =  hdr_oi_tmp['NAXIS2']
    hdr_oi['NAXIS3'][0] =  hdr_oi_tmp['NAXIS3']
    hdr_oi['CDELT1'][0] =  hdr_oi_tmp['CDELT1']
    hdr_oi['CDELT2'][0] =  hdr_oi_tmp['CDELT2']
    hdr_oi['CDELT3'][0] =  hdr_oi_tmp['CDELT3']
    hdr_oi['CRPIX1'][0] =  hdr_oi_tmp['CRPIX1']
    hdr_oi['CRPIX2'][0] =  hdr_oi_tmp['CRPIX2']
    hdr_oi['CRPIX2'][0] =  hdr_oi_tmp['CRPIX3']
    hdr_oi['CRVAL1'][0] =  hdr_oi_tmp['CRVAL1']
    hdr_oi['CRVAL2'][0] =  hdr_oi_tmp['CRVAL2']
    hdr_oi['CRVAL3'][0] =  hdr_oi_tmp['CRVAL3']
    hdr_oi['CUNIT1'][0] =  hdr_oi_tmp['CUNIT1']
    hdr_oi['CUNIT2'][0] =  hdr_oi_tmp['CUNIT2']
    hdr_oi['CUNIT3'][0] =  hdr_oi_tmp['CUNIT3']
    hdr_oi['exptimen'][:,0] = hdr_raster_vals[:,4]
    hdr_oi['exptimef'][:,0] = hdr_raster_vals[:,3]
    hdr_oi['xcen'][:,0] = hdr_raster_vals[:,13]
    hdr_oi['ycen'][:,0] = hdr_raster_vals[:,14]
   # hdr_oi['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']

    hdr_siiv1403['NAXIS'][0] =  hdr_siiv1403_tmp['NAXIS']
    hdr_siiv1403['NAXIS1'][0] =  hdr_siiv1403_tmp['NAXIS1']
    hdr_siiv1403['NAXIS2'][0] =  hdr_siiv1403_tmp['NAXIS2']
    hdr_siiv1403['NAXIS3'][0] =  hdr_siiv1403_tmp['NAXIS3']
    hdr_siiv1403['CDELT1'][0] =  hdr_siiv1403_tmp['CDELT1']
    hdr_siiv1403['CDELT2'][0] =  hdr_siiv1403_tmp['CDELT2']
    hdr_siiv1403['CDELT3'][0] =  hdr_siiv1403_tmp['CDELT3']
    hdr_siiv1403['CRPIX1'][0] =  hdr_siiv1403_tmp['CRPIX1']
    hdr_siiv1403['CRPIX2'][0] =  hdr_siiv1403_tmp['CRPIX2']
    hdr_siiv1403['CRPIX2'][0] =  hdr_siiv1403_tmp['CRPIX3']
    hdr_siiv1403['CRVAL1'][0] =  hdr_siiv1403_tmp['CRVAL1']
    hdr_siiv1403['CRVAL2'][0] =  hdr_siiv1403_tmp['CRVAL2']
    hdr_siiv1403['CRVAL3'][0] =  hdr_siiv1403_tmp['CRVAL3']
    hdr_siiv1403['CUNIT1'][0] =  hdr_siiv1403_tmp['CUNIT1']
    hdr_siiv1403['CUNIT2'][0] =  hdr_siiv1403_tmp['CUNIT2']
    hdr_siiv1403['CUNIT3'][0] =  hdr_siiv1403_tmp['CUNIT3']
    hdr_siiv1403['exptimen'][:,0] = hdr_raster_vals[:,4]
    hdr_siiv1403['exptimef'][:,0] = hdr_raster_vals[:,3]
    hdr_siiv1403['xcen'][:,0] = hdr_raster_vals[:,13]
    hdr_siiv1403['ycen'][:,0] = hdr_raster_vals[:,14]
   # hdr_siiv1403['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']

    minx_coord[:] = hdr_raster_vals[:,13]
    maxx_coord[:] = hdr_raster_vals[:,13]
    slitxpos = hdr_raster_vals[:,13]
    
    for tind in range(ntimes):
        miny_coord[tind] = np.min(hdr_mgiik['ycen'][tind,0]-hdr_mgiik['CDELT2'][0]*(hdr_mgiik['NAXIS2'][0]-1.e0)/2.e0)
        maxy_coord[tind]= np.max(hdr_mgiik['ycen'][tind,0]+hdr_mgiik['CDELT2'][0]*(hdr_mgiik['NAXIS2'][0]-1.e0)/2.e0)
        slitypos[:,tind] = np.linspace(miny_coord[tind],maxy_coord[tind],num=hdr_mgiik['NAXIS2'][0])

    del hdr_mgiik_tmp
    del hdr_oi_tmp
    del hdr_siiv1403_tmp
    del hdr_raster_vals
    del iris_sg
                                     
else:    
    ## Extension indices
    ind_mgii = (np.where(lines == 'Mg II k 2796'))[0][0]+1
    ind_oi = (np.where(lines == 'O I 1356'))[0][0]+1
    ind_siiv1403 = (np.where(lines == 'Si IV 1403'))[0][0]+1


    # Create header to hold bulk of info from extension 0
    hdr_primary = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfi], extension = 0)

    # Hold the slit positions
    rast_dims = iris_sg_init.raster['Mg II k 2796']['_raster__dim_data']
    slitxpos = np.zeros([rast_dims[1], nraster], dtype=np.float64)
    slitypos = np.zeros([rast_dims[0], nraster], dtype=np.float64)
    minx_coord = np.zeros(nraster, dtype=np.float64)
    maxx_coord = np.zeros(nraster, dtype=np.float64)
    miny_coord = np.zeros(nraster, dtype=np.float64)
    maxy_coord = np.zeros(nraster, dtype=np.float64)


    ## Create some dictionaries to hold the data. 
    ## There will be some repeated info, but I like having a header associated with each data cube

    hdr_mgiik = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(nraster,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }
    hdr_oi = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(nraster,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }
    hdr_siiv1403 = {'NAXIS':np.zeros(nraster,dtype=int),'NAXIS1':np.zeros(nraster,dtype=int),
                    'NAXIS2':np.zeros(nraster,dtype=int),'NAXIS3':np.zeros(nraster,dtype=int),
                    'CDELT1':np.zeros(nraster,dtype=float),'CDELT2':np.zeros(nraster,dtype=float),
                    'CDELT3':np.zeros(nraster,dtype=float), 'CRPIX1':np.zeros(nraster,dtype=float),
                    'CRPIX2':np.zeros(nraster,dtype=float), 'CRPIX3':np.zeros(nraster,dtype=float),
                    'CRVAL1':np.zeros(nraster,dtype=float),'CRVAL2':np.zeros(nraster,dtype=float),
                    'CRVAL3':np.zeros(nraster,dtype=float), 'CUNIT1':np.zeros(nraster,dtype=str),
                    'CUNIT2':np.zeros(nraster,dtype=str),'CUNIT3':np.zeros(nraster,dtype=str),
                    'exptimen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'exptimef':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'xcen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'ycen':np.zeros([hdr_primary['nexp'],nraster],dtype=np.float),
                    'dateobs':np.zeros([hdr_primary['nexp'],nraster],dtype='<U22')
                    }


    for i in range(nraster):

        hdr_mgiik_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[i]], extension = ind_mgii)
        hdr_oi_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[i]], extension = ind_oi)
        hdr_siiv1403_tmp = iris_lmsalpy.extract_irisL2data.only_header(dir1+rasterfiles[rfinds[i]], extension = ind_siiv1403)
        hdr_raster_vals = iris_lmsalpy.extract_irisL2data.only_data(dir1+rasterfiles[rfinds[i]], extension = -2)
        iris_sg = extract_irisL2data.load(dir1+rasterfiles[rfinds[i]], window_info = lines[lineid])

        # Define the reference pixels in python indexing (counting from zero)
        x_px_adj = hdr_mgiik_tmp['CRPIX3']- 1.0
        y_px_adj = hdr_mgiik_tmp['CRPIX2']- 1.0
        px_adj = [x_px_adj, y_px_adj]

        # Define the minimum X coordinate and the min/max y coordinates in arcsec
        # this is the box bounded by the raster. 
        minx_coord[i] = hdr_mgiik_tmp['CRVAL3']-px_adj[0]*hdr_mgiik_tmp['CDELT3']
        maxx_coord[i] = hdr_mgiik_tmp['CRVAL3']+px_adj[0]*hdr_mgiik_tmp['CDELT3']
        miny_coord[i] = hdr_mgiik_tmp['CRVAL2']-px_adj[1]*hdr_mgiik_tmp['CDELT2']
        maxy_coord[i] = hdr_mgiik_tmp['CRVAL2']+px_adj[1]*hdr_mgiik_tmp['CDELT2']

        # The pixel-centered slit positions are then 
        slitxpos[0,i] = minx_coord[i]#+0.33/2.0
        for j in range(1,rast_dims[1]):
            slitxpos[j,i] = slitxpos[j-1,i]+hdr_mgiik_tmp['CDELT3']
        slitypos[0,i] = miny_coord[i]#+hdr_mgiik_tmp['CDELT2']/2.0
        for j in range(1,rast_dims[0]):
            slitypos[j,i] = slitypos[j-1,i]+hdr_mgiik_tmp['CDELT2']


        hdr_mgiik['NAXIS'][i] =  hdr_mgiik_tmp['NAXIS']
        hdr_mgiik['NAXIS1'][i] =  hdr_mgiik_tmp['NAXIS1']
        hdr_mgiik['NAXIS2'][i] =  hdr_mgiik_tmp['NAXIS2']
        hdr_mgiik['NAXIS3'][i] =  hdr_mgiik_tmp['NAXIS3']
        hdr_mgiik['CDELT1'][i] =  hdr_mgiik_tmp['CDELT1']
        hdr_mgiik['CDELT2'][i] =  hdr_mgiik_tmp['CDELT2']
        hdr_mgiik['CDELT3'][i] =  hdr_mgiik_tmp['CDELT3']
        hdr_mgiik['CRPIX1'][i] =  hdr_mgiik_tmp['CRPIX1']
        hdr_mgiik['CRPIX2'][i] =  hdr_mgiik_tmp['CRPIX2']
        hdr_mgiik['CRPIX2'][i] =  hdr_mgiik_tmp['CRPIX3']
        hdr_mgiik['CRVAL1'][i] =  hdr_mgiik_tmp['CRVAL1']
        hdr_mgiik['CRVAL2'][i] =  hdr_mgiik_tmp['CRVAL2']
        hdr_mgiik['CRVAL3'][i] =  hdr_mgiik_tmp['CRVAL3']
        hdr_mgiik['CUNIT1'][i] =  hdr_mgiik_tmp['CUNIT1']
        hdr_mgiik['CUNIT2'][i] =  hdr_mgiik_tmp['CUNIT2']
        hdr_mgiik['CUNIT3'][i] =  hdr_mgiik_tmp['CUNIT3']
        hdr_mgiik['exptimen'][:,i] = hdr_raster_vals[:,4]
        hdr_mgiik['exptimef'][:,i] = hdr_raster_vals[:,3]
        hdr_mgiik['xcen'][:,i] = hdr_raster_vals[:,13]
        hdr_mgiik['ycen'][:,i] = hdr_raster_vals[:,14]
       # hdr_mgiik['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']

        hdr_oi['NAXIS'][i] =  hdr_oi_tmp['NAXIS']
        hdr_oi['NAXIS1'][i] =  hdr_oi_tmp['NAXIS1']
        hdr_oi['NAXIS2'][i] =  hdr_oi_tmp['NAXIS2']
        hdr_oi['NAXIS3'][i] =  hdr_oi_tmp['NAXIS3']
        hdr_oi['CDELT1'][i] =  hdr_oi_tmp['CDELT1']
        hdr_oi['CDELT2'][i] =  hdr_oi_tmp['CDELT2']
        hdr_oi['CDELT3'][i] =  hdr_oi_tmp['CDELT3']
        hdr_oi['CRPIX1'][i] =  hdr_oi_tmp['CRPIX1']
        hdr_oi['CRPIX2'][i] =  hdr_oi_tmp['CRPIX2']
        hdr_oi['CRPIX2'][i] =  hdr_oi_tmp['CRPIX3']
        hdr_oi['CRVAL1'][i] =  hdr_oi_tmp['CRVAL1']
        hdr_oi['CRVAL2'][i] =  hdr_oi_tmp['CRVAL2']
        hdr_oi['CRVAL3'][i] =  hdr_oi_tmp['CRVAL3']
        hdr_oi['CUNIT1'][i] =  hdr_oi_tmp['CUNIT1']
        hdr_oi['CUNIT2'][i] =  hdr_oi_tmp['CUNIT2']
        hdr_oi['CUNIT3'][i] =  hdr_oi_tmp['CUNIT3']
        hdr_oi['exptimen'][:,i] = hdr_raster_vals[:,4]
        hdr_oi['exptimef'][:,i] = hdr_raster_vals[:,3]
        hdr_oi['xcen'][:,i] = hdr_raster_vals[:,13]
        hdr_oi['ycen'][:,i] = hdr_raster_vals[:,14]
       # hdr_oi['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']

        hdr_siiv1403['NAXIS'][i] =  hdr_siiv1403_tmp['NAXIS']
        hdr_siiv1403['NAXIS1'][i] =  hdr_siiv1403_tmp['NAXIS1']
        hdr_siiv1403['NAXIS2'][i] =  hdr_siiv1403_tmp['NAXIS2']
        hdr_siiv1403['NAXIS3'][i] =  hdr_siiv1403_tmp['NAXIS3']
        hdr_siiv1403['CDELT1'][i] =  hdr_siiv1403_tmp['CDELT1']
        hdr_siiv1403['CDELT2'][i] =  hdr_siiv1403_tmp['CDELT2']
        hdr_siiv1403['CDELT3'][i] =  hdr_siiv1403_tmp['CDELT3']
        hdr_siiv1403['CRPIX1'][i] =  hdr_siiv1403_tmp['CRPIX1']
        hdr_siiv1403['CRPIX2'][i] =  hdr_siiv1403_tmp['CRPIX2']
        hdr_siiv1403['CRPIX2'][i] =  hdr_siiv1403_tmp['CRPIX3']
        hdr_siiv1403['CRVAL1'][i] =  hdr_siiv1403_tmp['CRVAL1']
        hdr_siiv1403['CRVAL2'][i] =  hdr_siiv1403_tmp['CRVAL2']
        hdr_siiv1403['CRVAL3'][i] =  hdr_siiv1403_tmp['CRVAL3']
        hdr_siiv1403['CUNIT1'][i] =  hdr_siiv1403_tmp['CUNIT1']
        hdr_siiv1403['CUNIT2'][i] =  hdr_siiv1403_tmp['CUNIT2']
        hdr_siiv1403['CUNIT3'][i] =  hdr_siiv1403_tmp['CUNIT3']
        hdr_siiv1403['exptimen'][:,i] = hdr_raster_vals[:,4]
        hdr_siiv1403['exptimef'][:,i] = hdr_raster_vals[:,3]
        hdr_siiv1403['xcen'][:,i] = hdr_raster_vals[:,13]
        hdr_siiv1403['ycen'][:,i] = hdr_raster_vals[:,14]
       # hdr_siiv1403['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok']


        del hdr_mgiik_tmp
        del hdr_oi_tmp
        del hdr_siiv1403_tmp
        del hdr_raster_vals
        del iris_sg

'''
Grab the data from the raster object
> - Similar to the above, but this time grabbing the actual spectra.
> - data = the spectra, in DN
> - wl = the wavelengths at which each spectra is defined

> - data = [Y,X,wave,R].
> - For sit-and-stare: X = time as well as solar-X (since the Sun rotates over time, so Solar-X changes; Y = solar-Y along slit; wave = wavelength; R = raster number.
> - For multiple slit positions: X = solar-X, one for each slit (really, this is kind of time also, since a raster takes time to produce, but its easier to think of as Solar-X); Y = solar-Y along slit; wave = wavelength; R = raster number, which is also time since it can represent the repeat cadence at each location.
'''
##
## Extract the data from the RoI object
##

wl_mgiik = np.zeros([hdr_mgiik['NAXIS1'][0],nraster],dtype=np.float)
wl_oi = np.zeros([hdr_oi['NAXIS1'][0],nraster],dtype=np.float)
wl_siiv1403 = np.zeros([hdr_siiv1403['NAXIS1'][0],nraster],dtype=np.float)



data_mgiik = np.zeros([hdr_mgiik['NAXIS2'][0],hdr_mgiik['NAXIS3'][0],hdr_mgiik['NAXIS1'][0],nraster],dtype=np.float)
data_oi = np.zeros([hdr_oi['NAXIS2'][0],hdr_oi['NAXIS3'][0],hdr_oi['NAXIS1'][0],nraster],dtype=np.float)
data_siiv1403 = np.zeros([hdr_siiv1403['NAXIS2'][0],hdr_siiv1403['NAXIS3'][0],hdr_siiv1403['NAXIS1'][0],nraster],dtype=np.float)


for i in range(nraster):
    iris_sg = extract_irisL2data.load(dir1+rasterfiles[rfinds[i]], window_info = lines[lineid])
    
    data_mgiik[:,:,:,i] = iris_sg.raster['Mg II k 2796']['data'].copy()
    data_oi[:,:,:,i] = iris_sg.raster['O I 1356']['data'].copy()
    data_siiv1403[:,:,:,i] = iris_sg.raster['Si IV 1403']['data'].copy()

    wl_mgiik[:,i] = iris_sg.raster['Mg II k 2796']['wl'].copy()
    wl_oi[:,i] = iris_sg.raster['O I 1356']['wl'].copy()
    wl_siiv1403[:,i] = iris_sg.raster['Si IV 1403']['wl'].copy()

    hdr_mgiik['dateobs'][:,i] = iris_sg.raster['Mg II k 2796']['date_time_acq_ok'].copy()
    hdr_oi['dateobs'][:,i] = iris_sg.raster['O I 1356']['date_time_acq_ok'].copy()
    hdr_siiv1403['dateobs'][:,i] = iris_sg.raster['Si IV 1403']['date_time_acq_ok'].copy()

    del iris_sg

'''
---
### Save the raw data 
> - Saves a different dictionary per passband
> - These are the spectra in DN
'''

##
## Create a dictionary in which to save the data
##

sg_mgiik_dict = {'data':data_mgiik, 
                'wl':wl_mgiik,
                'hdr':hdr_mgiik,
                'minx_coord':minx_coord,
                'maxx_coord':maxx_coord,
                'miny_coord':miny_coord,
                'maxy_coord':maxy_coord,
                'slitposx':slitxpos,
                'slitposy':slitypos,
                'readme':'SG data is in DN/px'}
file_sg_mgiik = './IRIS_SG_2014_Sept_10_mgiik.pkl'
with open(file_sg_mgiik, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_mgiik_dict, output, pickle.HIGHEST_PROTOCOL)

sg_oi_dict = {'data':data_oi, 
                'wl':wl_oi,
                'hdr':hdr_oi,
                'minx_coord':minx_coord,
                'maxx_coord':maxx_coord,
                'miny_coord':miny_coord,
                'maxy_coord':maxy_coord,
                'slitposx':slitxpos,
                'slitposy':slitypos,
                'readme':'SG data is in DN/px'}
file_sg_oi = './IRIS_SG_2014_Sept_10_oi.pkl'
with open(file_sg_oi, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_oi_dict, output, pickle.HIGHEST_PROTOCOL)

sg_siiv1403_dict = {'data':data_siiv1403, 
                'wl':wl_siiv1403,
                'hdr':hdr_siiv1403,
                'minx_coord':minx_coord,
                'maxx_coord':maxx_coord,
                'miny_coord':miny_coord,
                'maxy_coord':maxy_coord,
                'slitposx':slitxpos,
                'slitposy':slitypos,
                'readme':'SG data is in DN/px'}
file_sg_siiv1403 = './IRIS_SG_2014_Sept_10_siiv1403.pkl'
with open(file_sg_siiv1403, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_siiv1403_dict, output, pickle.HIGHEST_PROTOCOL)

###
### The code to reload the files
###

file_sg = './IRIS_SG_2014_Sept_10_mgiik.pkl'
with open(file_sg, 'rb') as output:  
    sg_mgiik_dict = pickle.load(output)
    
file_sg = './IRIS_SG_2014_Sept_10_oi.pkl'
with open(file_sg, 'rb') as output:  
    sg_oi_dict = pickle.load(output)
    
file_sg = './IRIS_SG_2014_Sept_10_siiv1403.pkl'
with open(file_sg, 'rb') as output:  
    sg_siiv1403_dict = pickle.load(output)

'''
### Correct for exposure time 
> - Corrects for exposure time.
> - Saves the new data
> - Note that the NUV and FUV channels have different exposure times
'''
##
## Correct for exposure time
##
for i in range(nraster):
    for j in range(rast_dims[1]):
        sg_mgiik_dict['data'][:,j,:,i]/=sg_mgiik_dict['hdr']['exptimen'][j,i]
        sg_oi_dict['data'][:,j,:,i]/=sg_oi_dict['hdr']['exptimef'][j,i]
        sg_siiv1403_dict['data'][:,j,:,i]/=sg_siiv1403_dict['hdr']['exptimef'][j,i]
    
sg_mgiik_dict['readme'] = 'SG data is has been corrected for exposure time, DN/s/px'
sg_oi_dict['readme'] = 'SG data is has been corrected for exposure time, DN/s/px'
sg_siiv1403_dict['readme'] = 'SG data is has been corrected for exposure time, DN/s/px'

file_sg_mgiik = './IRIS_SG_2014_Sept_10_mgiik_expcorr.pkl'
with open(file_sg_mgiik, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_mgiik_dict, output, pickle.HIGHEST_PROTOCOL)

file_sg_siiv1403 = './IRIS_SG_2014_Sept_10_siiv1403_expcorr.pkl'
with open(file_sg_siiv1403, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_siiv1403_dict, output, pickle.HIGHEST_PROTOCOL)

file_sg_oi = './IRIS_SG_2014_Sept_10_oi_expcorr.pkl'
with open(file_sg_oi, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sg_oi_dict, output, pickle.HIGHEST_PROTOCOL)


