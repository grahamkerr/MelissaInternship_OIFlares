'''
Graham Kerr
graham.s.kerr@nasa.gov ; kerrg@cua.edu
June 2022

O I in Flares

NAME: IRIS_SJI_initial.py

PURPOSE: To perform some intial analysis on IRIS SJI data, reading it into memory, 
         normalising exposure time, and performing quicklooks. 

INPUTS: None (probably best to open this up and run it line by line in the terminal/ipython)

OUTPUTS: Some .pkl files (these can get large, so maybe isolate the flare times before saving)

NOTES:

'''
###################
### SOME SET UP 
###################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
import re
from scipy import io
import pickle
import cmocean
import copy
import iris_lmsalpy
import astropy
import pandas as pd

'''
### Read in the data 
> - Just comment out if there if any of the filters are missing
'''
##
## The SJI filenames
## 

dir1 = '/Users/gskerr1/Documents/Research/Melissa_OI_IRIS/2014_09_10_1130/'

filename_1400 = dir1+'iris_l2_20140910_112825_3860259453_SJI_1400_t000.fits'
# filename_1330 = dir1+'iris_l2_20140910_112825_3860259453_SJI_1330_t000.fits'
# filename_2832 = dir1+'iris_l2_20140910_112825_3860259453_SJI_2832_t000.fits'
filename_2796 = dir1+'iris_l2_20140910_112825_3860259453_SJI_2796_t000.fits'

##
## Create SJI objects
##
sji_1400 = iris_lmsalpy.extract_irisL2data.load(filename_1400)
# sji_1330 = iris_lmsalpy.extract_irisL2data.load(filename_1330)
# sji_2832 = iris_lmsalpy.extract_irisL2data.load(filename_2832)
sji_2796 = iris_lmsalpy.extract_irisL2data.load(filename_2796)

'''
> - The data cubes [y,x,time] are saved as memory map objects, but if we want to modify them (e.g. normalise by exposure time etc.,) then we need to have them in the computer memory properly.
'''
##
## Bring them into local memory so we can modify them
##
sji_1400.flush()
# sji_1330.flush()
# sji_2832.flush()
sji_2796.flush()

'''
> - The SJI data are image cubes of [x,y,time], and are accessed from each object.
> - Each object is a dictionary, containing certain 'keys'
'''
## This tell us that the variable 'SJI_1400' is a key of sji_1400.SJI object
sji_1400.SJI.keys()

## held within that is another dictionary (I know...)
sji_1400.SJI['SJI_1400'].keys()

## data holds the actual images, the size of which is
sji_1400.SJI['SJI_1400'].data.shape

'''
The header files contain lots of metadata and useful info we need for the analysis***
> - There are different extensions to the fits headers that each contain different info, some of them are part of the time-series info, some are more general for the full observing run.
> - To be honest I find them rather confusing sometimes, but some experimenting, checking the code for iris_lmsalpy, and Alberto's docs for iris_lmsalpy helped me figure out what to grab 
'''
##
## Create header to hold bulk of info from extension 0
##
hdr_1400_ext0 = iris_lmsalpy.extract_irisL2data.only_header(filename_1400, extension = 0)
# hdr_1330_ext0 = iris_lmsalpy.extract_irisL2data.only_header(filename_1330, extension = 0)
# hdr_2832_ext0 = iris_lmsalpy.extract_irisL2data.only_header(filename_2832, extension = 0)
hdr_2796_ext0 = iris_lmsalpy.extract_irisL2data.only_header(filename_2796, extension = 0)

##
## Create header to hold bulk of info from extension 1
##
hdr_1400_ext1 = iris_lmsalpy.extract_irisL2data.only_header(filename_1400, extension = 1)
# hdr_1330_ext1 = iris_lmsalpy.extract_irisL2data.only_header(filename_1330, extension = 1)
# hdr_2832_ext1 = iris_lmsalpy.extract_irisL2data.only_header(filename_2832, extension = 1)
hdr_2796_ext1 = iris_lmsalpy.extract_irisL2data.only_header(filename_2796, extension = 1)

##
## Create header to hold bulk of info from extension 2
##
hdr_1400_ext2 = iris_lmsalpy.extract_irisL2data.only_header(filename_1400, extension = 2)
# hdr_1330_ext2 = iris_lmsalpy.extract_irisL2data.only_header(filename_1330, extension = 2)
# hdr_2832_ext2 = iris_lmsalpy.extract_irisL2data.only_header(filename_2832, extension = 2)
hdr_2796_ext2 = iris_lmsalpy.extract_irisL2data.only_header(filename_2796, extension = 2)

'''
> - Lets look at what is inside the header ext 0
> - Then, access one of the variables... we will choose the naxis 3 variable, which is the number of frames we have of 1400A images. 
'''
hdr_1400_ext0.keys

numframes_1400 = hdr_1400_ext0['naxis3']
numframes_2796 = hdr_2796_ext0['naxis3']


print("number of 1400A frames = %d" %(numframes_1400))
print("number of 2796A frames = %d" %(numframes_2796))

'''
> - Lets look at what is inside the header ext 1 and 2
> - These are a bit different, and tell us what index of the header file data (next set of cells) each parameter is held in
'''
hdr_1400_ext1.keys

hdr_1400_ext2.keys

'''
> - Now we grab the actual values from the header files
'''
##
## Grab the data for the header extension 1
##
data_1400_ext1 = iris_lmsalpy.extract_irisL2data.only_data(filename_1400, extension=1)
# data_1330_ext1 = iris_lmsalpy.extract_irisL2data.only_data(filename_1330, extension=1)
# data_2832_ext1 = iris_lmsalpy.extract_irisL2data.only_data(filename_2832, extension=1)
data_2796_ext1 = iris_lmsalpy.extract_irisL2data.only_data(filename_2796, extension=1)

##
## Grab the data for the header extension 2
##
data_1400_ext2 = iris_lmsalpy.extract_irisL2data.only_data(filename_1400, extension=2)
# data_1330_ext2 = iris_lmsalpy.extract_irisL2data.only_data(filename_1330, extension=2)
# data_2832_ext2 = iris_lmsalpy.extract_irisL2data.only_data(filename_2832, extension=2)
data_2796_ext2 = iris_lmsalpy.extract_irisL2data.only_data(filename_2796, extension=2)

'''
Now we assign some of those values to variables that we want to use later
> - number of exposures (frames) for each filter
'''
##
## Number of exposures
##
nexp_1400 = hdr_1400_ext0['naxis3']
# nexp_1330 = hdr_1330_ext0['naxis3']
# nexp_2832 = hdr_2832_ext0['naxis3']
nexp_2796 = hdr_2796_ext0['naxis3']

'''
> - number of pixels in the x- and y- direction (usually this is also x- and y- on the Sun, but sometimes the spacecraft is rotated by 45 or 90 degrees)
'''
##
## Number of x and y pixels
##
nx_1400 = hdr_1400_ext0['naxis1']
# nx_1330 = hdr_1330_ext0['naxis1']
# nx_2832 = hdr_2832_ext0['naxis1']
nx_2796 = hdr_2796_ext0['naxis1']

ny_1400 = hdr_1400_ext0['naxis2']
# ny_1330 = hdr_1330_ext0['naxis2']
# ny_2832 = hdr_2832_ext0['naxis2']
ny_2796 = hdr_2796_ext0['naxis2']

'''
> - The exposure time for each exposure (how long was the shutter on the camera open). This is sometimes constant, but usually there is some level of variation. This variation may be small or if the automatic flare detection flag was triggered the exposure time can drop suddenely (in flares the exposure time might drop in order to not saturate that camera).
> - We know what index this is held in by looking at the ext1 header info above. 
> - Lets plot the exposure times to see how much they vary.
'''
##
## Grab exposure times 
##
exptimes_1400 = data_1400_ext1[:,3]
# exptimes_1330 = data_1330_ext1[:,3]
# exptimes_2832 = data_2832_ext1[:,3]
exptimes_2796 = data_2796_ext1[:,3]

%matplotlib inline
plt.plot(exptimes_1400, color = 'black', label = '1400', linestyle = '-')
# plt.plot(exptimes_1330, color = 'tomato', label = '1330',linestyle = '--')
# plt.plot(exptimes_2832, color = 'forestgreen', label = '2832', linestyle = ':')
plt.plot(exptimes_2796, color = 'dodgerblue', label = '2796', linestyle = '-.')
plt.show()

'''
---
### Correct for Exposure time 
> - Since exposure time varies over time, we want to make sure that if we compare different frames we comparing apples-to-apples. 
> - The data are in Counts, often referred to as Data Numbers (DN).
> - Dividing by the exposure time (how long the telescope was looking for) we can get DN/s, which can be more directly compared from frame to frame.
'''
##
## Correct for Exposure Times
##
sji_1400.SJI['SJI_1400'].data/=exptimes_1400
# sji_1330.SJI['SJI_1330'].data/=exptimes_1330
# sji_2832.SJI['SJI_2832'].data/=exptimes_2832
sji_2796.SJI['SJI_2832'].data/=exptimes_2796

'''
### Lets use the in-built iris_lmsalpy methods for quicklook 

> - First set some limits (trial and error, you'll probably have to fiddle with the intensity ranges to clip).
> - Then run the quicklook method, which opens an interactive window (check Alberto's docs for some keystroke commands to manually change the intensity ranges etc.,)
'''
##
## Some stuff for plotting
##
sji_1400.SJI['SJI_1400'].show_slit = True
sji_1400.SJI['SJI_1400'].clip_ima = [0, 1000]

# sji_1330.SJI['SJI_1330'].show_slit = False
# sji_1330.SJI['SJI_1330'].clip_ima = [0, 200]

# sji_2832.SJI['SJI_2832'].show_slit = False
# sji_2832.SJI['SJI_2832'].clip_ima = [0, 2500]

sji_2796.SJI['SJI_2796'].show_slit = True
sji_2796.SJI['SJI_2796'].clip_ima = [0, 2500]

#
# Quick look plotting
#
%matplotlib tk
sji_1400.quick_look()

## You might have to run this cell twice when loading the 
#tk backend to make it plot in an interactive window

#
# Quick look plotting
#
#%matplotlib tk
sji_2796.quick_look()


## Bring it back to plotting within the notebook (might not be needed in ipython??)
%matplotlib inline


'''
### Grab the slit locations, in arcseconds 
> - Uses a few variables from the header files, we can discuss what those variables are, or check the header ext info
'''
##
## Slit location in arcsecs
## 

lcslit_1400 = ( (sji_1400.SJI['SJI_1400']['SLTPX1IX'] - np.repeat(hdr_1400_ext0['CRPIX1'],nexp_1400)) * 
          hdr_1400_ext0['CDELT1'] + 
          sji_1400.SJI['SJI_1400'].XCENIX)
# lcslit_1330 = ( (sji_1330.SJI['SJI_1330']['SLTPX1IX'] - np.repeat(hdr_1330_ext0['CRPIX1'],nexp_1330)) * 
#           hdr_1330_ext0['CDELT1'] + 
#           sji_1330.SJI['SJI_1330'].XCENIX)
# lcslit_2832 = ( (sji_2832.SJI['SJI_2832']['SLTPX1IX'] - np.repeat(hdr_2832_ext0['CRPIX1'],nexp_2832)) * 
#           hdr_2832_ext0['CDELT1'] + 
#           sji_2832.SJI['SJI_2832'].XCENIX)
lcslit_2796 = ( (sji_2796.SJI['SJI_2796']['SLTPX1IX'] - np.repeat(hdr_2796_ext0['CRPIX1'],nexp_2796)) * 
          hdr_2796_ext0['CDELT1'] + 
          sji_2796.SJI['SJI_2796'].XCENIX)

'''
### If you want to save the data and header info for use later 

> - Saves the data, and header info, in a dictionary that you can load in another notebook for analysis
> - Doing this a lot can make some pretty large files... some people like just running through the above each time they analyse the SJI, but I prefer to save them and load them in later for future analysis
> - This saves them all in one big file, so all the SJI wavelengths in a single place, but you can modify the code to save each filter separately if you wish.
'''
##
## Create a dictionary in which to save the data
##
sji_dict = {
#             'sji_2832':sji_2832.SJI['SJI_2832']['data'], 
#             'sji_1330':sji_1330.SJI['SJI_1330']['data'],
            'sji_1400':sji_1400.SJI['SJI_1400']['data'],
            'sji_2796':sji_2796.SJI['SJI_2796']['data'],
#             'xcen_2832':sji_2832.SJI['SJI_2832'].XCENIX,
#             'xcen_1330':sji_1330.SJI['SJI_1330'].XCENIX,
            'xcen_1400':sji_1400.SJI['SJI_1400'].XCENIX,
            'xcen_2796':sji_2796.SJI['SJI_2796'].XCENIX,
#             'ycen_2832':sji_2832.SJI['SJI_2832'].YCENIX,
#             'ycen_1330':sji_1330.SJI['SJI_1330'].YCENIX,
            'ycen_1400':sji_1400.SJI['SJI_1400'].YCENIX,
            'ycen_2796':sji_2796.SJI['SJI_2796'].YCENIX,
#             'SLTPX1IX_2832':sji_2832.SJI['SJI_2832']['SLTPX1IX'],
#             'SLTPX1IX_1330':sji_1330.SJI['SJI_1330']['SLTPX1IX'],
            'SLTPX1IX_1400':sji_1400.SJI['SJI_1400']['SLTPX1IX'],
            'SLTPX1IX_2796':sji_2796.SJI['SJI_2796']['SLTPX1IX'],
#             'lcslit_2832':lcslit_2832,
#             'lcslit_1330':lcslit_1330,
            'lcslit_1400':lcslit_1400,
            'lcslit_2796':lcslit_2796,
#             'hdr_2832_ext0':hdr_2832_ext0,
#             'hdr_1330_ext0':hdr_1330_ext0,
            'hdr_1400_ext0':hdr_1400_ext0,
            'hdr_2796_ext0':hdr_2796_ext0,
#             'time_2832':sji_2832.SJI['SJI_2832'].date_time_acq_ok,
#             'time_1330':sji_1330.SJI['SJI_1330'].date_time_acq_ok,
            'time_1400':sji_1400.SJI['SJI_1400'].date_time_acq_ok,
            'time_2796':sji_2796.SJI['SJI_2796'].date_time_acq_ok,
            'readme':'SJI data has been exposure corrected'}
file_sji = './IRIS_SJI_2014_Sept_10_exptimecorr.pkl'
with open(file_sji, 'wb') as output:  # Overwrites any existing file.
    pickle.dump(sji_dict, output, pickle.HIGHEST_PROTOCOL)


'''
> - This is the command to read the pickle (pkl) file back into memory
'''
file_sji = 'IRIS_SJI_2014_Sept_10_exptimecorr.pkl'
with open(file_sji, 'rb') as output:  
    sji_dict = pickle.load(output)
