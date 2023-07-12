# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 15:23:49 2022

@author: Jack Cribbin - 19328253

This file reads in raw FITS data files for bias, flat, scientific and star 
field exposures stored alongside this file, and combines
and normalises them into a single Master file
"""

# Code written to perform reduction of images taken at OHP.
import numpy as np
import ccdproc

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ZScaleInterval, ImageNormalize, AsinhStretch
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename

import os
import matplotlib.pyplot as plt
import re

from scipy.stats import mode 

import warnings
warnings.filterwarnings("ignore")

# Function to get a list of files matching a given regular expression in 
# a given directory
def file_grabber(directory, regex):
    '''
    

    Parameters
    ----------
    directory : String
        The directory name where the fiels being fetched from are stored
    regex : String
        A String that contains the regular experession that will fetch
        the correct files and only the correct files from the given folder

    Returns
    -------
    String []
        A list of Strings that are filenames for the fetched files

    '''
    
    # Move to the directory where the file is located 
    directorysuper = os.path.dirname(os.path.abspath(__file__))
    os.chdir(directorysuper)
    print('\nIn:',directorysuper)
    
    image_list = []
    for filename in os.scandir(directory):
        if re.search(regex+'.fits',filename.name) or re.search(regex+'.fit',filename.name):
            image_list.append(filename.name)
            
    if(len(image_list)>1):
        return image_list
    elif(len(image_list) == 1):
        return image_list[0]
    
  

# Move to the directory where the file is located 
directorysuper = os.path.dirname(os.path.abspath(__file__))
os.chdir(directorysuper)
print('\nIn:',directorysuper)

'''
############ Code for making Masterbias File ############

# Get a list of the bias filenames using file_grabber
biaslist = file_grabber('Files/Biases and Flats/18','.*Bias.*')
    
# Combine the bias files into a MasterBias file
os.chdir('Files/Biases and Flats/18')
MasterBias = ccdproc.combine(biaslist, None, 'median', unit='adu')

############ End of Code for making Masterbias File ############





############ Code for making Masterflat File ############

os.chdir(directorysuper)

# Get a list of the flat filenames using file_grabber
flatlist = file_grabber('Files/Biases and Flats/18','.*Flat.*')

# Combine the flat files into a MasterFlat file
os.chdir('Files/Biases and Flats/18')
MasterFlat = ccdproc.combine(flatlist, None, 'median', unit='adu')

############ End of Code for making Masterbflat File ############
'''




############ Code for making MasterSci File ############

os.chdir(directorysuper)

# Get a list of the scientific data filenames using file_grabber
Scilistnames = file_grabber('Files/Night3Data', '.*_SII')

# Convert the scientific data from filenames into CCDData objects
os.chdir('Files/Night3Data')
Scilist = []
for file in Scilistnames:
    Scilist.append(CCDData(fits.getdata(file), unit='adu'))
'''
# Subtract the bias from the flat
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias)

# Subtract the bias from each of the scientific data files
BiasSubScilist = []
#for file in Scilist:
 #   BiasSubScilist.append(ccdproc.subtract_bias(file, MasterBias))

# Divide the scientific data files by the flat file
FinalScilist = []
#for file in BiasSubScilist:
 #   FinalScilist.append(ccdproc.flat_correct(file, BiasSubFlat))

# Subtract the bias from each of the scientific data files
FinalScilistdel = []
#for file in FinalScilist:
 #   FinalScilistdel.append(CCDData(ccdproc.subtract_bias(file, MasterBias), unit='adu'))
'''
# Combine the reduced scientific data into a master file
MasterSci = ccdproc.combine(Scilist, None, 'median', unit='adu', overwrite=True)

############ End of Code for making MasterSci File ############




'''
############ Code for making Masterstar File ############

# Get a list of the broadband star data filenames using file_grabber
os.chdir(directorysuper)
Starlistnames = file_grabber('Files/Night3Data', '.*_R')

# Convert the star data from filenames into CCDData objects
os.chdir('Files/Night3Data')
Starlist = []
for file in Starlistnames:
    Starlist.append(CCDData(fits.getdata(file), unit='adu'))

# Subtract the bias from each of the star data files
BiasSubStarlist = []
for file in Starlist:
    BiasSubStarlist.append(ccdproc.subtract_bias(file, MasterBias))

# Divide the star data files by the flat file
FinalStarlist = []
for file in BiasSubStarlist:
    FinalStarlist.append(ccdproc.flat_correct(file, BiasSubFlat))

# Combine the reduced scientific data into a master file
MasterStar = ccdproc.combine(FinalStarlist, None, 'median', unit='adu', overwrite=True)

# Shift the star image to line up with the narrowband image
MasterStar = np.roll(np.asarray(MasterStar), 0, axis=0)
MasterStar = np.roll(np.asarray(MasterStar), 14, axis=1)
MasterStar = CCDData((MasterStar), unit='adu')

# Normalise the star image so that is can be subtracte form the narrowband image
ratio = (mode((np.asarray(MasterSci)), keepdims=True)[0][0][0]/
         mode((np.asarray(MasterStar)), keepdims=True)[0][0][0])
MasterStar = CCDData(np.asarray(MasterStar) * ratio, unit = 'adu')

############ End of Code for making Masterstar File ############
'''




# Subtract the star image from the narrowband image
#MasterSci = ccdproc.subtract_bias(MasterSci, MasterStar)

# Initialise normalisaton for the image 
interval = ZScaleInterval()
stretch = AsinhStretch(.85)
norm = ImageNormalize(MasterSci,interval=interval, stretch=stretch)

# Get WCS coordinates a file that has been fed into astrometry.net
directorysuper = os.path.dirname(os.path.abspath(__file__))
os.chdir(directorysuper)
directory = 'Files/ProcessedFilesNew'
fn = get_pkg_data_filename(directory+'/'+'CasA2022_WCS.fits')
f = fits.open(fn)
wcs = WCS(f[0].header)

# Normalise and stretch the image, and add the WCS
MasterSci = interval(MasterSci)
MasterSci = stretch(MasterSci)
MasterSci = CCDData(MasterSci,unit='adu',wcs=wcs)

# Flip the image about the x axis (due to telescope formatting)
MasterSci = CCDData(np.flip(np.asarray(MasterSci),0),unit='adu',wcs=wcs)

# Save the fits file to an external file
os.chdir(directorysuper)
os.chdir('Files/ProcessedFilesNew')
MasterSci.write('Sci-Reduction.fits', overwrite=True)
    
# Display the image with WCS as the coordinate system
plt.figure()
plt.subplot(projection=wcs)
plt.imshow(MasterSci)
plt.xlim(300,750)
plt.ylim(300,750)
plt.colorbar()
plt.show()

print('\nDone!')






