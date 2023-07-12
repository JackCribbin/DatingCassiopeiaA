# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 17:09:33 2022

@author: jackp
"""

# Code written to perform reduction of images taken at OHP.
import numpy as np
import ccdproc
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ZScaleInterval, ImageNormalize, AsinhStretch
from astropy.stats import sigma_clip

import os
import matplotlib.pyplot as plt
import re

from scipy.stats import mode 

#import twirl

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


# Get a list of the bias filenames using file_grabber
biaslist = file_grabber('Files/2007Data/Biases','.*bias.*')

# Combine the bias files into a MasterBias file
os.chdir('Files/2007Data/Biases')
MasterBias = ccdproc.combine(biaslist, None, 'median', unit='adu')

os.chdir(directorysuper)

# Get a list of the flat filenames using file_grabber
flatlist = file_grabber('Files/2007Data/Flats','.*flat.*')

flatlist=[flatlist]

# Combine the flat files into a MasterFlat file
os.chdir('Files/2007Data/Flats')

#f(type(flatlist) != 'str'):
MasterFlat = ccdproc.combine(flatlist, None, 'median', unit='adu')
#MasterFlat = CCDData(flatlist,unit='adu')


os.chdir(directorysuper)

# Get a list of the scientific data filenames using file_grabber
Scilistnames = file_grabber('Files/2007Data/SciFiles', '.*_SII')

# Convert the scientific data from filenames into CCDData objects
os.chdir('Files/2007Data/SciFiles')
Scilist = []
for file in Scilistnames:
    Scilist.append(CCDData(fits.getdata(file), unit='adu'))

# Subtract the bias from the flat
BiasSubFlat = ccdproc.subtract_bias(MasterFlat, MasterBias)

# Subtract the bias from each of the scientific data files
BiasSubScilist = []
for file in Scilist:
    BiasSubScilist.append(ccdproc.subtract_bias(file, MasterBias))

# Divide the scientific data files by the flat file
FinalScilist = []
for file in BiasSubScilist:
    FinalScilist.append(ccdproc.flat_correct(file, BiasSubFlat))

# Subtract the bias from each of the scientific data files
FinalScilistdel = []
for file in FinalScilist:
    FinalScilistdel.append(CCDData(ccdproc.subtract_bias(file, MasterBias), unit='adu'))

# Combine the reduced scientific data into a master file
MasterSci = ccdproc.combine(FinalScilist, None, 'median', unit='adu', overwrite=True)

# Initialise normalisaton for the image 




### Code to save file ###

from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename

directory = 'Files/ProcessedFiles'
fn = get_pkg_data_filename(directory+'/'+'CasA2007_WCS.fits')
f = fits.open(fn)
wcs = WCS(f[0].header)

#MasterSci = interval(MasterSci)
#MasterSci = stretch(MasterSci)
MasterSci = CCDData(MasterSci,unit='adu',wcs=wcs)


interval = ZScaleInterval()
stretch = AsinhStretch(.85)
norm = ImageNormalize(MasterSci,interval=interval, stretch=stretch)
MasterSci = interval(MasterSci)
MasterSci = stretch(MasterSci)
MasterSci = CCDData(MasterSci,unit='adu')

os.chdir(directorysuper)
os.chdir('Files/ProcessedFilesNew')
#MasterSci.write('CasA2007_WCS.fits', overwrite=True)

### End of code to save file ###
    
# Plot the image
plt.subplot(projection=wcs)
plt.imshow(MasterSci)#, norm=norm)
plt.colorbar()
plt.show()

print('\nDone!')






