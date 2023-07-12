# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 16:19:13 2022

@author: Jack Cribbin - 19328253

This file plots the knot locations and knot trajectories, as well as the
expansion centre calculated, on a processed and reduced image of Cas A
"""

# Code written to perform reduction of images taken at OHP.
import numpy as np

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ZScaleInterval, ImageNormalize, AsinhStretch
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename

import os
import matplotlib.pyplot as plt

import re


import warnings
warnings.filterwarnings("ignore")

# Function to get a list of files matching a given regular expression in 
# a given directory
def file_grabber(directory, regex):
    '''
    

    Parameters
    ----------
    directory : String
        The directory name where the files being fetched from are stored
    regex : String
        A String that contains the regular experession that will fetch
        the correct files and only the correct files from the given folder

    Returns
    -------
    If multiple files match the given regex in the given directory
    String []
        A list of Strings that are filenames for the fetched files
        
    If only one file matches the given regex in the given directory
    CCDData
        The file converted to fits data and then converted to a CCDData object
    '''
    
    # Move to the directory where the file is located 
    directorysuper = os.path.dirname(os.path.abspath(__file__))
    os.chdir(directorysuper)
    
    
    image_list = []
    for filename in os.scandir(directory):
        if re.search(regex+'.fits',filename.name) or re.search(regex+'.fit',filename.name):
            image_list.append(filename.name)
            
    if(len(image_list)>1):
        return image_list
    elif(len(image_list) == 1):
        os.chdir(directory)
        return CCDData(fits.getdata(image_list[0]), unit='adu')

# Function to retrieve WCS from the header of a given file
def get_WCS(directory, filename):
    '''
    

    Parameters
    ----------
    directory : String
        The directory name where the files being fetched from are stored
    filename : String
        The name of the fits file that contains the wcs data

    Returns
    -------
    WCS Object
        The WCS object extracted from the file header.

    '''
    directorysuper = 'C:\\Users\jackp\OneDrive\Documents\Year 4\EP FYP'
    os.chdir(directorysuper)
    
    # Change the working directory to the location of the python script
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    
    os.chdir(directory)
    fn = get_pkg_data_filename(directory+'\\'+filename)
    f = fits.open(fn)
    
    # Extract and return the WCS Object from the header data of the file
    return WCS(f[0].header)
    

# Move to the directory where the file is located 
directorysuper = os.path.dirname(os.path.abspath(__file__))
os.chdir(directorysuper)
print('\nIn:',directorysuper)

# Fetch the required files by using the file_grabber() function
directory = 'Files/ProcessedFilesNew'
regex = 'CasA[0-9]*2007_WCS'
WCS2007 = file_grabber(directory, regex)

directory = 'Files/ProcessedFilesNew'
regex = 'CasA[0-9]*2022_WCS'
WCS2022 = file_grabber(directory, regex)

directory = 'Files/2007Data/SciFiles'
regex = 'MasterSci'
Data2007 = file_grabber(directory, regex)

directory = 'Files/ProcessedFiles'
regex = 'Stars2022'
Stars2022 = file_grabber(directory, regex)


# Change the working directory to the location of the python script
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Add WCS to the file being displayed
wcs2022 = get_WCS('Files/ProcessedFilesNew','CasA2022_WCS.fits')
WCS2022 = CCDData(WCS2022,unit='adu',wcs=wcs2022)

# Normalise the image 
interval = ZScaleInterval()
stretch = AsinhStretch(.85)
norm = ImageNormalize(WCS2022,interval=interval, stretch=stretch)


# Plot the image
ax = plt.subplot(projection=wcs2022)
ax.imshow(WCS2022, cmap='gray')
plt.xlim(300,750)
plt.ylim(300,750)
plt.xlabel('Right Ascension [°]')
plt.ylabel('Declination [°]')



# Plot the knots found in 2022
x2022 = [350.8371466666666, 350.82079625, 350.82887708333334, 
         350.86690208333334, 350.93376374999, 350.93997]
y2022 = [58.84604444444443, 58.84571805555555, 58.78896138888888, 
         58.7871033333, 58.8051972222, 58.80908861]
plt.scatter(x2022,y2022, s=50, edgecolor='red', facecolor='red', 
            transform = ax.get_transform('fk5'))

# Plot the knots found in 2007
x2007 = [350.83848916666676, 350.82300208333334, 350.83071291666664, 
         350.86628874999, 350.93119625, 350.9373354]
y2007 = [58.84447388888889, 58.84432527777777, 58.790007777777774, 
         58.7874608333, 58.8055413888, 58.80904444]
plt.scatter(x2007,y2007, s=50, edgecolor='black', facecolor='black', 
            transform = ax.get_transform('fk5'))

# Plot the Expansion center calculated by the paper
xcenter = 350.86693625 
ycenter = 58.8130525
ax.scatter(xcenter, ycenter, s=100, edgecolor='black', 
            facecolor='yellow', label='Center found', 
            transform = ax.get_transform('fk5'))



# Plot the knot trajectories for all knot pairs
x = np.linspace(-400,400,10000)
f = (lambda x: slope*x+c)
invert = 1

# Plot the trajectory of each knot pair
for num in range(0,len(x2007)):
    
    # Calculate the slope and constant of each line
    slope = invert*(y2022[num]-y2007[num])/(x2022[num]-x2007[num])
    c = y2022[num]-slope*x2022[num]
    
    print('\ny = '+str(slope)+'x + '+str(c))
    
    # Plot the lines
    ax.plot(x,f(x),label='Path of Knot '+str(num+1), transform = ax.get_transform('fk5'))

# Add a legend to the graph
pos = ax.get_position()
ax.set_position([pos.x0, pos.y0, pos.width * .6, pos.height])
#ax.legend(loc='center right', bbox_to_anchor=(2, .5),ncol=1)
plt.show()


print('\nDone!')






