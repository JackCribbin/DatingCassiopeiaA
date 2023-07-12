# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:38:48 2022

@author: jackp
"""

# Code written to perform reduction of images taken at OHP.
import numpy as np
import ccdproc

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.visualization import ZScaleInterval, ImageNormalize, AsinhStretch
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename

import os
import matplotlib.pyplot as plt
import re

from scipy.stats import mode 

#import twirl

import warnings
warnings.filterwarnings("ignore")

# Function to get a list of files matching a given regular expression in 
# a given directory
def file_grabber(regex):
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
    
    directory = os.getcwd()
    image_list = []
    for filename in os.scandir(directory):
        if re.search(regex+'.fits',filename.name) or re.search(regex+'.fit',filename.name):
            image_list.append(filename.name)
            
    if(len(image_list)>1):
        return image_list
    elif(len(image_list) == 1):
        return image_list[0]
    
  
def array_padder(data, size):
    padded = np.empty((size,size))
    print(padded[0,0])
    dif = int((size-len(data))/2)
    count = 0

    for i in range(0,size):
        if(i < dif):
            count = count + 1
        elif(i < size-dif):
            count = 0
            for j in range(dif, len(data)+dif):                
                padded[i,j] = data[i-dif,j-dif]
    return padded

def array_shortener(ar, size):
    '''
    Function to shorten 2D array by 1/9

    Parameters
    ----------
    ar : Square numpy array
        The array to be shortened
    num : int
        The size the array is desired to be shortened to
        

    Returns
    -------
    ar : numpy array
        The shortened array

    '''
    
    num = int((1-((len(ar)-size)/len(ar)))*10)
    length = len(ar)-1
    for i in range(0,length):
        if(i<length):
            if(i%num == 0):
                ar = np.delete(ar, i-1,1)
                ar = np.delete(ar, i-1,0)
                length = len(ar)-1
        
    return ar
            
            
# Move to the directory where the file is located 
directorysuper = os.path.dirname(os.path.abspath(__file__))
os.chdir(directorysuper)
print('\nIn:',directorysuper)

directory = 'Files/Night3Data'
os.chdir(directory)
CasA2022 = file_grabber('MasterSci*')
CasA2022 = CCDData(fits.getdata(CasA2022), unit='adu')

os.chdir(directorysuper)
directory = 'Files/2007Data/SciFiles'
os.chdir(directory)
CasA2007 = file_grabber('MasterSci*')
CasA2007 = CCDData(fits.getdata(CasA2007), unit='adu')

os.chdir(directorysuper)
directory = 'Files/ProcessedFiles'
os.chdir(directory)
WCS2007 = file_grabber('CasA[0-9]*2007_WCS')
WCS2007 = CCDData(fits.getdata(WCS2007), unit='adu')

os.chdir(directorysuper)
directory = 'Files/ProcessedFiles'
os.chdir(directory)
WCS2022 = file_grabber('CasA[0-9]*2022_WCS')
WCS2022 = CCDData(fits.getdata(WCS2022), unit='adu')





### Code to focus on Cas A and zoom in ###########
from astropy.nddata import Cutout2D
position = (530,500)
size = (500,500)

print(CasA2007.shape,'to',Cutout2D(CasA2007, position, size).shape)

WCS2007 = np.asarray(Cutout2D(np.asarray(WCS2007), position, size).data)

WCS2007 = array_shortener(WCS2007, 450)

WCS2007 = np.roll(np.asarray(WCS2007), -6, axis=0)
WCS2007 = np.roll(np.asarray(WCS2007), -1, axis=1)

WCS2007 = CCDData(np.asarray(WCS2007), unit='adu')

position = (530,500)
size = (450,450)

print(WCS2022.shape,'to',Cutout2D(WCS2022, position, size).shape)
WCS2022 = CCDData(np.asarray(Cutout2D(np.asarray(WCS2022), position, size).data), unit='adu')

#WCS2007 = CCDData(array_padder(np.asarray(WCS2007), 500), unit='adu')

########### end of zoom and focus code ###########

'''
import cv2
import numpy as np

#WCS2022 = cv2.resize(np.asarray(WCS2022), dsize=(500,500), interpolation=cv2.INTER_LINEAR)

#WCS2022 = CCDData(np.asarray(WCS2022), unit='adu')


#WCSFiles = ccdproc.combine([WCS2007,WCS2022], None, 'median', unit='adu')
'''


print(WCS2022.shape)
print(WCS2007.shape)


CasA2007 = CCDData(np.flip(np.asarray(CasA2007), 0), unit='adu')
CasA2007 = CCDData(np.roll(np.asarray(CasA2007), -50, axis=0), unit='adu')
CasA2007 = CCDData(np.roll(np.asarray(CasA2007), 115, axis=1), unit='adu')



'''

from astropy.nddata import Cutout2D
position = (530,500)
size = (500,500)

print(CasA2007.shape,'to',Cutout2D(CasA2007, position, size).shape)

CasA2007 = np.asarray(Cutout2D(np.asarray(CasA2007), position, size).data)
CasA2007 = CCDData(np.asarray(CasA2007), unit='adu')

position = (530,490)
size = (450,450)
CasA2022 = CCDData(np.asarray(Cutout2D(np.asarray(CasA2022), position, size).data), unit='adu')


#CasA2022 = array_padder(CasA2022, 500)
#CasA2022 = CCDData(np.asarray(CasA2022), unit='adu')


print(CasA2022.shape,'vs',CasA2007.shape)

'''




#final_image = ccdproc.combine([CasA2022,CasA2007], None, 'median', unit='adu')

final_image = ccdproc.combine([WCS2022,WCS2007], None, 'median', unit='adu')






# Initialise normalisaton for the image 
interval = ZScaleInterval()
stretch = AsinhStretch(.85)
norm = ImageNormalize(CasA2022,interval=interval, stretch=stretch)



os.chdir(directorysuper)
directory = 'Files/ProcessedFiles'
os.chdir(directory)
#wcs = file_grabber('CasA2022_WCS')


fn = get_pkg_data_filename(directory+'/'+'CasA2022_WCS.fits')
f = fits.open(fn)
wcs = WCS(f[0].header)




# Plot the image


plt.subplot(projection=wcs)
plt.imshow(CasA2022, norm = norm)
plt.colorbar()
plt.show()


'''
stretch = AsinhStretch(.85)
norm = ImageNormalize(CasA2007,interval=interval, stretch=stretch)
plt.figure()
plt.imshow(CasA2007, origin='upper', norm = norm)
plt.colorbar()
#plt.show()

'''


############ Code to save files

os.chdir(directorysuper)
os.chdir('Files/ProcessedFiles')
#WCS2022.write('CasA2022mini_WCS.fits', overwrite=True)
#WCS2007.write('CasA2007mini_WCS.fits', overwrite=True)
#final_image.write('CasACombined.fits', overwrite=True)

############




norm = ImageNormalize(WCS2022,interval=interval, stretch=stretch)
plt.subplot(projection=wcs)
plt.imshow(WCS2022, norm = norm)
plt.colorbar()
plt.show()


norm = ImageNormalize(WCS2007,interval=interval, stretch=stretch)
plt.figure()
plt.imshow(WCS2007, origin='upper', norm = norm)
plt.colorbar()
plt.show()




norm = ImageNormalize(final_image,interval=interval, stretch=stretch)
plt.figure()
plt.imshow(final_image, origin='upper', norm = norm)
plt.colorbar()
plt.show()



'''
norm = ImageNormalize(WCSFiles,interval=interval, stretch=stretch)
plt.figure()
plt.imshow(WCSFiles, origin='upper', norm = norm)
plt.colorbar()
plt.show()
'''
print('\nDone!')






