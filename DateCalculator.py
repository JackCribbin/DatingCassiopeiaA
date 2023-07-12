# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 17:32:47 2022

@author: Jack Cribbin - 19328253

This file takes in the knot locations and expansion centre, and from these
calculates the date at which the supernova occurred
"""

from numpy import cos,sin,pi,arccos
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt

def RAtoRad(a):
    return ((a[0] + a[1]/60 + a[2]/3600) * 15) * pi/180

def DectoRad(d):
    return (d[0] + d[1]/60 + d[2]/3600) * pi/180

def RadtoDeg(x):
    return x * 180/pi

def DegtoRad(x):
    return x * pi/180

# The locations of the 6 knots in 2022 and 2007
x2022 = [350.8371466666666, 350.82079625, 350.82887708333334, 
         350.86690208333334, 350.93376374999, 350.93997]
y2022 = [58.84604444444443, 58.84571805555555, 58.78896138888888, 
         58.7871033333, 58.8051972222, 58.80908861]

x2007 = [350.83848916666676, 350.82300208333334, 350.83071291666664, 
         350.86628874999, 350.93119625, 350.9373354]
y2007 = [58.84447388888889, 58.84432527777777, 58.790007777777774, 
         58.7874608333, 58.8055413888, 58.80904444]
    
# Convert the knot coordinates from degrees to radians
for i in range(len(x2022)): 
    x2022[i] = DegtoRad(x2022[i])
    y2022[i] = DegtoRad(y2022[i])
    x2007[i] = DegtoRad(x2007[i])
    y2007[i] = DegtoRad(y2007[i])
    
# Calculate the distance between each knot pair in radians
dists = []
for i in range(len(x2022)): 
    dists.append(arccos(sin(y2022[i])*sin(y2007[i]) + 
                    cos(y2022[i])*cos(y2007[i])*cos(x2022[i]-x2007[i])))
    
    print('Distance is:','%.4E' % Decimal(str(dists[i])),'radians')
print()

# Distance from Earth to Cas A in meters
d = 1.0492e+20

# Find the distance travelled by the knots in meters
for i in range(len(x2022)): 
    dists[i] = dists[i] * d
    print('Distance is:','%.4E' % Decimal(str(dists[i])),'[m]')
print()


# Find the velocity of the knots in m/s
vels = []
for i in range(len(x2022)): 
    vels.append(((dists[i])/((2022.833-2007.75))))
    print('Velocity [m s-1]:', '%.4E' % Decimal(str(vels[i])))
print()

avgvel=[]
# Print the velocities to km/y
for i in range(len(x2022)): 
    if(i!=3):
        avgvel.append(vels[i]/1000/3.154e+7)
    print('Velocity [km y-1]:', '%.4E' % Decimal(str((vels[i]/1000/3.154e+7))))
print()

# Calculate the average velocity of the 5 knots
avgvel = np.average(avgvel)
print(avgvel)

# Plot the velocity for each knot
x = np.linspace(0,len(vels)-1,len(vels))+1
plt.figure()
plt.scatter(x,vels)
plt.xlabel('Knot Number')
plt.ylabel('Knot Velocity [m s\u207B\u00b9]')
plt.title('Plot of Knot Velocities')
plt.show()


# Convert the expansion center coordinates to radians
center = [350.86955555555556, 58.811]
center[0] = DegtoRad(center[0])
center[1] = DegtoRad(center[1])
print('\nCenter:',center)
print()

# Find the distance each knot travelled from the expansion center in radians
travel = []
for i in range(len(x2007)):
    a1 = DegtoRad(x2007[i])
    d1 = DegtoRad(y2007[i])
    
     
    travel.append(arccos(sin(center[1])*sin(y2007[i]) + 
                    cos(center[1])*cos(y2007[i])*cos(center[0]-x2007[i])))
    
    print('Distance is:','%.4E' % Decimal(str(travel[i])),'radians')
print()


#Convert the distance travelled to meters
for i in range(len(x2022)): 
    travel[i] = travel[i] * d 
    print('Distance is:','%.4E' % Decimal(str(travel[i])),'[m]')
print()

# Find the time each knot spent travelling from center to current location
time = []
for i in range(len(x2022)): 
    vels[i] = vels[i]
    time.append(travel[i] / vels[i])
    print('Time Travelled is:',time[i],'[y]')
print()


# Find the date that the supernova occured
date = []
for i in range(len(x2022)): 
    date.append(2007.75 - time[i])
    print('Date of Supernova is:',date[i],'[y]')
print()

# Calculate the average date the supernova occured
avg = 0
count = 0
for i in range(len(x2022)): 
    if(i != 3):
        avg = avg + date[i]
        count += 1
avg = avg/count
print('The average Supernova Date is:',avg)

# Plot the velocity for each knot
x = np.linspace(0,len(vels)-1,len(vels))+1
plt.figure()
plt.scatter(x,time)
plt.xlabel('Knot Number')
plt.ylabel('Knot Time Travelling [years]')
plt.title('Plot of Time Knots have been Travelling')
plt.show()









