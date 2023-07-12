# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:07:12 2022

@author: Jack Cribbin - 19328253

This file takes in the coordinates of any knots identified, and finds
the trajectory of each pair and from these estimates the expansion centre
of the supernova remnant
"""

import matplotlib.pyplot as plt
import numpy as np


# Plot the image
fig, ax = plt.subplots()
ax.set_xlabel('Right Ascension [°]')
ax.set_ylabel('Declination [°]')
ax.set_xlim(-350.95,-350.815)
ax.set_ylim(58.78,58.85)

# Plot the knots found in 2022
x2022 = [350.8371466666666, 350.82079625, 350.82887708333334, 
         350.86690208333334, 350.93376374999, 350.93997]
x2022 = [x * -1 for x in x2022]
y2022 = [58.84604444444443, 58.84571805555555, 58.78896138888888, 
         58.7871033333, 58.8051972222, 58.80908861]
ax.scatter(x2022,y2022, s=100, edgecolor='red', facecolor='red')

# Plot the knots found in 2007
x2007 = [350.83848916666676, 350.82300208333334, 350.83071291666664, 
         350.86628874999, 350.93119625, 350.9373354]
x2007 = [x * -1 for x in x2007]
y2007 = [58.84447388888889, 58.84432527777777, 58.790007777777774, 
         58.7874608333, 58.8055413888, 58.80904444]
ax.scatter(x2007,y2007, s=100, edgecolor='black', facecolor='black')


# Plot the knot trajectories for all knot pairs
x = np.linspace(-400,400,10000)
f = (lambda x: slope*x+c)
invert = 1
    
# Plot the trajectory of each knot pair
for num in range(0,len(x2007)):
    
    # Calculate the slope and constant of each line
    slope = invert*(y2022[num]-y2007[num])/(x2022[num]-x2007[num])
    c = y2022[num]-slope*x2022[num]
    
    # Plot the lines
    ax.plot(x,f(x),label='Path of Knot '+str(num+1))


# Get the intercepts of each line
intercepts = [[350.87187,58.81346],
              [350.87191,58.81349],
              [350.87182,58.81350],
              [350.86725,58.81083],
              [350.86974,58.80791],
              [350.86189,58.80778],
              [350.86408,58.81454],
              [350.85693,58.82291]]

# Get the average of the x and y position of each intercept to get
# the expansion center
sumx = 0
sumy = 0
for i in range(0, len(intercepts)):
    sumx = sumx + (intercepts[i][0])
    sumy = sumy + (intercepts[i][1])
    
     #ax.scatter(-intercepts[i][0], intercepts[i][1], s=100, edgecolor='black', 
     #           facecolor='yellow', label='Center found')
    
center = [sumx/(len(intercepts)), sumy/(len(intercepts))]

print('\nCenter is:',center)

'''
# Plot the expansion center found 
ax.scatter(-center[0], center[1], s=100, edgecolor='black', 
            facecolor='yellow', label='Center found')
'''
# Add a legend to the graph
pos = ax.get_position()
ax.set_position([pos.x0, pos.y0, pos.width * 0.6, pos.height])
ax.legend(loc='center right', bbox_to_anchor=(1.65, .5),ncol=1)
plt.show()








