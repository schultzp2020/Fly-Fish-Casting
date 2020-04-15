# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp
import math

numberOfPoints = int(input("How many points do you want to measure? "))
initialXPos = int(input("What is the rod's starting x position? "))
initialYPos = int(input("What is the rod's starting y position? "))

phiList = []
massList = []
lengthList = []

for i in range(numberOfPoints):
    phiList.append(int(input("What is your angle (in degrees) at point " + str(i + 1) + "? ")))
    
for i in range(numberOfPoints):
    massList.append(int(input("What is your mass (in grams) at point " + str(i + 1) + "? ")))
    
for i in range(numberOfPoints):
    lengthList.append(int(input("What is your length (in mm) at point " + str(i + 1) + "? ")))
    
xPosition = []
yPosition = []

for i in range(numberOfPoints):
    xPosition.append(initialXPos)
    yPosition.append(initialYPos)
    for j in range(i + 1):
        if(i == j):
            xPosition[i] = xPosition[i] + lengthList[j] * 0.5 * np.cos(math.radians(phiList[j]))
            yPosition[i] = yPosition[i] + lengthList[j] * 0.5 * np.sin(math.radians(phiList[j]))
        else:
            xPosition[i] = xPosition[i] + lengthList[j] * np.cos(math.radians(phiList[j]))
            yPosition[i] = yPosition[i] + lengthList[j] * np.sin(math.radians(phiList[j]))
        