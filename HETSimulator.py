# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 17:57:36 2024

@author: garan

Simulates the motion of an electron in the channel of a hall effect thruster
"""
import numpy as maths
from matplotlib import pyplot as plt 
import sys
from timeit import default_timer as timer
#fundamental constants
electronCharge = -1.6E-19
electronMass = 9.11E-31
mu_0 = 1.26E-6

def setup():
    #Magnetic Field
    print("First we are defining the magnetic field strength\n")
    temp = str(input("Do you want to define the electro magnets? [yes/no]\n"))
    B = [0,0,0]
    if temp.lower() == "yes":
        turns = int(input("Enter number of turns\n"))
        solenoidCurrent = float(input("Enter current through solenoid in amps\n"))
        solenoidLength = float(input("Enter length of the solenoid in cm\n")) * 10**-2
        B[2] = mu_0 * turns * solenoidCurrent / solenoidLength
        print("magnetic field strenght: " + str(B) + " T\n")
    else:
        B[2] = float(input("Enter the magnetic field strength in teslas\n"))
        
    #Electric Field 
    print("Now we define the electric field\n")
    temp = str(input("Do you wish to define the electrodes? [yes/no]\n"))
    E = [0,0,0]
    if temp.lower() == "yes":
        anodeVoltage = float(input("Enter the anode voltage in volts\n"))
        cathodeVoltage = float(input("Enter the cathode voltage in volts)\n"))
        electrodeSeperation = float(input("Enter the electrode seperation in cm\n")) * 10**-2
        E[1] = (anodeVoltage - cathodeVoltage)/electrodeSeperation
        print("electric field strenght: " + str(E) + " V/m\n")
    else:
        E[1] = float(input("Enter the electric field strength in V/m\n"))

    #HET parameters
    channelLength = 0
    channelMeanDiameter = 0
    electrodeSeperation = 0
    temp = str(input("Do you wish to define the HETs boundaries? [yes/no]\n"))
    if temp.lower() == "yes":
        channelLength = float(input("Enter the channel length in cm\n")) * 10**-2
        channelMeanDiameter = float(input("Enter the mean channel diameter in mm\n")) * 10**-3
        electrodeSeperation = float(input("Enter the cathode and anode seperation in cm\n")) * 10**-2
    
    #Initial conditions
    v = [0,0,0]
    print("Now you must define the initial velocity in m/s\n")
    v[0] = float(input("Enter the inital horozontal velcoity\n"))
    v[1] = float(input("Enter the initial vertical velocity\n"))

    #Simulation setup
    temp = str(input("Do you want to customise the running conditions of the simulation? [yes/no]\n"))
    simulationDuration = 10**-8
    timeStep = 10**-13
    maxRunTime = 600
    if temp.lower() == "yes":
        simulationDuration = float(input("Enter the duration that the electron will be simulated for (not real time) in seconds\n"))
        timeStep = float(input("Enter the simulation timeStep in seconds\n"))
        maxRunTime = float(input("Enter the duration that the simulation can run for (real time) in seconds\n"))
    
    run_simulation(E, B, v, simulationDuration, timeStep, maxRunTime, channelLength, channelMeanDiameter, electrodeSeperation)
    
def lorentz_force(E, B, v):
    F = maths.multiply(electronCharge, (E + maths.cross(v, B)))
    return F

def magnitude(x):
    return maths.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    
def run_simulation(E, B, v, duration, delta, maxTime, c, d, s):
    maxV = 0
    positions = [[], []]
    velocities = []
    pos = [0,0,0]
    a = [0,0,0]
    runTime = 0
    t0 = int(timer())
    tn = 0
    flag = False
    maxR = 0
    while runTime < duration and tn < maxTime and flag == False:
        positions[0].append(pos[0])
        positions[1].append(pos[1])
        velocities.append(v)
        Eu = E
        Bu = B
        
        if c != 0 or d != 0 or s != 0:
            if pos[1] < -c:
                flag = True
            if pos[1] > (c * 1.1):
                Bu = [0,0,0]
            if pos[1] > (s-c):
                Eu = [0,0,0]
         
        F = lorentz_force(Eu, Bu, v)
        a = F/electronMass

        #print("velocity magnitude: " + str(magnitude(v)) + " acceleration magnitude: " + str(magnitude(a)) + " dot: " + str(maths.dot(a, v)))
        v = v + maths.multiply(a, delta)
        pos = pos + maths.multiply(v, delta) - 0.5 * maths.multiply(a, delta**2)
        
        #print("accerleration: " + str(a) + " V: " + str(v) + " dot: " + str(maths.dot(v, F)))
        runTime = runTime + delta
        tn = int(timer()) - t0
        print(str(tn) + "s progress: " + str(maths.round(runTime/duration * 100, 3)) + "%")
        
        if magnitude(v) > maxV:
            maxV = magnitude(v)
        if maths.abs(pos[1])/2 > maxR:
            maxY = maths.abs(pos[1])/2
    tN = timer()

    print("Percent difference in calculations: " + str(calculate_error(velocities[0], B, maxY) * 100) + "%")
    print("Processing time in seconds: " + str(tN-t0))
    print("maximum velocity: " + str(maxV))
    plot(positions[0], positions[1])
    
def calculate_error(v, B, maxR):
    r = electronMass * magnitude(v)/(maths.abs(electronCharge) * magnitude(B))
    #print(r)
    err = maths.abs(maxR - r)/r
    return err

def plot(x, y):
    plt.title("Electron Motion through a HET")
    plt.xlabel("X position/m")
    plt.ylabel("Y position/m")
    gradient = maths.linspace(-1, 1, len(x))
    dotColour = maths.tan(gradient)
    plt.scatter(x, y, c=dotColour, s=(maths.exp(-len(x)/1000)+0.1))
    plt.show()

setup()


