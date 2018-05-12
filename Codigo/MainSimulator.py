# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 11:30:08 2017

@author: Nicholas de Almeida Pinto
"""

#libraries
from Inputs import *
from Method import *
import numpy as np
import time as tm

                       ##################################
                       #Define variables & start program#
                       ##################################

StartTime = tm.time() #save the time of the beggining of program

N, deltx, delty, deltt, mio, miw, time, injection, simulator = Initial_values()
Sw = np.zeros((N,N, time+1))
porosity = Porosity(N)

#Start running time
for j in range(time):
    Loading(time,j)
    Sw = Injection_well(Sw, j, N, injection)
    
    
    #############METHODS##############
    if simulator == 1:
        Sw = First_TVD(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 2:
        Sw = TVD(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 3:
        Sw = Original_NonStandard(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 4:
        Sw = NonStandard(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 5:
        Sw = Test_New_Simulator(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 6:
        Sw = NewTVD(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
    elif simulator == 7:
        Sw = TVD_NonStandard_Not_Both_Directions(Sw, N, mio, miw, j, porosity, deltt, deltx, delty)
        
volume  = sum(sum(Sw[:,:,time]))
Graphic(Sw,N,time)

#Show the total time to run the program
EndTime = tm.time()
print("Time: ",EndTime - StartTime)
print("Total volume: ", volume)
del EndTime, StartTime, j, injection, simulator
                       ##################################
                       ########End of the program########
                       ##################################
