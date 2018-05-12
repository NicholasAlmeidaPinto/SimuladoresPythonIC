# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 08:31:05 2017

@author: Nicholas de Almeida Pinto
"""
import numpy as np

#'Clear' the console for beauty
def clc(): print('\n'*20)#'Clear' IDLE
#-----------------------------------------------------------------------------#


#Function to input initial values
def Initial_values(): 
    """Function to user input values, return N, deltx, delty, deltt, mio, miw,
    time"""
    clc()
    N = 30
    deltx = 0.003#0.003#5
    delty = 0.003#0.003#5
    deltt = 0.0005#0.0005#1
    mio = 1.
    miw = 1.
    time = int(input("Time: "))
    print("Injection")
    print("1 - Well      2 -  Line")
    print("3 - FiveSpot  4 - Simulator test")
    injection = int(input("==> "))
    
    clc()
    print("1 - First TVD                2 - TVD complete")
    print("3 - Original NonStandard     4 - NonStandard")
    print("5 - NewSimulator             6 - NewTVD")
    print("7 - TVD_NonStandard")
    simulator = int(input("==> "))
    clc()
    print("Number of cells: ", N)
    print("Deltx: ", deltx, "          Delty: ", delty)
    print("Deltt: ", deltt, "         Oil visc.: ", mio)
    print("Water visc.: ", miw, "      Time: ", time)
    choice = 2#int(input("1 to change value above: "))
    while choice == 1:
        N = int(input("Number of cells: "))
        deltx = float(input("Deltx: "))
        delty = float(input("Delty: "))
        deltt = float(input("Deltt: "))
        mio = float(input("Oil viscosity: "))
        miw = float(input("Water viscosity: "))
        time = int(input("Time: "))
        
        print("Number of cells: ", N)
        print("Deltx: ", deltx, "          Delty: ", delty)
        print("Deltt: ", deltt, "         Oil visc.: ", mio)
        print("Water visc.: ", miw, "     Time: ", time)
        choice = int(input("1 to change value above: "))
    return N, deltx, delty, deltt, mio, miw, time, injection, simulator
#-----------------------------------------------------------------------------#
    
def Loading(time, j):#Only to show progress 
    if j%10 == 0:
        print(j*100/time, "%")
#-----------------------------------------------------------------------------#
        
def Porosity(N):
    """Function to define matrix of porosity, Porosity(N), which N is number of
    cells"""
    print("Value porosity, bigger than 1 to random")
    choice = 1#float(input("==> "))
    if choice > 1.:
        porosity = np.random.rand(N, N)
    else:
        porosity = np.empty((N,N))
        porosity[:,:] = choice
    return porosity
#-----------------------------------------------------------------------------#

def Injection_well(Sw, j, N, injection):
    if injection == 1:#Well
        Sw[int(N/2),int(N/2),:] = 1 
    elif injection == 2:#Line
        Sw[:,0,:] = 1
    elif injection == 3:#FiveSpot
        Sw[0,0,:] = 1
        Sw[0,N-1,:] = 1
        Sw[N-1,0,:] = 1
        Sw[N-1,N-1,:] = 1
    elif injection == 4:#Simulator's test
        Sw[int(N/2-5):int(N/2+5),int(N/2-5):int(N/2+5),0] = 1 
    return Sw
#-----------------------------------------------------------------------------#
    
def Graphic(Sw,N,time): #Function that show the graphic 
    from mpl_toolkits.mplot3d import axes3d
    import matplotlib.pyplot as plt
    from matplotlib import cm

    X = np.linspace(0, N, N)
    Y = np.linspace(0, N, N)
    X, Y = np.meshgrid(X, Y)
    
    fig = plt.figure(figsize=plt.figaspect(0.5))

    #SURFACE
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_surface(X, Y, Sw[:,:,time-1], cmap=cm.coolwarm, linewidth=1, antialiased=True)
    # Add a color bar which maps values to colors.
    ax.set_zlim(0, 1)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    #ANIMATION
    an = fig.add_subplot(1, 2, 2, projection='3d')
    # Set the z axis limits so they aren't recalculated each frame.
    an.set_zlim(0, 1)
    fig.colorbar(surf, shrink=0.8, aspect=5)
    # Begin plotting.
    wframe = None
    for j in range(time):
        # If a line collection is already remove it before drawing.
        if wframe:
            an.collections.remove(wframe)
        # Plot the new wireframe and pause briefly before continuing.
        wframe = an.plot_surface(X, Y, Sw[:,:,j],  cmap=cm.coolwarm, linewidth=1, antialiased=True)
        an.set_xlabel('Tempo: ' + str(j))
        plt.pause(.01)
