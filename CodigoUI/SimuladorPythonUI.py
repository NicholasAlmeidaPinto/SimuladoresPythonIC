from tkinter import *
import numpy as np

def mainSimulator(N, time, simulator, deltx = 0.003, delty = 0.003, 
                      deltt = 0.0005, mio = 1, miw = 1, injection = 4):
    Sw = np.zeros((N,N, time+1))
    porosity = Porosity(N)
    
    #Start running time
    for j in range(time):
#        Loading(time,j)
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
            
#    print('Volume:' + sum(sum(Sw[:,:,time])))
    Graphic(Sw,N,time)
#-----------------------------------------------------------------------------#
    
#def Loading(time, j):#Only to show progress 
#    if j%10 == 0:
#        has=str(j*100/time)
#        print(has, "%")
#-----------------------------------------------------------------------------#
        
def Porosity(N, choice=1):
    """Function to define matrix of porosity, Porosity(N), which N is number of
    cells"""
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
    surf = ax.plot_surface(X, Y, Sw[:,:,time-1], cmap=cm.inferno, linewidth=1, antialiased=True, alpha=0.7)
    # Add a color bar which maps values to colors.
    ax.set_zlim(0, 1)
    fig.colorbar(surf, shrink=0.5, aspect=5)

    #ANIMATION
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    # Set the z axis limits so they aren't recalculated each frame.
    ax.set_zlim(0, 1)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    # Begin plotting.
    wframe = None
    for j in range(time):
        # If a line collection is already remove it before drawing.
        if wframe:
            ax.collections.remove(wframe)
        # Plot the new wireframe and pause briefly before continuing.
        wframe = ax.plot_surface(X, Y, Sw[:,:,j],  cmap=cm.inferno, linewidth=1, antialiased=True, alpha=0.7)
        plt.pause(.05)

###############################################################################
####################################METODOS####################################
###############################################################################
#Function that defines and calculates the TVD's method
def First_TVD(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    for i in range(N):
        for k in range(N):
            if i == 0:
                Fax=Sw[i,k,j]
                Fbx=Sw[i,k,j]
            elif i == 1:
                Fax=Sw[i-1,k,j]
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
            elif i == N-1:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]
            else:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))

            if k == 0:
                Fay=Sw[i,k,j]
                Fby=Sw[i,k,j]
            elif k == 1:
                Fay=Sw[i,k-1,j]
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))
            elif k == N-1:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]
            else:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))
           
            Fx = Safe((Fax**2),(Fax**2 + (miw/mio)*((1-Fax)**2))) - Safe((Fbx**2),(Fbx**2 + (miw/mio)*((1-Fbx)**2)))
            Fy = Safe((Fay**2),(Fay**2 + (miw/mio)*((1-Fay)**2))) - Safe((Fby**2),(Fby**2 + (miw/mio)*((1-Fby)**2)))
            
            Sw[i,k,j+1] =  Sw[i,k,j] + ((Fx/deltx) + (Fy/delty))*deltt/porosity[i,k]
    return Sw
#-----------------------------------------------------------------------------#
    

#Function that defines and calculates the TVD's method
def TVD(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    
    """This simulator is different from the previous 'First_TVD' because the is
    a 'reverse copy' of all the math, this means that, while the for makes the
    simulator advance from 0 to N, there are the same math doing from N to 0, 
    using the oposite signal inside the index."""

    
    for i in range(N):
        for k in range(N):
            if i == 0:
                Fax=Sw[i,k,j]
                Fbx=Sw[i,k,j]
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]
            elif i == 1:
                Fax=Sw[i-1,k,j]
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
            elif i == N-2:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))                
            elif i == N-1:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]
                
                Faxn=Sw[i,k,j]
                Fbxn=Sw[i,k,j]
            else:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
    
    
            if k == 0:
                Fay=Sw[i,k,j]
                Fby=Sw[i,k,j]
                
                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]
            elif k == 1:
                Fay=Sw[i,k-1,j]
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))
                
                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
            elif k == N-2:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))

                Fayn=Sw[i,k+1,j]
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))                
            elif k == N-1:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]
                
                Fayn=Sw[i,k,j]
                Fbyn=Sw[i,k,j]
            else:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))

                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
                
                
            if i < N-1 and i > 0 and Sw[i,k,j] < Sw[i-1,k,j]:
                Faxn=0
                Fbxn=0
            elif i > 0 and i < N-1:
                Fax=0
                Fbx=0
                
            if k < N-1 and k > 0 and Sw[i,k,j] < Sw[i,k-1,j]:
                Fayn=0
                Fbyn=0
            elif k > 0 and k < N-1:
                Fay=0
                Fby=0 
                
                
            Fx = Conta(Fax,miw,mio) - Conta(Fbx,miw,mio) + Conta(Faxn,miw,mio) - Conta(Fbxn,miw,mio)
            Fy = Conta(Fay,miw,mio) - Conta(Fby,miw,mio) + Conta(Fayn,miw,mio) - Conta(Fbyn,miw,mio)
            Sw[i,k,j+1]=Sw[i,k,j]+((Fx/deltx)+(Fy/delty))*deltt/porosity[i,k]
    return Sw
#-----------------------------------------------------------------------------#

def Original_NonStandard(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    
    import math #uses to calculate exp
    alpha = 2
    phix=(1/(2*alpha))*(1-math.exp(-alpha*deltt/deltx))
    phiy=(1/(2*alpha))*(1-math.exp(-alpha*deltt/delty))
    for i in range(N):
        for k in range(N):

            if i==0:
                Fax=(Sw[i+1,k,j]**2)/((Sw[i+1,k,j]**2)+(miw/mio)*((1-Sw[i+1,k,j])**2)) 
                Gx=phix*(alpha*(Sw[i+1,k,j]-2*Sw[i,k,j]+Sw[i,k,j])-Fax)*.5/porosity[i,k]
            elif i==N-1:
                Fax=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2))
                Gx=phix*(alpha*(Sw[i,k,j]-2*Sw[i,k,j]+Sw[i-1,k,j])+Fax)*.5/porosity[i,k]
            else:
                Fax=(Sw[i+1,k,j]**2)/((Sw[i+1,k,j]**2)+(miw/mio)*((1-Sw[i+1,k,j])**2)) 
                Gx=phix*(alpha*(Sw[i+1,k,j]-2*Sw[i,k,j]+Sw[i-1,k,j])-Fax)*.5/porosity[i,k]

            if k==0:
                Fay=(Sw[i,k+1,j]**2)/((Sw[i,k+1,j]**2)+(miw/mio)*((1-Sw[i,k+1,j])**2)) 
                Gy=phiy*(alpha*(Sw[i,k+1,j]-2*Sw[i,k,j]+Sw[i,k,j])-Fay)*.5/porosity[i,k]
            elif k==N-1:
                Fay=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2))
                Gy=phiy*(alpha*(Sw[i,k,j]-2*Sw[i,k,j]+Sw[i,k-1,j])+Fay)*.5/porosity[i,k]
            else:
                Fay=(Sw[i,k+1,j]**2)/((Sw[i,k+1,j]**2)+(miw/mio)*((1-Sw[i,k+1,j])**2)) 
                Gy=phiy*(alpha*(Sw[i,k+1,j]-2*Sw[i,k,j]+Sw[i,k-1,j])-Fay)*.5/porosity[i,k]

            Sw[i,k,j+1]=Sw[i,k,j]+Gx+Gy
    return Sw
#-----------------------------------------------------------------------------#

def NonStandard(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    """This simulator is called NonStandard. The principal changes comparing 
    from TVD: for each value of saturation, is needed the previous and next value,
    this is easier to do than TVD because only uses one previous or next cell."""
    
    import math #uses to calculate exp
    alpha = 2
    phix=(1/(2*alpha))*(1-math.exp(-alpha*deltt/deltx))
    phiy=(1/(2*alpha))*(1-math.exp(-alpha*deltt/delty))
    for i in range(N):
        for k in range(N):

            if i==0:
                Fax=(Sw[i+1,k,j]**2)/((Sw[i+1,k,j]**2)+(miw/mio)*((1-Sw[i+1,k,j])**2)) 
                Fbx=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2)) #Pay attention to lines below with comments, probabily the correct value is 0
                Gx=phix*(alpha*(Sw[i+1,k,j]-2*Sw[i,k,j]+Sw[i,k,j])-Fax+Fbx)*.5/porosity[i,k]
            elif i==N-1:
                Fax=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2)) #this one
                Fbx=(Sw[i-1,k,j]**2)/((Sw[i-1,k,j]**2)+(miw/mio)*((1-Sw[i-1,k,j])**2))  
                Gx=phix*(alpha*(Sw[i,k,j]-2*Sw[i,k,j]+Sw[i-1,k,j])+Fax-Fbx)*.5/porosity[i,k]
            else:
                Fax=(Sw[i+1,k,j]**2)/((Sw[i+1,k,j]**2)+(miw/mio)*((1-Sw[i+1,k,j])**2)) 
                Fbx=(Sw[i-1,k,j]**2)/((Sw[i-1,k,j]**2)+(miw/mio)*((1-Sw[i-1,k,j])**2))  
                if  Sw[i,k,j]>Sw[i+1,k,j]:
                    Gx=phix*(alpha*(Sw[i+1,k,j]-2*Sw[i,k,j]+Sw[i-1,k,j])-Fax+Fbx)*.5/porosity[i,k]
                else:
                    Gx=phix*(alpha*(Sw[i+1,k,j]-2*Sw[i,k,j]+Sw[i-1,k,j])+Fax-Fbx)*.5/porosity[i,k]

            if k==0:
                Fay=(Sw[i,k+1,j]**2)/((Sw[i,k+1,j]**2)+(miw/mio)*((1-Sw[i,k+1,j])**2)) 
                Fby=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2)) # this one
                Gy=phiy*(alpha*(Sw[i,k+1,j]-2*Sw[i,k,j]+Sw[i,k,j])-Fay+Fby)*.5/porosity[i,k]
            elif k==N-1:
                Fay=(Sw[i,k,j]**2)/((Sw[i,k,j]**2)+(miw/mio)*((1-Sw[i,k,j])**2)) #this one too
                Fby=(Sw[i,k-1,j]**2)/((Sw[i,k-1,j]**2)+(miw/mio)*((1-Sw[i,k-1,j])**2)) 
                Gy=phiy*(alpha*(Sw[i,k,j]-2*Sw[i,k,j]+Sw[i,k-1,j])+Fay-Fby)*.5/porosity[i,k]
            else:
                Fay=(Sw[i,k+1,j]**2)/((Sw[i,k+1,j]**2)+(miw/mio)*((1-Sw[i,k+1,j])**2)) 
                Fby=(Sw[i,k-1,j]**2)/((Sw[i,k-1,j]**2)+(miw/mio)*((1-Sw[i,k-1,j])**2)) 
                if Sw[i,k,j]>Sw[i,k+1,j]:
                    Gy=phiy*(alpha*(Sw[i,k+1,j]-2*Sw[i,k,j]+Sw[i,k-1,j])-Fay+Fby)*.5/porosity[i,k]
                else:
                    Gy=phiy*(alpha*(Sw[i,k+1,j]-2*Sw[i,k,j]+Sw[i,k-1,j])+Fay-Fby)*.5/porosity[i,k]

            Sw[i,k,j+1]=Sw[i,k,j]+Gx+Gy
    return Sw
#-----------------------------------------------------------------------------#
    
def Test_New_Simulator(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    Fax = 0.
    Fbx = 0.
    Fay = 0.
    Fby = 0.
    Faxn = 0.
    Fbxn = 0.
    Fayn = 0.
    Fbyn = 0.
    Fy  = 0.
    Fx  = 0.
    import math #uses to calculate exp
    alpha = 2
    phix=(1/(2*alpha))*(1-math.exp(-alpha*deltt/deltx))
    phiy=(1/(2*alpha))*(1-math.exp(-alpha*deltt/delty))
    
    for i in range(N):
        for k in range(N):
            if i == 0:
                Fax=Sw[i,k,j]
                Fbx=Sw[i,k,j]
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]
            elif i == 1:
                Fax=Sw[i-1,k,j]
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
            elif i == N-2:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))                
            elif i == N-1:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]
                
                Faxn=Sw[i,k,j]
                Fbxn=Sw[i,k,j]
            else:
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))
                
                Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
    
    
            if k == 0:
                Fay=Sw[i,k,j]
                Fby=Sw[i,k,j]
                
                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]
            elif k == 1:
                Fay=Sw[i,k-1,j]
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))
                
                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
            elif k == N-2:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))

                Fayn=Sw[i,k+1,j]
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))                
            elif k == N-1:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]
                
                Fayn=Sw[i,k,j]
                Fbyn=Sw[i,k,j]
            else:
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))

                Fayn=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))      
                
            if i < N-1 and i > 0 and Sw[i,k,j] < Sw[i-1,k,j]:
                Faxn=0
                Fbxn=0
            elif i > 0 and i < N-1:
                Fax=0
                Fbx=0
                
            if k < N-1 and k > 0 and Sw[i,k,j] < Sw[i,k-1,j]:
                Fayn=0
                Fbyn=0
            elif k > 0 and k < N-1:
                Fay=0
                Fby=0 
                
            Fx=(Safe((Faxn**2),(Faxn**2+(miw/mio)*((1-Faxn)**2)))-Safe((Fbxn**2),(Fbxn**2+(miw/mio)*((1-Fbxn)**2)))) + (Safe((Fax**2),(Fax**2+(miw/mio)*((1-Fax)**2)))-Safe((Fbx**2),(Fbx**2+(miw/mio)*((1-Fbx)**2))))
            Fy=(Safe((Fayn**2),(Fayn**2+(miw/mio)*((1-Fayn)**2)))-Safe((Fbyn**2),(Fbyn**2+(miw/mio)*((1-Fbyn)**2)))) + (Safe((Fay**2),(Fay**2+(miw/mio)*((1-Fay)**2)))-Safe((Fby**2),(Fby**2+(miw/mio)*((1-Fby)**2))))
            
            Sw[i,k,j+1]=Sw[i,k,j]+(phix*Fx+phiy*Fy)/porosity[i,k]
    
    return Sw

#-----------------------------------------------------------------------------#
def NewTVD(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    '''Fa means the water is increasing, Fb means the water is going out '''
    for i in range(N):
        for k in range(N):
            
            Fax = 0.
            Fbx = 0.
            Faxn = 0.
            Fbxn = 0.
            Fay = 0.
            Fby = 0.
            Fayn = 0.
            Fbyn = 0.
            
            if i == 0:
                if Sw[i,k,j]>Sw[i+1,k,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fbx=Sw[i,k,j]#0=+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i+1,k,j]:
                    Fax=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
            elif i == 1:
                if Sw[i,k,j]>Sw[i+1,k,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i+1,k,j]:
                    Fax=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i-1,k,j]>Sw[i,k,j]:
                    Faxn=Sw[i-1,k,j]#0=+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-1,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i-1,k,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
            elif i == N-1:
                if Sw[i-1,k,j]>Sw[i,k,j]:
                    Faxn=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i-1,k,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbxn=Sw[i,k,j]#0=+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
            elif i == N-2:
                if Sw[i,k,j]>Sw[i+1,k,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i+1,k,j]:
                    Fax=Sw[i+1,k,j]#0=+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+1,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i-1,k,j]>Sw[i,k,j]:
                    Faxn=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i-1,k,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
            else:
                if Sw[i,k,j]>Sw[i+1,k,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i+1,k,j]:
                    Fax=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i-1,k,j]>Sw[i,k,j]:
                    Faxn=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i-1,k,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbxn=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                    

            if k == 0:
                if Sw[i,k,j]>Sw[i,k+1,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fby=Sw[i,k,j]#0=+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k,j],Sw[i,k+1,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i,k+1,j]:
                    Fay=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
            elif k == 1:
                if Sw[i,k,j]>Sw[i,k+1,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i,k+1,j]:
                    Fay=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i,k-1,j]>Sw[i,k,j]:
                    Fayn=Sw[i,k-1,j]#0=+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-1,j],Sw[i,k,j]-Sw[i,k-1,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k-1,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
            elif k == N-1:
                if Sw[i,k-1,j]>Sw[i,k,j]:
                    Fayn=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k-1,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbyn=Sw[i,k,j]#0=+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k,j],Sw[i,k-1,j]-Sw[i,k,j])))
            elif k == N-2:
                if Sw[i,k,j]>Sw[i,k+1,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i,k+1,j]:
                    Fay=Sw[i,k+1,j]#0=+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+1,j],Sw[i,k,j]-Sw[i,k+1,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i,k-1,j]>Sw[i,k,j]:
                    Fayn=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k-1,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
            else:
                if Sw[i,k,j]>Sw[i,k+1,j]:
                    #Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                    Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))   
                elif Sw[i,k,j]<Sw[i,k+1,j]:
                    Fay=Sw[i,k+1,j]+.5*(Sw[i,k,j]-Sw[i,k+1,j])*max(0,min(1,Safe(Sw[i,k+1,j]-Sw[i,k+2,j],Sw[i,k,j]-Sw[i,k+1,j])))
                    #Fbx=Sw[i,k,j]+.5*(Sw[i-1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i+1,k,j],Sw[i-1,k,j]-Sw[i,k,j])))
                
                if Sw[i,k-1,j]>Sw[i,k,j]:
                    Fayn=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                    #Fbxn=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))   
                elif Sw[i,k-1,j]<Sw[i,k,j]:
                    #Faxn=Sw[i+1,k,j]+.5*(Sw[i,k,j]-Sw[i+1,k,j])*max(0,min(1,Safe(Sw[i+1,k,j]-Sw[i+2,k,j],Sw[i,k,j]-Sw[i+1,k,j])))
                    Fbyn=Sw[i,k,j]+.5*(Sw[i,k-1,j]-Sw[i,k,j])*max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k+1,j],Sw[i,k-1,j]-Sw[i,k,j])))
           
            Fx = Conta(Fax,miw,mio) - Conta(Fbx,miw,mio) + Conta(Faxn,miw,mio) - Conta(Fbxn,miw,mio)
            Fy = Conta(Fay,miw,mio) - Conta(Fby,miw,mio) + Conta(Fayn,miw,mio) - Conta(Fbyn,miw,mio)
            
            Sw[i,k,j+1] =  Sw[i,k,j] + ((Fx/deltx) + (Fy/delty))*deltt/porosity[i,k]
    
    return Sw
#-----------------------------------------------------------------------------#
def Conta(F,miw,mio):
    return Safe((F**2),(F**2 + (miw/mio)*((1-F)**2)))
#-----------------------------------------------------------------------------#
#Safe division, prevines division by zero
def Safe(a,b):
    if b == 0:
        return 0
    else:
        return a/b
#-----------------------------------------------------------------------------#
###############################################################################
###############################################################################
###############################################################################
        
 
    
###################################Interface###################################
def evaluate():
#    global Nx 
#    global Te
#    global Si
#    Nx = int(celulas.get())
#    Te = int(tempo.get())
#    Si = int(opcao.get())
    mainSimulator(int(celulas.get()), int(tempo.get()), int(var.get()))
    
    
root = Tk()
root.title('Simulador Iniciacao Cientifica - version 0.2')
root.iconbitmap(r'simulador.ico')
Label(root,text='SIMULADOR - IC - UENF', font=('bold')).grid(columnspan=5, pady=10, row=0)

Label(root, text="Celulas:").grid(row=1,pady=10,sticky=E)
celulas = Entry(root, width='5',justify='center')
celulas.grid(row=1,column=1,sticky=W)
Label(root,text='', width='2').grid(row=2,column=2)
Label(root, text="Tempo:",width=0).grid(row=1,column=3,sticky=E)

tempo = Entry(root,width='4',justify='center')
tempo.grid(row=1,column=4,sticky=W)

#Label(root, text="1-TVD").grid(row=2,column=0,columnspan=2,sticky=W)
#Label(root, text="2-TVD Completo").grid(row=2,column=3,columnspan=2,sticky=W)
#
#Label(root, text="3-NonStandard").grid(row=3,columnspan=2,sticky=W)
#Label(root, text="4-NonStandard Completo").grid(row=3,column=3,columnspan=2,sticky=W)
#
#Label(root, text="5-NewSimulator").grid(row=4,columnspan=2,sticky=W)
#Label(root, text="6-New TVD").grid(row=4,column=3,columnspan=2,sticky=W)
var = IntVar()
Radiobutton(root, text="TVD", value=1, variable=var).grid(row=2,column=0,columnspan=2,sticky=W)
Radiobutton(root, text="TVD Completo", value=2, variable=var).grid(row=2,column=3,columnspan=2,sticky=W)

Radiobutton(root, text="NonStandard", value=3, variable=var).grid(row=3,columnspan=2,sticky=W)
Radiobutton(root, text="NonStandard Completo", value=4, variable=var).grid(row=3,column=3,columnspan=2,sticky=W)

Radiobutton(root, text="NewSimulator", value=5, variable=var).grid(row=4,columnspan=2,sticky=W)
Radiobutton(root, text="New TVD", value=6, variable=var).grid(row=4,column=3,columnspan=2,sticky=W)

#opcao = Entry(root, width='2',justify='center')
#opcao.grid(row=5,column=1,columnspan=3)

Button(root,text='Simular!!', command = evaluate).grid(row=6, column=1, columnspan=3,pady=10)

root.mainloop()
del E,N,S,W,TkVersion,TclVersion,X,Y,wantobjects
