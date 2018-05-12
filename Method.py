# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 09:03:01 2017

@author: Nicholas de Almeida Pinto
"""

#import numpy as np



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
def TVD_NonStandard_Not_Both_Directions(Sw,N, mio, miw, j, porosity, deltt, deltx, delty):
    for i in range(N):
        for k in range(N):
            Nax,Nbx,Nay,Nby = 0,0,0,0;
            if i == 0:
                Fax=Sw[i,k,j]
                Fbx=Sw[i,k,j]
            elif i == N-1:
                Nxa=(Sw[i-1,k,j]**2)/((Sw[i-1,k,j]**2)+(miw/mio)*((1-Sw[i-1,k,j])**2))
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*Nax#max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Fbx=Sw[i,k,j]
            else:
                Nax=(Sw[i-1,k,j]**2)/((Sw[i-1,k,j]**2)+(miw/mio)*((1-Sw[i-1,k,j])**2))
                Fax=Sw[i-1,k,j]+.5*(Sw[i,k,j]-Sw[i-1,k,j])*Nax#max(0,min(1,Safe(Sw[i-1,k,j]-Sw[i-2,k,j],Sw[i,k,j]-Sw[i-1,k,j])))
                Nbx=(Sw[i+1,k,j]**2)/((Sw[i+1,k,j]**2)+(miw/mio)*((1-Sw[i+1,k,j])**2))
                Fbx=Sw[i,k,j]+.5*(Sw[i+1,k,j]-Sw[i,k,j])*Nbx#max(0,min(1,Safe(Sw[i,k,j]-Sw[i-1,k,j],Sw[i+1,k,j]-Sw[i,k,j])))

            if k == 0:
                Fay=Sw[i,k,j]
                Fby=Sw[i,k,j]
            elif k == N-1:
                Nay=(Sw[i,k-1,j]**2)/((Sw[i,k-1,j]**2)+(miw/mio)*((1-Sw[i,k-1,j])**2))
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*Nay#max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Fby=Sw[i,k,j]
            else:
                Nay=(Sw[i,k-1,j]**2)/((Sw[i,k-1,j]**2)+(miw/mio)*((1-Sw[i,k-1,j])**2))
                Fay=Sw[i,k-1,j]+.5*(Sw[i,k,j]-Sw[i,k-1,j])*Nay#max(0,min(1,Safe(Sw[i,k-1,j]-Sw[i,k-2,j],Sw[i,k,j]-Sw[i,k-1,j])))
                Nby=(Sw[i,k+1,j]**2)/((Sw[i,k+1,j]**2)+(miw/mio)*((1-Sw[i,k+1,j])**2))
                Fby=Sw[i,k,j]+.5*(Sw[i,k+1,j]-Sw[i,k,j])*Nby#max(0,min(1,Safe(Sw[i,k,j]-Sw[i,k-1,j],Sw[i,k+1,j]-Sw[i,k,j])))
           
            Fx = Safe((Fax**2),(Fax**2 + (miw/mio)*((1-Fax)**2))) - Safe((Fbx**2),(Fbx**2 + (miw/mio)*((1-Fbx)**2)))
            Fy = Safe((Fay**2),(Fay**2 + (miw/mio)*((1-Fay)**2))) - Safe((Fby**2),(Fby**2 + (miw/mio)*((1-Fby)**2)))
            
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