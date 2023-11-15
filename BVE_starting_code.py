#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:50:18 2023

@author: Paolo
"""
###################################################
##############PARAMETERS###########################
################################################### 
DAYLEN  = 1   #  Forecast length in days.
DtHours = .5   #  Timestep in hours.
###################################################
##############IMPORT LIBRARIES#####################
################################################### 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

###################################################
##############DEFINE FUNCTIONS#####################
###################################################     


def make_Laplacian(Z):
    #  Compute the Laplacian 
    #of the geopotential height
    # within the boundary
    #(boundary excluded)
    M       = Z.shape[0]
    N       = Z.shape[1]
    Zxx  = np.zeros([M,N])     #  second x derivative of Z
    Zyy  = np.zeros([M,N])     #  second y derivative of Z
    L0in     = np.zeros([M-1,N-1])     #  Laplacian of Z
    # Compute within the domain (no boundaries) 
    # Second x-derivative of Z
    Zxx[1:M-1,:] = (Z[2:M,:]+Z[0:M-2,:]-2*Z[1:M-1,:])/(736e+3**2)
    # Second y-derivative of Z
    Zyy[:,1:N-1] = (Z[:,2:N]+Z[:,0:N-2]-2*Z[:,1:N-1])/(736e+3**2)    
    ##  Laplacian of height (or vorticity)
    L0in = Zxx[1:M-1,1:N-1]+Zyy[1:M-1,1:N-1]
    return L0in

def make_Jacobian(Z,ABS_VOR):
    M       = Z.shape[0]
    N       = Z.shape[1]
    Zx    = np.zeros([M,N])     #  x derivative of Z
    Zy    = np.zeros([M,N])     #  y derivative of Z
    ABS_VORx  = np.zeros([M,N])     #  x derivative of ABS_VOR
    ABS_VORy  = np.zeros([M,N])     #  y derivative of ABS_VOR
    # Compute within 
    #the domain (boundary excluded) 
    # x-derivative of Z
    Zx[1:M-1,:] = (Z[2:M,:]-Z[0:M-2,:])/(2*736e+3)
    # y-derivative of Z
    Zy[:,1:N-1] = (Z[:,2:N]-Z[:,0:N-2])/(2*736e+3)
    # x-derivative of the absolute vorticity 
    ABS_VORx[1:M-1,:] = (ABS_VOR[2:M,:]-ABS_VOR[0:M-2,:])/(2*736e+3)
    # y-derivative of the absolute vorticity 
    ABS_VORy[:,1:N-1] = (ABS_VOR[:,2:N]-ABS_VOR[:,0:N-2])/(2*736e+3)
    ##  Compute the Jacobian J(ABS_VOR,Z)
    Jacobi = ABS_VORx * Zy - ABS_VORy * Zx
    return Jacobi

def Poisson_solver(Jacobi):
    M       = Jacobi.shape[0]
    N       = Jacobi.shape[1]
    SM=np.zeros([M-2,M-2])
    SN=np.zeros([N-2,N-2])
    EIGEN=np.zeros([M-2,N-2])    
    ##  Coefficients for x-transformation
    for m1 in range(0,M-2):
     for m2 in range(0,M-2):
      SM[m1,m2] = np.sin(np.pi*(m1+1)*(m2+1)/(M-1))       
    ##  Coefficients for y-transformation
    for n1 in range(0,N-2):
     for n2 in range(0,N-2):
      SN[n1,n2] = np.sin(np.pi*(n1+1)*(n2+1)/(N-1))        
    ##  Eigenvalues of Laplacian operator
    for mm in range(0,M-2):
     for nn in range(0,N-2):
      eigen = (np.sin(np.pi*(mm+1)/(2*(M-1))))**2 +(np.sin(np.pi*(nn+1)/(2*(N-1))))**2
      EIGEN[mm,nn] = (-4/736e+3**2) * eigen
    #  Tendency values in interior.
    Ldot = Jacobi[1:M-1,1:N-1] 
    #  Compute the transform of the solution
    LDOT = np.dot(SM,np.dot(Ldot,SN))
    #  Convert transform of d(xi)/dt to transform of d(Z)/dt
    ZDOT = LDOT / EIGEN 
    #  Compute inverse transform to get the height tendency.
    Zdot = (4/((M-1)*(N-1))) *np.dot(SM,np.dot(ZDOT,SN))
    return Zdot
      
 
def make_f_and_h(N,M,Xp,Yp):
    FCOR=np.zeros([M,N])
    h=np.zeros([M,N])    
    a = (4*10**7)/(2*np.pi)      #  Radius of the Earth
    grav = 9.80665           #  Gravitational acceleration
    Omega = 2*np.pi/(24*60*60)  #  Angular velocity of Earth.
    ##  Compute Coriolis Parameter and Map Factor
    ##  and parameter h = g*m**2/f used in the BVE
    for ny in range(0,N):
     for nx in range(0,M):
      xx = (nx-Xp)*736e+3
      yy = (ny-Yp)*736e+3
      rr = np.sqrt(xx**2+yy**2)
      phi = 2*((np.pi/4)-np.arctan(rr/(2*a)))
      mapPS = 2 / (1+np.sin(phi))
      f = 2*Omega*np.sin(phi)
      FCOR[nx,ny] = f
      h[nx,ny] = grav * mapPS**2 / f
    return FCOR,h

###################################################
##############FIXED PARAMETERS###########################
###################################################   
M  = 19 # Points in x direction
N  = 16 # Points in y direction
Xp = 8 # Coord. of North Pole
Yp = 12 # Coord. of North Pole

###################################################
###############COORDINATES AND TIME################
###################################################
daylen = DAYLEN                   #  Integration time (in days)
seclen = int(daylen*24*60*60)     #  Integration time (in secon736e+3)
Dt = DtHours*60*60                #  Timestep in secon736e+3
nt = int(seclen//Dt )             #  Total number of time-steps.
# Define the (X,Y) grid (for plotting)
X, Y  = np.meshgrid(np.linspace(1,M,M),np.linspace(1,N,N))
X = np.transpose(X)
Y = np.transpose(Y)
#Coriolis and map factor
FCOR,h=make_f_and_h(N,M,Xp,Yp)

###################################################
###############DEFINE WORKING ARRAYS###############
###################################################
Zout=np.zeros([nt+1,M,N])   
L0=np.zeros([nt+1,M,N])
###################################################
###############READ INPUT DATA##########################
###################################################
###   Read and plot the initial and verification height data
## Input files
File1 = './Case1-1949010503.z00' #The initial value
File2 = './Case1-1949010603.z00' #The final value
Z0  = np.genfromtxt(File1)
Z24 = np.genfromtxt(File2)

Zout[0,:,:]  = Z0      #  Copy initial height field
##  Compute the Laplacian of the height/psi
#################################################
#################################################
############      MAIN LOOP     #################
#################################################
#################################################
###                                           ###
###      Integrate the BVE in time            ###
###                                           ###
#################################################
#################################################
######## Start of the time-stepping loop  #######

# step 0

# linear extrapolation of boundary
def extrapolate(A_core):
  m,n = np.shape(A_core)
  A = np.zeros((m+2,n+2))
  A[1:-1,1:-1] = A_core

  # rows
  A[1:-1,0] = 2*A_core[:,0] - A_core[:,1]
  A[1:-1,-1] = 2*A_core[:,-1] - A_core[:,-2]

  # cols
  A[0,1:-1] = 2*A_core[0,:] - A_core[1,:]
  A[-1,1:-1] = 2*A_core[-1,:] - A_core[-2,:]

  # corners
  A[0,0] = 0.5*(A[1,0]+A[0,1])
  A[0,-1] = 0.5*(A[1,-1]+A[0,-2])
  A[-1,0] = 0.5*(A[-1,1]+A[-2,0])
  A[-1,-1] = 0.5*(A[-1,-2]+A[-2,-1])

  return A

L0[0,:,:] = extrapolate(make_Laplacian(Zout[0,:,:]))
print(L0[0,0,:])
J = make_Jacobian(Zout[0,:,:],h*L0[0,:,:]+FCOR)

# Euler
Zout[1,:,:] = Zout[0,:,:]
Zout[1,1:-1,1:-1] = Zout[0,1:-1,1:-1] + Dt*Poisson_solver(J)

Zout_unfiltered = Zout

# step n > 0

for i in range(1,nt):
  L0[i,:,:] = L0[i-1,:,:] # save boundary
  #print(L0[i,0,:])
  L0[i,1:-1,1:-1] = make_Laplacian(Zout[i,:,:])
  J = make_Jacobian(Zout[i,:,:],h*L0[i,:,:]+FCOR)
  #print(Poisson_solver(J)[5,5])
  # Leapfrog
  Zout[i+1,:,:] = Zout[i,:,:] # save boundary
  Zout_unfiltered[i+1,1:-1,1:-1] = Zout[i-1,1:-1,1:-1] + (Poisson_solver(J) * Dt*2)
  # Robert-Asselin
  Zout[i+1,1:-1,1:-1] = 0.1*Zout[i+1,1:-1,1:-1] + 0.9*Zout_unfiltered[i+1,1:-1,1:-1]

fig, ax = plt.subplots(nrows=1)

cntr0 = ax.contour(X, Y, Zout[0,:,:], levels=8, linewidths=0.5, colors='k')
#cntr1 = ax.contourf(X, Y, Zout[-1,:,:], levels=14, cmap="RdBu_r")
cntr1 = ax.contour(X, Y, Zout[-1,:,:], levels=8, linewidths=0.5, colors='r')
ax.clabel(cntr1, inline=True, fontsize=10)
ax.clabel(cntr0, inline=True, fontsize=10)
labels = ['IC', 'forecast']

artist = cntr0.legend_elements()[0]
artist.append(cntr1.legend_elements()[0])
pers_label = mpatches.Patch(color='k', label='persistence', fill=False, linewidth=0.5)
forecast_label = mpatches.Patch(color='r', label='forecast', fill=False, linewidth=0.5)
plt.legend(handles=[pers_label, forecast_label], loc='lower left')
plt.tight_layout()
fig.suptitle("Figure 1: persistence and forecast", fontsize=8)

plt.savefig("forecast.png")

fig, ax = plt.subplots(nrows=1)

cntr0 = ax.contour(X, Y, Z24[:,:], levels=8, linewidths=0.5, colors='k')
#cntr1 = ax.contourf(X, Y, Zout[-1,:,:], levels=14, cmap="RdBu_r")
cntr1 = ax.contour(X, Y, Zout[-1,:,:], levels=8, linewidths=0.5, colors='r')
ax.clabel(cntr1, inline=True, fontsize=10)
ax.clabel(cntr0, inline=True, fontsize=10)
obs_label = mpatches.Patch(color='k', label='observation', fill=False, linewidth=0.5)
plt.legend(handles=[obs_label, forecast_label], loc='lower left')
plt.tight_layout()
fig.suptitle("Figure 2: observation and forecast", fontsize=8)

plt.savefig("analysis.png")

fig, ax = plt.subplots(nrows=1)

cntr0 = ax.contour(X, Y, Zout[-1,:,:]-Zout[0,:,:], levels=6, linewidths=0.5, colors='k')
#cntr1 = ax.contourf(X, Y, Zout[-1,:,:], levels=14, cmap="RdBu_r")
cntr1 = ax.contour(X, Y, Z24[:,:] - Zout[0,:,:], levels=6, linewidths=0.5, colors='r')
ax.clabel(cntr1, inline=True, fontsize=10)
ax.clabel(cntr0, inline=True, fontsize=10)
obs_tend_label = mpatches.Patch(color='k', label='observed tendency', fill=False, linewidth=0.5)
forecast_tend_label = mpatches.Patch(color='r', label='forecast tendency', fill=False, linewidth=0.5)
plt.legend(handles=[obs_tend_label, forecast_tend_label], loc='lower right')
plt.tight_layout()
fig.suptitle("Figure 3: observation and forecast tendencies", fontsize=8)

plt.savefig("tendency.png")

def RMSE(x_true, x):
    return np.sqrt((x_true - x)**2)

rmse_comp = RMSE(Zout[-1,:,:], Z24[:,:])
rmse_pers = RMSE(Zout[0,:,:], Z24[:,:])

print(np.mean(rmse_comp))
print(np.mean(rmse_pers))

fig, ax = plt.subplots(nrows=1)
cntr1 = ax.contourf(X, Y, rmse_comp, levels=14, cmap="RdBu_r")
plt.savefig("rmse.png")

