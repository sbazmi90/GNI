#!/usr/bin/env python
# coding: utf-8

# # Heun Method for solar system

# ### import libraries




import numpy as np
import matplotlib.pyplot as plt


# ### Reading data file




data = np.loadtxt('data.txt')


# ### Importing data to variables




m = data[:,0]                             # Mass
xi = data[:,1]                            # X position
yi = data[:,2]                            # Y Position
zi = data[:,3]                            # Z Position
vxi = data[:,4]                           # Velocity in X direction
vyi = data[:,5]                           # Velocity in Y direction
vzi = data[:,6]                           # Velocity in Z direction
t_max = 365000                            # Total time of integration
G = 0.000295912208286                     # Gravitational constants
delta =1                                  # Time steps that is 1 day
N = 365000                                # Define the number of elements in our matrices





#====================== Define arrays for storing answers ========================


#====================== Positions ======================
x = np.zeros((N,10),float)
y = np.zeros((N,10),float)
z = np.zeros((N,10),float)

#==================== Momentums and energies ========================
px = np.zeros((N,10),float)
py = np.zeros((N,10),float)
pz = np.zeros((N,10),float)
UE = np.zeros((N-1),float) 
KE = np.zeros((N-1),float)

#=================== Temporary arrays for storing measurement in different time steps ===================
dx =np.zeros((10,10),float)
dy = np.zeros((10,10),float)
dz = np.zeros((10,10),float)
ux = np.zeros((10),float)
uy = np.zeros((10),float)
uz = np.zeros((10),float)
Ux = np.zeros((10),float)
Uy = np.zeros((10),float)
Uz = np.zeros((10),float)
r = np.zeros((10,10),float)

x1 = np.zeros((10),float)
x2 = np.zeros((10),float)
y1 = np.zeros((10),float)
y2 = np.zeros((10),float)
z1 = np.zeros((10),float)
z2 = np.zeros((10),float)

px1 = np.zeros((10),float)
px2 = np.zeros((10),float)
py1 = np.zeros((10),float)
py2 = np.zeros((10),float)
pz1 = np.zeros((10),float)
pz2 = np.zeros((10),float)

ep = 0                        # In order to prevent any singularities in denominator


# ### Define initial conditions



for i in range(10):
    x[0,i] = xi[i]
    y[0,i] = yi[i]
    z[0,i] = zi[i]
    px[0,i] = m[i]*vxi[i]
    py[0,i] = m[i]*vyi[i]
    pz[0,i] = m[i]*vzi[i]


# ### Main algorithm



for t in range(1,N):
    for k in range(10):
        ux[k] = uy[k] = uz[k] = UE[k] = KE[k] = 0
    for i in range(10):
        KE[t-1] = KE[t-1] + 0.5 * (px[t-1,i] * px[t-1,i] + py[t-1,i] * py[t-1,i] + pz[t-1,i] * pz[t-1,i])/m[i]
        for j in range(10):
            if(i != j):
                dx[i,j] = x[t-1,i]-x[t-1,j]
                dy[i,j] = y[t-1,i]-y[t-1,j]
                dz[i,j] = z[t-1,i]-z[t-1,j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2 + dz[i,j]**2)
                UE[t-1]=UE[t-1] - G* m[i]*m[j] / (r[i,j]+ep)
    for i in range(10):
        x1[i] = x[t-1,i]
        y1[i] = y[t-1,i]
        z1[i] = z[t-1,i]
        px1[i] = px[t-1,i]
        py1[i] = py[t-1,i]
        pz1[i] = pz[t-1,i]
    for i in range(10):
        for j in range(10):
            if(i != j):
                dx[i,j] = x1[i]-x1[j]
                dy[i,j] = y1[i]-y1[j]
                dz[i,j] = z1[i]-z1[j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2 + dz[i,j]**2)
                ux[i] = ux[i]-G*m[i]*m[j]*dx[i,j]/((r[i,j]+ep)**3)
                uy[i] = uy[i]-G*m[i]*m[j]*dy[i,j]/((r[i,j]+ep)**3)
                uz[i] = uz[i]-G*m[i]*m[j]*dz[i,j]/((r[i,j]+ep)**3)
    for i in range(10):
        x2[i] = x1[i] + delta * px1[i]/m[i]
        y2[i] = y1[i] + delta * py1[i]/m[i]
        z2[i] = z1[i] + delta * pz1[i]/m[i]
        px2[i] = px1[i] + delta*ux[i]
        py2[i] = py1[i] + delta*uy[i]
        pz2[i] = pz1[i] + delta*uz[i]
    for k in range(10):
        Ux[k] = Uy[k] = Uz[k] = 0
    for i in range(10):
        for j in range(10):
            if(i != j):
                dx[i,j] = x2[i]-x2[j]
                dy[i,j] = y2[i]-y2[j]
                dz[i,j] = z2[i]-z2[j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2 + dz[i,j]**2)
                Ux[i] = Ux[i]-G*m[i]*m[j]*dx[i,j]/((r[i,j]+ep)**3)
                Uy[i] = Uy[i]-G*m[i]*m[j]*dy[i,j]/((r[i,j]+ep)**3)
                Uz[i] = Uz[i]-G*m[i]*m[j]*dz[i,j]/((r[i,j]+ep)**3)   
    for l in range(10):
            x[t,l] = x[t-1,l] + 0.5*delta*(px1[l]+px2[l])/m[l]
            y[t,l] = y[t-1,l] + 0.5*delta*(py1[l]+py2[l])/m[l]
            z[t,l] = z[t-1,l] + 0.5*delta*(pz1[l]+pz2[l])/m[l]
            px[t,l] = px[t-1,l] + 0.5*delta*(ux[l]+Ux[l])
            py[t,l] = py[t-1,l] + 0.5*delta*(uy[l]+Uy[l])
            pz[t,l] = pz[t-1,l] + 0.5*delta*(uz[l]+Uz[l])


# ### Writing answers to new arrays

# In[7]:


xsun = x[:,0];ysun = y[:,0];zsun = z[:,0];
xmer = x[:,1];ymer = y[:,1];zmer = z[:,1];
xven = x[:,2];yven = y[:,2];zven = z[:,2];
xear = x[:,3];year = y[:,3];zear = z[:,3];
xmar = x[:,4];ymar = y[:,4];zmar = z[:,4];
xjup = x[:,5];yjup = y[:,5];zjup = z[:,5];
xsat = x[:,6];ysat = y[:,6];zsat = z[:,6];
xura = x[:,7];yura = y[:,7];zura = z[:,7];  
xnep = x[:,8];ynep = y[:,8];znep = z[:,8];                                     
xplo = x[:,9];yplo = y[:,9];zplo = z[:,9];


# ### Define time and total energy of the system



time = np.arange(N-1)
H_Heun = np.zeros((N-1),float)




for i in range(0,N-1):
    H_Heun[i] = UE[i]/2 + KE[i]





datafile_path = "./H-Heun.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, H_Heun)
















