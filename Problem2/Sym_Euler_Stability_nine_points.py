#!/usr/bin/env python
# coding: utf-8

# # Symplectic Euler - Polar Vortex - Check the stability

# ### Import libraries




import numpy as np
import matplotlib.pyplot as plt
import math as mt
n = 9





N = 180*24*60*60
delta = 1

####################### initializtion ########################
x = np.zeros((N,n),float)
y = np.zeros((N,n),float)
dx =np.zeros((n,n),float)
dy = np.zeros((n,n),float)

###################### Store answers ########################
x1 = np.zeros(N,float)
y1 = np.zeros(N,float)
x2 = np.zeros(N,float)
y2 = np.zeros(N,float)
x3 = np.zeros(N,float)
y3 = np.zeros(N,float)
x4 = np.zeros(N,float)
y4 = np.zeros(N,float)
x5 = np.zeros(N,float)
y5 = np.zeros(N,float)
x6 = np.zeros(N,float)
y6 = np.zeros(N,float)
x7 = np.zeros(N,float)
y7 = np.zeros(N,float)
x8 = np.zeros(N,float)
y8 = np.zeros(N,float)
x9 = np.zeros(N,float)
y9 = np.zeros(N,float)


r = np.zeros((n,n),float)
ux = np.zeros((n),float)
uy = np.zeros((n),float)
######################## Constants #########################
H = np.zeros((N-1),float)
px = np.zeros((N-1),float)
py = np.zeros((N-1),float)
lz = np.zeros((N-1),float)


gamma = np.array([6*10**8,0.7*10**8,0.7*10**8,0.7*10**8,0.7*10**8,0.7*10**8,0.7*10**8,0.7*10**8,0.7*10**8])





gamma[0]


# ### Initial values 




x[0,0] = 0
y[0,0] = 250000

x[0,1] = 10000000
y[0,1] = 0

x[0,2] = (np.sqrt(2)/2)*10000000
y[0,2] = (np.sqrt(2)/2)*10000000

x[0,3] = 0
y[0,3] = 10000000

x[0,4] = -(np.sqrt(2)/2)*10000000
y[0,4] = (np.sqrt(2)/2)*10000000

x[0,5] = -10000000
y[0,5] = 0

x[0,6] = -(np.sqrt(2)/2)*10000000
y[0,6] = -(np.sqrt(2)/2)*10000000

x[0,7] = 0
y[0,7] = -10000000

x[0,8] = (np.sqrt(2)/2)*10000000
y[0,8] = -(np.sqrt(2)/2)*10000000


# ### Main algorithms




for t in range(1,N):
    for k in range(n):
         ux[k] = uy[k] = 0
    for i in range(n):
        lz[t-1] = lz[t-1] - 0.5 * (gamma[i]) * (x[t-1,i] * x[t-1,i] + y[t-1,i] * y[t-1,i])
        px[t-1] = px[t-1] + (gamma[i]) * y[t-1,i]
        py[t-1] = py[t-1] + (gamma[i]) * x[t-1,i]
        for j in range(n):
            if(i != j):
                dx[i,j] = x[t-1,i]-x[t-1,j]
                dy[i,j] = y[t-1,i]-y[t-1,j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2)
                H[t-1] = H[t-1] - (gamma[i]*gamma[j])*(1/mt.pi)*np.log(r[i,j])               
    for i in range(n):
        for j in range(n):
            if(i != j):
                dx[i,j] = x[t-1,i]-x[t-1,j]
                dy[i,j] = y[t-1,i]-y[t-1,j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2)
                ux[i] = ux[i] - (gamma[j]) * (1/4*mt.pi)*dy[i,j]/(r[i,j]**2)
    for l in range(n):
            x[t,l] = x[t-1,l] + delta*ux[l]           
    for i in range(n):
        for j in range(n):
            if(i != j):
                dx[i,j] = x[t,i]-x[t,j]
                dy[i,j] = y[t-1,i]-y[t-1,j]
                r[i,j] = np.sqrt(dx[i,j]**2 + dy[i,j]**2)
                uy[i] = uy[i] + (gamma[j]) * (1/4*mt.pi)*dx[i,j]/(r[i,j]**2)
    for l in range(n):
            y[t,l] = y[t-1,l] + delta*uy[l]


# ### Answers to arrays




x1 = x[:,0]
x2 = x[:,1]
x3 = x[:,2]
x4 = x[:,3]
x5 = x[:,4]
x6 = x[:,5]
x7 = x[:,6]
x8 = x[:,7]
x9 = x[:,8]

y1 = y[:,0]
y2 = y[:,1]
y3 = y[:,2]
y4 = y[:,3]
y5 = y[:,4]
y6 = y[:,5]
y7 = y[:,6]
y8 = y[:,7]
y9 = y[:,8]


# ### plot answers




plt.figure(figsize = (20,20))
plt.plot(x1,y1,color = 'green')
plt.plot(x2,y2,color = 'red')
plt.plot(x3,y3,color = 'blue')
plt.plot(x4,y4,color = 'pink')
plt.plot(x5,y5,color = 'brown')
plt.plot(x6,y6,color = 'violet')
plt.plot(x7,y7,color = 'black')
plt.plot(x8,y8,color = 'orange')
plt.plot(x9,y9,color = 'gray')
#plt.title("Polar Vortex (N = 3) with symplectic Euler Method", fontsize = 18)
plt.xlabel('x',fontsize = 18)
plt.ylabel('y', fontsize = 18)
plt.savefig('polar votex-sym_Euler.png')


# ### Define time and constants




time = np.arange(N-1)





datafile_path = "./H-Sym_Euler.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, H)
    
datafile_path = "./L-Sym_Euler.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, lz)
    
datafile_path = "./Px-Sym_Euler.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, px)
    
datafile_path = "./Py-Sym_Euler.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, py)


# ### Plots




plt.figure(figsize = (20,10))
plt.plot(time,H)
plt.xlabel("Time", fontsize = 18)
plt.ylabel("Energy", fontsize = 18)
plt.title("Energy for polar vortex - Symplectic Euler", fontsize = 18)
plt.xlim(0,N-1)
plt.savefig("H-Sym-Euler.png")





plt.figure(figsize = (20,10))
plt.plot(time,lz)
plt.xlabel("Time", fontsize = 18)
plt.ylabel("Angular Momentum", fontsize = 18)
plt.title("Angular Momentum for polar vortex - Symplectic Euler", fontsize = 18)
plt.xlim(0,N-1)
plt.savefig("L-Sym-Euler.png")





plt.figure(figsize = (20,10))
plt.plot(time,px)
plt.xlabel("Time", fontsize = 18)
plt.ylabel("$P_{x}$", fontsize = 18)
plt.title("Momentum for polar vortex (x-direction) - Symplectic Euler", fontsize = 18)
plt.xlim(0,N-1)
plt.savefig("px-Sym-Euler.png")





plt.figure(figsize = (20,10))
plt.plot(time,py)
plt.xlabel("Time", fontsize = 18)
plt.ylabel("$P_{y}$", fontsize = 18)
plt.title("Momentum for polar vortex (y-direction) - Symplectic Euler", fontsize = 18)
plt.xlim(0,N-1)
plt.savefig("py-Sym-Euler.png")







