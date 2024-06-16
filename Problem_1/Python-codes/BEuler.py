#!/usr/bin/env python
# coding: utf-8



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import Axes3D






# This function implements the right-hand side of the q-equation
def nPointSystem(m, q, G, n):
    f = np.zeros((n, 3))

    # Loop over all bodies
    for ii in range(0, n):

        # Loop over the three coordinates
        for jj in range(0, 3):
            temp = -m * (q[ii, jj] - q[:, jj]) / np.sqrt((q[ii, 0] - q[:, 0]) ** 2 + (q[ii, 1] - q[:, 1]) ** 2 + (q[ii, 2] - q[:, 2]) ** 2) ** 3
            not_nan_values = [not np.isnan(k) for k in temp]
            f[ii, jj] = G * m[ii] * np.sum(temp[not_nan_values])

    return f





def forBackwardEuler(zn, data):
    # Unpack the arguments
    z = data[0]
    m = data[1]
    m_repeated = data[2]
    G = data[3]
    dt = data[4]

    F = np.zeros(6 * n)

    # Unpack the candidate new values
    qnx = zn[0:n];
    qny = zn[n:2 * n];
    qnx = zn[2 * n:3 * n]
    pnx = zn[3 * n:4 * n];
    pny = zn[4 * n:5 * n];
    pnz = zn[5 * n:6 * n]

    # Unpack the old values
    qx = z[0:n];
    qy = z[n:2 * n];
    qx = z[2 * n:3 * n]
    px = z[3 * n:4 * n];
    py = z[4 * n:5 * n];
    pz = z[5 * n:6 * n]

    # Implicit form of the Hamiltonian equations for backward Euler
    Tq = (qn - q) / dt - pn / m_repeated
    Tp = (pn - p) / dt - nPointSystem(m, qn, G, n)

    # Stack the variables in a 6n dimensional vector
    F[0:n] = Tq[:, 0]
    F[n:2 * n] = Tq[:, 1]
    F[2 * n:3 * n] = Tq[:, 2]
    F[3 * n:4 * n] = Tp[:, 0]
    F[4 * n:5 * n] = Tp[:, 1]
    F[5 * n:6 * n] = Tp[:, 2]

    return F





data = np.loadtxt("data.txt")





m = data[:,0]
q = data[:,1:4]
qdot = data[:,4:]
G = 2.95912208286e-4





# Integration time in years
tyear = 1000

# Integration time in days and time step in days
t = tyear * 365
dt = 1
steps = np.int(t / dt)

# Number of planets
n = 10

# Initialize the solution vectors
solx = np.zeros((n, steps))
soly = np.zeros((n, steps))
solz = np.zeros((n, steps))

# Initial momentum
m_repeated = np.repeat(np.expand_dims(m, 1), 3, axis=1)  # This just creates the numpy array [m,m,m]
p = qdot * m_repeated


# Save initial solution vector
solx[:, 0] = q[:, 0]
soly[:, 0] = q[:, 1]
solz[:, 0] = q[:, 2]


# Energies
UE = np.zeros((steps),float)
KE = np.zeros((steps),float)
H = np.zeros((steps),float)

for i in range(10):
  KE[0] = KE[0] + 0.5 * (p[i,0]*p[i,0] + p[i,1]*p[i,1] + p[i,2]*p[i,2] )/m[i]    # Kinetic Energy
  for j in range(10):
    if(i != j):      
      dx = q[i,0] - q[j,0] 
      dy = q[i,1] - q[j,1] 
      dz = q[i,2] - q[j,2] 
      r = np.sqrt(dx**2 + dy**2 + dz**2)
      UE[0] = UE[0] - G* m[i]*m[j] / (r)                                      # Potential Energy





q.shape                                                                       # In order to be sure from the size of solutions







# Main time iiping loop
for ii in range(1, steps):

  if np.mod(ii,10)==0:
    print('Computing step {} out of {} steps.'.format(ii, steps))



  # Euler forward method (starting step for implicit solver)
  qn = q + dt*p/m_repeated
  pn = p + dt*nPointSystem(m, q, G, n)

  # Convert to state vector
  zn = np.concatenate((qn[:,0],qn[:,1],qn[:,2],pn[:,0],pn[:,1],pn[:,2]))
  z = np.concatenate((q[:,0],q[:,1],q[:,2],p[:,0],p[:,1],p[:,2]))

  # Auxiliary data needed to solve implicit equations
  data = np.array([z,m,m_repeated,G,dt])

  # Fsolve
  zn = fsolve(forBackwardEuler,zn,args=data,xtol=1e-9)

  # Unpack the values and store in q and p
  qn[:,0] = zn[0:n];  qn[:,1] = zn[n:2 * n]; qn[:,2] = zn[2 * n:3 * n]
  pn[:,0] = zn[3 * n:4 * n];  pn[:,1] = zn[4 * n:5 * n];  pn[:,2] = zn[5 * n:6 * n]

  # Move values
  p = pn
  q = qn



  # Save solution vector
  solx[:, ii] = q[:, 0]; soly[:, ii] = q[:, 1]; solz[:, ii] = q[:, 2]


# energy
  for i in range(10):
    KE[ii] = KE[ii] + 0.5 * (p[i,0]*p[i,0] + p[i,1]*p[i,1] + p[i,2]*p[i,2]  )/m[i]
    for j in range(10):
      if(i != j):
              
        dx = q[i,0] - q[j,0] 
        dy = q[i,1] - q[j,1] 
        dz = q[i,2] - q[j,2] 
        r = np.sqrt(dx**2 + dy**2 + dz**2)
        UE[ii]=UE[ii] - G* m[i]*m[j] / (r)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for ii in range(0, n):
    ax.plot(solx[ii, :], soly[ii, :], solz[ii, :])


## Write the Hemilotonian


H = UE/2 + KE





time = np.arange(t)


########## Plot the Hamiltonian ##############


plt.figure(figsize = (20,10))
plt.plot(time,H)
plt.xlim(0,365000)
plt.xlabel("Time",fontsize = 18)
plt.ylabel("Energy", fontsize = 18)
plt.title("Energy for Backward Euler method",fontsize = 18)
plt.savefig("Energy-Backward.png")



######## Plot the planets ##############

fig = plt.figure(figsize=(30,15))
ax = plt.axes(projection='3d')
for ii in range(0, n):
    ax.plot(solx[ii, :], soly[ii, :], solz[ii, :])
ax.set_title("Solar system by Backward Euler",fontsize = 25)
ax.set_xlabel('x',fontsize = 25)
ax.set_ylabel('y',fontsize = 25)
ax.set_zlabel('z',fontsize = 25)
fig.savefig("Backward_1.png")



###### Write Hamiltonian into a file ##########

datafile_path = "./H-BEuler.txt"
with open(datafile_path, 'w+') as datafile_id:
    np.savetxt(datafile_id, H)







