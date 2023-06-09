import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import solve


num_dt = 1000 # number of timesteps
dt = 0.001 # difference between timesteps
Length = 1 # Length of porous medium
num_dx = 10

def solve_pressure_1D(num_dx, L, num_dt, dt, dx=Length/num_dx):
    pres_line = np.zeros(shape=[num_dx+1, num_dt+1])
    pres_line[0,:] = 1 # left boundary condition
    pres_line[-1,:] = 0 # right boundary condition
    dtdx2 = dt/(dx*dx)
    for time in range(num_dt-1):
        pres_line[1:-1, time+1] = pres_line[1:-1, time] + dtdx2*(pres_line[:-2, time] + pres_line[2:, time] - 2*pres_line[1:-1, time])
    return pres_line

def solve_pressure_1D_implicit(num_dx, Length, num_dt, dt):
    pres_line = np.zeros(shape=[num_dx+1, num_dt])
    pres_line[0,:] = 1 # left boundary condition
    pres_line[-1,:] = 0 # right boundary condition
    dx = Length/num_dx
    dx2dt = np.power(dx, 2)/dt
    for time in np.arange(num_dt-1):
        A = np.zeros([num_dx-1, num_dx-1])
        mid = -(2+np.power(dx, 2)/dt)
        for n in range(num_dx-2):
            A[n, n] = mid
            
            A[n+1, n] = 1
            
            A[n, n+1] = 1
            
            if n == num_dx-3:
                A[n+1,n+1] = mid
        b = -(dx2dt)*pres_line[1:-1, time]
        # first and last elements
        b[0] = -(dx2dt)*pres_line[1, time] - pres_line[0, time]
        b[-1] = -(dx2dt)*pres_line[-2, time] - pres_line[-1, time]
        
        pres_line[1:-1, time+1] = solve(A, b)
    return pres_line

pressure_explicit = solve_pressure_1D(10, Length, num_dt, dt) # explicit
pressure_implicit = solve_pressure_1D_implicit(10, Length, num_dt, dt)



plt.figure()
plt.xlabel('x')
plt.ylabel('pressure')
time = 200

x = np.linspace(0, Length, num_dx+1)
plt.plot(x, pressure_explicit[:, time])
plt.plot(x, pressure_implicit[:, time])
