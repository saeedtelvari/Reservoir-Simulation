import matplotlib.pyplot as plt
import numpy as np


num_dt = 1000 # number of time steps
dt = 0.001 # difference between time steps
Length = 1 # Length of porous medium
num_dx = 10 # number of grids

def solve_pressure_1D(num_dx, L, num_dt, dt, dx=Length/num_dx):
    pres_line = np.zeros(shape=[num_dx+1, num_dt+1])
    pres_line[0,:] = 1 # left boundary condition
    pres_line[-1,:] = 0 # right boundary condition
    dtdx2 = dt/(dx*dx)
    for time in range(num_dt-1):
        pres_line[1:-1, time+1] = pres_line[1:-1, time] + dtdx2*(pres_line[:-2, time] + pres_line[2:, time] - 2*pres_line[1:-1, time])
    return pres_line


pressure_ndx_3 = solve_pressure_1D(3, Length, num_dt, dt)
pressure_ndx_6 = solve_pressure_1D(6, Length, num_dt, dt)
pressure_ndx_12 = solve_pressure_1D(12, Length, num_dt, dt)

plt.figure()
plt.title('Different number of grids')
plt.xlabel('x')
plt.ylabel('pressure')
time = 80

plt.plot(np.linspace(0, Length, 3+1), pressure_ndx_3[:, time])
plt.plot(np.linspace(0, Length, 6+1), pressure_ndx_6[:, time])
plt.plot(np.linspace(0, Length, 12+1), pressure_ndx_12[:, time])

###### different time steps #######


pressure_dt_1 = solve_pressure_1D(12, Length, num_dt, 0.001)
pressure_dt_2 = solve_pressure_1D(12, Length, num_dt, 0.005)
pressure_dt_3 = solve_pressure_1D(12, Length, num_dt, 0.0001)

plt.figure()
plt.title('Different number of grids')
plt.xlabel('x')
plt.ylabel('pressure')
time = 50

x = np.linspace(0, Length, 12+1)
plt.plot(x, pressure_dt_1[:, time])
plt.plot(x, pressure_dt_2[:, time])
plt.plot(x, pressure_dt_3[:, time])
