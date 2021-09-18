from matplotlib import pyplot as plt
from matplotlib.pyplot import MultipleLocator
import matplotlib.ticker as ticker
from matplotlib import animation
import numpy as np
import math
import threading
from SpringSimulator.Simulator import Spring
from SpringSimulator.Simulator import Particle
from SpringSimulator.Simulator import coordinate

f = 0.8  # friction coefficient
g = -9.81  # gravity
dt = 0.05  # differential time
stop_time = 30  # tot running time
time = 0  # time mark
t = np.arange(0, stop_time, dt)

mass = np.arange(1, 5, 0.1)
mass = np.around(mass, decimals=2)
spring_coefficient = np.arange(1, 4, 0.1)  # spring_coefficient
spring_coefficient = np.around(spring_coefficient, decimals=1)
damping_coefficient = np.arange(0.1, 3, 0.05)  # damping_coefficient
damping_coefficient = np.around(damping_coefficient, decimals=1)


mass_plot_x = []
mass_plot_y = []
k_plot_x = []
k_plot_y = []
b_plot_x = []
b_plot_y = []
t_plot_x = []
t_plot_y = []


V_ms = []
x1y1_mv = []
av_x = []
for m in mass:
    for sc_v in spring_coefficient:
        dc_v = 0.3  # damping_coefficient
        P1 = Particle(0, 20, m)  # X, Y, mass
        P2 = Particle(20, 50, m)
        P3 = Particle(40, 20, m)

        S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
        S23 = Spring(sc_v, dc_v, P2, P3)
        S13 = Spring(sc_v, dc_v, P1, P3)

        for i in t:
            # additional force
            # S12.input_force = 80 * math.sin(i)
            if i == 0:
                x1y1_mv.clear()
            coordinate(S12, S13, P1)
            x1y1_mv.append(P1.velocity_x)
        av_x.append(np.mean(x1y1_mv))

av_x = np.array(av_x)
#print(av_x)

#matrix = av_x.reshape(len(spring_coefficient),len(mass))
matrix = av_x.reshape(len(mass),len(spring_coefficient))


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(matrix, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
fig.colorbar(cax)

labels_x = mass
loc_x = np.array(range(0,len(mass)))
labels_y = spring_coefficient
loc_y = np.array(range(0,len(spring_coefficient)))
plt.xticks(loc_y, labels_y)
plt.yticks(loc_x, labels_x)
plt.xlabel("spring_cof")
plt.ylabel("mass")


ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))



####################################################################################################



x1y1_md = []
bv_x = []
for m in mass:
    for dc_v in damping_coefficient:
        sc_v = 2  # damping_coefficient
        P1 = Particle(0, 20, m)  # X, Y, mass
        P2 = Particle(20, 50, m)
        P3 = Particle(40, 20, m)

        S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
        S23 = Spring(sc_v, dc_v, P2, P3)
        S13 = Spring(sc_v, dc_v, P1, P3)

        for i in t:
            # additional force
            # S12.input_force = 80 * math.sin(i)
            if i == 0:
                x1y1_md.clear()
            coordinate(S12, S13, P1)
            x1y1_md.append(P1.velocity_x)
        bv_x.append(np.mean(x1y1_md))

bv_x = np.array(bv_x)
#print(av_x)

#matrix = av_x.reshape(len(spring_coefficient),len(mass))
matrix1 = bv_x.reshape(len(mass),len(damping_coefficient))


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
cax1 = ax1.matshow(matrix1, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
fig1.colorbar(cax1)

labels_x = mass
loc_x = np.array(range(0,len(mass)))
labels_y = damping_coefficient
loc_y = np.array(range(0,len(damping_coefficient)))
plt.xticks(loc_y, labels_y)
plt.yticks(loc_x, labels_x)
plt.xlabel("damping_cof")
plt.ylabel("mass")


ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
##############################################################################


x1y1_sd = []
cv_x = []
for sc_v in spring_coefficient:
    for dc_v in damping_coefficient:
        m = 2  # damping_coefficient
        P1 = Particle(0, 20, m)  # X, Y, mass
        P2 = Particle(20, 50, m)
        P3 = Particle(40, 20, m)

        S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
        S23 = Spring(sc_v, dc_v, P2, P3)
        S13 = Spring(sc_v, dc_v, P1, P3)

        for i in t:
            # additional force
            # S12.input_force = 80 * math.sin(i)
            if i == 0:
                x1y1_sd.clear()
            coordinate(S12, S13, P1)
            x1y1_sd.append(P1.velocity_x)
        cv_x.append(np.mean(x1y1_sd))

cv_x = np.array(cv_x)
#print(av_x)

#matrix = av_x.reshape(len(spring_coefficient),len(mass))
matrix1 = cv_x.reshape(len(spring_coefficient),len(damping_coefficient))


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
cax1 = ax1.matshow(matrix1, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
fig1.colorbar(cax1)

labels_x = spring_coefficient
loc_x = np.array(range(0,len(spring_coefficient)))
labels_y = damping_coefficient
loc_y = np.array(range(0,len(damping_coefficient)))
plt.xticks(loc_y, labels_y)
plt.yticks(loc_x, labels_x)
plt.xlabel("damping_cof")
plt.ylabel("spring_coefficient")


ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))

plt.show()