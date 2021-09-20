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

def ad_force(var, fun, a1, a2, a3):
    force = a1 * fun(a2 * var) + a3
    return force


def initpos(x, y, l, theta):
    l = math.sqrt(3) / 3 * l

    x1 = x - l * math.sin(theta + 1 / 3 * math.pi)
    y1 = y - l * math.cos(theta + 1 / 3 * math.pi)

    x2 = x + l * math.sin(theta)
    y2 = y + l * math.cos(theta)

    x3 = x + l * math.cos(theta + 1 / 6 * math.pi)
    y3 = y - l * math.sin(theta + 1 / 6 * math.pi)

    return x1, y1, x2, y2, x3, y3






mass = np.arange(0.3, 3, 0.3)
mass = np.around(mass, decimals=2)

spring_coefficient = np.arange(0.5, 5, 0.5)  # spring_coefficient
spring_coefficient = np.around(spring_coefficient, decimals=1)

damping_coefficient = np.arange(0.1, 1, 0.1)  # damping_coefficient
damping_coefficient = np.around(damping_coefficient, decimals=1)

hight = np.arange(50, 100, 5)
hight = np.around(hight, decimals=2)

angle = np.arange(-math.pi / 3, 5*math.pi / 3, math.pi / 6)
angle = np.around(angle, decimals=3)

length = np.arange(10, 40, 3)
length = np.around(length, decimals=1)

a1R = np.arange(1, 30, 3)
a1R = np.around(a1R, decimals=3)

a2R = np.arange(1, 30, 3)
a2R = np.around(a2R, decimals=3)

a3R = np.arange(1, 30, 3)
a3R = np.around(a3R, decimals=3)

mass_plot_x = []
mass_plot_y = []
k_plot_x = []
k_plot_y = []
b_plot_x = []
b_plot_y = []
t_plot_x = []
t_plot_y = []

centre_of_mass = []
V_ms = []
x1y1_mv = []
av_x = []
for sc_v in spring_coefficient:
    for a3 in a3R:
        m = 1
        dc_v = 0.8
        # sc_v = 2.8
        an = -math.pi / 3
        le = 40
        h = 50
        a2 = 1
        a1 = 20
        POS = initpos(0, h, le, -math.pi / 3)  # x, y, edge_length, rot angle(rad)
        P1 = Particle(POS[0], POS[1], m)  # X, Y, mass
        P2 = Particle(POS[2], POS[3], m)
        P3 = Particle(POS[4], POS[5], m)

        S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
        S23 = Spring(sc_v, dc_v, P2, P3)
        S13 = Spring(sc_v, dc_v, P1, P3)

        for i in t:
            if i == 0:
                x1y1_mv.clear()
            S12.input_force = ad_force(i, math.sin, a1, a2, a3)
            coordinate(S12, S13, P1)
            coordinate(S12, S23, P2)
            coordinate(S13, S23, P3)
            B = 1 / 3 * (P1.x + P2.x + P3.x)
            centre_of_mass.append(B)
        av_x.append((centre_of_mass[-1] - centre_of_mass[0]) / stop_time)


av_x = np.array(av_x)
#print(av_x)

#matrix = av_x.reshape(len(spring_coefficient),len(mass))
matrix = av_x.reshape(len(a3R),len(spring_coefficient))


fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(matrix, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
cb = fig.colorbar(cax)
cb.set_label('Average velocity along X-axis (mm/s)', size=16)

# labels_x = mass
# loc_x = np.array(range(0,len(mass)))
labels_x = a3R
loc_x = np.array(range(0,len(a3R)))



labels_y = spring_coefficient
loc_y = np.array(range(0,len(spring_coefficient)))




plt.xticks(loc_y, labels_y, size=16)
plt.yticks(loc_x, labels_x, size=16)
# plt.xlabel("k (N/mm)")
plt.xlabel("k (N/mm)", size=16)


plt.ylabel("a3", size=16)


ax.xaxis.set_major_locator(ticker.MultipleLocator(2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(2))

plt.show()

####################################################################################################



# x1y1_md = []
# bv_x = []
# for m in mass:
#     for dc_v in damping_coefficient:
#         sc_v = 2  # damping_coefficient
#         P1 = Particle(0, 20, m)  # X, Y, mass
#         P2 = Particle(20, 50, m)
#         P3 = Particle(40, 20, m)
#
#         S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
#         S23 = Spring(sc_v, dc_v, P2, P3)
#         S13 = Spring(sc_v, dc_v, P1, P3)
#
#         for i in t:
#             # additional force
#             # S12.input_force = 80 * math.sin(i)
#             if i == 0:
#                 x1y1_md.clear()
#             coordinate(S12, S13, P1)
#             x1y1_md.append(P1.velocity_x)
#         bv_x.append(np.mean(x1y1_md))
#
# bv_x = np.array(bv_x)
# #print(av_x)
#
# #matrix = av_x.reshape(len(spring_coefficient),len(mass))
# matrix1 = bv_x.reshape(len(mass),len(damping_coefficient))
#
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# cax1 = ax1.matshow(matrix1, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
# fig1.colorbar(cax1)
#
# labels_x = mass
# loc_x = np.array(range(0,len(mass)))
# labels_y = damping_coefficient
# loc_y = np.array(range(0,len(damping_coefficient)))
# plt.xticks(loc_y, labels_y)
# plt.yticks(loc_x, labels_x)
# plt.xlabel("damping_cof")
# plt.ylabel("mass")
#
#
# ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
# ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
##############################################################################


# x1y1_sd = []
# cv_x = []
# for sc_v in spring_coefficient:
#     for dc_v in damping_coefficient:
#         m = 2  # damping_coefficient
#         P1 = Particle(0, 20, m)  # X, Y, mass
#         P2 = Particle(20, 50, m)
#         P3 = Particle(40, 20, m)
#
#         S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
#         S23 = Spring(sc_v, dc_v, P2, P3)
#         S13 = Spring(sc_v, dc_v, P1, P3)
#
#         for i in t:
#             # additional force
#             # S12.input_force = 80 * math.sin(i)
#             if i == 0:
#                 x1y1_sd.clear()
#             coordinate(S12, S13, P1)
#             x1y1_sd.append(P1.velocity_x)
#         cv_x.append(np.mean(x1y1_sd))
#
# cv_x = np.array(cv_x)
# #print(av_x)
#
# #matrix = av_x.reshape(len(spring_coefficient),len(mass))
# matrix1 = cv_x.reshape(len(spring_coefficient),len(damping_coefficient))
#
#
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# cax1 = ax1.matshow(matrix1, cmap=plt.get_cmap('Greens'), alpha=0.5, origin='upper')  # , alpha=0.3
# fig1.colorbar(cax1)
#
# labels_x = spring_coefficient
# loc_x = np.array(range(0,len(spring_coefficient)))
# labels_y = damping_coefficient
# loc_y = np.array(range(0,len(damping_coefficient)))
# plt.xticks(loc_y, labels_y)
# plt.yticks(loc_x, labels_x)
# plt.xlabel("damping_cof")
# plt.ylabel("spring_coefficient")
#
#
# ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))
# ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
#
# plt.show()