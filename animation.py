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

# class Particle:
#     def __init__(self, x, y, m):
#         self.x = x
#         self.y = y
#         self.m = m
#         self.velocity_x = 0
#         self.velocity_y = 0
#
#
# class Spring:
#     def __init__(self, spring_coefficient, damping_coefficient, particle_1, particle_2):
#         self.spring_coefficient = spring_coefficient
#         self.damping_coefficient = damping_coefficient
#         self.particle_1 = particle_1
#         self.particle_2 = particle_2
#
#         self.original_distance = math.sqrt((particle_2.x - particle_1.x) ** 2 + (particle_2.y - particle_1.y) ** 2)
#         self.original_distance_x = abs(particle_2.x - particle_1.x)
#         self.original_distance_y = abs(particle_2.y - particle_1.y)
#
#     input_force = 0
#
#
# def spring_displacement(spring):
#     # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
#     displacement = math.sqrt((spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2) - spring.original_distance
#     # print(spring.original_distance)
#     return displacement
#
#
# def spring_displacement_x(spring):
#     # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
#     displacement = abs(spring.particle_2.x - spring.particle_1.x) - spring.original_distance_x
#     # print(spring.original_distance)
#     return displacement
#
#
# def spring_displacement_y(spring):
#     # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
#     displacement = abs(spring.particle_2.y - spring.particle_1.y) - spring.original_distance_y
#     # print(displacement)
#     return displacement
#
#
# def acceleration(spring, particle):
#
#     # if particle.x == spring.particle_1.x and particle.y == spring.particle_1.y:
#     #     sin_angle = (spring.particle_2.y - spring.particle_1.y) / math.sqrt(
#     #         (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#     #     cos_angle = (spring.particle_2.x - spring.particle_1.x) / math.sqrt(
#     #         (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#     #
#     # else:
#     #     sin_angle = (spring.particle_1.y - spring.particle_2.y) / math.sqrt(
#     #         (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)
#     #     cos_angle = (spring.particle_1.x - spring.particle_2.x) / math.sqrt(
#     #         (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)
#     #
#     # # sin_angle = (spring.particle_2.y - spring.particle_1.y) / math.sqrt((spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#     # # cos_angle = (spring.particle_2.x - spring.particle_1.x) / math.sqrt((spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#     #
#     # spring_force = spring.spring_coefficient * spring_displacement(spring)  # Fs = k * x
#     # damper_force_x = - spring.damping_coefficient * particle.velocity_x + spring.input_force * cos_angle    # Fb = b * x'
#     # damper_force_y = - spring.damping_coefficient * particle.velocity_y + spring.input_force * sin_angle
#     # # print(spring.input_force)
#     # # print(sin_angle)
#     #
#     # acceleration_x = (spring_force * cos_angle + damper_force_x) / particle.m
#     # acceleration_y = (spring_force * sin_angle + damper_force_y) / particle.m
#     #
#     # return acceleration_x, acceleration_y
#     if particle == spring.particle_1:
#         sin_angle = (spring.particle_2.y - spring.particle_1.y) / math.sqrt(
#             (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#         cos_angle = (spring.particle_2.x - spring.particle_1.x) / math.sqrt(
#             (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
#
#         velocity_x = spring.particle_1.velocity_x - spring.particle_2.velocity_x
#         velocity_y = spring.particle_1.velocity_y - spring.particle_2.velocity_y
#
#     else:
#         sin_angle = (spring.particle_1.y - spring.particle_2.y) / math.sqrt(
#             (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)
#         cos_angle = (spring.particle_1.x - spring.particle_2.x) / math.sqrt(
#             (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)
#
#         velocity_x = spring.particle_2.velocity_x - spring.particle_1.velocity_x
#         velocity_y = spring.particle_2.velocity_y - spring.particle_1.velocity_y
#
#     # if abs(velocity_x) < 1:
#     #     velocity_x = 0
#
#     spring_force_x = spring.spring_coefficient * spring_displacement(spring) * cos_angle  # Fs = k * x
#     spring_force_y = spring.spring_coefficient * spring_displacement(spring) * sin_angle  # Fs = k * x
#
#     damper_force_x = - spring.damping_coefficient * velocity_x  # Fb = b * x'
#     damper_force_y = - spring.damping_coefficient * velocity_y
#
#     Ft_x = spring.input_force * cos_angle
#     Ft_y = spring.input_force * sin_angle
#
#     F_tot_x = spring_force_x + damper_force_x + Ft_x
#     F_tot_y = spring_force_y + damper_force_y + Ft_y + particle.m * g
#
#     if particle.y < 0:
#         if F_tot_x < 0 and abs(F_tot_x) > abs(F_tot_y) * f:
#             F_tot_x = F_tot_x + abs(F_tot_y) * f
#         if F_tot_x > 0 and abs(F_tot_x) > abs(F_tot_y) * f:
#             F_tot_x = F_tot_x - abs(F_tot_y) * f
#         else:
#             F_tot_x = 0
#
#     # print(velocity_y)
#     # acceleration_x = (spring_force_x + damper_force_x + Ft_x) / particle.m
#     acceleration_x = F_tot_x / particle.m
#     acceleration_y = (spring_force_y + damper_force_y + Ft_y) / particle.m
#
#     return acceleration_x, acceleration_y
#
#
# def coordinate(spring_1, spring_2, particle):
#     (X_1, Y_1) = acceleration(spring_1, particle)
#     (X_2, Y_2) = acceleration(spring_2, particle)
#
#     acceleration_x = X_1 + X_2
#     acceleration_y = Y_1 + Y_2 + g
#
#     particle.velocity_x += (acceleration_x * dt)  # Integral(a) = v
#     particle.velocity_y += (acceleration_y * dt)
#
#     particle.x += (particle.velocity_x * dt)  # Integral(v) = x
#     particle.y += (particle.velocity_y * dt)
#
#     if particle.y < 0 and particle.velocity_y < 0:
#         particle.velocity_y = - particle.velocity_y * 0.9
#         particle.y = -particle.y
#
#     return particle.x, particle.y
#
#
# def energy_tot(spring_1, spring_2, spring_3, particle_1, particle_2, particle_3):
#
#      tot_energy = 0.5 * spring_1.spring_coefficient * (spring_displacement(spring_1) * 0.01)**2 \
#              + 0.5 * spring_2.spring_coefficient * (spring_displacement(spring_2) * 0.01)**2 \
#              + 0.5 * spring_3.spring_coefficient * (spring_displacement(spring_3) * 0.01)**2 \
#              + particle_1.m * particle_1.y * 0.01 * -g \
#              + particle_2.m * particle_2.y * 0.01 * -g \
#              + particle_3.m * particle_3.y * 0.01 * -g \
#              + 0.5 * particle_1.m * (particle_1.velocity_x**2 + particle_1.velocity_y**2)\
#              + 0.5 * particle_2.m * (particle_2.velocity_x**2 + particle_2.velocity_y**2)\
#              + 0.5 * particle_3.m * (particle_3.velocity_x**2 + particle_3.velocity_y**2)
#
#      return tot_energy
#
#
# def energy_kn(particle_1, particle_2, particle_3):
#     kn_energy = 0.5 * particle_1.m * (particle_1.velocity_x ** 2 + particle_1.velocity_y ** 2) \
#                  + 0.5 * particle_2.m * (particle_2.velocity_x ** 2 + particle_2.velocity_y ** 2) \
#                  + 0.5 * particle_3.m * (particle_3.velocity_x ** 2 + particle_3.velocity_y ** 2)
#
#     return kn_energy
#
#
# def energy_pt(spring_1, spring_2, spring_3, particle_1, particle_2, particle_3):
#     pt_energy = 0.5 * spring_1.spring_coefficient * (spring_displacement(spring_1) * 0.01) ** 2 \
#                  + 0.5 * spring_2.spring_coefficient * (spring_displacement(spring_2) * 0.01) ** 2 \
#                  + 0.5 * spring_3.spring_coefficient * (spring_displacement(spring_3) * 0.01) ** 2 \
#                  + particle_1.m * particle_1.y * 0.01 * -g \
#                  + particle_2.m * particle_2.y * 0.01 * -g \
#                  + particle_3.m * particle_3.y * 0.01 * -g
#     return pt_energy


'''
#########################################################################
Main
########################################################################
'''


def initpos(x, y, l, theta):

    l = math.sqrt(3) / 3 * l

    x1 = x - l * math.sin(theta + 1 / 3 * math.pi)
    y1 = y - l * math.cos(theta + 1 / 3 * math.pi)

    x2 = x + l * math.sin(theta)
    y2 = y + l * math.cos(theta)

    x3 = x + l * math.cos(theta + 1 / 6 * math.pi)
    y3 = y - l * math.sin(theta + 1 / 6 * math.pi)

    return x1, y1, x2, y2, x3, y3


def ad_force(var, fun, a1, a2, a3):
    force = a1 * fun(a2 * var) + a3
    return force


mass = 1
k = 0.5   # spring_coefficient
b = 0.8  # damping_coefficient
a1 = 20
a2 = 1
a3 = 4
an = -math.pi / 3
le = 40
h = 50

# P1 = Particle(0, 50, mass)  # X, Y, mass
# P2 = Particle(20, 20, mass)
# P3 = Particle(40, 50, mass)
POS = initpos(0, h, le, an)  # x, y, edge_length, rot angle(rad)
P1 = Particle(POS[0], POS[1], mass)  # X, Y, mass
P2 = Particle(POS[2], POS[3], mass)
P3 = Particle(POS[4], POS[5], mass)


S12 = Spring(k, b, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
S23 = Spring(k, b, P2, P3)
S13 = Spring(k, b, P1, P3)

f = 0.8  # friction coefficient
g = -9.81  # gravity
dt = 0.05  # differential time
stop_time = 30  # tot running time
time = 0  # time mark

t = np.arange(0, stop_time, dt)

x1y1 = []
x2y2 = []
x3y3 = []
energy = []
energy_k = []
energy_p = []
real_time = []
centre_of_mass = []

for i in t:
    # additional force
    S12.input_force = ad_force(i, math.sin, a1, a2, a3)
    # S12.input_force = ad_force(i, math.sin, 50, 2, 0)
    x1y1.append(coordinate(S12, S13, P1))
    x2y2.append(coordinate(S12, S23, P2))
    x3y3.append(coordinate(S13, S23, P3))
    # energy.append(energy_tot(S12, S23, S13, P1, P2, P3))
    # energy_k.append(energy_kn(P1, P2, P3))
    # energy_p.append(energy_pt(S12, S23, S13, P1, P2, P3))
    A = (P1.x + P2.x) / 2
    B = 2 / 3 * A + 1 / 3 * (P3.x)
    centre_of_mass.append(B)
    real_time.append(i)


x1 = []; x2 = []; x3 = []
y1 = []; y2 = []; y3 = []

for i in range(len(x1y1)):
    x1.append(x1y1[i][0])
    x2.append(x2y2[i][0])
    x3.append(x3y3[i][0])
    y1.append(x1y1[i][1])
    y2.append(x2y2[i][1])
    y3.append(x3y3[i][1])

# print(x1)

'''
########################################################
                    DRAWING
########################################################
'''
# plt.figure()

# plt.plot(real_time, energy, label = 'Total Energy')
# plt.plot(real_time, energy_k, label = 'Kinetic Energy')
# plt.plot(real_time, energy_p, label = 'Potential Energy')
#
# plt.legend(loc='upper right', fontsize=10)


#plt.figure()
fig, ax = plt.subplots()
# passive spring
line, = ax.plot([], [], 'o-', lw=2)
# active spring
line1, = ax.plot([], [], 'ro-', lw=2)

time_template = 'Distance = %.2fcm'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    ax.set_xlim(-200, 200)
    ax.set_ylim(0, 100)
    return line,


def update(i):

    thisx = [x2[i], x3[i], x1[i]]
    thisy = [y2[i], y3[i], y1[i]]

    thisx1 = [x1[i], x2[i]]
    thisy1 = [y1[i], y2[i]]

    line1.set_data(thisx1, thisy1)
    line.set_data(thisx, thisy)

    # if x1[i] > x_last:
    #     x_max = x1[i]
    #     x_last = x1[i]
    # else:

    time_text.set_text(time_template % (centre_of_mass[i]))

    return line, line1, time_text


ani = animation.FuncAnimation(
    fig =fig,
    func=update,
    init_func=init,
    #frames = np.linspace(0, stop_time / dt, stop_time / dt),
    frames=len(x1),
    interval = dt*1000,
    repeat = False,
    blit = False
)

plt.axis("scaled")
plt.show()
# ani.save('motion1.gif', writer='pillow', dpi=40) #save
###################################################################################
