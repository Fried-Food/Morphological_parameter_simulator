from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import math
import threading


f = 0.5  # friction coefficient
g = -9.81  # gravity
dt = 0.05  # differential time
stop_time = 30  # tot running time
time = 0  # time mark
t = np.arange(0, stop_time, dt)


class Parameters:
    def __init__(self, f, g, dt, stop_time):
        f = 0.5  # friction coefficient
        g = -9.81  # gravity
        dt = 0.05  # differential time


class Particle:
    def __init__(self, x, y, m):
        self.x = x
        self.y = y
        self.m = m
        self.velocity_x = 0
        self.velocity_y = 0


class Spring:
    def __init__(self, spring_coefficient, damping_coefficient, particle_1, particle_2):
        self.spring_coefficient = spring_coefficient
        self.damping_coefficient = damping_coefficient
        self.particle_1 = particle_1
        self.particle_2 = particle_2

        self.original_distance = math.sqrt((particle_2.x - particle_1.x) ** 2 + (particle_2.y - particle_1.y) ** 2)
        self.original_distance_x = abs(particle_2.x - particle_1.x)
        self.original_distance_y = abs(particle_2.y - particle_1.y)

    input_force = 0


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


def spring_displacement(spring):
    # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
    displacement = math.sqrt((spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2) - spring.original_distance
    # print(spring.original_distance)
    return displacement


def spring_displacement_x(spring):
    # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
    displacement = abs(spring.particle_2.x - spring.particle_1.x) - spring.original_distance_x
    # print(spring.original_distance)
    return displacement


def spring_displacement_y(spring):
    # sqrt((x2-x1)^2+(y2-y1)^2) - original distance. The distance between two particles
    displacement = abs(spring.particle_2.y - spring.particle_1.y) - spring.original_distance_y
    # print(displacement)
    return displacement


def acceleration(spring, particle):

    if particle == spring.particle_1:
        sin_angle = (spring.particle_2.y - spring.particle_1.y) / math.sqrt(
            (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)
        cos_angle = (spring.particle_2.x - spring.particle_1.x) / math.sqrt(
            (spring.particle_2.x - spring.particle_1.x) ** 2 + (spring.particle_2.y - spring.particle_1.y) ** 2)

        velocity_x = spring.particle_1.velocity_x - spring.particle_2.velocity_x
        velocity_y = spring.particle_1.velocity_y - spring.particle_2.velocity_y

    else:
        sin_angle = (spring.particle_1.y - spring.particle_2.y) / math.sqrt(
            (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)
        cos_angle = (spring.particle_1.x - spring.particle_2.x) / math.sqrt(
            (spring.particle_1.x - spring.particle_2.x) ** 2 + (spring.particle_1.y - spring.particle_2.y) ** 2)

        velocity_x = spring.particle_2.velocity_x - spring.particle_1.velocity_x
        velocity_y = spring.particle_2.velocity_y - spring.particle_1.velocity_y

    spring_force_x = spring.spring_coefficient * spring_displacement(spring) * cos_angle  # Fs = k * x
    spring_force_y = spring.spring_coefficient * spring_displacement(spring) * sin_angle  # Fs = k * x

    damper_force_x = - spring.damping_coefficient * velocity_x  # Fb = b * x'
    damper_force_y = - spring.damping_coefficient * velocity_y

    Ft_x = spring.input_force * cos_angle
    Ft_y = spring.input_force * sin_angle

    F_tot_x = spring_force_x + damper_force_x + Ft_x
    F_tot_y = spring_force_y + damper_force_y + Ft_y + particle.m * g

    # if particle.y < 0 and particle.velocity_y < 0:
    if -0.1 < particle.y < 0.1 and F_tot_y < 0:
        if F_tot_x < 0 and abs(F_tot_x) > abs(F_tot_y) * f:
            F_tot_x = F_tot_x - abs(F_tot_y) * f
        if F_tot_x > 0 and abs(F_tot_x) > abs(F_tot_y) * f:
            F_tot_x = F_tot_x + abs(F_tot_y) * f
        else:
            F_tot_x = 0
            particle.velocity_x = 0
        # if F_tot_x > 0 and abs(F_tot_x) < abs(F_tot_y) * f:
        #     F_tot_x = F_tot_x + abs(F_tot_y) * f
        # if F_tot_x < 0 and abs(F_tot_x) < abs(F_tot_y) * f:
        #     F_tot_x = F_tot_x - abs(F_tot_y) * f

    acceleration_x = F_tot_x / particle.m
    acceleration_y = (spring_force_y + damper_force_y + Ft_y) / particle.m

    return acceleration_x, acceleration_y


def coordinate(spring_1, spring_2, particle):
    (X_1, Y_1) = acceleration(spring_1, particle)
    (X_2, Y_2) = acceleration(spring_2, particle)

    acceleration_x = X_1 + X_2
    acceleration_y = Y_1 + Y_2 + g

    particle.velocity_x += (acceleration_x * dt)  # Integral(a) = v
    particle.velocity_y += (acceleration_y * dt)

    particle.x += (particle.velocity_x * dt)  # Integral(v) = x
    particle.y += (particle.velocity_y * dt)

    # if (particle.y < 0 and particle.velocity_y).all() < 0:
    # if all([particle.y < 0, particle.velocity_y < 0]):
    if all([particle.y < 0, particle.velocity_y < 0]):
        particle.velocity_y = - particle.velocity_y * 0.9
        particle.y = -particle.y
        # particle.y = 0

    return particle.x, particle.y


def energy_tot(spring_1, spring_2, spring_3, particle_1, particle_2, particle_3):

     tot_energy = 0.5 * spring_1.spring_coefficient * (spring_displacement(spring_1) * 0.01)**2 \
             + 0.5 * spring_2.spring_coefficient * (spring_displacement(spring_2) * 0.01)**2 \
             + 0.5 * spring_3.spring_coefficient * (spring_displacement(spring_3) * 0.01)**2 \
             + particle_1.m * particle_1.y * 0.01 * -g \
             + particle_2.m * particle_2.y * 0.01 * -g \
             + particle_3.m * particle_3.y * 0.01 * -g \
             + 0.5 * particle_1.m * (particle_1.velocity_x**2 + particle_1.velocity_y**2)\
             + 0.5 * particle_2.m * (particle_2.velocity_x**2 + particle_2.velocity_y**2)\
             + 0.5 * particle_3.m * (particle_3.velocity_x**2 + particle_3.velocity_y**2)

     return tot_energy


def energy_kn(particle_1, particle_2, particle_3):
    kn_energy = 0.5 * particle_1.m * (particle_1.velocity_x ** 2 + particle_1.velocity_y ** 2) \
                 + 0.5 * particle_2.m * (particle_2.velocity_x ** 2 + particle_2.velocity_y ** 2) \
                 + 0.5 * particle_3.m * (particle_3.velocity_x ** 2 + particle_3.velocity_y ** 2)

    return kn_energy


def energy_pt(spring_1, spring_2, spring_3, particle_1, particle_2, particle_3):
    pt_energy = 0.5 * spring_1.spring_coefficient * (spring_displacement(spring_1) * 0.01) ** 2 \
                 + 0.5 * spring_2.spring_coefficient * (spring_displacement(spring_2) * 0.01) ** 2 \
                 + 0.5 * spring_3.spring_coefficient * (spring_displacement(spring_3) * 0.01) ** 2 \
                 + particle_1.m * particle_1.y * 0.01 * -g \
                 + particle_2.m * particle_2.y * 0.01 * -g \
                 + particle_3.m * particle_3.y * 0.01 * -g
    return pt_energy