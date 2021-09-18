import numpy as np
from global_var import parameters as pa
from simulator import position
from simulator.position import Particle
from simulator.position import Spring
from simulator.position import AddForce
from simulator.position import cal
from simulator.position import ft


from scipy.integrate import odeint, solve_bvp, solve_ivp
from numpy import sin, cos, linspace, pi
import matplotlib.pyplot as plt
import math


##### PARAMETERS SETTING #######################################################################
pa.f = 0.8  # friction coefficient
pa.g = -9.81  # gravity
pa.dt = 0.02  # differential time
pa.stop_time = 10  # tot running time
pa.time = 0  # time mark
pa.cel = 0.9 # Coefficient of collision energy loss


force_p1 = 0
force_p2 = 1
force_p3 = 0

####   FUNC DEFINE ########################################################################


def ad_force(var, fun, a1, a2, a3):
    force = a1 * fun(a2 * var) + a3
    return force


def sina(a, b, c, d): # x1, y1, x2, y2
    fx = (b - d)/math.sqrt((a - c)**2 + (b - d)**2)
    return fx


def cosa(a, b, c, d): # x1, y1, x2, y2
    fx = (a - c)/math.sqrt((a - c)**2 + (b - d)**2)
    return fx


def ysecd(py1, py2, py3, vy1, spring1, spring2, part, sin1, sin2, t):
    fx = (- spring1.damping_coefficient * vy1
            - spring1.spring_coefficient * ((abs(py1- py2) - spring1.original_distance_y))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3) * sin1
            - spring2.damping_coefficient * vy1
            - spring2.spring_coefficient * ((abs(py1 - py3) - spring2.original_distance_y))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3) * sin2
            + part.m * pa.g
            ) / part.m
    return fx


def xsecd(px1, px2, px3, vx1, spring1, spring2, part, cos1, cos2, t):
    fx = (- spring1.damping_coefficient * vx1
            - spring1.spring_coefficient * ((abs(px1- px2) - spring1.original_distance_x))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3) * cos1
            - spring2.damping_coefficient * vx1
            - spring2.spring_coefficient * ((abs(px1 - px3) - spring2.original_distance_x))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3) * cos2
            ) / part.m
    return fx


def initpos(x, y, l, theta):

    l = math.sqrt(3) / 3 * l

    x1 = x - l * math.sin(theta + 1 / 3 * math.pi)
    y1 = y - l * math.cos(theta + 1 / 3 * math.pi)

    x2 = x + l * math.sin(theta)
    y2 = y + l * math.cos(theta)

    x3 = x + l * math.cos(theta + 1 / 6 * math.pi)
    y3 = y - l * math.sin(theta + 1 / 6 * math.pi)

    return x1, y1, x2, y2, x3, y3


mass = 1.8
k = 4   # spring_coefficient
b = 0.2

POS = initpos(0, 50, 40, -math.pi/3)  # x, y, edge_length, rot angle(rad)
P1 = Particle(POS[0], POS[1], mass)  # X, Y, mass
P2 = Particle(POS[2], POS[3], mass)
P3 = Particle(POS[4], POS[5], mass)


S12 = Spring(k, b, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
S23 = Spring(k, b, P2, P3)
S13 = Spring(k, b, P1, P3)


def func(t,r):
    px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3 = r  # position and volocity

    ############ P1 ##############################
    sin12 = sina(px1, py1, px2, py2)
    sin13 = sina(px1, py1, px3, py3)
    cos12 = cosa(px1, py1, px2, py2)
    cos13 = cosa(px1, py1, px3, py3)

    if py1 < 0 and vy1 < 0:  # bounce
        vy1 = - vy1 * pa.cel
        py1 = 0

    # first order P1y, velocity
    fy1f = vy1
    # second order P1y, accle
    fy1s = (#S12
            (- S12.damping_coefficient * vy1
            + S12.spring_coefficient * ((abs(py1- py2) - S12.original_distance_y))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3)) * sin12  #var, fun, a1, a2, a3

            # S13
            +(- S13.damping_coefficient * vy1
            + S13.spring_coefficient * ((abs(py1 - py3) - S13.original_distance_y)))*sin13

           #mg
            + P1.m * pa.g
           ) / P1.m

    fx1f = vx1  # first order P1x
    fx1s = (#S12
            (- S12.damping_coefficient * vx1
            + S12.spring_coefficient * ((abs(px1- px2) - S12.original_distance_x))
            + ad_force(t, np.sin, force_p1, force_p2, force_p3)) * cos12
            #S13
            +(- S13.damping_coefficient * vx1
            + S13.spring_coefficient * ((abs(px1 - px3) - S13.original_distance_x))) * cos13
            ) / P1.m


    if 0 <= py1 < 0.1:  # frinction
        if fx1s < 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
            fx1s = fx1s + abs(fy1s * pa.f)
        if fx1s > 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
            fx1s = fx1s - abs(fy1s * pa.f)
        else:
            fx1f = 0
            fx1s = 0

    ############ P2 ##############################
    sin21 = sina(px2, py2, px1, py1)
    sin23 = sina(px2, py2, px3, py3)
    cos21 = cosa(px2, py2, px1, py1)
    cos23 = cosa(px2, py2, px3, py3)

    if py2 < 0 and vy2 < 0:  # bounce
        vy2 = - vy2 * pa.cel
        py2 = 0

    # first order P2y, velocity
    fy2f = vy2
    # second order P2y, accleration
    fy2s = (       # S12
                   (- S12.damping_coefficient * vy2
                    + S12.spring_coefficient * (abs(py1 - py2) - S12.original_distance_y)
                    + ad_force(t, np.sin, force_p1, force_p2, force_p3)) * sin21  # var, fun, a1, a2, a3

                   # S23
                   + (- S23.damping_coefficient * vy1
                      + S23.spring_coefficient * (abs(py2 - py3) - S23.original_distance_y)) * sin23

                   # mg
                   + P2.m * pa.g
           ) / P2.m

    fx2f = vx2  # first order P2x
    fx2s = (  # S12
                   (- S12.damping_coefficient * vx2
                    + S12.spring_coefficient * (abs(px1 - px2) - S12.original_distance_x)
                    + ad_force(t, np.sin, force_p1, force_p2, force_p3)) * cos21
                   # S23
                   + (- S23.damping_coefficient * vx2
                      + S23.spring_coefficient * (abs(px2 - px3) - S23.original_distance_x)) * cos23
           ) / P2.m

    if 0 <= py2 < 0.1:  # frinction
        if fx2s < 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
            fx2s = fx2s + abs(fy2s * pa.f)
        if fx2s > 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
            fx2s = fx2s - abs(fy2s * pa.f)
        else:
            fx2f = 0
            fx2s = 0

    ############ P1 ##############################
    sin31 = sina(px3, py3, px1, py1)
    sin32 = sina(px3, py3, px2, py2)
    cos31 = cosa(px3, py3, px1, py1)
    cos32 = cosa(px3, py3, px2, py2)

    if py3 < 0 and vy3 < 0:  # bounce
        vy3 = - vy3 * pa.cel
        py3 = 0

    # first order P1y, velocity
    fy3f = vy3
    # second order P1y, accle
    fy3s = (  # S13
                   (- S13.damping_coefficient * vy3
                    + S13.spring_coefficient * (abs(py1 - py3) - S13.original_distance_y)) * sin31  # var, fun, a1, a2, a3

                   # S23
                   + (- S23.damping_coefficient * vy3
                      + S23.spring_coefficient * (abs(py2 - py3) - S13.original_distance_y)) * sin32

                   # mg
                   + P3.m * pa.g
           ) / P3.m

    fx3f = vx3  # first order P1x
    fx3s = (  # S12
                   (- S13.damping_coefficient * vx3
                    + S13.spring_coefficient * (
                                abs(px1 - px3) - S13.original_distance_x)) * cos31  # var, fun, a1, a2, a3

                   # S23
                   + (- S23.damping_coefficient * vx3
                      + S23.spring_coefficient * (abs(px2 - px3) - S13.original_distance_x)) * cos32
           ) / P3.m

    if 0 <= py3 < 0.1:  # frinction
        if fx3s < 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
            fx3s = fx3s + abs(fy3s * pa.f)
        if fx3s > 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
            fx3s = fx3s - abs(fy3s * pa.f)
        else:
            fx3f = 0
            fx3s = 0

    return fx1f, fx1s, fy1f, fy1s, fx2f, fx2s, fy2f, fy2s, fx3f, fx3s, fy3f, fy3s

# def func(t, r):
#
#     px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3 = r  # position and volocity
#
#     sin12 = sina(px1, py1, px2, py2)
#     sin13 = sina(px1, py1, px3, py3)
#     cos12 = cosa(px1, py1, px2, py2)
#     cos13 = cosa(px1, py1, px3, py3)
#
#     if py1 < 0 and vy1 < 0:  # bounce
#         vy1 = - vy1 * pa.cel
#         py1 = - py1
#
#     fy1f = vy1  # first order P1y
#     fy1s = ysecd(py1, py2, py3, vy1, S12, S13, P1, sin12, sin13, t)
#
#     fx1f = vx1  # first order P1x
#     fx1s = xsecd(px1, px2, px3, vx1, S12, S13, P1, cos12, cos13, t)
#
#     if 0 <= py1 < 0.1:  # frinction
#         if fx1s < 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
#             fx1s = fx1s + abs(fy1s * pa.f)
#         if fx1s > 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
#             fx1s = fx1s - abs(fy1s * pa.f)
#         else:
#             fx1s = 0
#
#
#     ########################################## x2 y2
#     sin21 = sina(px2, py2, px1, py1)
#     sin23 = sina(px2, py2, px3, py3)
#     cos21 = cosa(px2, py2, px3, py3)
#     cos23 = cosa(px2, py2, px3, py3)
#
#     if py2 < 0 and vy2 < 0:  # bounce
#         vy1 = - vy2 * pa.cel
#     py2 = - py2
#
#     fy2f = vy2  # first order P1y
#     fy2s = ysecd(py2, py1, py3, vy2, S12, S23, P2, sin21, sin23,t)
#
#     fx2f = vx2  # first order P1x
#     fx2s = xsecd(px2, px1, px3, vx2, S12, S23, P2, cos21, cos23,t)
#
#     if py2 < 0.1 and py2 >= 0:  # frinction
#         if fx2s < 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
#             fx2s = fx2s + abs(fy2s * pa.f)
#         if fx2s > 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
#             fx2s = fx2s - abs(fy2s * pa.f)
#         else:
#             fx2s = 0
#     #######################################  x3 y3
#     sin31 = sina(px3, py3, px1, py1)
#     sin32 = sina(px3, py3, px1, py1)
#     cos31 = cosa(px3, py3, px1, py1)
#     cos32 = cosa(px3, py3, px1, py1)
#
#     if py3 < 0 and vy3 < 0:  # bounce
#         vy3 = - vy3 * pa.cel
#         py3 = - py3
#
#     fy3f = vy3  # first order P1y
#     fy3s = ysecd(py3, py1, py2, vy3, S23, S13, P3, sin31, sin32, t)
#
#     fx3f = vx3  # first order P1x
#     fx3s = xsecd(px3, px1, px2, vx3, S23, S13, P3, cos31, cos32, t)
#
#     if py3 < 0.1 and py3 >= 0:  # frinction
#         if fx3s < 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
#             fx3s = fx3s + abs(fy3s * pa.f)
#         if fx3s > 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
#             fx3s = fx3s - abs(fy3s * pa.f)
#         else:
#             fx3s = 0
#
#     P1.x = px1
#     P1.y = py1
#     P2.x = px2
#     P2.y = py2
#     P3.x = px3
#     P3.y = py3
#
#     return fx1f, fx1s, fy1f, fy1s, fx2f, fx2s, fy2f, fy2s, fx3f, fx3s, fy3f, fy3s


px1 = P1.x
vx1 = 0
py1 = P1.y
vy1 = 0
px2 = P2.x
vx2 = 0
py2 = P2.y
vy2 = 0
px3 = P3.x
vx3 = 0
py3 = P3.y
vy3 = 0

t1 = linspace(0, pa.stop_time, int((pa.stop_time/pa.dt)))
y0 = [px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3]

sol = solve_ivp(func, (0, pa.stop_time), y0, method='RK23', t_eval=t1)
X1, VX1, Y1, VY1, X2, VX2, Y2, VY2, X3, VX3, Y3, VY3 = sol.y

plt.plot(sol.y[0,:], sol.y[2,:],'g--',label='px1')
plt.plot(sol.y[5,:], sol.y[7,:],'r-',label='vx1')
plt.plot(sol.y[9,:], sol.y[11,:],'r-',label='vx1')
plt.show()


