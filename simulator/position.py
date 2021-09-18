import numpy as np
import math
import global_var
from global_var import parameters as pa


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


class AddForce:
    def __init__(self, f, c1, c2):
        self.f = f
        self.c1 = c1
        self.c2 = c2


def ft(t, f, c1, c2):
    fx = c1 * f(t) + c2
    return fx


def sina(a, b, c, d): # x1, y1, x2, y2
    fx = (b - d)/math.sqrt((a - c)**2 + (b - d)**2)
    return fx


def cosa(a, b, c, d): # x1, y1, x2, y2
    fx = (a - c)/math.sqrt((a - c)**2 + (b - d)**2)
    return fx


def ysecd(py1, py2, py3, vy1, spring1, spring2, part, force1, force2, sin1, sin2, t):
    fx = (- spring1.damping_coefficient * vy1\
            - spring1.spring_coefficient * ((abs(py1- py2) - spring1.original_distance_y))\
            + ft(t, force1.f,force1.c1, force1.c2) * sin1 \
            - spring2.damping_coefficient * vy1 \
            - spring2.spring_coefficient * ((abs(py1 - py3) - spring2.original_distance_y)) \
            + ft(t, force2.f, force2.c1, force2.c2) * sin2 \
            + part.m * pa.g
            ) / part.m
    return fx


def xsecd(px1, px2, px3, vx1, spring1, spring2, part, force1, force2, cos1, cos2, t):
    fx = (- spring1.damping_coefficient * vx1\
            - spring1.spring_coefficient * ((abs(px1- px2) - spring1.original_distance_x))\
            + ft(t, force1.f,force1.c1, force1.c2) * cos1 \
            - spring2.damping_coefficient * vx1 \
            - spring2.spring_coefficient * ((abs(px1 - px3) - spring2.original_distance_x)) \
            + ft(t, force2.f, force2.c1, force2.c2) * cos2
            ) / part.m
    return fx


# def calculator(spring_1, spring_2, spring_3, particle1, particle2, particle3, ad_fo1, ad_fo2): # S12 S23 S13

    # def ft(f, c1, c2):
    #     fx = c1 * f(t) + c2
    #     return fx
    #
    #
    # def sina(a, b, c, d): # x1, y1, x2, y2
    #     fx = (b - d)/math.sqrt((a - c)**2 + (b - d)**2)
    #     return fx
    #
    #
    # def cosa(a, b, c, d): # x1, y1, x2, y2
    #     fx = (a - c)/math.sqrt((a - c)**2 + (b - d)**2)
    #     return fx
    #
    #
    # def ysecd(py1, py2, py3, vy1, spring1, spring2, part, force1, force2, sin1, sin2, t):
    #     fx = (- spring1.damping_coefficient * vy1\
    #             - spring1.spring_coefficient * ((abs(py1- py2) - spring_1.original_distance_y))\
    #             + ft(force1.f,force1.c1, force1.c2) * sin1 \
    #             - spring2.damping_coefficient * vy1 \
    #             - spring2.spring_coefficient * ((abs(py1 - py3) - spring2.original_distance_y)) \
    #             + ft(force2.f, force2.c1, force2.c2) * sin2 \
    #             + part.m * pa.g
    #             ) / part.m
    #     return fx
    #
    #
    # def xsecd(px1, px2, px3, vx1, spring1, spring2, part, force1, force2, cos1, cos2, t):
    #     fx = (- spring1.damping_coefficient * vx1\
    #             - spring1.spring_coefficient * ((abs(px1- px2) - spring_1.original_distance_x))\
    #             + ft(force1.f,force1.c1, force1.c2) * cos1 \
    #             - spring2.damping_coefficient * vx1 \
    #             - spring2.spring_coefficient * ((abs(px1 - px3) - spring2.original_distance_x)) \
    #             + ft(force2.f, force2.c1, force2.c2) * cos2
    #             ) / part.m
    #     return fx


class cal():
    def __init__(self, spring_1, spring_2, spring_3, particle1, particle2, particle3, ad_fo1, ad_fo2):
        self.spring_1 = spring_1
        self.spring_2 = spring_2
        self.spring_3 = spring_3
        self.particle1 = particle1
        self.particle2 =  particle2
        self.particle3 = particle3
        self.ad_fo1 = ad_fo1
        self.ad_fo2 = ad_fo2


def calculator(spring_1, spring_2, spring_3, particle1, particle2, particle3, ad_fo1, ad_fo2):
        def func(t, r):
            px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3 = r  # position and volocity

    #############################  x1 y1
            sin12 = sina(px1, py1, px2, py2)
            sin13 = sina(px1, py1, px3, py3)
            cos12 = cosa(px1, py1, px2, py2)
            cos13 = cosa(px1, py1, px3, py3)

            if py1 < 0 and vy1 < 0:  # bounce
                vy1 = - vy1 * pa.cel
                py1 = - py1

            fy1f = vy1  # first order P1y
            fy1s = ysecd(py1, py2, py3, vy1, spring_1, spring_3, particle1, ad_fo1, ad_fo2, sin12, sin13, t)

            fx1f = vx1  # first order P1x
            fx1s = xsecd(px1, px2, px3, vx1, spring_1, spring_3, particle1, ad_fo1, ad_fo2, cos12, cos13, t)

            if py1 < 0.1 and py1 >= 0: # frinction
                if fx1s < 0 and abs(fx1s) > abs(fy1s * pa.f * particle1.m):
                    fx1s = fx1s + abs(fy1s * pa.f)
                if fx1s > 0 and abs(fx1s) > abs(fy1s * pa.f * particle1.m):
                    fx1s = fx1s - abs(fy1s * pa.f)
                else:
                    fx1s = 0


    ########################################## x2 y2
            sin21 = sina(px2, py2, px1, py1)
            sin23 = sina(px2, py2, px3, py3)
            cos21 = cosa(px2, py2, px3, py3)
            cos23 = cosa(px2, py2, px3, py3)

            if py2 < 0 and vy2 < 0:  # bounce
                vy1 = - vy2 * pa.cel
            py2 = - py2

            fy2f = vy2  # first order P1y
            fy2s = ysecd(py1, py2, py3, vy1, spring_1, spring_2, particle2, ad_fo1, ad_fo2, sin21, sin23)

            fx2f = vx2  # first order P1x
            fx2s = xsecd(px1, px2, px3, vx1, spring_1, spring_2, particle2, ad_fo1, ad_fo2, cos21, cos23)

            if py2 < 0.1 and py2 >= 0:  # frinction
                if fx2s < 0 and abs(fx2s) > abs(fy2s * pa.f * particle2.m):
                    fx2s = fx2s + abs(fy2s  * pa.f)
                if fx2s > 0 and abs(fx2s) > abs(fy2s * pa.f * particle2.m):
                    fx2s = fx2s - abs(fy2s * pa.f)
                else:
                    fx2s = 0
     #######################################  x3 y3
            sin31 = sina(px3, py3, px1, py1)
            sin32 = sina(px3, py3, px1, py1)
            cos31 = cosa(px3, py3, px1, py1)
            cos32 = cosa(px3, py3, px1, py1)

            if py3 < 0 and vy3 < 0:  # bounce
                vy3 = - vy3 * pa.cel
                py3 = - py3

            fy3f = vy3 # first order P1y
            fy3s = ysecd(py1, py2, py3, vy1, spring_2, spring_3, particle3, ad_fo2, ad_fo2, sin31, sin32)

            fx3f = vx3  # first order P1x
            fx3s = xsecd(px1, px2, px3, vx1, spring_2, spring_3, particle3, ad_fo2, ad_fo2, cos31, cos32)

            if py3 < 0.1 and py3 >= 0:  # frinction
                if fx3s < 0 and abs(fx3s) > abs(fy3s * pa.f * particle3.m):
                    fx3s = fx3s + abs(fy3s * pa.f)
                if fx3s > 0 and abs(fx3s) > abs(fy3s * pa.f * particle3.m):
                    fx3s = fx3s - abs(fy3s * pa.f)
                else:
                    fx3s = 0

            return fx1f, fx1s, fy1f, fy1s, fx2f, fx2s, fy2f, fy2s, fx3f, fx3s, fy3f, fy3s
