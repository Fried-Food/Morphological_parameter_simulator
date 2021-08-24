import numpy as np
from global_var import parameters as pa
from simulator import position
from simulator.position import Particle
from simulator.position import Spring
from simulator.position import AddForce
from simulator.position import cal
from simulator.position import ft
from simulator.position import sina
from simulator.position import cosa
from simulator.position import ysecd
from simulator.position import xsecd
from scipy.integrate import odeint, solve_bvp, solve_ivp
from numpy import sin, cos, linspace, pi
import numpy as np
import matplotlib.pyplot as plt


pa.f = 0.5  # friction coefficient
pa.g = -9.81  # gravity
pa.dt = 0.05  # differential time
pa.stop_time = 30  # tot running time
pa.time = 0  # time mark
pa.cel = 0.9 # Coefficient of collision energy loss



P1 = Particle(0, 20, 1)  # X, Y, mass
P2 = Particle(20, 50, 1)
P3 = Particle(40, 20, 1)

S12 = Spring(6, 0.3, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
S23 = Spring(6, 0.3, P2, P3)
S13 = Spring(6, 0.3, P1, P3)

active = AddForce(sin, 1, 0)
passive = AddForce(sin, 0, 0)

cal(S12, S23, S13, P1, P2, P3, active, passive)

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

t1 = linspace(0, 0.01, 10)
y0 = [px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3]

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
    fy1s = ysecd(py1, py2, py3, vy1, S12, S13, P1, active, passive, sin12, sin13, t)

    fx1f = vx1  # first order P1x
    fx1s = xsecd(px1, px2, px3, vx1, S12, S13, P1, active, passive, cos12, cos13, t)

    if py1 < 0.1 and py1 >= 0: # frinction
        if fx1s < 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
            fx1s = fx1s + abs(fy1s * pa.f)
        if fx1s > 0 and abs(fx1s) > abs(fy1s * pa.f * P1.m):
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
    fy2s = ysecd(py1, py2, py3, vy1, S12, S23, P2, active, passive, sin21, sin23, t)

    fx2f = vx2  # first order P1x
    fx2s = xsecd(px1, px2, px3, vx1, S12, S23, P2, active, passive, cos21, cos23, t)

    if py2 < 0.1 and py2 >= 0:  # frinction
        if fx2s < 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
            fx2s = fx2s + abs(fy2s  * pa.f)
        if fx2s > 0 and abs(fx2s) > abs(fy2s * pa.f * P2.m):
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
    fy3s = ysecd(py1, py2, py3, vy1, S23, S13, P3, active, passive, sin31, sin32, t)

    fx3f = vx3  # first order P1x
    fx3s = xsecd(px1, px2, px3, vx1, S23, S13, P3, active, passive, cos31, cos32, t)

    if py3 < 0.1 and py3 >= 0:  # frinction
        if fx3s < 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
            fx3s = fx3s + abs(fy3s * pa.f)
        if fx3s > 0 and abs(fx3s) > abs(fy3s * pa.f * P3.m):
            fx3s = fx3s - abs(fy3s * pa.f)
        else:
            fx3s = 0

    return fx1f, fx1s, fy1f, fy1s, fx2f, fx2s, fy2f, fy2s, fx3f, fx3s, fy3f, fy3s



from matplotlib import animation

sol = solve_ivp(func, (0, 1), y0, method='RK45', t_eval=t1)
px1, vx1, py1, vy1, px2, vx2, py2, vy2, px3, vx3, py3, vy3 = sol.y
#############################################################################

plt.figure()
fig, ax = plt.subplots()
# passive spring
line, = ax.plot([], [], 'o-', lw=2)
# active spring
line1, = ax.plot([], [], 'ro-', lw=2)

time_template = 'Distance = %.2fcm'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    ax.set_xlim(-100, 100)
    ax.set_ylim(0, 100)
    return line,


def update(i):

    thisx = [px2[i], px3[i], px1[i]]
    thisy = [py2[i], py3[i], py1[i]]

    thisx1 = [px1[i], px2[i]]
    thisy1 = [py1[i], py2[i]]

    line1.set_data(thisx1, thisy1)
    line.set_data(thisx, thisy)

    # if x1[i] > x_last:
    #     x_max = x1[i]
    #     x_last = x1[i]
    # else:

    time_text.set_text(time_template % (px1[i]))

    return line, line1, time_text


ani = animation.FuncAnimation(
    fig =fig,
    func=update,
    init_func=init,
    #frames = np.linspace(0, stop_time / dt, stop_time / dt),
    frames=len(px1),
    interval = 20,
    repeat = False,
    blit = False
)


plt.show()