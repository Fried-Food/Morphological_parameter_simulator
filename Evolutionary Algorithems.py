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
from SpringSimulator.Simulator import ad_force
from SpringSimulator.Simulator import initpos


f = 0.8  # friction coefficient
g = -9.81  # gravity
dt = 0.05  # differential time
stop_time = 30  # tot running time
time = 0  # time mark
t = np.arange(0, stop_time, dt)


DNA_SIZE = 10            # DNA length
POP_SIZE = 100           # population size
CROSS_RATE = 0.8         # mating probability (DNA crossover)
MUTATION_RATE = 0.003    # mutation probability
N_GENERATIONS = 10      # generation number
X_BOUND = [0, 5]         # k,b upper and lower bounds



x1y1_mv = []
av_x = []
centre_of_mass = []


def F(sc_v, dc_v):
    POS = initpos(0, 50, 40, -math.pi / 3)  # x, y, edge_length, rot angle(rad)
    P1 = Particle(POS[0], POS[1], 1)  # X, Y, mass
    P2 = Particle(POS[2], POS[3], 1)
    P3 = Particle(POS[4], POS[5], 1)

    S12 = Spring(sc_v, dc_v, P1, P2)  # spring_coefficient, damping_coefficient, particle_1, particle_2
    S23 = Spring(sc_v, dc_v, P2, P3)
    S13 = Spring(sc_v, dc_v, P1, P3)

    for i in t:
        if i == 0:
            x1y1_mv.clear()
        # S12.input_force = ad_force(i, math.sin, 50, 2, 0)
        coordinate(S12, S13, P1)
        coordinate(S12, S23, P2)
        coordinate(S13, S23, P3)
        B = 1 / 3 * (P1.x + P2.x + P3.x)
        centre_of_mass.append(B)

    result = (centre_of_mass[-1] - centre_of_mass[0]) / stop_time  # average speed in x-axis

    return result    # to find the maximum of this function


# find non-zero fitness for selection
def get_fitness(pred):
    return pred + 1e-3 - np.min(pred)


# convert binary DNA to decimal and normalize it to a range(0, 5)
def translateDNA(pop):
    return pop.dot(2 ** np.arange(DNA_SIZE)[::-1]) / float(2**DNA_SIZE-1) * X_BOUND[1]
    # return int(pop, 2)*0.1


def select(pop, fitness):    # nature selection wrt pop's fitness
    idx = np.random.choice(np.arange(POP_SIZE), size=POP_SIZE, replace=True,
                           p=fitness/fitness.sum())
    return pop[idx]


def crossover(parent, pop):     # mating process (genes crossover)
    if np.random.rand() < CROSS_RATE:
        i_ = np.random.randint(0, POP_SIZE, size=1)                             # select another individual from pop
        cross_points = np.random.randint(0, 2, size=2*DNA_SIZE).astype(np.bool)   # choose crossover points
        parent[cross_points] = pop[i_, cross_points]                            # mating and produce one child
    return parent


def mutate(child):
    for point in range(DNA_SIZE):
        if np.random.rand() < MUTATION_RATE:
            child[point] = 1 if child[point] == 0 else 0
    return child


pop = np.random.randint(2, size=(POP_SIZE, 2*DNA_SIZE))   # initialize the pop DNA
pop_k = pop[:, :DNA_SIZE]
pop_b = pop[:, DNA_SIZE:]



HAHA = []
HAHA = np.array(HAHA)
# TL = np.zeros((1, 1))
# FA = np.zeros((1, 1))
k_drwaing = []
b_drawing = []
fit_drawing = []


for _ in range(N_GENERATIONS):
    # F_values = F(translateDNA(pop_k)[i], translateDNA(pop_b)[i])    # compute function value by extracting DNA
    for i in range(POP_SIZE):
        F_values = F(translateDNA(pop_k)[i], translateDNA(pop_b)[i])
        if i == 0:
            HAHA = []
            HAHA = np.array(HAHA)
        HAHA = np.append(HAHA, F_values)
        # print(HAHA)
        # HAHA = np.append(F_values, i)
    print("GENERATIONS: %s" %_, 'Progress:{:.2%}'.format(_/N_GENERATIONS))

    # GA part (evolution)
    fitness = get_fitness(HAHA)
    print("Most fitted DNA: ", pop[np.argmax(fitness), :])
    print("K: ", translateDNA(pop[np.argmax(fitness), :][:10]))
    print("B: ", translateDNA(pop[np.argmax(fitness), :][10:]))
    pop = select(pop, fitness)
    pop_copy = pop.copy()
    for parent in pop:
        child = crossover(parent, pop_copy)
        child = mutate(child)
        parent[:] = child       # parent is replaced by its child

    k_drwaing.append(translateDNA(pop[np.argmax(fitness), :][:10]))
    b_drawing.append(translateDNA(pop[np.argmax(fitness), :][10:]))
    fit_drawing.append(np.argmax(fitness))


# saving
np.save('k_drwaing.npy',k_drwaing)
np.save('b_drawing.npy',b_drawing)
np.save('fit_drawing.npy',fit_drawing)
# loading
list1=np.load('k_drwaing.npy')
list2=np.load('b_drawing.npy')
list3=np.load('fit_drawing.npy')

############  DRAWING  ###############
fig, ax = plt.subplots(1, 1)
# set same x-axis
ax_sub = ax.twinx()
# x-axis values
x = np.arange(0,len(list1))+1
x[0] = 1
l1, = ax.plot(x,list1,label='k',marker='o')
l2, = ax.plot(x,list2,label='b',marker='o')
l3, = ax_sub.plot(x,list3,label='fitness',marker='o',color='g')
# add legend
plt.legend(handles=[l1, l2, l3], labels=['k', 'b','fitness'], loc=0)
# set y-axis lable
ax.set_ylabel('Values of k and b')
ax_sub.set_ylabel('Fitness')
# set x-axis title
ax.set_xlabel('Generation')
# set figure title
ax.set_title('Evolutionary Algorithms')
plt.show()