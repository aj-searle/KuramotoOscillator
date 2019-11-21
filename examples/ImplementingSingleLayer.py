import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import KuramotoOscillator as ko
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sympy import *
import numpy as np
import random
import pickle
import cmath
import scipy.signal as sg


def random_initial_conditions(num, stop=np.pi/2):
    X0 = []

    for i in range(num):
        X0.append(random.uniform(0, stop))

    return X0


if __name__ == '__main__' :
    # Define the parameters of the network

    # oscillator frequency
    om = 1
    # phase frustration parameter
    al = 0.3

    alp_array = []
    # Adjacency matrix
    A1 = Matrix(
        [[0, 0, 1, 1, 0, 0, 1],
         [0, 0, 1, 0, 0, 0, 0],
         [1, 1, 0, 0, 0, 0, 0],
         [1, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [1, 0, 0, 0, 0, 1, 0]])

    init_params1 = {'om': om, 'A': A1, 'X0': random_initial_conditions(4), 'al': al}

    my_kuramoto1 = ko.Kuramoto(init_params1)

    my_kuramoto1.reset_initial_conditions()

    # link_array = [1 * float(i) / (100 - 1) for i in range(100)]
    # position_array = [[0,1], [0,2], [0,3], [0,6], [1,2], [3,4], [4,5], [5,6]]

    alp_array = [np.pi / 2 * float(i) / (50 - 1) for i in range(50)]

    t, y = my_kuramoto1.solve_ode(stop_time=300, num_points=10000)

    plt.plot(t, np.sin(y[:, 0]), color='aqua')
    plt.plot(t, np.sin(y[:, 1]), color='orange')
    plt.plot(t, np.sin(y[:, 2]), '--', color='red')
    plt.plot(t, np.sin(y[:, 3]), color='blue')
    plt.plot(t, np.sin(y[:, 4]), color='coral')
    plt.plot(t, np.sin(y[:, 5]), '--', color='crimson')
    plt.plot(t, np.sin(y[:, 6]), '--', color='darkgreen')



    plt.legend((r' $\phi_{1}$', r'$\phi_{2}$', r'$\phi_{3}$', r'$\phi_{4}$', r'$\phi_{5}$',
                  r'$\phi_{6}$', r'$\phi_{7}$'), prop=FontProperties(size=10))

    # plt.text(t[-1], (y[:, 0])[-1], r' $\phi_{2,1}$')
    # plt.text(t[-1], (y[:, 1])[-1], r' $\phi_{2,2}$')
    # plt.text(t[-1]+12.5, (y[:, 1])[-1], r' $\phi_{2,3}$')
    # plt.text(t[-1], (y[:, 3])[-1], r' $\phi_{2,4}$')
    # plt.text(t[-1], (y[:, 4])[-1], r' $\phi_{2,5}$')
    # plt.text(t[-1]+12.5, (y[:, 4])[-1], r' $\phi_{2,6}$')
    # plt.text(t[-1]+12.5, (y[:, 3])[-1], r' $\phi_{2,7}$')

    plt.xlim(0, 340)


    # plt.title(
    #     r'$\alpha=%f$' % (al))
    plt.xlabel('time, $t$')
    plt.ylabel('displacement $sin(\phi_{i})$ of $i_{th}$ oscillator')

    plt.show()


    # loop through alpha (phase frustration) and sigma (interlayer coupling) and put the solved ode time series into file



