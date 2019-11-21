import ../KuramotoOscillator as ko
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from sympy import *
import numpy as np
import random



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
    # interlayer coupling strength
    lm = 1
    A1 = Matrix(
        [[0, 1, 1, 1, 0, 0, 1],
         [1, 0, 1, 0, 0, 0, 0],
         [1, 1, 0, 0, 0, 0, 0],
         [1, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 1, 0],
         [0, 0, 0, 0, 1, 0, 1],
         [1, 0, 0, 0, 0, 1, 0]])

    A2 = Matrix(
        [[0, 1, 1, 1, 0, 0, 1],
         [1, 0, 1, 0, 0, 0, 0],
         [1, 1, 0, 0, 0, 0, 0],
         [1, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [1, 0, 0, 0, 0, 1, 0]])


    # A1 = Matrix(
    #     [[0, 1, 1],
    #      [1, 0, 1],
    #      [1, 1, 0]])
    #
    # A2 = Matrix(
    #     [[0, 1, 1],
    #      [1, 0, 0],
    #      [1, 0, 0]])


    # create a dictionary to hold the parameters of the network

    Z0 = []

    init_params1 = {'om': om, 'A': A1, 'X0': random_initial_conditions(7), 'al': al}

    my_kuramoto1 = ko.Kuramoto(init_params1)

    # my_kuramoto1.reset_initial_conditions()
    #
    # t1, y1 = my_kuramoto1.solve_ode(stop_time=500, num_points=500)
    #
    # for i in y1[-1]:
    #     Z0.append(i)
    #

    init_params2 = {'om': om, 'A': A2, 'X0': random_initial_conditions(7), 'al': al}

    my_kuramoto2 = ko.Kuramoto(init_params2)

    # my_kuramoto2.reset_initial_conditions()

    # t2, y2 = my_kuramoto2.solve_ode(stop_time=500, num_points=500)
    #
    # for j in y2[-1]:
    #     Z0.append(j)
    #
    # print(Z0)


    multiplex_params = {'om': om, 'al': al, 'X0': random_initial_conditions(14), 'lm': lm, 'sig':1}

    my_network = ko.Multiplex([my_kuramoto1, my_kuramoto2], multiplex_params)

    t, y = my_network.solve_multiplex_ode(num_points=1000, stop_time=1000)

    plt.plot(t, np.sin(y[:, 0]), color='aqua')
    plt.plot(t, np.sin(y[:, 1]), color='orange')
    plt.plot(t, np.sin(y[:, 2]), '--', color='red')
    plt.plot(t, np.sin(y[:, 3]), color='blue')
    plt.plot(t, np.sin(y[:, 4]), color='coral')
    plt.plot(t, np.sin(y[:, 5]), '--', color='crimson')
    plt.plot(t, np.sin(y[:, 6]), '--', color='darkgreen')
    plt.plot(t, np.sin(y[:, 7]), color='gold')
    plt.plot(t, np.sin(y[:, 8]), color='indigo')
    plt.plot(t, np.sin(y[:, 9]), '--', color='lightblue')
    plt.plot(t, np.sin(y[:, 10]), color='navy')
    plt.plot(t, np.sin(y[:, 11]), color='olive')
    plt.plot(t, np.sin(y[:, 12]), '--', color='salmon')
    plt.plot(t, np.sin(y[:, 13]), '--', color='tomato')

    # plt.plot(t, (y[:, 0]), color='aqua')
    # plt.plot(t, (y[:, 1]), color='orange')
    # plt.plot(t, (y[:, 2]), '--', color='red')
    # plt.plot(t, (y[:, 3]), color='blue')
    # plt.plot(t, (y[:, 4]), color='coral')
    # plt.plot(t, (y[:, 5]), '--', color='crimson')
    # plt.plot(t, (y[:, 6]), '--', color='darkgreen')
    # plt.plot(t, (y[:, 7]), color='gold')
    # plt.plot(t, (y[:, 8]), color='indigo')
    # plt.plot(t, (y[:, 9]), '--', color='lightblue')
    # plt.plot(t, (y[:, 10]), color='navy')
    # plt.plot(t, (y[:, 11]), color='olive')
    # plt.plot(t, (y[:, 12]), '--', color='salmon')
    # plt.plot(t, (y[:, 13]), '--', color='tomato')


    plt.text(t[-1], (y[:, 0])[-1]+100, r' $\phi_{1, 1}$')
    plt.text(t[-1], (y[:, 0])[-1]+50, r' $\phi_{1, 2}$')
    plt.text(t[-1] + 45, (y[:, 0])[-1]+50, r' $\phi_{1, 3}$')
    plt.text(t[-1], (y[:, 3])[-1]-50, r' $\phi_{1, 4}$')
    plt.text(t[-1], (y[:, 4])[-1]+50, r' $\phi_{1, 5}$')
    plt.text(t[-1] + 45, (y[:, 4])[-1]+50, r' $\phi_{1, 6}$')
    plt.text(t[-1] + 45, (y[:, 3])[-1]-50, r' $\phi_{1, 7}$')

    plt.text(t[-1], (y[:, 7])[-1]-50, r' $\phi_{2, 1}$')
    plt.text(t[-1], (y[:, 7])[-1]-95, r' $\phi_{2, 2}$')
    plt.text(t[-1] + 45, (y[:, 7])[-1]-95, r' $\phi_{2, 3}$')
    plt.text(t[-1], (y[:, 10])[-1]-50, r' $\phi_{2, 4}$')
    plt.text(t[-1], (y[:, 11])[-1]+50, r' $\phi_{2, 5}$')
    plt.text(t[-1] + 45, (y[:, 11])[-1]+50, r' $\phi_{2, 6}$')
    plt.text(t[-1] + 45, (y[:, 10])[-1]-50, r' $\phi_{2, 7}$')

    plt.xlim(0, 1116)

    plt.legend((r' $\phi_1$', r'$\phi_2$', r'$\phi_{3}$', r'$\phi_{4}$', r'$\phi_{5}$',
                 r'$\phi_{6}$', r'$\phi_{7}$', r'$\phi_{8}$', r'$\phi_{9}$',
                 r'$\phi_{10}$', r'$\phi_{11}$', r'$\phi_{12}$', r'$\phi_{13}$',
                 r'$\phi_{14}$'
                 ), prop=FontProperties(size=10))
    plt.title(
        r'Oscillator phases for multiplex network $G_{1,\gamma}$ with phase frustration $\alpha$= %f and $\lambda$ = %f' % (al, lm))
    plt.xlabel('time, $t$')
    plt.ylabel('phase $sin(\phi_i)$ of the $i_{th}$ oscillator')

    plt.show()







