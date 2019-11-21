import numpy as np
from scipy.integrate import odeint
from sympy import *
import cmath
import random

class Kuramoto:

    def __init__(self, p):
        """
                Initialises the n-oscillator system as a single object

                Arguments:
                    p :  dictionary containing the intrinsic frequencies of the oscillators om (degenerate), the initial phases init_phase,
                        the adjacency matrix A, and a phase frustration parameter al (=alpha). The number of oscillators in the network
                        and the connectivity of each oscillator is stored in variables n_osc and m_osc respectively.

                """

        self.init_phase = np.array(p['X0'])
        self.om = np.array(p['om'])
        self.A = p['A']
        self.al = p['al']
        self.n_osc = self.A.shape[1]
        self.m_order = self.A.shape[0]




    def define_ode(self, w, t):
        """
        Defines the differential equations for the coupled spring-mass system.

        Arguments:
            w :  vector of the state variables:
                      w = [theta_1, theta_2, theta_3, theta_4, theta_5, theta_6, theta_7]
            t :  time
            p :  vector of the parameters:
                      p = [al, A]

            A is the adjacency matrix
        """

        f = []  # create an array to hold n differential equations for variables X_1-X_7

        # Create f = (x1',x2',..., x7'):

        #  loop through the adjacency array to build f.
        for i in range(0, self.n_osc):

            temp_var = []  # array to hold the contributions to the differential equation for the nth oscillator.

            for j in range(0, self.m_order):
                temp_var.append(self.A[i,j] * np.sin(w[i] - w[j] + self.al))

            # use this to build elements of f
            f.append(self.om - np.sum(temp_var))

        return f

    # solves a set of differential equations using the scipy odeint method.
    def solve_ode(self, stop_time = 1000, num_points = 10000000):
        """
        Takes the set of differential equations encoded in f and solves them.
        This uses the odeint method in the scipy package.

        """
        t = [stop_time * float(i) / (num_points - 1) for i in range(num_points)]

        w_sol = odeint(self.define_ode, self.init_phase, t)

        return t, w_sol

    def reset_initial_conditions(self):
        X0 = []
        num = self.n_osc

        j=0

        for i in range(num):
            if j == 2:
                X0.append(X0[1])
            elif j == 6:
                X0.append(X0[3])
            elif j == 5:
                X0.append(X0[4])
            else:
                X0.append(random.uniform(0, np.pi / 2))

            j += 1

        self.init_phase = X0

    def change_link_strength(self, link, strength):
        J = self.A
        J[link[0]][link[1]]= strength
        J[link[1]][link[0]] = strength
        self. A =J

class Multiplex(Kuramoto):

    def __init__(self, kur_obj, p):

        self.layers=kur_obj
        self.layer_number = len(kur_obj)
        self.n_osc = kur_obj[0].n_osc

        self.om = (p['om'])
        self.al = p['al']
        self.init_phase = np.array(p['X0'])
        self.lm = p['lm']
        self.sigma = p['sig']

        self.subnetworks = kur_obj

        K_1 = self.sigma * (self.subnetworks[0].A).row_join(self.lm * eye(self.subnetworks[0].n_osc))
        K_2 = self.sigma * (self.lm * eye(self.subnetworks[1].n_osc)).row_join(self.subnetworks[1].A)
        self.M = K_1.col_join(K_2)


    def define_multiplex_ode(self, w, t):
        """
                Defines the differential equations for the coupled spring-mass system.

                Arguments:
                    w :  vector of the state variables:
                              w = [theta_11, theta_12,..., theta_17, theta_21, theta_22, ..., theta_27]
                    t :  time
                    p :  vector of the parameters:
                              p = [al1, al2, A1, A2]

                    A1 and A2 are the adjacency matrices of the networks in layer 1 and layer 2 respectively.
                """

        f = []  # create an array to hold n differential equations for variables X_11-X_27

        #  loop through the adjacency array to build f.
        for i in range(0, self.n_osc * self.layer_number):

            temp_var = []  # array to hold the contributions to the differential equation for the nth oscillator.

            for j in range(0, self.n_osc * self.layer_number):
                temp_var.append(self.M[i, j] * np.sin(w[i] - w[j] + self.al))

            # use this to build elements of f
            f.append(self.om - np.sum(temp_var))

        return f

    def update_m(self, lm, sig):
        self.lm = lm
        self.sigma = sig

        K_1 = self.sigma * (self.subnetworks[0].A).row_join(self.lm * eye(self.subnetworks[0].n_osc))
        K_2 = self.sigma * (self.lm * eye(self.subnetworks[1].n_osc)).row_join(self.subnetworks[1].A)
        self.M = K_1.col_join(K_2)


    def update_sigma(self, s):
        self.sigma = s

        self.update_m()

    def update_lambda(self, lam):
        self.lam = lam

        self.update_m()

    def reset_initial_conditions(self):
        X0 = []
        num = self.n_osc * self.layer_number
        j=0

        j = 0

        for i in range(num):
            if i==2:
                X0.append(X0[1])
            elif i==6:
                X0.append(X0[3])
            elif i==5:
                X0.append(X0[4])
            elif i==7:
                X0.append(X0[0])
            elif i==9:
                X0.append(X0[8])
            elif i==13:
                X0.append(X0[10])
            elif i==12:
                X0.append(X0[11])
            else:
                X0.append(random.uniform(0, np.pi/2))

            j+=1

        self.init_phase=X0

    def solve_multiplex_ode(self, stop_time = 500, num_points = 10000):
        """
            Takes the set of differential equations encoded in f and solves them.
            This uses the odeint method in the scipy package.

        """
        t = [stop_time * float(i) / (num_points - 1) for i in range(num_points)]

        w_sol = odeint(self.define_multiplex_ode, self.init_phase, t)

        return t, w_sol


