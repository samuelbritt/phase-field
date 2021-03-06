"""
Reformulation of the 1-D equations in Kobayashi, but this time an using
an implicit finite difference method
"""

import kobayashi as kb

from scipy import *
import scipy.linalg as linalg
import numpy as np
import matplotlib.pyplot as plt

def boundary(t):
    """ values at the boundaries, x = 0 and x = n-1 """
    return 1,0

def fprime(phi, a, Df):
    return 2 * a * a * phi * (1 - phi) * (phi - .5 + 15 * Df / (a * a) * phi *
                                          (1 - phi))

def phi_np1(phi_n, x, dx, t, dt, M, eps, a, Df, verbose=False):
    """ returns phi(x, t=n+1) given phi(x, t=n) using finite difference method.

    Uses *backward* finite difference in time (implicit), central finite
    difference in x
    """

    mu = M * dt / (dx * dx)
    eps2 = eps * eps

    # Set up matrix from 1 to n-2, and handle x = 0 and x = n-1 via boundary
    # conditions
    alpha, beta = boundary(t)
    boundary_vec = zeros_like(phi_n[1:-1])
    boundary_vec[0] = - mu * eps2 * alpha
    boundary_vec[-1] = - mu * eps2 * beta

    fp = fprime(phi_n[1:-1], a, Df)

    b = phi_n[1:-1] - boundary_vec + M * dt * fp

    d = ones_like(phi_n[1:-1]) * (1 + 2 * mu * eps2)
    u = ones_like(phi_n[1:-1]) * (- mu * eps2)
    u[0] = 0
    l = ones_like(phi_n[1:-1]) * (- mu * eps2)
    l[-1] = 0

    ab = np.matrix([u, d, l])

    res = zeros_like(phi_n)
    res[1:-1] = linalg.solve_banded((1, 1), ab, b)
    res[0] = alpha
    res[-1] = beta
    return res

if __name__ == '__main__':
    # require | Df | < a^2 / 6
    # require eps < a

    x0 = -5.
    x_max = 150.
    dx = 0.06
    x_all = arange(x0, x_max, dx)

    t0 = 0.
    t_max = 40
    dt = 0.02
    t_all = arange(t0, t_max, dt)
    plot_count = 10
    plot_every = int(t_max / dt / plot_count) 

    a = 10.
    eps = 2.
    Df = 3.
    M = .5

    phi0 = 0.5 * (1 - tanh(3 * x_all))
    phi0 = 0.5 * (1 - tanh(a * x_all / (2 * eps)))

    f = plt.figure()
    ax = f.add_subplot(111)

    # phi_old = zeros_like(x_all)
    phi_old = phi0
    t = t0
    i = 0
    while (t < t_max):
        if i % plot_every == 0:
            ax.plot(x_all, phi_old, label="$t={}$".format(t))
        i += 1
        t = t + dt
        phi_new = phi_np1(phi_old, x_all, dx, t, dt, M, eps, a, Df)
        phi_old = phi_new

    ax.legend(loc='best')
    plt.show()
