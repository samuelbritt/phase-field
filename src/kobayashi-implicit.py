"""
Reformulation of the 1-D equations in Kobayashi, but this time an using
an implicit finite difference method
"""

import kobayashi as kb

from scipy import *
import scipy.linalg as linalg
import numpy as np
import matplotlib.pyplot as plt

def phi_np1(phi_n, x, dx, t, dt, M, eps, Df, a, verbose=False):
    """ returns phi(x, t=n+1) given phi(x, t=n) using finite difference method.

    Uses *backward* finite difference in time (implicit), central finite
    difference in x
    """
    mu = M * dt / (dx * dx)
    eps2 = eps * eps
    b = phi_n - M * dt * kb.fprime(phi_n, Df, a)

    d = ones_like(phi_n) * (1 + 2 * mu * eps2)
    u = ones_like(phi_n) * (- mu * eps2)
    u[0] = 0
    l = ones_like(phi_n) * (- mu * eps2)
    l[-1] = 0

    ab = np.matrix([u, d, l])

    return linalg.solve_banded((1, 1), ab, b)

if __name__ == '__main__':
    # require | Df | < a^2 / 6
    # require eps < a

    x0 = -10.
    x_max = 50.
    dx = 0.06
    x_all = arange(x0, x_max, dx)

    t0 = 0.
    t_max = 20
    dt = 0.01
    t_all = arange(t0, t_max, dt)
    plot_count = 10
    plot_every = int(t_max / dt / plot_count) 

    eps = 1.
    a = 10. * eps
    Df = -10.
    tau = 2.

    phi0 = 0.5 * (1 - tanh(a * x_all / (2 * eps)))

    f = plt.figure()
    ax = f.add_subplot(111)

    # phi_old = zeros_like(x_all)
    phi_old = phi0
    M = 1. / tau
    t = t0
    i = 0
    while (t < t_max):
        if i % plot_every == 0:
            ax.plot(x_all, phi_old, label="$t={}$".format(t))
        i += 1
        t = t + dt
        phi_new = phi_np1(phi_old, x_all, dx, t, dt, M, eps, Df, a) 
        phi_old = phi_new

    ax.legend(loc='best')
    plt.show()
