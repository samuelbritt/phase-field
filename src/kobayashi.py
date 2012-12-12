#!/usr/bin/env python
# encoding: utf-8

r""" 1D solution of the Cahn Hilliard eqn asusming constant driving force
    \Delta f:

        \tau \partial\phi / \partial t =
            \epsilon^2 \grad^2\phi
            + 2 a^2 \phi (1-\phi) * (\phi - 1/2 + 3 \Delta f / a^2)
"""

from scipy import *
import matplotlib.pyplot as plt

def h1(phi):
    return phi * phi * phi * (10 - 15 * phi + 6 * phi * phi)

def h2(phi):
    return phi * phi * (3 - 2 * phi)

def g(phi):
    return phi * phi * (1 - phi) * (1 - phi)

def f(phi, Df, a, h=h1):
    return a * a / 2. * g(phi) + Df * h(phi)

def fprime(phi, Df, a):
    res = 2 * a * a * phi * (1-phi) * (.5 - phi + 3 * Df / (a * a))
    # if res != 0:
    #     print ("fprime({}, {}, {}) = {}".format(phi, Df, a, res))
    return res

def central_diff(phi, x, dx):
    d2phi_dx2 = zeros_like(phi)
    for i in range(1, len(phi) - 1):
        d2phi_dx2[i] = (phi[i+1] - 2 * phi[i] + phi[i-1])
    d2phi_dx2[0] = d2phi_dx2[1]
    d2phi_dx2[-1] = d2phi_dx2[-2]
    return d2phi_dx2 / (dx * dx)

def phi_np1(phi_n, x, dx, t, dt, M, eps, Df, a, verbose=False):
    """ returns phi(x, t=n+1) given phi(x, t=n) using finite difference method.

    Uses forward finite difference in time, central finite difference in x
    """
    c = central_diff(phi_n, x, dx)
    fp = fprime(phi_n, Df, a)
    res = zeros_like(phi_n)

    if verbose:
        print a, Df, eps, M, dx, dt
        f = plt.figure()
        ax = f.add_subplot(111)
        ax.set_title("t={}".format(t))
        ax.plot(x, c, label="central diff")
        ax.plot(x, phi_n, label="phi n")
        ax.plot(x, fp, label="f'")
        plt.show()

    for i, phi in enumerate(phi_n[:-1]):
        res[i] = phi_n[i] + M * dt * (eps * eps * c[i] - fp[i])
    res[-1] = res[-2]
    return res

def vel_1D(tau, eps, Df, a):
    return - 6 * eps * Df / (a * tau)

def phi_1D_explicit(x, t, tau, eps, Df, a):
    """ explicit solution to the diff eq
    """
    return .5 * (1 - tanh((a*(x - vel_1D(tau, eps, Df, a) * t) / 2 * eps)))

if __name__ == '__main__':
    # require | Df | < a^2 / 6
    # require eps < a

    x0 = -10.
    x_max = 50.
    dx = 0.06
    x_all = arange(x0, x_max, dx)

    t0 = 0.
    t_max = 1
    dt = 0.00036
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

    analytical = False
    if analytical:
        for t in t_all:
            ax.plot(x_all, phi_1D_explicit(x_all, t, tau, eps, Df, a),
                    label="$t={}$".format(t))
    else:
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
