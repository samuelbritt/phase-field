#!/usr/bin/env python
# encoding: utf-8

import matplotlib.pyplot as plt
import numpy as np

import os
import sys

def parse_t(fname):
    for line in open(fname):
        if line.startswith("#"):
            continue
        else:
            t = float(line.strip())
            break
    return t

def plot(fname, ax):
    t = parse_t(fname)
    x, phi = np.loadtxt(fname, skiprows=2, unpack=True)
    ax.plot(x, phi, label = "$t = {}$".format(t))


if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111)

    data_dir = sys.argv[1]
    data_files = os.listdir(data_dir)

    for f in data_files:
        plot(os.path.join(data_dir, f), ax)

    ax.legend(loc='best')
    plt.show()
