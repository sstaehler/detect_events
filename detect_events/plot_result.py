#!/usr/bin/env python
"""

"""
__author__ = "Simon Stähler"

import numpy as np
import matplotlib.pyplot as plt

with open('result.txt', 'r') as f:
    for line in f.readlines():
        dist, mag, det = line.split(',')
        plt.plot(float(dist), float(mag), 'o', c='C%d' % (int(det)))

plt.show()
