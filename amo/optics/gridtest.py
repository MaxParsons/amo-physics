'''
Created on Feb 9, 2015

@author: Max
'''
# What's the fastest way to evaluate a function over a grid'

import numpy as np
import datetime
import matplotlib.pyplot as plt


def gaussian(x, y, z):
    return np.exp(-x ** 2 - y ** 2 - z ** 2)

# evaluate with empty arrays
xaxis = np.linspace(0, 1, 50)
yaxis = np.linspace(0, 1, 50)
zaxis = np.linspace(0, 1, 300)
start = datetime.datetime.now()
result = gaussian(xaxis[:, None, None], yaxis[None, :, None], zaxis[None, None, :])
end = datetime.datetime.now()
interval = end - start
print "Time using empty arrays:" + str(interval.microseconds)

# evaluate with meshgrid
xaxis = np.linspace(-1, 1, 50)
yaxis = np.linspace(-1, 1, 50)
zaxis = np.linspace(-1, 1, 50)
start = datetime.datetime.now()
x, y, z = np.meshgrid(xaxis, yaxis, zaxis)
result = gaussian(x, y, 0.0)[:, :, 0]
plt.contour(xaxis, yaxis, result)
end = datetime.datetime.now()
interval = end - start
print "Time using meshgrid:" + str(interval.microseconds)
plt.show()
