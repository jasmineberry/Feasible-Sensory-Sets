import numpy as np 
from scipy.io import loadmat
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

x = -0.3231091
y =  0.4434951

c = (math.pow(x,2) + math.pow(y,2) - math.pow(.35,2) - math.pow(.27,2)) / (2 * .35 * .27)
s = np.sqrt(1 - math.pow(c,2))
a = np.arcsin((y * (.35 + .27 * c) - x * .27 * s) / (math.pow(x,2) + math.pow(y,2)))
b = np.arccos((math.pow(x,2) + math.pow(y,2) - math.pow(.35,2) - math.pow(.27,2)) / (2 * .35 * .27)) 

print(a *180/np.pi, b*180/np.pi)