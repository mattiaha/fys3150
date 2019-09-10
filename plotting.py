# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:36:26 2019

@author: matti
"""

import numpy as np
import matplotlib.pyplot as plt

f = open("1000.txt","r")
a =np.loadtxt(f,dtype=float)

n = 1000
v = np.zeros(n)
b = np.zeros(n)
for i in range(n):
    v[i] = a[i]
    b[i] =a[1000+i]
x = np.linspace(0,1,n)

plt.plot(x,v,'r')
plt.plot(x,b,'b')
plt.title('n = 1000')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.show()