# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 19:44:18 2022

@author: ardau
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import *
init_printing(use_unicode=True)
    
ma = 1/2
mb = 1/4

w = symbols('w')

# Ma = np.array([[-ma*(w**2)+2, -1], [1, 0]])
# Mb = np.array([[-mb*(w**2)+2, -1], [1, 0]])
Ma = np.array([[-(w**2)+1/ma, -1], [1, 0]])
Mb = np.array([[-(w**2)+1/mb, -1], [1, 0]])

#For case  ABA
ABA = Ma.dot(Mb.dot(Ma))
eqn1 = np.trace(ABA)-2
result1 = nroots(eqn1, n = 15, maxsteps=100)

#For case  BAABA
BAABA = Mb.dot(Ma.dot(Ma.dot(Mb.dot(Ma))))
eqn2 = np.trace(BAABA)-2
result2 = nroots(eqn2, n = 15, maxsteps=100)

#For case  ABABAABA
ABABAABA = Ma.dot(Mb.dot(Ma.dot(Mb.dot(Ma.dot(Ma.dot(Mb.dot(Ma)))))))
eqn3 = np.trace(ABABAABA)-2
result3 = nroots(eqn3, n = 15, maxsteps=100)

#For case  BAABAABABAABA
BAABAABABAABA = Mb.dot(Ma.dot(Ma.dot(Mb.dot(Ma.dot(Ma.dot(Mb.dot(Ma.dot(Mb.dot(Ma.dot(Ma.dot(Mb.dot(Ma))))))))))))
eqn4 = np.trace(BAABAABABAABA)-2
result4 = nroots(eqn4, n = 15, maxsteps=100)

plt.figure(figsize = (10, 8))

def y_tick_create(result):
    
    y_list = []
    i = 0

    while True:
        if i+1 <= len(result):
            y_list.append(list(np.linspace(float(result[i]), float(result[i+1]), 1000)))
        else:
            break
        i += 2
    return y_list


def y_plot(y_list):
    
    color_list = ["blue", "red", "green", "black"]
    label_name_list = ["ABA", "BAABA", "ABABAABA", "BAABAABABAABA"]
    x_ticks_list = [[1 for i in range(1000)], [2 for i in range(1000)], [3 for i in range(1000)], [4 for i in range(1000)]]
    
    if len(y_list) == 3:
        label_name = label_name_list[0]
        colorr = color_list[0]
        x_ticks = x_ticks_list[0]
    elif len(y_list) == 5:
        label_name = label_name_list[1]
        colorr = color_list[1]
        x_ticks = x_ticks_list[1]
    elif len(y_list) == 8:
        label_name = label_name_list[2]
        colorr = color_list[2]
        x_ticks = x_ticks_list[2]
    else:
        label_name = label_name_list[3]
        colorr = color_list[3]
        x_ticks = x_ticks_list[3]
    
    for root in range(len(y_list)):
        if root == 1:
            plt.scatter(x_ticks, y_list[root], color = colorr, label = label_name, s=4)
        else:
            plt.scatter(x_ticks, y_list[root], color = colorr, s=4)
            
    plt.xlabel("Patterns")
    plt.ylabel("w")
    plt.title("w Spectrum of Patterns")
    plt.legend()
    plt.grid("on")
       
        
y_plot(y_tick_create(result1))
y_plot(y_tick_create(result2))
y_plot(y_tick_create(result3))
y_plot(y_tick_create(result4))
plt.show()