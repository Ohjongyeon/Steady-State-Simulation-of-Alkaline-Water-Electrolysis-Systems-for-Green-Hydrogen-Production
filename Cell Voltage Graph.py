# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:37:44 2021

@author: user
"""

from matplotlib import pyplot as plt
from AWECV import Cell_Voltage as CV

Current = []
CV40 = []
CV60 = []
CV80 = []

for i in range(105):
    
    CV1 = CV(313, 1.013, 30, i).CV_Abdin()
    CV2 = CV(333, 1.013, 30, i).CV_Abdin()
    CV3 = CV(353, 1.013, 30, i).CV_Abdin()

    Current.append(i*1000/300)
    CV40.append(CV1)
    CV60.append(CV2)
    CV80.append(CV3)
    
plt.plot(Current, CV40, label='40°C')
plt.plot(Current, CV60, label='60°C')
plt.plot(Current, CV80, label='80°C')
plt.xlim(0, 350)
plt.ylim(1, 2.5)
plt.xlabel('current density[mA/cm2]')
plt.ylabel('Voltage[V]')
plt.legend()
plt.grid()
plt.show()