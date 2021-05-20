# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:37:44 2021

@author: user
"""

from matplotlib import pyplot as plt
from AWECV import Hydrogen_Production as HV

Current = []
HV40 = []
HV60 = []
HV80 = []

for i in range(105):
    
    HV1 = HV(313, 1.013, 30, i).HP_Ulleberg()
    HV2 = HV(333, 1.013, 30, i).HP_Ulleberg()
    HV3 = HV(353, 1.013, 30, i).HP_Ulleberg()
    
    Current.append(i*1000/300)
    HV40.append(HV1)
    HV60.append(HV2)
    HV80.append(HV3)
    
plt.plot(Current, HV40, label='40°C')
plt.plot(Current, HV60, label='60°C')
plt.plot(Current, HV80, label='80°C')
plt.xlim(0, 350)
plt.ylim(0, 0.05)
plt.xlabel('current_density[mA/cm2]')
plt.ylabel('Voltage[V]')
plt.legend()
plt.grid()
plt.show()