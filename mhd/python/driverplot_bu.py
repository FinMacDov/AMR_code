# -*- coding: utf-8 -*-
"""
Created on Fri May 18 11:21:03 2018

@author: fionnlagh
"""
import numpy as np
import matplotlib.pyplot as plt

unit_time = 8.5873143947605932  # s
cm_2_km = 1e-5
A = 6e6*cm_2_km  # km

phase = 2*np.pi/10
phaseqt = 2*np.pi/(10*unit_time/unit_time)

endpt = 200  # s
endptqt = endpt/unit_time  # no dim

t_change = endpt/2  # s
qt_change = t_change/unit_time  # no dim


t0 = 0
t1 = 500

nbpts = 1000
time = np.linspace(t0, t1, nbpts)
qtime = time/unit_time
v2 = []
for qt in range(np.size(qtime)):
    if qtime[qt] >= qt_change:
        if qtime[qt]-endptqt > 0:
            v2.append(0.0)
        else:
            v2.append(A*-np.tanh(phaseqt*(qtime[qt]-endptqt)))
    else:
        v2.append(A*np.tanh(phaseqt*qtime[qt]))

plt.xlabel('time [s]')
plt.ylabel('v2 [km s-1]')
plt.ylim(np.min(v2)-5, np.max(v2)+5)
plt.plot(time, v2)
plt.grid()
plt.show()
