#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 18:36:38 2020

@author: tusnin
"""

import numpy as np
from sys import argv

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors

mpl.rcParams['figure.dpi'] = 300
#path = str(argv[1])
path='/home/tusnin/Documents/C_LLE/data/'
Nphi = 2**9
Ndet = 400
dir_name = '../'#'det_1.88_ndet_3000_f2_9_T_0.1_1e-4_1_2_pi2_045_J-9/'
path+=dir_name
Field = np.fromfile(path+'Field.bin', dtype=complex)
Field = np.reshape(Field,(Ndet,Nphi))
detuning = np.fromfile(path+'Detuning.bin')
phi = np.fromfile(path+'Phi.bin')
#print x, type(x)

fig = plt.figure(figsize=[3.6*2/3,2.2],frameon=False)
ax = fig.add_subplot(1,1,1)
#ax.set_yticks([-25,-20,-15,-10,-5,0,5,10,15,20,25])

ax.pcolormesh(phi,detuning,np.log(abs(Field)))
plt.show()

Power = np.zeros(detuning.size)
i_det = np.arange(0,detuning.size)

Power[i_det] = np.linalg.norm(Field[:,i_det],axis=1)


fig2 = plt.figure(figsize=[3.6*2/3,2.2],frameon=False)
ax = fig2.add_subplot(1,1,1)
#ax.set_xticks([-25,-20,-15,-10,-5,0,5,10,15,20,25])

ax.plot(detuning,Power)
plt.show()