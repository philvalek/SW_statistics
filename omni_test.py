# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:31:27 2019

@author: Phil
"""



import heliopy.data.omni as omni
import matplotlib.pyplot as plt
from datetime import datetime

starttime = datetime(1970, 1, 1, 0, 0, 0)
endtime = datetime(1971, 1, 3, 0, 0, 0)

omni_data = omni.low(starttime, endtime)

fig, axs = plt.subplots(3, 1, sharex=True)
axs[0].plot(omni_data.data['Bx GSE, GSM'])
axs[1].plot(omni_data.data['By GSE'])
axs[1].plot(omni_data.data['Bz GSE'])
axs[2].plot(omni_data.data['By GSM'])
axs[2].plot(omni_data.data['Bz GSM'])

for ax in axs:
    ax.legend()
plt.show()
