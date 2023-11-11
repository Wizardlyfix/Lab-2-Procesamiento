# -*- coding: utf-8 -*-
"""
Created on Tue May  9 21:46:33 2023

@author: Andres_Studio
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftfreq

fsamp=10e3
fr=100
Ts=1/fsamp
ti=-1
tf=1

t=np.arange(ti,tf,Ts)
x=signal.square(2*np.pi*fr*t, 0.5)

N=4096
X=fft(x,N)/len(x)
freq=fftfreq(N,Ts)

Sx= (np.abs(X))**2

fcorte=[250,350]
fcorte_rad=[2*np.pi*250,2*np.pi*350]
Ganancia=1
Atenuacion=40
Orden=10
frechazo=[10,190]



#N,Wn=scipy.signal.buttord(fcorte,frechazo,Ganancia,Atenuacion,True)
sos=scipy.signal.butter(Orden,fcorte, 'bandpass',fs=fsamp,output='sos')
y=signal.sosfilt(sos,x)

b,a = scipy.signal.butter(Orden,fcorte_rad, 'bandpass',True,output='ba')
w,H =signal.freqs(b,a)
ws=w/(2*np.pi)


plt.plot(t,x)
plt.xlim(0,4/fr)
plt.ylim(-1.5,1.5)
plt.grid()
