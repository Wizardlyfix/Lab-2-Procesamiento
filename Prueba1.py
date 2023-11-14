import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
#import mir_eval.sonify
import sounddevice as sd
import soundfile as sf

y, sr = librosa.load('/Users/wizard/Library/CloudStorage/OneDrive-INSTITUTOTECNOLOGICOMETROPOLITANO-ITM/Procesamiento de Señales/Laboratorio_2/Grabacion.wav')

if len(y.shape)==2:
    # if audio is stereo
    channel1 = y[:,0]
    channel2 = y[:,1]
    y = (channel1 + channel2)/2

btype = 'bandpass' #N = 5 o 6
N = 2
Wn = [4e3, 7e3] 
Fs = sr

b, a = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs, output='ba')

w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
y_filtrada = signal.lfilter(b, a, y)

f = w*Fs/(2*np.pi)

plt.title('Digital filter frequency response')
plt.plot(f,np.abs(H),'b') # Blue color line
plt.ylabel('Magnitude', color='b')
plt.xlabel('Frequency [Hz]')

plt.xlim((0,20))
plt.show()

#Periodo para gráficar los audios
T = 1/sr

tam = np.size(y)
t = np.arange(0, tam*T,T)

tam1 = np.size(y_filtrada)
t1 = np.arange(0, tam1*T,T)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))

ax1.plot(t, y, 'b')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Amplitud')
ax1.set_title('Audio sin filtrar')
ax1.grid(True)

ax2.plot(t1, y_filtrada, 'r')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Amplitud')
ax2.set_title('Audio filtrado')
ax2.grid(True)

plt.subplots_adjust(wspace=0.5)
plt.show()
