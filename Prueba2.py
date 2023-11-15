import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
#import mir_eval.sonify
import sounddevice as sd
import soundfile as sf

y, Fs = librosa.load('/Users/wizard/Library/CloudStorage/OneDrive-INSTITUTOTECNOLOGICOMETROPOLITANO-ITM/Procesamiento de Señales/Laboratorio_2/Grabacion.wav')

if len(y.shape)==2:
    # if audio is stereo
    channel1 = y[:,0]
    channel2 = y[:,1]
    y = (channel1 + channel2)/2

w, H = signal.freqz(b, a)

y_filtrada = signal.lfilter(b, a, y)

plt.plot(w*Fs/(2*np.pi), 20*np.log10(np.abs(H)))
plt.show()

angles = np.unwrap(np.angle(H))
plt.plot(w*Fs/(2*np.pi), angles*180/np.pi)
plt.xlabel("Frecuencia [Hz]")
plt.ylabel('Fase [°]')
plt.title("Respuesta en frecuencia de Fase")
plt.show()

#####################################MAPA DE CALOR##########################################

f1, t1, Sxx1 = signal.spectrogram(y, Fs, scaling = 'density')
plt.pcolormesh(t1, f1, 10*np.log10(Sxx1), shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

f2, t2, Sxx2 = signal.spectrogram(y_filtrada, Fs, scaling = 'density')
plt.pcolormesh(t2, f2, 10*np.log10(Sxx2), shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()

#Periodo para gráficar los audios
T = 1/Fs

tam = np.size(y)
t = np.arange(0, tam*T,T)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))

ax1.plot(t, y, 'b')
ax1.set_xlabel('Tiempo (s)')
ax1.set_ylabel('Amplitud')
ax1.set_title('Audio sin filtrar')
ax1.grid(True)

ax2.plot(t, y_filtrada, 'r')
ax2.set_xlabel('Tiempo (s)')
ax2.set_ylabel('Amplitud')
ax2.set_title('Audio filtrado')
ax2.grid(True)

plt.subplots_adjust(wspace=0.5)
plt.show() 