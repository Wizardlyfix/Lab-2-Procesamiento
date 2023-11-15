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


b,a = signal.cheby2(5, 30, [1e3, 3e3], btype='bandpass', analog=False, output='ba', fs=Fs)
#w = np.linspace(0, np.pi, int(10e3))
w, H = signal.freqz(b,a)
plt.plot(w*Fs/(2*np.pi),np.abs(H))
plt.show()

angles = np.unwrap(np.angle(H))
plt.plot(w*Fs/(2*np.pi), angles*180/np.pi)
plt.xlabel("Frecuencia [Hz]")
plt.ylabel('Fase [°]')
plt.title("Respuesta en frecuencia de Fase")
plt.show()

#####################################MAPA DE CALOR##########################################

f1, t1, Sxx = signal.spectrogram(y, Fs)
plt.pcolormesh(t1, f1, 20*np.log10(Sxx), shading='gouraud')
plt.ylabel('Frequency [Hz]')
plt.xlabel('Time [sec]')
plt.show()