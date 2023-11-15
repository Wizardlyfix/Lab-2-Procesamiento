import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
#import mir_eval.sonify
import sounddevice as sd
import soundfile as sf

y, Fs = librosa.load('\\Users\wizar\OneDrive - INSTITUTO TECNOLOGICO METROPOLITANO - ITM\Procesamiento de Señales\Laboratorio_2\Grabacion.wav')

if len(y.shape)==2:
    # if audio is stereo
    channel1 = y[:,0]
    channel2 = y[:,1]
    y = (channel1 + channel2)/2


b,a = signal.ellip(5, 3, 30, [1e3, 3e3], btype='bandpass', analog=False, output='ba', fs=Fs)
#w = np.linspace(0, np.pi, int(10e3))
w, H = signal.freqz(b,a)
plt.plot(w*Fs/(2*np.pi),np.abs(H))
plt.show()
