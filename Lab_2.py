import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
#import mir_eval.sonify
import sounddevice as sd
import soundfile as sf
"""
# Preguntar al usuario si quiere que la interacción sea en inglés o en español
idioma = input("¿Quieres que la interacción sea en inglés o en español? (inglés/español): ")

# Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
ruta_audio = input("Ingresa la ruta del archivo de audio: ")
"""
#data, fs = sf.read('signal_disturbed.wav')

#sd.play(data,fs)

# Procesar el archivo de audio
y, sr = librosa.load(ruta_audio)


sd.play(y, sr)
sd.wait()


Audio(data=y, rate=sr)

# switch audio to mono

if len(y.shape)==2:
    # if audio is stereo
    channel1 = y[:,0]
    channel2 = y[:,1]
    y = (channel1 + channel2)/2

Qtn=input("¿Desea escuchar el audio? (S/N)")
if Qtn=='S':
    
    sd.play(y, sr)
    sd.wait()
else: 
    print("Entendido.")