import numpy as no
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
import mir_eval.sonify

# Preguntar al usuario si quiere que la interacción sea en inglés o en español
idioma = input("¿Quieres que la interacción sea en inglés o en español? (inglés/español): ")

# Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
ruta_audio = input("Ingresa la ruta del archivo de audio: ")

# Procesar el archivo de audio
y, sr = librosa.load(ruta_audio)

Audio(ruta_audio)


