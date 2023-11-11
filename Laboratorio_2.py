import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
#import mir_eval.sonify
import sounddevice as sd
import soundfile as sf

# Preguntar al usuario si quiere que la interacción sea en inglés o en español
idioma = input("¿Quieres que la interacción sea en inglés o en español? (inglés/español): ")

# Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
ruta_audio = input("Ingresa la ruta del archivo de audio: ")

#data, fs = sf.read('signal_disturbed.wav')

#sd.play(data,fs)

# Procesar el archivo de audio
y, sr = librosa.load(ruta_audio)

#sd.play(y, sr)
#sd.wait()
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

# Solicitamos al usuario que seleccione el tipo de filtro FIR o IIR
seleccion = input("Qué tipo de filtro desea aplicar, ¿FIR o IIR?: ")

if seleccion == 'FIR':

################################################################

    # Lista de filtros FIR disponibles
    senales = ["Enventanado", "Muestreo en frecuencia", "Parks-McClellan"]

    #print(señales)
    # Mostramos el menú
    print("Tipos de filtros FIR disponibles:")
    for i, x in enumerate(senales):
        print(f"{i + 1}) {x}")

    seleccion_1 = int(input("Ingrese su selección: "))

################################################################
    btype_filtros = ["Pasa-bajas", "Pasa-altas", "Pasa-banda", "Rechaza-banda", "Arbitrario"]

    print("Filtros disponibles:")
    for i, x in enumerate(btype_filtros):
        print(f"{i + 1}) {x}")

    btype_S = int(input("Ingrese el tipo de filtro que desea: "))

    if btype_S == 1:
        btype = 'lowpass'
    elif btype_S == 2:
        btype = 'highpass'
    elif btype_S == 3:
        btype = 'bandpass'
    elif btype_S == 4:
        btype = 'bandstop'
    elif btype_S == 5:
        valores = ['lowpass', 'highpass', 'bandpass', 'bandstop']
        btype = np.random.choice(valores)

    print(btype)

elif seleccion == 'IIR':

    ################################################################
    btype_filtros = ["Pasa-bajas", "Pasa-altas", "Pasa-banda", "Rechaza-banda", "Arbitrario"]

    print("Filtros disponibles:")
    for i, x in enumerate(btype_filtros):
        print(f"{i + 1}) {x}")

    btype_S = int(input("Ingrese el tipo de filtro que desea: "))

    if btype_S == 1:
        btype = 'lowpass'
    elif btype_S == 2:
        btype = 'highpass'
    elif btype_S == 3:
        btype = 'bandpass'
    elif btype_S == 4:
        btype = 'bandstop'
    elif btype_S == 5:
        valores = ['lowpass', 'highpass', 'bandpass', 'bandstop']
        btype = np.random.choice(valores)

    print(btype)

    # Lista de filtros IIR disponibles
    senales_1 = ["Butterworth", "Chebyshov I", "Chebyshov II", "Elíptico"]

    #print(señales)
    # Mostramos el menú
    print("Tipos de filtros IIR disponibles:")
    for i, x in enumerate(senales_1):
        print(f"{i + 1}) {x}")

    seleccion_2 = int(input("Ingrese su selección: "))

    N = int(input("Ingrese el orden del filtro: "))          # Order
    Wn = int(input("Ingrese la frecuencia de corte: "))      # Cutoff frequency in Hz
    
    #Este solo aplica para el cheby2 y ellipt    
    rs = int(input("Ingrese el orden del filtro: "))        # Stopband ripple

    #Este solo aplica para el cheby1 y ellipt
    rp = int(input("Ingrese el : "))         # Bandpass ripple
    
    if seleccion_2==1:
        
        print("Eligio Butter")