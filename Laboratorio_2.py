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

Qtn=input("¿Desea escuchar el audio? (S/N): ")
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
        #Meterle A(k) con sus correspondientes frecuencias y no esto de los "valores"
        valores = ['lowpass', 'highpass', 'bandpass', 'bandstop']
        btype = np.random.choice(valores)

    print(btype)
    
    senales_1_1 = ["Enventanado", "Muestreo en Frecuencia", "Parks", "Elíptico"]
    #print(señales)
    # Mostramos el menú
    print("Tipos de filtros IIR disponibles:")
    for i, x in enumerate(senales_1_1):
        print(f"{i + 1}) {x}")

    seleccion_2 = int(input("Ingrese su selección: "))

    """    N = int(input("Ingrese el orden del filtro: "))          # Order
    Wn = int(input("Ingrese la frecuencia de corte: "))      # Cutoff frequency in Hz
    
    #Este solo aplica para el cheby2 y ellipt    
    rs = int(input("Ingrese el orden del filtro: "))        # Stopband ripple

    #Este solo aplica para el cheby1 y ellipt
    rp = int(input("Ingrese el : "))       """  # Bandpass ripple
    
##################################BUTTER SET#######################################################  
    rp = int(input("Ingrese el : "))         # Bandpass ripple
    
    if seleccion_2 == 1:

        print("Eligió Butter")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        Fs = 4*Wn      

        b, a = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs)


#####################################CHEBYSHOV_1###################################################
            
    elif seleccion_2 == 2:
        
        print("Eligió Chebyshov I")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      

        b, a = signal.cheby1(N, rp, Wn, btype=btype, analog=False, fs=Fs)

################################CHEBYSHOV II################################
    elif seleccion_2 == 3:
        
        print("Eligió Chebyshov II")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        Fs = 4*Wn      

        b, a = signal.cheby2(N, rs, Wn, btype=btype, analog=False, fs=Fs)

######################################ELIPTICO#######################################################            

    elif seleccion_2 == 4:
        
        print("Eligió Eliptico: ")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      

        b, a = signal.ellip(N, rp, rs, Wn, btype=btype, analog=False, fs=Fs)


elif seleccion == 'IIR':

    ################################################################
    btype_filtros = ["Pasa-bajas", "Pasa-altas", "Pasa-banda", "Rechaza-banda"]

    print("Filtros disponibles:")
    for i, x in enumerate(btype_filtros):
        print(f"{i + 1}) {x}")

    btype_S = int(input("Ingrese el tipo de filtro que desea: "))

    if btype_S == 1:
        btype = 'lowpass'  #N = 10
        N = 10
    elif btype_S == 2:
        btype = 'highpass'  #N = 10
        N = 10
    elif btype_S == 3:
        btype = 'bandpass' #N = 5 o 6
        N = 5
    elif btype_S == 4:
        btype = 'bandstop' #N = 5 o 6
        N = 5

    rs = 3        
    rp = 30
    Fs = sr

    #print(btype)

    # Lista de filtros IIR disponibles
    senales_1 = ["Butterworth", "Chebyshov I", "Chebyshov II", "Elíptico"]
    # Mostramos el menú
    print("Tipos de filtros IIR disponibles:")
    for i, x in enumerate(senales_1):
        print(f"{i + 1}) {x}")

    seleccion_2 = int(input("Ingrese su selección: "))


##################################BUTTER SET#######################################################  

    if seleccion_2 == 1:
        
        print("Eligió Butter")
        #Nosotros definimos N
        #N = int(input("Ingrese el orden del filtro: "))          
        #Wn = int(input("Ingrese la frecuencia de corte: ")) 
        #De momento 
        #Ingrese la frecuencia de corte 1 y la f de corte 2 hay que preguntarlas al usuario
        #Si es pasa altas o pasabajas solo se pide un Wn sino se piden dos Wn
        W_n = []  # Default initialization


        if btype_S == 1:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
            print(W_n)
        elif btype_S == 2:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
        elif btype_S == 3:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))
            print(W_n)
        elif btype_S == 4:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))

        # Other parameters
        #Fs = sr

        # Design the Butterworth filter
        #b, a = signal.butter(N, W_n, btype=btype, analog=False, fs=Fs, output='ba')
        #Wn = [1e3, 3e3]
        b, a = signal.butter(N, W_n, btype, analog=False, output='ba', fs=Fs)

#####################################CHEBYSHOV_1###################################################
            
    elif seleccion_2 == 2:
        
        print("Eligió Chebyshov I")

        W_n = []  # Default initialization


        if btype_S == 1:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
            print(W_n)
        elif btype_S == 2:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
        elif btype_S == 3:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))
            print(W_n)
        elif btype_S == 4:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))

        #Wn = [1e3, 3e3]

        b, a = signal.cheby1(N, rp, W_n, btype, analog=False, output='ba', fs=Fs)

################################CHEBYSHOV II################################
    elif seleccion_2 == 3:
        
        print("Eligió Chebyshov II")

        #Wn = [1e3, 3e3]         
        W_n = []  # Default initialization


        if btype_S == 1:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
            print(W_n)
        elif btype_S == 2:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
        elif btype_S == 3:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))
            print(W_n)
        elif btype_S == 4:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))

        b,a = signal.cheby2(N, rs, W_n, btype, analog=False, output='ba', fs=Fs)

######################################ELIPTICO#######################################################            
        
    elif seleccion_2 == 4:
        
        print("Eligió Eliptico: ")

        #Wn = [1e3, 3e3] 
        W_n = []  # Default initialization


        if btype_S == 1:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
            print(W_n)
        elif btype_S == 2:
            W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
        elif btype_S == 3:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))
            print(W_n)
        elif btype_S == 4:
            W_n = list(map(float, input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: ").split()))

        b, a = signal.ellip(N, rs, rp, W_n, btype, analog=False, output='ba', fs=Fs)

####################################GRÁFICOS#######################################################

w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
w, H = signal.freqz(b,a) #Sacar la respuesta en frecuencia
y_filtrada = signal.lfilter(b, a, y)

Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N): ")
if Qtn == 'S':
    sd.play(y, sr)
    sd.wait()
    sd.play(y_filtrada, sr)
    sd.wait()
else:
    print("Entendido.")

print(y)
print(y_filtrada)

#f = w*Fs/(2*np.pi)
plt.plot(w*Fs/(2*np.pi),np.abs(H))
plt.xlabel("Frecuencia")
plt.ylabel("Amplitud (dB)")
plt.title("Respuesta en frecuencia de Magnitud")
plt.show()


#Periodo para gráficar los audios
T = 1/sr

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