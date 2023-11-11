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
    
    senales_1_1 = ["Enventanado", "Muestreo en Frecuencia", "Parks", "Elíptico"]
    #print(señales)
    # Mostramos el menú
    print("Tipos de filtros IIR disponibles:")
    for i, x in enumerate(senales_1):
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
    
    if seleccion_2==1:
        
        print("Eligió Butter")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        Fs = 4*Wn      


        b, a = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs)

        y_filtrada = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N): ")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada, sr)
            sd.wait()
        else:
            print("Entendido.")
           
 #####################################CHEBYSHOV_1###################################################
            
    if seleccion_2==2:
        
        print("Eligió Chebyshov I")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      


        b, a = signal.cheby1(N, rp, Wn, btype=btype, analog=False, fs=Fs)

        y_filtrada_Chb1 = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_Chb1, sr)
            sd.wait()
        else:
            print("Entendido.")

################################CHEBYSHOV II################################
    if seleccion_2==3:
        
        print("Eligió Chebyshov II")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        Fs = 4*Wn      


        b, a = signal.cheby2(N, rs, Wn, btype=btype, analog=False, fs=Fs)

        y_filtrada_Chb2 = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_Chb2, sr)
            sd.wait()
        else:
            print("Entendido.")
            
######################################ELIPTICO#######################################################            
            

    if seleccion_2==4:
        
        print("Eligió Eliptico: ")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      


        b, a = signal.ellip(N, rp, rs, Wn, btype=btype, analog=False, fs=Fs)

        y_filtrada_Chb2 = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_Chb2, sr)
            sd.wait()
        else:
            print("Entendido.")

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


##################################BUTTER SET#######################################################  

    
    if seleccion_2==1:
        
        print("Eligió Butter")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        Fs = 4*Wn      


        sos = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs, output='sos')
        b, a = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs, output='ba')
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
        y_filtrada = signal.sosfilt(sos, y)

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
        
        f = w*Fs/(2*np.pi)

        fig, ax1 = plt.subplots()
        ax1.set_title('Digital filter frequency response')
        ax1.plot(f,np.abs(H),'b') # Blue color line
        ax1.set_ylabel('Magnitude', color='b')
        ax1.set_xlabel('Frequency [Hz]')

        ax2 = ax1.twinx()

        angles = np.unwrap(np.angle(H))
        ax2.plot(f, angles*180/np.pi, 'g') # Phase converted to degrees, and green color line
        ax2.set_ylabel('Phase [°]', color='g')
        ax2.grid()
        ax2.set_xlim((0, 20))
        ax2.axis('tight')

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

#####################################CHEBYSHOV_1###################################################
            
    if seleccion_2==2:
        
        print("Eligió Chebyshov I")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      


        b, a = signal.cheby1(N, rp, Wn, btype=btype, analog=False, fs=Fs)
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
        y_filtrada_Chb1 = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_Chb1, sr)
            sd.wait()
        else:
            print("Entendido.")

################################CHEBYSHOV II################################
    if seleccion_2==3:
        
        print("Eligió Chebyshov II")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        Fs = 4*Wn      


        b, a = signal.butter(N, Wn, btype=btype, analog=False, fs=Fs)
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
        y_filtrada_Chb2 = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_Chb2, sr)
            sd.wait()
        else:
            print("Entendido.___")
            print("Entendido.")
            
######################################ELIPTICO#######################################################            
        
    if seleccion_2==4:
        
        print("Eligió Eliptico: ")
        N = int(input("Ingrese el orden del filtro: "))          
        Wn = int(input("Ingrese la frecuencia de corte: "))  
        rs = int(input("Ingrese el Stopband ripple : "))         
        rp = int(input("Ingrese el Bandpass Ripple : "))         
        Fs = 4*Wn      


        b, a = signal.ellip(N, rp, rs, Wn, btype=btype, analog=False, fs=Fs)
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia
        y_filtrada_ellip = signal.lfilter(b, a, y)

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N)")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada_ellip, sr)
            sd.wait()
        else:
            print("Entendido.")

