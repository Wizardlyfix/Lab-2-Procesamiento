import numpy as np
import matplotlib.pyplot as plt
import librosa
from scipy import signal
from IPython.display import Audio
import sounddevice as sd
import soundfile as sf
from matplotlib.gridspec import GridSpec

# Preguntar al usuario si quiere que la interacción sea en inglés o en español
idioma = input("¿Quieres que la interacción sea en inglés o en español? (ingles/español): ").lower()

if idioma == 'español':

    # Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
    ruta_audio = input("Ingresa la ruta del archivo de audio: ")

    # Procesar el archivo de audio
    y, sr = librosa.load(ruta_audio)

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
        Fs = sr

        if btype_S == 1:
            btype = 'lowpass'
            bandera = 1
            while True:
                print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                f_cutoff = float(input(f"Ingrese una única frecuencia de corte para el filtro {btype} elegido: "))
                if 0 <= f_cutoff <= Fs / 2:
                    break
                else:
                    print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 2:
            btype = 'highpass' 
            bandera = 1
            while True:
                    print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                    f_cutoff = float(input(f"Ingrese una única frecuencia de corte para el filtro {btype} elegido: "))
                    if 0 <= f_cutoff <= Fs / 2:
                        break
                    else:
                        print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 3:
            btype = 'bandpass'
            bandera = 1
            while True:
                    print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                    f_cutoff = list(map(float, input(f"Ingrese las frecuencias de corte para el filtro {btype} elegido (separadas por espacio): ").split()))
                    if all(0 <= w <= Fs / 2 for w in f_cutoff):
                        break
                    else:
                        print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 4:
            btype = 'bandstop'
            bandera = 1
            while True:
                    print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                    f_cutoff = list(map(float, input(f"Ingrese las frecuencias de corte para el filtro {btype} elegido (separadas por espacio): ").split()))
                    if all(0 <= w <= Fs / 2 for w in f_cutoff):
                        break
                    else:
                        print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 5:
            try:
                bandera = 2
                while True:
                        print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                        #f_cutoff = list(map(float, input(f"Ingrese las frecuencias de corte para el filtro {btype} elegido (separadas por espacio): ").split()))
                        N=1001
                        fre_hz = list(map(float, input(f"Ingrese el vector de frecuencias en Hz (separadas por espacio): ").split()))
                        magn_lin = list(map(float, input(f"Ingrese el vector de magnitudes lineales (números entre 0 y 1 - separadas por espacio): ").split()))
                        
                        if all(0 <= w <= Fs / 2 for w in fre_hz):
                            break
                        else:
                            print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
            except NameError: 
                print("WARNING")

    ############################################### ENVENTANADO #######################################################
        if seleccion_1 == 1:

            print("Eligió Enventanado")
        
            if bandera == 1:
                N=1001        
                h= signal.firwin(N, f_cutoff , window='hann', fs=sr, pass_zero=btype) 
                print(h)
            elif bandera == 2:
                print("Eligió filtro arbitrario")
                N=1001
                h = signal.firwin2(N, fre_hz, magn_lin, window='hann', fs=sr)
                #print(h)

    ##################################### MUESTREO EN FRECUENCIA ###################################################
                
        elif seleccion_1 == 2:
            
            print("Eligió Muestreo en frecuencia")

            cos, pi, flip = np.cos, np.pi, np.flip

            # Orden del filtro
            M=1001

            if M%2==0:
                len_A = int(M/2)
                flag_inv = 0
            else:
                len_A = int((M+1)/2)
                flag_inv = -1
            
            if btype_S ==1:
                print("Filtro pasabajas")
                K_=int((f_cutoff/sr)*M)
                A = np.concatenate((np.ones(K_), np.zeros(len_A - K_)))
                #print(A)
            elif btype_S==2:
                print("Filtro pasa-altas")
                K_=int((f_cutoff/sr)*M)
                A = np.concatenate((np.zeros(K_), np.ones(len_A - K_)))
            
            elif btype_S==3:
                print("Filtro Pasabanda")
                fc = np.array([f_cutoff])
                k = (fc / Fs) * M
                K_ = k.astype(int).ravel()
                A = np.concatenate((np.zeros(K_[0]), np.ones(K_[1]-K_[0]+1), np.zeros(len_A-K_[1]-1)))
            elif btype_S==4:
                print("Filtro rechazabanda")
                fc = np.array([f_cutoff])
                k = (fc / Fs) * M
                K_ = k.astype(int).ravel()
                A = np.concatenate((np.ones(K_[0]), np.zeros(K_[1]-K_[0]+1), np.ones(len_A-K_[1]-1)))
            #[0, fc,    ,fc2,0]

            #print(A)

            # Inicializamos h(n)
            h = np.zeros(len(A))
            # Aplicamos la fórmula
            for n in range(len(h)):
                sum_k = 0
                for k in range(1,len(h)):
                    sum_k += A[k]*(-1)**k*cos(pi*k/M*(2*n+1))
                h[n] = 1/M*(A[0]+2*sum_k)
            # Aplicamos espejo por la simetría
            h_inv = flip(h[:len(h)+flag_inv])
            # Concatenamos para crear el h(n) completo
            h = np.concatenate((h, h_inv))
            
            #PRUEBA A_K
            """A_k=abs(fft(h))
            k=np.arange(M)
            plt.stem(k, A_k)
            plt.title('prueba')
            print(h)"""

################################ Parks - McClellan ################################

        elif seleccion_1 == 3:
            
            print("Eligió Parks-McClellan")
            N=1001 
            trans_width = 50  #Transición de banda
            numtaps = 1001  
            if btype_S==1:
                #Pasa-Bajas
                desired=[1,0]
                h = signal.remez(numtaps, [0, f_cutoff, f_cutoff + trans_width, 0.5*sr],desired, Hz=sr)         
            elif btype_S==2:
                #Pasa-Altas 
                desired=[0,1]
                h = signal.remez(numtaps, [0, f_cutoff, f_cutoff + trans_width, 0.5*sr],desired, Hz=sr)
                print(h)         
            elif btype_S==3:
                #Pasabanda
                desired=[0,1,0]
                edges = [0, f_cutoff[0] - trans_width, f_cutoff[0], f_cutoff[1], f_cutoff[1] + trans_width, 0.5*sr]
                h = signal.remez(numtaps, edges, desired, fs=sr)
                desired=[0,1,0]
                edges = [0, f_cutoff[0] - trans_width, f_cutoff[0], f_cutoff[1], f_cutoff[1] + trans_width, 0.5*sr]
                h = signal.remez(numtaps, edges, desired, fs=sr)
            elif btype_S==4:
                #Rechazabanda
                desired=[1,0,1]
                edges = [0, f_cutoff[0] - trans_width, f_cutoff[0], f_cutoff[1], f_cutoff[1] + trans_width, 0.5*sr]
                h = signal.remez(numtaps, edges, desired, fs=sr)

    ##################################IIR FILTERS####################################################################
    elif seleccion == 'IIR':

        ################################################################
        btype_filtros = ["Pasa-bajas", "Pasa-altas", "Pasa-banda", "Rechaza-banda"]

        print("Filtros disponibles:")
        for i, x in enumerate(btype_filtros):
            print(f"{i + 1}) {x}")

        btype_S = int(input("Ingrese el tipo de filtro que desea: "))
        W_n = []  # Default initialization
    #####################################################################################
        rs = 50 
        rp = 3.01
        Fs = sr
        if btype_S == 1:
            btype = 'lowpass'
            N = 10
            while True:
                print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
                if 0 <= W_n <= Fs / 2:
                    break
                else:
                    print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 2:
            btype = 'highpass'
            N = 10
            while True:
                print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                W_n = float(input(f"Ingrese la frecuencia de corte para el filtro {btype} elegido: "))
                if 0 <= W_n <= Fs / 2:
                    break
                else:
                    print("El valor ingresado está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 3:
            btype = 'bandpass'
            N = 5
            while True:
                print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                W_n = list(map(float, input(f"Ingrese las frecuencias de corte para el filtro {btype} elegido (separadas por espacio): ").split()))
                if all(0 <= w <= Fs / 2 for w in W_n):
                    break
                else:
                    print("Al menos uno de los valores ingresados está fuera del rango permitido. Por favor, inténtelo nuevamente.")
        elif btype_S == 4:
            btype = 'bandstop'
            N = 5
            while True:
                print(f"Tenga en cuenta que el rango de ingreso está en el rango 0--{int(Fs/2)}")
                W_n = list(map(float, input(f"Ingrese las frecuencias de corte para el filtro {btype} elegido (separadas por espacio): ").split()))
                if all(0 <= w <= Fs / 2 for w in W_n):
                    break
                else:
                    print("Al menos uno de los valores ingresados está fuera del rango permitido. Por favor, inténtelo nuevamente.")


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

            b, a = signal.butter(N, W_n, btype, analog=False, output='ba', fs=Fs)

    #####################################CHEBYSHOV_1###################################################

        elif seleccion_2 == 2:

            print("Eligió Chebyshov I")

            b, a = signal.cheby1(N, rp, W_n, btype, analog=False, output='ba', fs=Fs)

    ####################################CHEBYSHOV II###################################################
        elif seleccion_2 == 3:
            
            print("Eligió Chebyshov II")

            b,a = signal.cheby2(N, rs, W_n, btype, analog=False, output='ba', fs=Fs)

    ######################################ELIPTICO#######################################################            
            
        elif seleccion_2 == 4:
            
            print("Eligió Eliptico: ")

            b, a = signal.ellip(N, rp, rs, W_n, btype, analog=False, output='ba', fs=Fs)

    ####################################GRÁFICOS#######################################################

    if seleccion == 'FIR':
        w1, H1 = signal.freqz(h,1) #Sacar la respuesta en frecuencia
        y_filtrada1 = signal.lfilter(h, 1, y)
        
        T = 1/sr

        tam = np.size(y)
        t = np.arange(0, tam*T,T)

        angles1 = np.unwrap(np.angle(H1))

        f1_, t1_, Sxx_ = signal.spectrogram(y, sr, scaling = 'density')
        f2_, t2_, Sxx2_ = signal.spectrogram(y_filtrada1, sr, scaling = 'density')

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N): ")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada1, sr)
            sd.wait()
        else:
            print("Entendido.")
        
        f1 = w1*sr/(2*np.pi)

        gs1 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs2 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])
        gs3 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs4 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])

        fig = plt.figure(figsize=(12, 7))

        ax1 = fig.add_subplot(gs1[0, 0])
        ax1.plot(t, y, 'b', label = 'Señal Original')
        ax1.plot(t, y_filtrada1, 'r', label = 'Señal Filtrada')
        ax1.set_xlabel('Tiempo [s]')
        ax1.set_ylabel('Amplitud')
        ax1.set_title('Señal original y filtrada')
        ax1.legend()
        ax1.grid(True)

        ax2 = fig.add_subplot(gs2[0, 1])
        ax2.pcolormesh(t1_, f1_, 10*np.log10(Sxx_), shading='gouraud')
        ax2.set_ylabel('Frecuencia [Hz]')
        ax2.set_xlabel('Tiempo [s]')
        ax2.set_title('Espectrograma señal original')

        ax3 = fig.add_subplot(gs3[1, 0])
        ax3_fase = ax3.twinx()
        ax3.plot(f1, 20*np.log10(np.abs(H1)), 'b')
        ax3_fase.plot(f1, angles1*180/np.pi, 'r')
        ax3.set_ylabel('Magnitud [dB]', color = 'b')
        ax3_fase.set_ylabel('Fase [°]', color = 'r')
        ax3.set_xlabel('Frecuencia [Hz]')
        ax3.set_title('Respuesta en frecuencia')

        ax4 = fig.add_subplot(gs4[1, 1])
        cmap = plt.colormaps["seismic"]
        ax4.pcolormesh(t2_, f2_, 10*np.log10(Sxx2_), shading='gouraud', cmap=cmap)
        ax4.set_ylabel('Frecuencia [Hz]')
        ax4.set_xlabel('Tiempo [s]')
        ax4.set_title('Espectrograma señal filtrada')

        plt.tight_layout()
        plt.show()

    elif seleccion == 'IIR':
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia

        y_filtrada = signal.lfilter(b, a, y)
        T = 1/sr

        tam = np.size(y)
        t = np.arange(0, tam*T,T)

        angles = np.unwrap(np.angle(H))

        f1, t1, Sxx1 = signal.spectrogram(y, Fs, scaling = 'density')
        f2, t2, Sxx2 = signal.spectrogram(y_filtrada, Fs, scaling = 'density')

        Qtn = input("¿Desea escuchar la señal original y la señal filtrada? (S/N): ")
        if Qtn == 'S':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada, sr)
            sd.wait()
        else:
            print("Entendido.")

        f = w*Fs/(2*np.pi)

        gs1 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs2 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])
        gs3 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs4 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])

        fig = plt.figure(figsize=(12, 7))

        ax1 = fig.add_subplot(gs1[0, 0])
        ax1.plot(t, y, 'b', label = 'Señal Original')
        ax1.plot(t, y_filtrada, 'r', label = 'Señal Filtrada')
        ax1.set_xlabel('Tiempo [s]')
        ax1.set_ylabel('Amplitud')
        ax1.set_title('Señal original y filtrada')
        ax1.legend()
        ax1.grid(True)

        ax2 = fig.add_subplot(gs2[0, 1])
        ax2.pcolormesh(t1, f1, 10*np.log10(Sxx1), shading='gouraud')
        ax2.set_ylabel('Frecuencia [Hz]')
        ax2.set_xlabel('Tiempo [s]')
        ax2.set_title('Espectrograma señal original')
        
        ax3 = fig.add_subplot(gs3[1, 0])
        ax3_fase = ax3.twinx()
        ax3.plot(f, 20*np.log10(np.abs(H)), 'b')
        ax3_fase.plot(f, angles*180/np.pi, 'r')
        ax3.set_ylabel('Magnitud [dB]', color = 'b')
        ax3_fase.set_ylabel('Fase [°]', color = 'r')
        ax3.set_xlabel('Frecuencia [Hz]')
        ax3.set_title('Respuesta en frecuencia')

        ax4 = fig.add_subplot(gs4[1, 1])
        ax4.pcolormesh(t2, f2, 10*np.log10(Sxx2), shading='gouraud')
        ax4.set_ylabel('Frecuencia [Hz]')
        ax4.set_xlabel('Tiempo [s]')
        ax4.set_title('Espectrograma señal filtrada')

        plt.tight_layout()
        plt.show()

elif idioma == 'ingles':

    # Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
    ruta_audio = input("Enter the path of the audio: ")

    # Procesar el archivo de audio
    y, sr = librosa.load(ruta_audio)

    # switch audio to mono
    if len(y.shape)==2:
        # if audio is stereo
        channel1 = y[:,0]
        channel2 = y[:,1]
        y = (channel1 + channel2)/2

    Qtn=input("Do you want to listen the audio? (Y/N): ")
    if Qtn=='Y':
        sd.play(y, sr)
        sd.wait()
    else:
        print("Understood.")

    # Solicitamos al usuario que seleccione el tipo de filtro FIR o IIR
    seleccion = input("What type of filter do you want to apply, FIR or IIR?: ")

    if seleccion == 'FIR':

    ################################################################

        # Lista de filtros FIR disponibles
        senales = ["Windowing", "Frequency Sampling", "Parks-McClellan"]

        #print(señales)
        # Mostramos el menú
        print("Type of FIR filter available: ")
        for i, x in enumerate(senales):
            print(f"{i + 1}) {x}")

        seleccion_1 = int(input("Enter your selection: "))

    ################################################################
        btype_filtros = ["Lowpass", "Highpass", "Bandpass", "Bandstop", "Arbitrary"]

        print("Available filters: ")
        for i, x in enumerate(btype_filtros):
            print(f"{i + 1}) {x}")

        btype_S = int(input("Enter the type of filter that you want: "))
        Fs = sr

        if btype_S == 1:
            btype = 'lowpass'
            bandera = 1
            while True:
                print(f"Please note that the income range is in the range 0--.{int(Fs/2)}")
                f_cutoff = float(input(f"Enter a single cutoff frequency for the filter {btype} choosed: "))
                if 0 <= f_cutoff <= Fs / 2:
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 2:
            btype = 'highpass' 
            bandera = 1
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                f_cutoff = float(input(f"Enter a single cutoff frequency for the filter {btype} choosed: "))
                if 0 <= f_cutoff <= Fs / 2:
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 3:
            btype = 'bandpass'
            bandera = 1
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                f_cutoff = list(map(float, input(f"Enter the cutoff frequencies for the chosen {btype} filter (separated by space): ").split()))
                if all(0 <= w <= Fs / 2 for w in f_cutoff):
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 4:
            btype = 'bandstop'
            bandera = 1
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                f_cutoff = list(map(float, input(f"Enter the cutoff frequencies for the chosen {btype} filter (separated by space): ").split()))
                if all(0 <= w <= Fs / 2 for w in f_cutoff):
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 5:
            try:
                bandera = 2
                while True:
                        print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                        N=1001
                        fre_hz = list(map(float, input(f"Enter the vector of frequencies in Hz (space separated): ").split()))
                        magn_lin = list(map(float, input(f"Enter the vector of linear quantities (numbers between 0 and 1 - separated by space): ").split()))
                        
                        if all(0 <= w <= Fs / 2 for w in fre_hz):
                            break
                        else:
                            print("The value entered is outside the allowed range. Please try again.")
            except NameError: 
                print("WARNING")

    ###################################################################################################
        if seleccion_1 == 1:

            print("Chose Windowing")
            
            if bandera == 1:
                N=1001        
                h = signal.firwin(N, f_cutoff , window='hann', fs=sr, pass_zero=btype) 
                #print(h)
            elif bandera == 2:
                print("Chose arbitrary")
                N=1001
                h = signal.firwin2(N, fre_hz, magn_lin, window='hann', fs=sr)
                #print(h)
    ##################################### FREQUENCY SAMPLING ###################################################
                
        elif seleccion_1 == 2:
            
            print("Chose Frequency Sampling")
            cos, pi, flip = np.cos, np.pi, np.flip
            # Orden del filtro
            M=1001

            if M%2==0:
                len_A = int(M/2)
                flag_inv = 0
            else:
                len_A = int((M+1)/2)
                flag_inv = -1
            
            if btype_S ==1:
                print("Lowpass Filter")
                K_=int((f_cutoff/sr)*M)
                A = np.concatenate((np.ones(K_), np.zeros(len_A - K_)))
            elif btype_S==2:
                print("Highpass Filter")
                K_=int((f_cutoff/sr)*M)
                A = np.concatenate((np.zeros(K_), np.ones(len_A - K_)))
            elif btype_S==3:
                print("Bandpass Filter")
                fc = np.array([f_cutoff])
                k = (fc / Fs) * M
                K_ = k.astype(int).ravel()
                A = np.concatenate((np.zeros(K_[0]), np.ones(K_[1]-K_[0]+1), np.zeros(len_A-K_[1]-1)))
            elif btype_S==4:
                print("Bandstop Filter")
                fc = np.array([f_cutoff])
                k = (fc / Fs) * M
                K_ = k.astype(int).ravel()
                A = np.concatenate((np.ones(K_[0]), np.zeros(K_[1]-K_[0]+1), np.ones(len_A-K_[1]-1)))
    #[0, fc,    ,fc2,0]
            #print(A)

            # Inicializamos h(n)
            h = np.zeros(len(A))
            # Aplicamos la fórmula
            for n in range(len(h)):
                sum_k = 0
                for k in range(1,len(h)):
                    sum_k += A[k]*(-1)**k*cos(pi*k/M*(2*n+1))
                h[n] = 1/M*(A[0]+2*sum_k)
            # Aplicamos espejo por la simetría
            h_inv = flip(h[:len(h)+flag_inv])
            # Concatenamos para crear el h(n) completo
            h = np.concatenate((h, h_inv))
            
            #PRUEBA A_K
            """A_k=abs(fft(h))
            k=np.arange(M)
            plt.stem(k, A_k)
            plt.title('prueba')
            print(h)"""

        elif seleccion_1 == 3:
            
            print("Chose Parks-McClellan")
            N=1001 
            trans_width = 50  #Transición de banda
            numtaps = 1001  
            if btype_S==1:
                #Pasa-Bajas
                desired=[1,0]
                h = signal.remez(numtaps, [0, f_cutoff, f_cutoff + trans_width, 0.5*sr],desired, Hz=sr)         
            elif btype_S==2:
                #Pasa-Altas
                desired=[0,1]
                h = signal.remez(numtaps, [0, f_cutoff, f_cutoff + trans_width, 0.5*sr],desired, Hz=sr)
                print(h)         
            elif btype_S==3:
                #Pasabandas
                desired=[0,1,0]
                edges = [0, f_cutoff[0] - trans_width, f_cutoff[0], f_cutoff[1], f_cutoff[1] + trans_width, 0.5*sr]
                h = signal.remez(numtaps, edges, desired, fs=sr)
            elif btype_S==4:
                #Rechazabandas
                desired=[1,0,1]
                edges = [0, f_cutoff[0] - trans_width, f_cutoff[0], f_cutoff[1], f_cutoff[1] + trans_width, 0.5*sr]
                h = signal.remez(numtaps, edges, desired, fs=sr)

    ##################################IIR FILTERS####################################################################
    elif seleccion == 'IIR':

    #####################################################################################
        btype_filtros = ["Lowpass", "Highpass", "Bandpass", "Bandstop", "Arbitrary"]

        print("Available filters: ")
        for i, x in enumerate(btype_filtros):
            print(f"{i + 1}) {x}")

        btype_S = int(input("Enter the type of filter that you want: "))
        W_n = []  # Default initialization
    #####################################################################################
        rs = 50 
        rp = 3.01
        Fs = sr
        if btype_S == 1:
            btype = 'lowpass'
            N = 10
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                W_n = float(input(f"Enter a single cutoff frequency for the filter {btype} choosed: "))
                if 0 <= W_n <= Fs / 2:
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 2:
            btype = 'highpass' 
            N = 10
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                W_n = float(input(f"Enter a single cutoff frequency for the filter {btype} choosed: "))
                if 0 <= W_n <= Fs / 2:
                    break
                else:
                    print("The value entered is outside the allowed range. Please try again.")
        elif btype_S == 3:
            btype = 'bandpass'
            N = 5
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                W_n = list(map(float, input(f"Enter the cutoff frequencies for the chosen {btype} filter (separated by space): ").split()))
                if all(0 <= w <= Fs / 2 for w in W_n):
                    break
                else:
                    print("At least one of the entered values is out of the allowed range. Please try again.")
        elif btype_S == 4:
            btype = 'bandstop'
            N = 5
            while True:
                print(f"Please note that the income range is in the range 0--{int(Fs/2)}")
                W_n = list(map(float, input(f"Enter the cutoff frequencies for the chosen {btype} filter (separated by space): ").split()))
                if all(0 <= w <= Fs / 2 for w in W_n):
                    break
                else:
                    print("At least one of the entered values is out of the allowed range. Please try again.")


        # Lista de filtros IIR disponibles
        senales_1 = ["Butterworth", "Chebyshov I", "Chebyshov II", "Elliptical"]
        # Mostramos el menú
        print("Types of IIR filters available: ")
        for i, x in enumerate(senales_1):
            print(f"{i + 1}) {x}")

        seleccion_2 = int(input("Enter your selection: "))

    ##################################BUTTER SET#######################################################  

        if seleccion_2 == 1:

            print("Chose Butter")

            b, a = signal.butter(N, W_n, btype, analog=False, output='ba', fs=Fs)

    #####################################CHEBYSHOV_1###################################################

        elif seleccion_2 == 2:

            print("Chose Chebyshov I")

            b, a = signal.cheby1(N, rp, W_n, btype, analog=False, output='ba', fs=Fs)

    ####################################CHEBYSHOV II###################################################
        elif seleccion_2 == 3:
            
            print("Chose Chebyshov II")

            b,a = signal.cheby2(N, rs, W_n, btype, analog=False, output='ba', fs=Fs)

    ######################################ELIPTICO#######################################################            
            
        elif seleccion_2 == 4:
            
            print("Chose Elliptical")

            b, a = signal.ellip(N, rp, rs, W_n, btype, analog=False, output='ba', fs=Fs)

    ####################################GRÁFICOS#######################################################

    if seleccion == 'FIR':
        w1, H1 = signal.freqz(h,1) #Sacar la respuesta en frecuencia
        y_filtrada1 = signal.lfilter(h, 1, y)
        
        T = 1/sr

        tam = np.size(y)
        t = np.arange(0, tam*T,T)

        angles1 = np.unwrap(np.angle(H1))

        f1_, t1_, Sxx_ = signal.spectrogram(y, sr, scaling = 'density')
        f2_, t2_, Sxx2_ = signal.spectrogram(y_filtrada1, sr, scaling = 'density')
        
        f1 = w1*sr/(2*np.pi)

        Qtn = input("Do you want to listen to the original signal and the filtered signal (Y/N): ")
        if Qtn == 'Y':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada1, sr)
            sd.wait()
        else:
            print("Understood.")

        gs1 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs2 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])
        gs3 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs4 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])

        fig = plt.figure(figsize=(12, 7))

        ax1 = fig.add_subplot(gs1[0, 0])
        ax1.plot(t, y, 'b', label = 'Original Signal',)
        ax1.plot(t, y_filtrada1, 'r', label = 'Filtered Signal')
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('Amplitude')
        ax1.set_title('Original and filtered signal')
        ax1.legend()
        ax1.grid(True)

        ax2 = fig.add_subplot(gs2[0, 1])
        ax2.pcolormesh(t1_, f1_, 10*np.log10(Sxx_), shading='gouraud')
        ax2.set_ylabel('Frequency [Hz]')
        ax2.set_xlabel('Time [s]')
        ax2.set_title('Original signal spectrogram')

        ax3 = fig.add_subplot(gs3[1, 0])
        ax3_fase = ax3.twinx()
        ax3.plot(f1, 20*np.log10(np.abs(H1)), 'b')
        ax3_fase.plot(f1, angles1*180/np.pi, 'r')
        ax3.set_ylabel('Magnitude [dB]', color = 'b')
        ax3_fase.set_ylabel('Phase [°]', color = 'r')
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_title('Frequency Response')

        ax4 = fig.add_subplot(gs4[1, 1])
        #cmap = plt.colormaps["seismic"]
        ax4.pcolormesh(t2_, f2_, 10*np.log10(Sxx2_), shading='gouraud') #, cmap=cmap)
        ax4.set_ylabel('Frequency [Hz]')
        ax4.set_xlabel('Time [s]')
        ax4.set_title('Spectrogram filtered signal')

        plt.tight_layout()
        plt.show()

    elif seleccion == 'IIR':
        w, H = signal.freqz(b, a) #Sacar la respuesta en frecuencia

        y_filtrada = signal.lfilter(b, a, y)
        T = 1/sr

        tam = np.size(y)
        t = np.arange(0, tam*T,T)

        angles = np.unwrap(np.angle(H))

        f1, t1, Sxx1 = signal.spectrogram(y, Fs, scaling = 'density')
        f2, t2, Sxx2 = signal.spectrogram(y_filtrada, Fs, scaling = 'density')

        f = w*Fs/(2*np.pi)

        Qtn = input("Do you want to listen to the original signal and the filtered signal (Y/N): ")
        if Qtn == 'Y':
            sd.play(y, sr)
            sd.wait()
            sd.play(y_filtrada, sr)
            sd.wait()
        else:
            print("Understood.")

        gs1 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs2 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])
        gs3 = GridSpec(2, 2, height_ratios=[2/3, 1], width_ratios=[1,1])
        gs4 = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1,1])

        fig = plt.figure(figsize=(12, 7))

        ax1 = fig.add_subplot(gs1[0, 0])
        ax1.plot(t, y, 'b', label = 'Original Signal')
        ax1.plot(t, y_filtrada, 'r', label = 'Filtered Signal')
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('Amplitude')
        ax1.set_title('Original and filtered signal')
        ax1.legend()
        ax1.grid(True)

        ax2 = fig.add_subplot(gs2[0, 1])
        ax2.pcolormesh(t1, f1, 10*np.log10(Sxx1), shading='gouraud')
        ax2.set_ylabel('Frequency [Hz]')
        ax2.set_xlabel('Time [s]')
        ax2.set_title('Original signal spectrogram')
        
        ax3 = fig.add_subplot(gs3[1, 0])
        ax3_fase = ax3.twinx()
        ax3.plot(f, 20*np.log10(np.abs(H)), 'b')
        ax3_fase.plot(f, angles*180/np.pi, 'r')
        ax3.set_ylabel('Magnitude [dB]', color = 'b')
        ax3_fase.set_ylabel('Phase [°]', color = 'r')
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_title('Frequency Response')

        ax4 = fig.add_subplot(gs4[1, 1])
        ax4.pcolormesh(t2, f2, 10*np.log10(Sxx2), shading='gouraud')
        ax4.set_ylabel('Frequency [Hz]')
        ax4.set_xlabel('Time [s]')
        ax4.set_title('Spectrogram filtered signal')

        plt.tight_layout()
        plt.show()
