import librosa
from IPython.display import Audio
import sounddevice as sd

# Preguntar al usuario si quiere que la interacción sea en inglés o en español
idioma = input("¿Quieres que la interacción sea en inglés o en español? (inglés/español): ")

# Permitir al usuario cargar un archivo de audio cualquiera de formato .wav
ruta_audio = input("Ingresa la ruta del archivo de audio: ")

# Procesar el archivo de audio
y, sr = librosa.load(ruta_audio)

# Reproducir el audio
sd.play(y, sr)
sd.wait()

# Convertir a mono si es estéreo
if len(y.shape) == 2:
    channel1 = y[:, 0]
    channel2 = y[:, 1]
    y = (channel1 + channel2) / 2

    # Opcional: Guardar el audio mono
    guardar = input("¿Quieres guardar el audio mono? (s/n): ").lower()
    if guardar == 's':
        nueva_ruta = input("Ingresa la ruta para guardar el audio mono (.wav): ")
        librosa.output.write_wav(nueva_ruta, y, sr)

# Visualizar el audio
Audio(data=y, rate=sr)
