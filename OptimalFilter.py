import numpy as np
import matplotlib.pyplot as plt
import scipy

def GetExp(x,a,b):
    return a*np.exp(-b*x)

def OptimalFilter(Time, Noise, Signal, Data):
    Time = np.array(Time)
    SamplePoints = Time.size
    SampleSpace = 0.5/1E6
    TimeFFT = np.fft.fftfreq(SamplePoints, SampleSpace)[0:SamplePoints//2]
    print SamplePoints, SampleSpace, SamplePoints*SampleSpace

    NoiseFFT = GetFFT(Noise, SamplePoints)
    PowerSpectrumDensity = GetPowerSpectrumDensity(Noise, 1/SampleSpace,SamplePoints)

    SignalFFT = GetFFT(Signal, SamplePoints)
    SignalFFTConj = np.conj(SignalFFT)

    DataFFT = GetFFT(Data, SamplePoints)
    AmplitudeEstimator = GetAmplitudeEstimator(Data, DataFFT, SignalFFT, SignalFFTConj, PowerSpectrumDensity)

    return AmplitudeEstimator

def GetPowerSpectrumDensity(Noise, Sampling, SamplePoints):
    NoiseFFT = np.fft.rfft(Noise)
    PowerSpectrumDensity = np.abs(NoiseFFT)**2
    PowerSpectrumDensity = np.mean(PowerSpectrumDensity)
    return PowerSpectrumDensity

def GetNoiseWaveforms(Data):
    Template = np.mean(Data, axis=0)
    Noise = np.array([x-Template/(np.max(Template)/np.max(x)) for x in Data])
    return Noise

def GetFunc(x, A, Template):
    func = A * Template
    return func

def GetAmplitudeEstimator(Data, DataFFT, TemplateFFT, TemplateFFTConj, PowerSpectrumDensity):
    Amplitudes = []
    NrPoints = np.shape(Data)[0]/2
    print NrPoints
    j = complex(0, 1)
    best = []
    omega = np.linspace(-2*np.pi/NrPoints, 2*np.pi/NrPoints, NrPoints)
    datafft = DataFFT
    new = []
    print(Data, DataFFT, TemplateFFT, TemplateFFTConj, PowerSpectrumDensity)
    for t0 in range(0,1):
        Numerator = np.sum( np.exp(j*omega*t0) * TemplateFFTConj[1:] * datafft[1:] / PowerSpectrumDensity ).real
        Denominator = np.sum((np.abs(TemplateFFT[1:])**2)/PowerSpectrumDensity)
        AmplitudeEstimator = Numerator/Denominator
        new.append(AmplitudeEstimator)

    return AmplitudeEstimator

def GetFFT(Data, SamplePoints):
    DataFFT = np.fft.rfft(Data)
    return np.array(DataFFT)

if __name__ == '__main__':
    sig = 10*np.sin(np.arange(180)* np.pi / 180. )
    noise = np.random.random(180)
    datawf = sig + noise
    time = np.arange(180)
    print OptimalFilter(time, noise, sig/10, datawf)
