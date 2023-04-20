import numpy as np
import scipy.signal as signal
from scipy.fft import fft
from scipy import stats
import pandas as pd


def openibis(eeg):
    '''
    This is the main function of the algorithm. 
    It takes in the EEG signal and returns the depth of anesthesia.
    Fs: EEG sampling rate must be 128 Hz
    stride: the stride of the sliding window in seconds
    BSR: Burst Suppression Ratio
    BSRmap: Burst Suppression Ratio map to identify segments
    components: Three components of the depth of anesthesia
    '''
    Fs, stride = 128.00, 0.50 
    BSRmap, BSR = suppression(eeg, Fs, stride)
    components = logPowerRatios(eeg, Fs, stride, BSRmap)
    depthOfAnesthesia = mixer(components, BSR)
    return depthOfAnesthesia


def suppression(eeg, Fs, stride):
    '''
    This function calculates the Burst Suppression Ratio (BSR) and the BSR map.
    '''
    N, nStride = nEpochs(eeg, Fs, stride)
    BSRmap = np.zeros((int(N), 1))
    for n in range(N):
        x = segment(eeg, n + 6.5, 2, nStride)
        BSRmap[n] = np.all(np.abs(x-baseline(x)) <= 5)
    BSR = 100 * np.convolve(BSRmap.flatten()[::-1], np.ones((int(63/stride),)), mode='full')[::-1][(int(63/stride)-1):-(int(63/stride)-1)].reshape(-1, 1)
    return BSRmap, BSR



def logPowerRatios(eeg, Fs, stride, BSRmap):
    N, nStride = nEpochs(eeg, Fs, stride)
    B, A = signal.butter(2, 0.65 / (Fs / 2), 'high')
    eegHiPassFiltered = signal.lfilter(B, A, eeg)
    psd = np.empty((N, 4 * int(nStride // 2)))
    assert nStride // 2 == int(nStride // 2)
    psd[:] = np.nan
    suppressionFilter = np.piecewise(np.arange(0, 64, 0.5), [np.arange(0, 64, 0.5) < 3, np.arange(0, 64, 0.5) < 6], [0, 0.25, 1]) ** 2
    components = np.empty((N, 3))
    components[:] = np.nan

    for n in range(N):
        if isNotBurstSuppressed(BSRmap, n, 4):
            psd[n, :] = powerSpectralDensity(segment(eegHiPassFiltered, n + 4, 4, nStride))
            if not(sawtoothDetector(segment(eeg, n + 4, 4, nStride), nStride)[0]):
                psd[n, :] = suppressionFilter * psd[n, :]

        thirtySec = timeRange(30, n, stride)
        vhighPowerConc = np.sqrt(np.mean(psd[thirtySec, :][:, bandRange(39.5, 46.5, 0.5)] * psd[thirtySec, :][:, bandRange(40, 47, 0.5)], axis=1))
        wholePowerConc = np.sqrt(np.mean(psd[thirtySec, :][:, bandRange(0.5, 46.5, 0.5)] * psd[thirtySec, :][:, bandRange(1, 47, 0.5)], axis=1))
        midBandPower = prctmean(np.nanmean(10 * np.log10(psd[thirtySec, :][:, bandRange(11, 20, 0.5)]), axis=1), 50, 100)
        
        components[n, 0] = meanBandPower(psd[thirtySec, :], 30, 47, 0.5) - midBandPower
        components[n, 1] = stats.trim_mean( 10 * np.log10( vhighPowerConc / wholePowerConc ), 0.5)
        components[n, 2] = meanBandPower(psd[thirtySec, :], 0.5, 4, 0.5) - midBandPower

    return components


def powerSpectralDensity(x): 

    """
    This function computes the power spectral density of a signal 
    Implementing the method: Blackman-windowed Power Spectral Density
    
    x: input signal
    output: power spectral density of the signal
    
    """

    blackman_signal = signal.windows.blackman(len(x))
    baseline_signal = (x - baseline(x))

    multipl_signals = blackman_signal * baseline_signal

    f = np.fft.fft( multipl_signals) # Perform Blackman windowing FFT, removing the baseline drift

    y = 2 * np.abs(f[0:len(x)//2])**2 // (len(x) * np.sum(signal.windows.blackman(len(x))**2))

    return np.array(y)


def sawtoothDetector(eeg, n_stride):
    """

    This function detects sawtooth waves in an EEG signal

    """
    assert int(n_stride) == n_stride
    n_stride = int(n_stride)
    
    
    saw = np.hstack((np.zeros(n_stride - 5), np.arange(1, 6)))

    saw = (saw - np.mean(saw)) / np.std(saw, ddof=1)

    r = np.arange(len(eeg) - len(saw))
    v = np.apply_along_axis(lambda x: np.var(x), 0, np.lib.stride_tricks.as_strided(eeg, (len(eeg) - len(saw) + 1, len(saw)))) # moving window variance
    m = (np.array([np.convolve(eeg, np.flipud(saw), mode='valid'), np.convolve(eeg, saw, mode='valid')]) / len(saw))**2
    m = m.T # Transpose 
    
    
    y = np.max([(v > 10) * m[:len(v), 0] / v, (v > 10) * m[:len(v), 1] / v], axis=0) > 0.63

    return np.array(y)



def mixer(components, BSR): 

    """
    Mixer: Generate the output of the depth of anesthesia model by converting and weighting the components.

    :param components: number of seconds to cover
    :param BSR: Burst suppression ratio


    :varr sedationScore: Map the 1st component to the sedation score on a logistic curve
    :varr generalScore:  Map the 2nd component to the general score on a linear and logistic curve
    :varr bsrScore: Convert the BSR to a score on a linear function
    :varr generalWeight: Convert the 3rd component to a weight

    :output y: Calculate the depth of anesthesia score by a compressed and weighted sum of the components with BSR
    """

    sedationScore = scurve(components[:, 0], 104.4, 49.4, -13.9, 5.29)
    
    # Applies the piecewise function to each element of the components[:, 1] array. 
    # If an element is less than -30, the value in generalScore will be -40. 
    # Otherwise, the value will be 42 (since 43-1 = 42).
    generalScore = np.piecewise(components[:, 1], 
                            [components[:, 1] < -60.89, components[:, 1] >= -30],
                            [-40, lambda x: 43-1])  
    
    generalScore += scurve(components[:, 1], 61.3, 72.6, -24.0, 3.55) * (components[:, 1] >= -30)
    
    #Origianl: bsrScore = picewise(BSR, [0, 100], [50, 0])
    bsrScore = np.piecewise(BSR, [BSR<0, BSR>=100], [lambda x: 50, lambda x: 50 - (x / 2)]) #  same shape as BSR, with values ranging from 0 to 50.

    generalWeight = np.piecewise( components[:, 2], 
                                 [components[:, 2] < 0, components[:, 2]>= 5], 
                                 [0.5, 1] ) * (generalScore < sedationScore) 
    
    bsrWeight = np.piecewise(BSR,
                         [BSR < 10, BSR >= 50],
                         [0, 1])
    
    x = sedationScore * (1 - generalWeight) + generalScore * generalWeight
    
    # The first argument of np.piecewise() is the input array x. 
    # The second argument is a list of conditions that are used to define the regions over which the function is defined. 
    # In this case, we have 5 regions defined by the conditions: x < -40, -40 <= x < 10, 10 <= x < 97, 97 <= x < 110, and x >= 110.
    # The third argument is a list of functions that correspond to the regions defined by the conditions. 
    # For example, if x is less than -40, the function value will be 0. 
    # If x is between -40 and 10, the function value will also be 0. 
    # If x is between 10 and 97, the function value will be 10, and so on.
    y = np.piecewise(x, 
                      [x < -40, (x >= -40) & (x < 10), (x >= 10) & (x < 97), (x >= 97) & (x < 110), x >= 110], 
                      [0, 0, 10, 97, 100]) * (1 - bsrWeight) + bsrScore * bsrWeight
    
    return y


#### Helper functions

def nEpochs(eeg, Fs, stride):
    '''
    Calculates the number of epochs and the stride length based on the sampling frequency and a desired stride
    '''
    nStride = Fs * stride
    N = (len(eeg) - Fs) // nStride - 10

    return int(N), nStride

# meanBandPower: calculates the mean power in a frequency band of a power spectral density
def meanBandPower(psd, from_, to, bins):

    v = psd[:, bandRange(from_, to, bins)]

    return np.mean(10 * np.log10 (v[~np.isnan(v)]),  dtype=int)

def bandRange(from_, to, bins):
    """
    Generatesthe indices of bins for a given frequency band
    """

    return np.arange(from_ / bins, (to / bins) + 1, dtype=int)

def baseline(x):
    """
    calculates the baseline of a signal
    """

    v = np.hstack((np.arange(len(x))[:, np.newaxis], np.ones((len(x), 1))))

    return np.linalg.lstsq(v, x, rcond=None)[0][1]


def bound(x, lower_bound, upper_bound):
    """
    Constrains a value between lowerBound and upperBound
    """

    return np.minimum(np.maximum(x, lower_bound), upper_bound)


def segment(eeg, from_, number, nStride):
    """
    Extracts segments of the EEG signal
    
    :eeg: signal
    :from_: start time
    :number: duration
    :nStride: sliding window

    :returns: segment of the EEG signal

    """

    array_indices = int(from_ * nStride) + np.arange(int(number * nStride)) 
    signal_array = [eeg[i] for i in array_indices]

    return np.array(signal_array)



def isNotBurstSuppressed(BSRmap, n, p):
    '''
    isNotBurstSuppressed: checks if a given index is NOT in a burst-suppressed region
    n: number of epochs
    p: 4; fixed parameter   
    '''
    cond1 = n < p #

    list_index = n + np.arange(-(p-1), 1)

    cond2 = any(BSRmap[list_index] )

    return not(cond1 or cond2 )


def timeRange(seconds, n, stride):
    """
    timeRange: Generate indices for the most recent timepoints of the signal
    :param seconds: number of seconds to cover
    :param n: number of epochs
    :param stride: the stride of the signal
    :return: a numpy array containing the indices that cover the specified number of seconds
    """

    start = max(0, int((n -(seconds / stride))) +1)
    end = n + 1
    return np.arange(start, end)

# prctmean: calculates the mean of the values between the lower and upper percentiles of a vector
def prctmean(x, lo, hi):

    v = np.percentile(x, [lo, hi])

    return np.mean(x[(x>=v[0]) & (x<=v[1])])

# scurve: calculates an S-shaped curve
def scurve(x, Eo, Emax, x50, xwidth): 

    return Eo - Emax / (1 + np.exp((x - x50) / xwidth))


# HERE IS NOT USED 
def piecewise(x, xp, yp):

    '''
    Matlab code adaptation
    picewise : performs piecewise linear interpolation
    '''

    return np.interp(np.clip(x, xp[0], xp[-1]), xp, yp)


if __name__ == '__main__':
    import os; os.system('clear')
    import mat73
    eeg = mat73.loadmat('data/case1.mat')['EEG']
    
    print(eeg.shape, type(eeg))
    
    openibis(eeg)

