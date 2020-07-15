#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 11:02:13 2020

@author: Loic Viens

This code reads 1-h data files (20 Hz sampling frequency) from MeSO-net stations, computes the autocorrelation function (ACF), and stack the ACFs using the Phase Weighted Stack (PWS) method.
The data have been filtered between 0.05 and 5 Hz using a 4-pole 2-pass Butterworth bandpass filter and their instrument response has been corrected.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This is just an example with 2 days of data for one station, more data are used in the manuscript to obtain stable ACFs.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""

import os
from zipfile import ZipFile
import sys
import numpy as np
from scipy.signal import butter, filtfilt, hilbert, hann
from scipy.fftpack import next_fast_len
from obspy import read as readobs
import matplotlib.pyplot as plt


#%%

def moving_ave(A,N):
    '''
    Function does running smooth average for an array.
    PARAMETERS:
    ---------------------
    A: 1-D array of data to be smoothed
    N: integer, it defines the half window length to smooth
    
    RETURNS:
    ---------------------
    B: 1-D array with smoothed data
    '''
    A = np.concatenate((A[:N],A,A[-N:]),axis=0)
    B = np.zeros(A.shape,A.dtype)
    
    tmp=0.
    for pos in range(N,A.size-N):
        # Sum only once
        if pos==N:
            for i in range(-N,N+1):
                tmp+=A[pos+i]
        else:
            tmp=tmp-A[pos-N-1]+A[pos+N]
        B[pos]=tmp/(2*N+1)
        if B[pos]==0:
            B[pos]=1
    return B[N:-N]

def pws(arr, sampling_freq, power=2, pws_timegate=.1):
    '''
    Phase-weighted stack of an array of ACFs.
    Follows methods of Schimmel and Paulssen, 1997. 
    If s(t) is time series data (seismogram, or cross-correlation),
    S(t) = s(t) + i*H(s(t)), where H(s(t)) is Hilbert transform of s(t)
    S(t) = s(t) + i*H(s(t)) = A(t)*exp(i*phi(t)), where
    A(t) is envelope of s(t) and phi(t) is phase of s(t)
    Phase-weighted stack, g(t), is then:
    g(t) = 1/N sum j = 1:N s_j(t) * | 1/N sum k = 1:N exp[i * phi_k(t)]|^v
    where N is number of traces used, v is sharpness of phase-weighted stack
    
    PARAMETERS:
    ---------------------
    arr: Array with the time series data (numpy.ndarray)
    sampling_freq: sampling frequency of the time series in Hz (int)
    power: exponent of the phase weighted stack (int)
    pws_timegate: Smoothing of the  phase weighted stack in seconds (float)
    
    RETURNS:
    ---------------------
    weighted: Phase weighted stack of the data (numpy.ndarray)
    '''
    if arr.ndim == 1:
        return arr
    N,M = arr.shape
    analytic = hilbert(arr,axis=1, N=next_fast_len(M))[:,:M]
    phase = np.angle(analytic)
    phase_stack = np.mean(np.exp(1j*phase),axis=0)
    phase_stack = np.abs(phase_stack)**(power)

    # smoothing 
    timegate_samples = int(pws_timegate * sampling_freq)
    phase_stack = moving_ave(phase_stack,timegate_samples)
    weighted = np.multiply(arr,phase_stack)
    return np.mean(weighted,axis=0)



def butter_bandpass_filter(data, lowcut, highcut, fs, order =  4):
    '''
    Bandpass filter the waveforms with a Butterworth bandpass filter
    
    PARAMETERS:
    ---------------------
    data: time series to filter (numpy.ndarray)
    lowcut: Lower cutoff frequency (in Hz)
    highcut: Upper cutoff frequency (in Hz)
    fs: Sampling frequency of the data (in Hz)
    order: order of the filter
    RETURNS:
    ---------------------
    y: Bandpass filtered waveform
    
    '''
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y


def ACF_func(data, dt, timesave):
    '''
    Compute autocorrelation functions in the frequency domain
    
    PARAMETERS:
    ---------------------
    dat: time series (numpy.ndarray)
    dt: sampling rate of time series in second (float)
    timesave: lag time to save for the ACF
    RETURNS:
    ---------------------
    corr: autocorrelation function of time series data in the time domain (numpy.ndarray)
    '''
    n = len(data)  # define the length of the data
    delta = 1/dt # Sampling frequency of the data (in Hz)
    fft_r = np.fft.fft(data, n*4) # FFT with zero-padding
    cc_t1 = np.real(np.fft.ifft( (fft_r * np.conj(fft_r)) )) # Compute Autocorrelation
    corr2 = np.concatenate((cc_t1[int(len(cc_t1)/2):], cc_t1[:int(len(cc_t1)/2)])) # Shift the causal and anticausal ACFs to their right positions (does not really matter for ACFs as they are symmetric)
    corr  = corr2[int(len(corr2)/2-timesave*delta):int(len(corr2)/2+timesave*delta)] # Only select a few seconds (set by the timesave variable) for both the anticausal and causal parts
    return corr

#%%
dir_ini = './'
sys.path.append(dir_ini)
directory ='../Data/Raw/' # data directory (The full folder + file name is defined Line 177 )
if not os.path.exists(directory):
    zf = ZipFile('../Data/Raw.zip', 'r')
    zf.extractall('../Data/')
    zf.close()


sta = ['E.AYHM'] # Station name

year = '2019'
dy = ["%.3d" % i for i in range(182,184)] # Compute 20-min ACFs over 2 days of data (July 1st to 2nd, 2019)

cmpo1 = 'U' # Vertical component
timesave = 20 # save ACFs from -timesave to +timesave (in seconds)
len_min = 20 # Duration of the raw record to consider (in minutes)
len_time_win = (3600) / (len_min*60) # To split each 1-h time series into three 20-min time series

# Filter  parameters
filtlow = 0.1 #  Lower cutoff frequency (in Hz)
filthigh = 1 # Upper cutoff frequency (in Hz)

#%% Parameter of the PWS
pwer = 2 # Power of the PWS
pws_timegate = .1 # Time gate to perfom some smoothing (in seconds)

#%% Taper parameter for the ACFs before filtering 
hann2 = hann(41) # 40 sample -> 2 s given the 20 Hz sampling frequency
half_len =round(len(hann2)/2-.001)+1 # half taper: 1 s

#%% Loop to read all the 1-h MeSO-net data (20 Hz sampling frequency) and compute ACFs
for rec in sta: # Loop over the stations
    ZZcorr = []
    for day in dy:
        print(rec,day)
        for hr1 in np.arange(24):
            hr = str(hr1).zfill(2)
            name_rec = directory + rec + os.sep + year + os.sep + day + os.sep + rec + '_' + cmpo1 + '_' + year + '_' + day + '_' + str(hr) + '.sac'
            dat_r1 = readobs(str(name_rec)) # Read sac data with obspy
            dat_r1.detrend(type='constant') # Remove mean
            dat_r1.detrend(type='linear') # Remove trend
            dt = round(dat_r1[0].stats.delta,3) # Get dt sampling rate in seconds
            dat_r = dat_r1[0].data # get data
            
            if len(dat_r)==round(1/dt*3600): # Check if there is no problem with the 1-H time series
                a = np.split(dat_r1[0].data, len_time_win) # Split the 1-h time series into a 20-min time series
                for ab in a: # Loop over each of the 20-min time series
                    ab = ab - np.mean(ab) # Remove mean
                    corr = ACF_func(ab, dt, timesave) # Compute ACFs
                    ZZcorr.append(corr) # Save the ACF in a list

 #%% Compute PWS for the ACFs in ZZcorr
    delta = 1/dt # Sampling frequency from sampling rate
    ZZcorr = np.array(ZZcorr) # List to NP array
    mxv = np.empty(len(ZZcorr)) # Initiate array for maximum values
    for i in range(ZZcorr.shape[0]):
        ZZcorr[i] = ZZcorr[i] - np.mean(ZZcorr[i]) # Remove mean of the ACF
        ZZcorr[i,:half_len] = ZZcorr[i,:half_len] * hann2[:half_len] # Taper end of anticausal ACFs
        ZZcorr[i,-half_len:] = ZZcorr[i,-half_len:] * hann2[-half_len:] # Taper end of causal ACFs
        ZZcorr[i] = butter_bandpass_filter(ZZcorr[i], lowcut= filtlow, highcut = filthigh, fs = delta, order = 4) # Filter
        ZZcorr[i] = ZZcorr[i] - np.mean(ZZcorr[i])  # Remove mean of the ACF
        mxv[i] = np.max(ZZcorr[i]) # Compute the maximum value of each ACF

    valup = np.mean(mxv) + np.std(mxv) # Compute the metric used to remove ACFs with strong amplitudes
    data_F = ZZcorr[mxv<valup] # Select the data with maximum amplitudes lower than the metric
    ZZcorrf = pws(arr = data_F, sampling_freq = delta, power = pwer, pws_timegate = pws_timegate ) # Compute the PWS


#%% Plot the data
    t = np.linspace(-len(ZZcorrf)/2*1/delta,len(ZZcorrf)/2*1/delta, len(ZZcorrf)) # Time vector for the ACFs
    
    plt.figure(figsize =(6,7))
    plt.subplot(211)
    plt.imshow(data_F, aspect = 'auto', extent = (t[0], t[-1], data_F.shape[0] ,0) )
    plt.clim((-np.max(data_F)/100, np.max(data_F)/100)) # Amplitude clipping
    plt.xlabel('Time (s)')
    plt.ylabel('ACF #')
    plt.title(rec[2:] +' station\n20-min vertical ACFs (Filter: ' + str(filtlow) + '-' + str(filthigh) + ' Hz)'  )
    plt.grid(linewidth = 0.5)
    
    plt.subplot(212)
    plt.plot(t, ZZcorrf, linewidth =2)
    plt.xlim(-timesave,timesave)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude')
    plt.title( 'Stacked ACF')
    plt.grid()
    
    plt.tight_layout()
    plt.show()
