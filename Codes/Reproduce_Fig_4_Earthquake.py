#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 10:19:11 2021

@author: loic
"""

import os
import geopy
import numpy as np
import geopy.distance
import scipy.io as sio
from obspy import read
from scipy import signal
from zipfile import ZipFile
import matplotlib.pyplot as plt
from geopy.distance import  geodesic
from scipy.interpolate import griddata
from obspy.signal.filter import bandpass
from scipy.signal import butter, filtfilt
from scipy.fftpack import fft, ifft, next_fast_len


#%%
def butter_bandpass_filter(data, lowcut, highcut, fs, order = 4):
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

def moving_ave(A,N):
    '''
    this Numba compiled function does running smooth average for an array.
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
        # do summing only once
        if pos==N:
            for i in range(-N,N+1):
                tmp+=A[pos+i]
        else:
            tmp=tmp-A[pos-N-1]+A[pos+N]
        B[pos]=tmp/(2*N+1)
        if B[pos]==0:
            B[pos]=1
    return B[N:-N]

def cross_cor_Chengxin(dat_r, delta, smoothN, flilow, flihigh, ctaper, timesave):
    n = len(dat_r) 
    Nfft = int(next_fast_len(4*n))
    spec = fft(dat_r, Nfft) # return FFT
    temp  = moving_ave(np.abs(spec),smoothN)
    sfft1 = spec*np.conj(spec)/(temp*temp)
    tdata = np.real(np.fft.ifftshift(ifft(sfft1, Nfft, axis=0)))
    tdata2 = tdata[int(len(tdata)/2):int(len(tdata)/2+timesave*delta)]
    tdata2[:len(ctaper)] *=ctaper
    tdata2[-len(ctaper):] *=ctaper[::-1]
    corr = bandpass(tdata2,freqmin=flilow,freqmax= flihigh, df=delta,corners=4, zerophase=True)
    return corr

#%% Directory with all the data and metadata
#data_dir = '/Users/loic/Documents/Kanto_auto_corr/Data_Final/'
data_dir = '../Data/'
#%% Extract noise and earthquake ACF from zip file if needed
directory = '../Data/ACFs/'
if not os.path.exists(directory):
    zf = ZipFile('../Data/ACFs.zip', 'r')
    zf.extractall('../Data/')
    zf.close()
    
#%% Extract noise and earthquake ACF from zip file if needed
directoryAxi = '../Data/Axitra/'
if not os.path.exists(directoryAxi):
    zf = ZipFile('../Data/Axitra.zip', 'r')
    zf.extractall('../Data/')
    zf.close()
#%% Main parameter to plot the data along different Lines
Line_nb = '1' # Possible values are 1, 2, 3(str)

#%% Some more parameters
period1 = 1 # filter parameter for the Noise ACFs -> can be changed but will only impact Noise ACF results
period2 = 10 # filter parameter for the Noise ACFs-> can be changed but will only impact Noise ACF results
order = 4 # order of the filter -> can be changed but will only impact Noise ACF results
cut = 1/period2
cut2= 1/period1


#%% Load metadata for MeSOnet stations (station name, latitude, and longitude)
crs = open(data_dir + "Other_data/MeSOnet_channel_file", "r")
nm_sta = []
lat = []
lon = []
for columns in ( raw.strip().split() for raw in crs ):
     if (columns[4]) == 'U':
         nm_sta.append(columns[3])
         lat.append(columns[13])
         lon.append(columns[14])
crs.close()

#%% Load JIVSM
crs = open(data_dir +  "Other_data/Kanto_basin_shape.txt", "r")
topo = []
latg = []
long =[]
lay1 =[]
lay2 =[]
lay3 =[]
for columns in ( raw.strip().split() for raw in crs ):  
    long.append(columns[0])
    latg.append(columns[1])
    topo.append(columns[2])
    lay2.append(columns[8])
    lay1.append(columns[5])
    lay3.append(columns[15])
crs.close()   
    

#%%
# Parameters for each line
if Line_nb == '1':
    station = ['E.THCM','E.DICM','E.TYHM','E.TBKM','E.SSHM','E.NKMM','E.TK2M','E.HSDM','E.SRTM','E.NKNM','E.TSCM','E.KUDM', 'E.HRGM', 'E.RKGM' ,'E.MNKM','E.FUNM', 'E.KMRM', 'E.NNTM', 'E.STHM','E.TYNM','E.KSCM','E.MDHM', 'E.HKBM','E.TKNM','E.SDMM', 'E.YNMM', 'E.NSMM', 'E.HGCM','E.HYHM','E.SRCM','E.DAIM','E.OKCM', 'E.FMHM','E.FMNM', 'E.NISM', 'E.FKCM', 'E.OKDM', 'E.KSRM', 'E.IIDM', 'E.SYPM' , 'E.YSSM', 'E.SBCM', 'E.TACM'] 
    leftt = 'E'
    rightt = 'W'
    az_mean = 100
    mxdist = 128
    
elif Line_nb == '3':
    station = ['E.SYOM','E.FJSM','E.DGRM','E.KMHM','E.SMGM','E.HTTM','E.KSOM','OK.HRDM','E.HNOM', 'E.BKKM','E.MZPM','E.JDJM', 'E.SNHM', 'E.UNMM','E.KMKM','E.MKSM','E.GKSM','E.KDKM','E.ENZM', 'E.TKMM','E.SBAM','E.GNZM','E.YKKM','E.MKJM', 'E.KCBM','E.MZMM', 'E.TKNM','E.MBSM','E.YKSM','E.KGKM', 'E.NGSM','E.KKSM', 'E.KWHM', 'E.TNKM','E.KUYM', 'E.IN3M','E.INAM', 'E.KBRM','E.GSJM', 'E.YTBM', 'E.HNPM','E.YKBM','E.TSRM', 'E.TKZM','E.RMSM','E.RYGM', 'E.RYNM']
    az_mean = 34.2
    mxdist = 113

    
indv = nm_sta.index(station[0])
sourlat = float(lat[indv])
sourlon = float(lon[indv])
#%%
deltaeq = 20 # Sampling frequency of the data (in Hz)
lag_selct = 15*deltaeq # select 15 s of the positive lag of the ACFs
teq = np.linspace(0,lag_selct/deltaeq,lag_selct) # Create a time vector
#%%
datsaviEQ = np.empty((len(station),len(teq)) )
stalatsave = []
stalonsave = []
dist4 = []
for i in np.arange(len(station)):
    indr = nm_sta.index(station[i])
    stalat = float(lat[indr])
    stalon = float(lon[indr])
    stalatsave.append( float(lat[indr]) )
    stalonsave.append( float(lon[indr]) )
    dist4.append(geopy.distance.distance((sourlat, sourlon), (stalat, stalon)).km)

    fileEQ = data_dir +'/ACFs/EQ_ACFs/' + station[i]+ '.SAC'
    st = read(fileEQ) # Read sac file
    dat1 = st[0].data # get data
    dat1 -= np.mean(dat1) # Remove the mean
    datsaviEQ[i,:] = dat1 /np.max(np.abs(dat1)) # normalize the trace
datsavEQ = np.array(datsaviEQ)  # To array


#%% Assign Noise and EQ ACFs to bins for the plot
minkm = 5 # add 5 km to the begining of the line for visibility
nnb = 1 # width of the bin in km -> (can be changed: but might impact the final results)
dist_plot = np.linspace(-minkm,mxdist-minkm, int(mxdist/nnb)+1 ) # Create an array of distances along the line
dat_fin2EQ = np.zeros( (len(dist_plot), datsavEQ.shape[1]) ) # set up array for EQ ACFs
sta_name = list(str(' ') * len(dist_plot)) # np.empty( (len(dist_plot), datsav.shape[1]) ) 
sta_depth =np.empty(len(dist_plot))
sta_depth[:] = np.nan
for i in np.arange(len(dist_plot)):
    nbdat = 0
    for ii in np.arange(len(station)):
        if dist4[ii] >= dist_plot[i] and dist4[ii] < dist_plot[i] + nnb:
            dat_fin2EQ[i] = datsavEQ[ii] # Same fo the EQ ACF
            sta_name[i] = station[ii]
            nbdat += 1 # To check if there is more than one station in the bin
    if nbdat>1:
        print('This bin contains ' + str(nbdat) + ' stations')
    
#%% Get layers from the JIVSM along the line
dist_t = np.linspace(0, mxdist, int((mxdist)/nnb)+1 )
destination_ini = geodesic(kilometers = minkm).destination(geopy.Point(sourlat,sourlon), az_mean+180)
lat2, lon2 = destination_ini.latitude, destination_ini.longitude
latout = []
lonout = []
for dis in dist_t:
    destination2 = geodesic(kilometers=dis).destination(geopy.Point(lat2, lon2), az_mean)
    latouti,lonouti = destination2.latitude, destination2.longitude
    latout.append(latouti)
    lonout.append(lonouti)

dist_t = dist_t-minkm
grid_z0 = griddata((latg,long),topo, (latout,lonout), method='cubic')/1000
grid_z1 = griddata((latg,long),lay1, (latout,lonout), method='cubic')/1000
grid_z2 = griddata((latg,long),lay2, (latout,lonout), method='cubic')/1000
grid_z3 = griddata((latg,long),lay3, (latout,lonout), method='cubic')/1000


#%%
unit = 3
utext = 'acceleration'
t0 = .5
fname = data_dir + '/Axitra/Axitra_Earthquake_Line_' + Line_nb +'_stf_' + str(t0) + '_type_' + str(unit) +'.mat'
isfile = os.path.isfile(fname) 
dat = sio.loadmat(fname)
autosave = np.squeeze(dat['data'])
t = np.squeeze(dat['t'])


#%% Compute simulated earthquake ACFs
delta = 1/(t[10]-t[9])
smoothN = 30
flilow = .1
flihigh = 1
hann2 = signal.hann(int(delta*.5)*2+1)
half_len =round(len(hann2)/2-.001)+1
ctaper = hann2[:half_len]

ACFini = []
ACFiniall = []
for i in range(len(autosave)):    
    datACF = autosave[i] 
    corr = cross_cor_Chengxin(datACF, delta, smoothN, flilow, flihigh, ctaper, timesave = 15)
    corrst = corr.tolist()
    ACFini.append(corrst)
dfin = np.array(ACFini)
 

stalim = [5, -4]
T_P = np.zeros(len(grid_z0))
T_S = np.zeros(len(grid_z0))
for i in np.arange(len(grid_z0)):
    T_P[i] = (0 - grid_z1[i])/1.8 + (grid_z1[i] - grid_z2[i])/2.3  +(grid_z2[i] - grid_z3[i])/3.0
    T_S[i] = (0 - grid_z1[i])/.5 + (grid_z1[i] - grid_z2[i])/.9  +(grid_z2[i] - grid_z3[i])/1.5 
T_P =T_P*2 # multiply by two to get the two-way travel time


   
#%% Plot results
fig = plt.figure(figsize = (10,8) )
fnt  = 13

axs1 = plt.subplot(211)
plt.plot(dist_t, grid_z0, 'k') # Surface
plt.plot(dist_t, grid_z1, 'grey') # Layer 1
plt.plot(dist_t, grid_z2, 'grey') # Layer 2
plt.plot(dist_t, grid_z3, 'k' , linewidth = 4) #  Bedrock
plt.plot(dist_t, grid_z0+ .07, 'v', markersize = 5, mec = 'k') # Plot station location
plt.xticks(fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.title('Line '+Line_nb +  ': JIVSM P-wave velocity profile', fontsize=fnt +4)
if  Line_nb== '1':
    plt.text(58,-0.32, '1.8')
    plt.text(58,-.9, '2.3')
    plt.text(58,-1.5, '3.0')
    plt.text(49.25,-2.55, '$\mathrm{V_p}$ = 5.5 km/s')
    axs1.annotate("", xy=(-1.5, -2.6), xytext=(10.5, -2.6), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text( 2.5, -2.5, 'W')
    axs1.annotate("", xy=(119, -2.6), xytext=(107, -2.6), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(111, -2.5, 'E')
    plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)
    plt.ylim(-3,.3)
    plt.text(-14, .3, '(a)', fontsize = 18)
elif  Line_nb== '3':
        plt.text(41, -0.34, '1.8', fontsize =fnt)
        plt.text(41, -1, '2.3', fontsize =fnt)
        plt.text(41, -2.4, '3.0', fontsize =fnt)
        plt.text(41-7, -3.7, '$\mathrm{V_p}$ = 5.5 km/s', fontsize =fnt)
        axs1.annotate("", xy=(-2, -3.4), xytext=(10, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
        plt.text(1, -3.9, 'SW')
        axs1.annotate("", xy=(105, -3.4), xytext=(93, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
        plt.text(97, -3.9, 'NE')
        plt.ylim(-4,0.3)
        plt.text(-12, .3, '(d)', fontsize = 18)
   
plt.ylabel('Depth (km)',fontsize=fnt)
plt.xlabel('Distance along Line '  +Line_nb +' (km)',fontsize=fnt )
plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)
plt.xlim(dist_t[3],dist_t[-3])
plt.grid()


axs2 = plt.subplot(212)
clim = np.max(ACFini)/10
factval=1/np.max(ACFini)*2
tpl = np.linspace(0,dfin.shape[1]/delta,dfin.shape[1]) # Create a time vector
for i in np.arange(stalim[0],len(dfin)+stalim[-1]):
    y = dist_t[i]+ dfin[i]*factval 
    plt.plot(y ,tpl, 'k', linewidth = .5)
    axs2.fill_betweenx(tpl , dist_t[i], y, where=(y < dist_t[i]), color='gray', alpha=0.8)
    axs2.fill_betweenx(tpl , dist_t[i], y, where=(y > dist_t[i]), color='k', alpha=0.8)
plt.xlim(dist_t[3],dist_t[-3])
plt.xlabel('Distance along Line '  + Line_nb +' (km)',fontsize=fnt)
plt.xticks(fontsize=fnt)
plt.yticks(fontsize=fnt)
for i in np.arange(len(dat_fin2EQ)):
    factampeq = 2
    if np.mean(dat_fin2EQ[i,:50])==0:
        dat_fin2EQ[i] = np.nan      
    y = dist_t[i]+ dat_fin2EQ[i]*factampeq
    if not np.isnan(np.sum(y)):
        plt.plot(y ,teq,'k', linewidth =.5)
        axs2.fill_betweenx(teq , dist_t[i], y , where=(y< dist_t[i]) ,color='green')
        axs2.fill_betweenx(teq , dist_t[i], y , where=(y> dist_t[i]) ,color='blue')
  
if Line_nb =='1':
    plt.text(124, .4, '2p',  weight="bold", color = 'dodgerblue', fontsize = fnt)
    plt.text(124, .9, '$\mathbf{2p^2}$',  weight = "bold", color = 'mediumturquoise', fontsize = fnt)
    plt.text(-14, -.5, '(b)', fontsize = 18)
    plt.plot(dist_t,T_P,'--', color =  'dodgerblue', linewidth = 2)
    plt.plot(dist_t,T_P*2,'--', color = 'mediumturquoise', linewidth = 2)
    plt.ylim(5,0)  
elif Line_nb =='3':
    plt.text(109, .3, '2p',  weight="bold", color = 'dodgerblue', fontsize = 12)
    plt.text(109, .8, '$\mathbf{2p^2}$',  weight = "bold", color = 'mediumturquoise', fontsize = 12)
    plt.text(-12, -.5, '(e)', fontsize = 18)
    plt.plot(dist_t,T_P,'--', color =  'dodgerblue', linewidth = 2)
    plt.plot(dist_t,T_P*2,'--', color = 'mediumturquoise', linewidth = 2)
    plt.ylim(5,0)

plt.grid()
plt.title('Earthquake vs deep source Axitra ACFs \nFilter: ' + str(period1)+'-'+str(period2)+ ' s'  + ', STF duration: ' + str(t0) +' s ',fontsize=fnt +4 )
plt.ylabel('Time (s)',fontsize=fnt )
plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)



pos1 = axs1.get_position() # get the original position 
pos1 = [pos1.x0, pos1.y0+.17 , pos1.width, pos1.height-.1] 
axs1.set_position(pos1)
pos2 = axs2.get_position() # get the original position 
pos2 = [pos2.x0, pos2.y0-.03 , pos2.width, pos2.height+.1] 
axs2.set_position(pos2)
fig.savefig('../Figures/Fig_4_Line_' + Line_nb  + '_Earthquake.png', dpi = 100)

