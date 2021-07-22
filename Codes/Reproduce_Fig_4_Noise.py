#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 16:03:24 2021

@author: loic
"""


import os
import geopy
import numpy as np
import geopy.distance
import scipy.io as sio
from scipy import signal
from zipfile import ZipFile
import matplotlib.pyplot as plt
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
    this function runs a smoothing average of an array.
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
    tdata2 =  tdata[int(len(tdata)/2):int(len(tdata)/2+timesave*delta)]
    tdata2[:len(ctaper)] *=ctaper
    tdata2[-len(ctaper):] *=ctaper[::-1]
    corr = bandpass(tdata2,freqmin=flilow,freqmax= flihigh, df=delta,corners=4, zerophase=True)
    return corr
#%% Directory with all the data and metadata
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
Line_nb = '1' # Possible values are 1 or 3(str)

#%% Some more parameters
distlimm = 20  # radius for the station selection to remove the average trace (in km) -> can be changed but will only impact Noise ACF results
maxstation = 5 # Remove noise ACF if less than 5 stations are within a 20 km radius. -> can be changed but will only impact Noise ACF results
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
    
#%% MeSOnet station list 
all_sta = ['E.DSCM', 'E.AM1M', 'E.ANSM', 'E.AZMM', 'E.BKKM', 'E.BTSM', 'E.CF3M', 'E.DAIM', 'E.DGRM', 'E.DICM', 'E.ENZM', 'E.F3KM', 'E.FJSM', 'E.FKCM', 'E.FKMM', 'E.FKOM', 'E.FMNM', 'E.FRAM', 'E.FSMM', 'E.FTKM', 'E.FTPM', 'E.FUNM', 'E.GHGM', 'E.GKSM', 'E.GNZM', 'E.GSJM', 'E.HGCM', 'E.HKBM', 'E.HNOM', 'E.HNPM', 'E.HNSM', 'E.HONM', 'E.HRGM', 'E.HRMM', 'E.HSBM', 'E.HSDM', 'E.HTTM', 'E.HYDM', 'E.HYHM', 'E.HYSM', 'E.ICEM', 'E.IDUM', 'E.IIDM', 'E.IKBM', 'E.IKCM', 'E.IMIM', 'E.IN3M', 'E.INAM', 'E.INOM', 'E.IWNM', 'E.JDJM', 'E.JKPM', 'E.JNUM', 'E.JUWM', 'E.JYHM', 'E.KBAM', 'E.KBRM', 'E.KBTM', 'E.KCBM', 'E.KD9M', 'E.KDKM', 'E.KGKM', 'E.KH2M', 'E.KHBM', 'E.KHDM', 'E.KIMM', 'E.KK2M', 'E.KKHM', 'E.KKSM', 'E.KM5M', 'E.KMGM', 'E.KMHM', 'E.KMKM', 'E.KMNM', 'E.KMRM', 'E.KMTM', 'E.KNDM', 'E.KNJM', 'E.KNMM', 'E.KNZM', 'E.KOHM', 'E.KOYM', 'E.KRCM', 'E.KRPM', 'E.KSGM', 'E.KSHM', 'E.KSNM', 'E.KSOM', 'E.KSRM', 'E.KTOM', 'E.KUDM', 'E.KUYM', 'E.KW8M', 'E.KWHM', 'E.KWWM', 'E.KYNM', 'E.KYTM', 'E.KZ2M', 'E.KZMM', 'E.KZTM', 'E.MBSM', 'E.MBYM', 'E.MD1M', 'E.MDHM', 'E.MGRM', 'E.MHKM', 'E.MKBM', 'E.MKJM', 'E.MKSM', 'E.MNAM', 'E.MNHM', 'E.MNKM', 'E.MNMM', 'E.MNYM', 'E.MOKM', 'E.MRTM', 'E.MSKM', 'E.MSOM', 'E.MSRM', 'E.MYHM', 'E.MYMM', 'E.MZMM', 'E.MZPM', 'E.MZUM', 'E.NARM', 'E.NBKM', 'E.NDOM', 'E.NGOM', 'E.NGSM', 'E.NISM', 'E.NKDM', 'E.NKGM', 'E.NKMM', 'E.NKNM', 'E.NKSM', 'E.NMDM', 'E.NNBM', 'E.NNGM', 'E.NNSM', 'E.NNTM', 'E.NOBM', 'E.NS7M', 'E.NSHM', 'E.NSJM', 'E.NSKM', 'E.NSMM', 'E.NSOM', 'E.NSUM', 'E.NTNM', 'E.OA5M', 'E.OANM', 'E.OBRM', 'E.OBTM', 'E.OG6M', 'E.OHSM', 'E.OJCM', 'E.OKCM', 'E.OKDM', 'E.OMKM', 'E.OMNM', 'E.OMSM', 'E.ONMM', 'E.OOYM', 'E.OYMM', 'E.OYOM', 'E.OYTM', 'E.RKGM', 'E.RMSM', 'E.RYGM', 'E.RYNM', 'E.SBAM', 'E.SDMM', 'E.SEBM', 'E.SECM', 'E.SFHM', 'E.SGOM', 'E.SGWM', 'E.SICM', 'E.SJSM', 'E.SKCM', 'E.SKEM', 'E.SKHM', 'E.SKPM', 'E.SKYM', 'E.SMGM', 'E.SMNM', 'E.SNHM', 'E.SNJM', 'E.SNSM', 'E.SONM', 'E.SR2M', 'E.SRCM', 'E.SRGM', 'E.SRIM', 'E.SRKM', 'E.SRSM', 'E.SRTM', 'E.SSDM', 'E.SSHM', 'E.SSMM', 'E.SSNM', 'E.SSPM', 'E.SSSM', 'E.ST2M', 'E.STHM', 'E.STKM', 'E.SW2M', 'E.SYOM', 'E.SYPM', 'E.SYSM', 'E.TAKM', 'E.TBKM', 'E.TGNM', 'E.TK2M', 'E.TKKM', 'E.TKMM', 'E.TKNM', 'E.TKSM', 'E.TKWM', 'E.TKZM', 'E.TMHM', 'E.TMOM', 'E.TNKM', 'E.TNRM', 'E.TOCM', 'E.TOKM', 'E.TREM', 'E.TRSM', 'E.TSCM', 'E.TSRM', 'E.TTNM', 'E.TTOM', 'E.TWDM', 'E.TYHM', 'E.TYNM', 'E.TYPM', 'E.TYUM', 'E.U1CM', 'E.UHRM', 'E.UNHM', 'E.UNMM', 'E.URMM', 'E.USCM', 'E.UTKM', 'E.YGHM', 'E.YKBM', 'E.YKDM', 'E.YKKM', 'E.YKSM', 'E.YMKM', 'E.YMMM', 'E.YMPM', 'E.YNDM', 'E.YNMM', 'E.YROM', 'E.YRPM', 'E.YSKM', 'E.YSOM', 'E.YSPM', 'E.YSSM', 'E.YT2M', 'E.YTBM', 'E.YUKM', 'E.YYIM', 'E.KYDM', 'E.AYHM', 'E.KSCM' , 'E.MRJM' ,'E.THCM', 'E.FMHM', 'E.OHYM', 'E.SBCM' ,'E.TACM' ,'E.HSMM' ,'E.SIBM', 'OK.AOCM', 'OK.AONM', 'OK.ARMM', 'OK.HRDM', 'OK.KRHM', 'OK.KTGM', 'OK.NHMM', 'OK.NKYM', 'OK.NRAM', 'OK.TKCM' ] 


#%% Parameters for each line
if Line_nb == '1':
    station = ['E.THCM','E.DICM','E.TYHM','E.TBKM','E.SSHM','E.NKMM','E.TK2M','E.HSDM','E.SRTM','E.NKNM','E.TSCM','E.KUDM', 'E.HRGM', 'E.RKGM' ,'E.MNKM','E.FUNM', 'E.KMRM', 'E.NNTM', 'E.STHM','E.TYNM','E.KSCM','E.MDHM', 'E.HKBM','E.TKNM','E.SDMM', 'E.YNMM', 'E.NSMM', 'E.HGCM','E.HYHM','E.SRCM','E.DAIM','E.OKCM', 'E.FMHM','E.FMNM', 'E.NISM', 'E.FKCM', 'E.OKDM', 'E.KSRM', 'E.IIDM', 'E.SYPM' , 'E.YSSM', 'E.SBCM', 'E.TACM'] 
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
delta = 20 # Sampling frequency of the data (in Hz)
lag_selct = 15*delta # select 15 s of the positive lag of the ACFs
tnoise = np.linspace(0,lag_selct/delta,lag_selct) # Create a time vector

#%% Load Noise ACFs for stations along the line
dist4 = []
datsavi = np.empty( (len(station), len(tnoise)) )
indv = nm_sta.index(station[0])
sourlat = float(lat[indv])
sourlon = float(lon[indv])
stalonsave = []
stalatsave = []
for  i in np.arange(len(station)):
    indr = nm_sta.index(station[i])
    stalat = float(lat[indr])
    stalon = float(lon[indr])
    stalatsave.append( float(lat[indr]) )
    stalonsave.append( float(lon[indr]) )
    dist4.append(geopy.distance.distance((sourlat, sourlon), (stalat, stalon)).km)
    
    fileSta = data_dir +'ACFs/Noise_ACFs/' + station[i]+ '_ACF_ALL_pws_30.mat'
    data_mat = sio.loadmat(fileSta)  
    dat1 = np.squeeze(data_mat['AC'])
    dat1 -= np.mean(dat1)
    dat1 = butter_bandpass_filter(dat1, lowcut = cut, highcut = cut2, fs = delta, order = order) # Bandpass filter
    datsavi[i,:] = np.array( dat1[:lag_selct] )
datsavi = np.array(datsavi) # To array

toposta = griddata((latg,long), topo, (stalatsave,stalonsave), method='cubic')/1000 # Get topography at the stations

#%% Load Noise ACFs at all the MeSOnet stations
datALL = np.empty( (len(all_sta),len(tnoise)) )
for  i in np.arange(len(all_sta)):
    fileAll = data_dir +'ACFs/Noise_ACFs/' + all_sta[i] + '_ACF_ALL_pws_30.mat'
    data_mat = sio.loadmat(fileAll)  # Load Noise data
    dat1 = np.squeeze(data_mat['AC'])
    dat1 -= np.mean(dat1)
    dat1 = butter_bandpass_filter(dat1, lowcut = cut, highcut = cut2, fs = delta, order = order) # bandpass filter
    datALL[i,:] = np.array( dat1[:lag_selct] )  # Get the causal part
datALL = np.array(datALL) # To array

#%%  Average trace removal step
datsav = np.empty(datsavi.shape) # Create empty array
for i in np.arange(len(station)):
    val = []
    for j in  np.arange( len(datALL) ):
        spec = nm_sta.index(station[i])
        oth = nm_sta.index(all_sta[j])
        dist_inter = geopy.distance.distance( (float(lat[spec]),float(lon[spec]) ), (float(lat[oth]),float(lon[oth])) ).km
        if  dist_inter > .1 and dist_inter < distlimm: 
            val.append(j) # Get station indice if within 20 km
    if len(val) <= maxstation: 
        datsav[i,:] = np.NaN # If less than 5 stations around, dont consider the ACF 
    else: 
        mean_dat = np.nanmean(datALL[val,:],axis =0) # Compute the average trace #nanmedian
        inte = datsavi[i,:] - mean_dat # Remove the average trace
        datsav[i,:] = inte/np.max(np.abs(inte)) # Normalized trace

#%% Assign Noise and EQ ACFs to bins for the plot
minkm = 5 # add 5 km to the begining of the line for visibility
nnb = 1 # width of the bin in km -> (can be changed: but might impact the final results)
dist_plot = np.linspace(-minkm,mxdist-minkm, int(mxdist/nnb)+1 ) # Create an array of distances along the line
dat_fin2 = np.zeros( (len(dist_plot), datsav.shape[1]) ) # set up array for Noise ACFs
for i in np.arange(len(dist_plot)):
    nbdat = 0
    for ii in np.arange(len(station)):
        if dist4[ii] >= dist_plot[i] and dist4[ii] < dist_plot[i] + nnb:
            dat_fin2[i] = datsav[ii] # if there is a station at the within the distance, fill the array with the noise ACF
            nbdat += 1 # To check if there is more than one station in the bin
    if nbdat>1:
        print('This bin contains ' + str(nbdat) + ' stations')
    
#%% Load Axitra noise simulations
unit = 3
utext = 'acceleration'
t0 = 1.6
fname = data_dir + 'Axitra/Axitra_Surface_Simulations_Line_' + Line_nb +'_stf_' + str(t0) + '_type_' + str(unit) +'.mat'
fnametot = data_dir + 'Axitra/Axitra_Surface_Simulations_Line_MeSOnet_stations_stf_' + str(t0) + '_type_' + str(unit) +'.mat'

isfile = os.path.isfile(fname) 
isfile2 = os.path.isfile(fnametot) 
dat = sio.loadmat(fname)
autosave = np.squeeze(dat['data'])
tsim = np.squeeze(dat['t'])
latout = np.squeeze(dat['lat'])
lonout = np.squeeze(dat['lon'])

dattot = sio.loadmat(fnametot)
autosavetot = np.squeeze(dattot['data'])
tsimt = np.squeeze(dattot['t'])
latall = np.squeeze(dattot['lat'])
lonall = np.squeeze(dattot['lon'])

#%% Create JIVSM profiles and compute theoretical 2p travel time
grid_z0 = griddata((latg,long),topo, (latout,lonout), method='cubic')/1000
grid_z1 = griddata((latg,long),lay1, (latout,lonout), method='cubic')/1000
grid_z2 = griddata((latg,long),lay2, (latout,lonout), method='cubic')/1000
grid_z3 = griddata((latg,long),lay3, (latout,lonout), method='cubic')/1000
T_P = np.zeros(len(grid_z0))
for i in np.arange(len(grid_z0)):
    T_P[i] = (0 - grid_z1[i])/1.8 + (grid_z1[i] - grid_z2[i])/2.3  +(grid_z2[i] - grid_z3[i])/3.0
T_P =T_P*2 # multiply by two to get the two-way travel time


#%% Compute simulated waveforms ACFs 
deltasim = 1/(tsim[11]-tsim[10])
deltasimtot= 1/(tsimt[11]-tsimt[10])

datoi = 22 # Skip the first 2 s as there is a 2 s delay in the simulated waveforms
smoothN = 30 # Parameter for the denominator term of Equation 1
hann2 = signal.hann(int(delta*.5)*2+1) # Taper
half_len =round(len(hann2)/2-.001)+1
ctaper = hann2[:half_len]

timesave = 15 # in sec

# Compute ACFs at all the virtual stations along the line
ACFini = []
for i in range(len(autosave)):
    datACF = autosave[i] 
    corr = cross_cor_Chengxin(datACF[datoi:400], deltasim, smoothN, cut, cut2, ctaper, timesave)
    ACFini.append(corr)
ACFini = np.array(ACFini)

# Compute ACFs at all the MeSO-net stations
ACFtot = []
for i in range(len(autosavetot)):
    datACF = autosavetot[i]
    corr = cross_cor_Chengxin(datACF[datoi:400], deltasimtot, smoothN, cut, cut2, ctaper, timesave)
    ACFtot.append(corr)
ACFtot = np.array(ACFtot)

#%% Average trace removal for the simulated noise ACFs
datsim = np.empty(ACFini.shape) # Create empty array
for i in np.arange(len(ACFini)):
    val = []
    for j in  np.arange( len(ACFtot) ):
        dist_inter = geopy.distance.distance( (float(latout[i]),float(lonout[i]) ), (float(latall[j]),float(lonall[j])) ).km
        if  dist_inter > .1 and dist_inter < distlimm: 
            val.append(j) # Get station indice if within 25 km
    if len(val) <= maxstation: 
        datsim[i,:] = np.NaN # If less than 10 stations around, dont consider the ACF 
    else: 
        mean_dat = np.nanmean(ACFtot[val,:], axis = 0) # Compute the average trace #nanmedian
        inte = ACFini[i,:] - mean_dat # Remove the average trace
        datsim[i,:] = inte/np.max(np.abs(inte))
dfin = np.array(datsim)


#%%
tpl = np.linspace(0,dfin.shape[1]/deltasimtot,dfin.shape[1]) # Create a time vector for the simulated ACFs
dist_t = np.linspace(0, mxdist, int((mxdist)/nnb)+1 )# Create a distance vector for the simulated ACFs
dist_t = dist_t-minkm

#%% Plot simulated and observed Noise ACF profile

fig = plt.figure(figsize = (10,8) )
fnt = 13
# Plot the velocity profile for the line
axs1 = plt.subplot(211)
plt.plot(dist_t, grid_z0, 'k') # Surface
plt.plot(dist_t, grid_z1, 'grey') # Layer 1
plt.plot(dist_t, grid_z2, 'grey') # Layer 2
plt.plot(dist_t, grid_z3, 'k' , linewidth = 4) #  Bedrock
plt.plot(dist_t, grid_z0+ .07, 'v', markersize = 5, mec = 'k') # Plot station location

plt.title('Line '+Line_nb+  ': JIVSM P-wave velocity profile', fontsize =fnt +4)
if Line_nb == '1':
    plt.text(-14, .3, '(a)', fontsize = 18)
elif Line_nb == '3':
    plt.text(-12, .3, '(d)', fontsize = 18)
plt.ylabel('Depth (km)', fontsize =fnt)
plt.xlabel('Distance along Line '  +Line_nb +' (km)', fontsize =fnt)

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

plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)
plt.xlim(dist_t[3],dist_t[-3])
plt.xticks(fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.grid()




axs2 = plt.subplot(212)
# Plot simulated ACFs
factval= 2
stalim = [5, -4]
for i in np.arange(stalim[0],len(dfin)+stalim[-1]):
    y = dist_t[i]+ dfin[i]*factval 
    plt.plot(y ,tpl, 'k', linewidth = .5)
    axs2.fill_betweenx(tpl , dist_t[i], y, where=(y < dist_t[i]), color='gray')
    axs2.fill_betweenx(tpl , dist_t[i], y, where=(y > dist_t[i]), color='k')
plt.xlim(dist_t[3],dist_t[-3])
plt.xlabel('Distance along Line '  +Line_nb +' (km)', fontsize =fnt)
plt.xticks(fontsize=fnt)
plt.yticks(fontsize=fnt)
plt.ylim(10 ,0)
plt.grid()
plt.title('Noise vs surface source Axitra ACFs\nFilter: ' + str(period1)+'-'+str(period2)+ ' s'  + ', STF duration: ' + str(t0) +' s ' , fontsize = fnt +4) 
plt.ylabel('Time (s)', fontsize =fnt )
plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)


# Plot observed ACFs
for i in np.arange(len(dat_fin2)):
    factampeq = 2
    if np.mean(dat_fin2[i,:50])==0:
        dat_fin2[i] = np.nan      
    y = dist_t[i]+ dat_fin2[i]*factampeq
    if not np.isnan(np.sum(y)):
        plt.plot(y ,tnoise,'k', linewidth =.5)
        axs2.fill_betweenx(tnoise , dist_t[i], y , where=(y< dist_t[i]) ,color='green')
        axs2.fill_betweenx(tnoise , dist_t[i], y , where=(y> dist_t[i]) ,color='blue')
 
# Plot the parameters specific to a line
if Line_nb =='1':
        plt.text(122, .4, '2p',  weight="bold", color = 'dodgerblue', fontsize = 12)
        plt.text(122, 1, '$\mathbf{2p^2}$',  weight = "bold", color = 'mediumturquoise', fontsize = 12)
        plt.text(122, 1.6, '$\mathbf{2p^3}$', weight = "bold", color = 'orange', fontsize = 12)   
        plt.plot(dist_t,T_P,'--', color =  'dodgerblue', linewidth = 3)
        plt.plot(dist_t,T_P*2,'--', color = 'mediumturquoise', linewidth = 3)
        plt.plot(dist_t,T_P*3,'--', color = 'orange', linewidth = 3)
        plt.text(-14, -1, '(c)', fontsize =18)            
elif Line_nb =='3':
        plt.text(107, .3, '2p',  weight="bold", color = 'dodgerblue', fontsize = 12)
        plt.text(107, .8, '$\mathbf{2p^2}$',  weight = "bold", color = 'mediumturquoise', fontsize = 12)
        plt.text(107, 1.4, '$\mathbf{2p^3}$', weight = "bold", color = 'orange', fontsize = 12)   
        plt.text(-12, -.5, '(f)', fontsize = 18)
        plt.plot(dist_t,T_P,'--', color =  'dodgerblue', linewidth = 3)
        plt.plot(dist_t,T_P*2,'--', color = 'mediumturquoise', linewidth = 3)
        plt.plot(dist_t,T_P*3,'--', color = 'orange', linewidth = 3)

pos1 = axs1.get_position() # get the original position 
pos1 = [pos1.x0, pos1.y0+.17 , pos1.width, pos1.height-.1] 
axs1.set_position(pos1)
pos2 = axs2.get_position() # get the original position 
pos2 = [pos2.x0, pos2.y0-.03 , pos2.width, pos2.height+.1] 
axs2.set_position(pos2)

fig.savefig('../Figures/Fig_4_Line_' + Line_nb  + '_Noise.png', dpi = 100)
plt.show()
plt.close()

