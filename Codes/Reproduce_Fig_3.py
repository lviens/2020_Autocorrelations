#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 16:14:44 2021

@author: Loic Viens

Code to reproduce the Figure 3 of Viens et al. (2022, GJI) 
The station Line number (1-3) can be changed line 36. 
Some more parameters can be changed lines 70-77, but will only affect the noise ACF results.
More parameters in the code can be changed and have a "can be changed" comment next to them.
"""

from obspy import read
import numpy as np
import scipy.io as sio
import matplotlib
from scipy.signal import butter, filtfilt, find_peaks,hann
import matplotlib.pyplot as plt
import geopy.distance
import geopy
from geopy.distance import  geodesic
from scipy.interpolate import griddata
from matplotlib.ticker import  AutoMinorLocator
import os
from zipfile import ZipFile


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

#%% Extract noise and earthquake ACF from zip file if needed
directory = '../Data/ACFs/'
if not os.path.exists(directory):
    zf = ZipFile('../Data/ACFs.zip', 'r')
    zf.extractall('../Data/')
    zf.close()

#%% Directory with all the data and metadata
data_dir = '../Data/'

#%% Main parameter to plot the data along different Lines
Line_nb = '1' # Possible values are 1, 2, 3 (str)

#%% Some more parameters
distlimm = 20  # radius for the station selection to remove the average trace (in km) -> can be changed but will only impact Noise ACF results
maxstation = 5 # Remove noise ACF if less than 5 stations are within a 20 km radius. -> can be changed but will only impact Noise ACF results
period1 = 1 # filter parameter for the Noise ACFs -> can be changed but will only impact Noise ACF results
period2 = 10 # filter parameter for the Noise ACFs-> can be changed but will only impact Noise ACF results
order = 4 # order of the filter -> can be changed but will only impact Noise ACF results
cut = 1/period2
cut2= 1/period1
delta = 20 # Sampling frequency of the data (in Hz)
lag_selct = 15*delta # select 15 s of the positive lag of the ACFs
t = np.linspace(0,lag_selct/delta,lag_selct) # Create a time vector
print('Plotting ACFs for stations along Line ' + Line_nb)
print('The ACFs are bandpass filtered between ' + str(period1) + ' and ' + str(period2) + ' s')
print('For each noise ACF, the average trace is computed by averaging all the surrounding ACFs within a ' + str(distlimm) + ' km radius' )

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
all_sta = ['E.DSCM', 'E.AM1M', 'E.ANSM', 'E.AZMM', 'E.BKKM', 'E.BTSM', 'E.CF3M', 'E.DAIM', 'E.DGRM', 'E.DICM', 'E.ENZM', 'E.F3KM', 'E.FJSM', 'E.FKCM', 'E.FKMM', 'E.FKOM', 'E.FMNM', 'E.FRAM', 'E.FSMM', 'E.FTKM', 'E.FTPM', 'E.FUNM', 'E.GHGM', 'E.GKSM', 'E.GNZM', 'E.GSJM', 'E.HGCM', 'E.HKBM', 'E.HNOM', 'E.HNPM', 'E.HNSM', 'E.HONM', 'E.HRGM', 'E.HRMM', 'E.HSBM', 'E.HSDM', 'E.HTTM', 'E.HYDM', 'E.HYHM', 'E.HYSM', 'E.ICEM', 'E.IDUM', 'E.IIDM', 'E.IKBM', 'E.IKCM', 'E.IMIM', 'E.IN3M', 'E.INAM', 'E.INOM', 'E.IWNM', 'E.JDJM', 'E.JKPM', 'E.JNUM', 'E.JUWM', 'E.JYHM', 'E.KBAM', 'E.KBRM', 'E.KBTM', 'E.KCBM', 'E.KD9M', 'E.KDKM', 'E.KGKM', 'E.KH2M', 'E.KHBM', 'E.KHDM', 'E.KIMM', 'E.KK2M', 'E.KKHM', 'E.KKSM', 'E.KM5M', 'E.KMGM', 'E.KMHM', 'E.KMKM', 'E.KMNM', 'E.KMRM', 'E.KMTM', 'E.KNDM', 'E.KNJM', 'E.KNMM', 'E.KNZM', 'E.KOHM', 'E.KOYM', 'E.KRCM', 'E.KRPM', 'E.KSGM', 'E.KSHM', 'E.KSNM', 'E.KSOM', 'E.KSRM', 'E.KTOM', 'E.KUDM', 'E.KUYM', 'E.KW8M', 'E.KWHM', 'E.KWWM', 'E.KYNM', 'E.KYTM', 'E.KZ2M', 'E.KZMM', 'E.KZTM', 'E.MBSM', 'E.MBYM', 'E.MD1M', 'E.MDHM', 'E.MGRM', 'E.MHKM', 'E.MKBM', 'E.MKJM', 'E.MKSM', 'E.MNAM', 'E.MNHM', 'E.MNKM', 'E.MNMM', 'E.MNYM', 'E.MOKM', 'E.MRTM', 'E.MSKM', 'E.MSOM', 'E.MSRM', 'E.MYHM', 'E.MYMM', 'E.MZMM', 'E.MZPM', 'E.MZUM', 'E.NARM', 'E.NBKM', 'E.NDOM', 'E.NGOM', 'E.NGSM', 'E.NISM', 'E.NKDM', 'E.NKGM', 'E.NKMM', 'E.NKNM', 'E.NKSM', 'E.NMDM', 'E.NNBM', 'E.NNGM', 'E.NNSM', 'E.NNTM', 'E.NOBM', 'E.NS7M', 'E.NSHM', 'E.NSJM', 'E.NSKM', 'E.NSMM', 'E.NSOM', 'E.NSUM', 'E.NTNM', 'E.OA5M', 'E.OANM', 'E.OBRM', 'E.OBTM', 'E.OG6M', 'E.OHSM', 'E.OJCM', 'E.OKCM', 'E.OKDM', 'E.OMKM', 'E.OMNM', 'E.OMSM', 'E.ONMM', 'E.OOYM', 'E.OYMM', 'E.OYOM', 'E.OYTM', 'E.RKGM', 'E.RMSM', 'E.RYGM', 'E.RYNM', 'E.SBAM', 'E.SDMM', 'E.SEBM', 'E.SECM', 'E.SFHM', 'E.SGOM', 'E.SGWM', 'E.SICM', 'E.SJSM', 'E.SKCM', 'E.SKEM', 'E.SKHM', 'E.SKPM', 'E.SKYM', 'E.SMGM', 'E.SMNM', 'E.SNHM', 'E.SNJM', 'E.SNSM', 'E.SONM', 'E.SR2M', 'E.SRCM', 'E.SRGM', 'E.SRIM', 'E.SRKM', 'E.SRSM', 'E.SRTM', 'E.SSDM', 'E.SSHM', 'E.SSMM', 'E.SSNM', 'E.SSPM', 'E.SSSM', 'E.ST2M', 'E.STHM', 'E.STKM', 'E.SW2M', 'E.SYOM', 'E.SYPM', 'E.SYSM', 'E.TAKM', 'E.TBKM', 'E.TGNM', 'E.TK2M', 'E.TKKM', 'E.TKMM', 'E.TKNM', 'E.TKSM', 'E.TKWM', 'E.TKZM', 'E.TMHM', 'E.TMOM', 'E.TNKM', 'E.TNRM', 'E.TOCM', 'E.TOKM', 'E.TREM', 'E.TRSM', 'E.TSCM', 'E.TSRM', 'E.TTNM', 'E.TTOM', 'E.TWDM', 'E.TYHM', 'E.TYNM', 'E.TYPM', 'E.TYUM', 'E.U1CM', 'E.UHRM', 'E.UNHM', 'E.UNMM', 'E.URMM', 'E.USCM', 'E.UTKM', 'E.YGHM', 'E.YKBM', 'E.YKDM', 'E.YKKM', 'E.YKSM', 'E.YMKM', 'E.YMMM', 'E.YMPM', 'E.YNDM', 'E.YNMM', 'E.YROM', 'E.YRPM', 'E.YSKM', 'E.YSOM', 'E.YSPM', 'E.YSSM', 'E.YT2M', 'E.YTBM', 'E.YUKM', 'E.YYIM', 'E.KYDM', 'E.AYHM', 'E.KSCM' ,'E.ABHM', 'E.MRJM' ,'E.THCM', 'E.FMHM', 'E.OHYM', 'E.SBCM' ,'E.TACM' ,'E.HSMM' ,'E.SIBM', 'OK.AOCM', 'OK.AONM', 'OK.ARMM', 'OK.HRDM', 'OK.KRHM', 'OK.KTGM', 'OK.NHMM', 'OK.NKYM', 'OK.NRAM', 'OK.TKCM' ]

#%% Parameters for the 3 lines
 
if Line_nb == '1':
    station = ['E.THCM','E.DICM','E.TYHM','E.TBKM','E.SSHM','E.NKMM','E.TK2M','E.HSDM','E.SRTM','E.NKNM','E.TSCM','E.KUDM', 'E.HRGM', 'E.RKGM' ,'E.MNKM','E.FUNM', 'E.KMRM', 'E.NNTM', 'E.STHM','E.TYNM','E.KSCM','E.MDHM', 'E.HKBM','E.TKNM','E.SDMM', 'E.YNMM', 'E.NSMM', 'E.HGCM','E.HYHM','E.SRCM','E.DAIM','E.OKCM', 'E.FMHM','E.FMNM', 'E.NISM', 'E.FKCM', 'E.OKDM', 'E.KSRM', 'E.IIDM', 'E.SYPM' , 'E.YSSM', 'E.SBCM', 'E.TACM'] 
    leftt = 'E'
    rightt = 'W'
    az_mean = 100
    mxdist = 128   
elif Line_nb == '2':
    station = ['E.SICM','E.MSKM','E.SSNM','E.YSOM','E.MSRM','E.YGHM','E.SSDM','E.INOM','E.YMMM','E.SYPM','E.KMRM', 'E.SR2M', 'E.TGNM','E.NDOM','E.KKHM','E.TNKM','E.KUYM','E.MRJM','E.MYHM','E.SGOM','E.NNGM', 'E.NKGM','E.YSPM','E.NSUM', 'E.MRTM', 'E.SSMM', 'E.NOBM', 'E.SRSM', 'E.FKOM', 'E.MHKM', 'E.HYDM']
    az_mean = 142
    mxdist = 100
elif Line_nb == '3':
    station = ['E.SYOM','E.FJSM','E.DGRM','E.KMHM','E.SMGM','E.HTTM','E.KSOM','OK.HRDM','E.HNOM', 'E.BKKM','E.MZPM', 'E.JDJM', 'E.SNHM', 'E.UNMM','E.KMKM','E.MKSM','E.GKSM','E.KDKM','E.ENZM', 'E.TKMM','E.SBAM','E.GNZM','E.YKKM','E.MKJM', 'E.KCBM','E.MZMM', 'E.TKNM','E.MBSM','E.YKSM','E.KGKM', 'E.NGSM','E.KKSM', 'E.KWHM', 'E.TNKM','E.KUYM', 'E.IN3M','E.INAM', 'E.KBRM','E.GSJM', 'E.YTBM', 'E.HNPM','E.YKBM','E.TSRM', 'E.TKZM','E.RMSM','E.RYGM', 'E.RYNM']
    az_mean = 34.2
    mxdist = 113
    
#%% Load earthquake ACFs for stations along the line
datsaviEQ = np.empty((len(station),len(t)) )
for i in np.arange(len(station)):
    fileEQ = data_dir +'ACFs/EQ_ACFs/' + station[i]+ '.SAC'
    st = read(fileEQ) # Read sac file
    dat1 = st[0].data # get data
    dat1 -= np.mean(dat1) # Remove the mean
    datsaviEQ[i,:] = dat1 /np.max(np.abs(dat1)) # normalize the trace
datsavEQ = np.array(datsaviEQ)  # To array

#%% Load Noise ACFs for stations along the line
dist4 = []
datsavi = np.empty( (len(station), len(t)) )

indv = nm_sta.index(station[0])
sourlat = float(lat[indv])
sourlon = float(lon[indv])
stalonsave = []
stalatsave = []
for  i in np.arange(len(station)):
    dat2 = []
    indr = nm_sta.index(station[i])
    stalat = float(lat[indr])
    stalon = float(lon[indr])
    stalatsave.append( float(lat[indr]) )
    stalonsave.append( float(lon[indr]) )
    dist4.append(geopy.distance.distance((sourlat, sourlon), (stalat, stalon)).km)
    
    fileSta = data_dir +'/ACFs/Noise_ACFs/' + station[i]+ '_ACF_ALL_pws_30.mat'
    data_mat = sio.loadmat(fileSta)   # Load mat
    dat1 = np.squeeze(data_mat['AC']) # Get data
    dat1 -= np.mean(dat1) # remove mean
    dat1 = butter_bandpass_filter(dat1, lowcut = cut, highcut = cut2, fs = delta, order = order) # Bandpass filter
    dat2 = np.array(dat1[:lag_selct] ) # Get the causal part
    datsavi[i,:] = dat2 #/np.max(np.abs(dat2)) # Normalize the ACF
datsavi = np.array(datsavi) # To array

toposta = griddata((latg,long), topo, (stalatsave,stalonsave), method='cubic')/1000 # Get topography at the stations


#%% Load Noise ACFs at all the MeSOnet stations
datALL = np.empty( (len(all_sta),len(t)) )
for  i in np.arange(len(all_sta)):
    dat2=[]
    fileAll = data_dir +'ACFs/Noise_ACFs/' + all_sta[i]+ '_ACF_ALL_pws_30.mat'
    data_mat = sio.loadmat(fileAll)  # Load Noise data
    dat1 = np.squeeze(data_mat['AC'])
    dat1 -= np.mean(dat1)
    dat1 = butter_bandpass_filter(dat1, lowcut = cut, highcut = cut2, fs = delta, order = order) # bandpass filter
    dat2 = np.array( dat1[:lag_selct] )  # Get the causal part
    datALL[i,:] = dat2 #/np.max(np.abs(dat2)) # Normalize the ACF
datALL = np.array(datALL) # To array


#%% Compute and remove an average trace to each station. The average trace is computed by averaging the ACFs from the stations within a 25-km radius
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
        mean_dat = np.nanmean(datALL[val,:],axis = 0) # Compute the average trace
        inte = datsavi[i,:] - mean_dat # Remove the average trace
        datsav[i,:] = inte/np.max(np.abs(inte)) # Normalize ACF

#%% Assign Noise and EQ ACFs to bins for the plot
minkm = 5 # add 5 km to the begining of the line for visibility
nnb = .25 # width of the bin in km -> (can be changed: but might impact the final results)

dist_plot = np.linspace(-minkm,mxdist-minkm, int(mxdist/nnb)+1 ) # Create an array of distances along the line

dat_fin2 = np.zeros( (len(dist_plot), datsav.shape[1]) ) # set up array for Noise ACFs
dat_fin2EQ = np.zeros( (len(dist_plot), datsavEQ.shape[1]) ) # set up array for EQ ACFs
for i in np.arange(len(dist_plot)):
    nbdat = 0
    for ii in np.arange(len(station)):
        if dist4[ii] >= dist_plot[i] and dist4[ii] < dist_plot[i] + nnb:
            dat_fin2[i] = datsav[ii] # if there is a station at the within the distance, fill the array with the noise ACF
            dat_fin2EQ[i] = datsavEQ[ii] # Same fo the EQ ACF
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

#%% Get the JIVSM two-way travel time at the bin locations (The Kanto Basin has 3 layers)
# Layer 1: Vp = 1.8 km/s
# Layer 2: Vp = 2.3 km/s
# Layer 3: Vp = 3.0 km/s
T_PJIVSM = np.zeros(len(grid_z0))
for i in np.arange(len(grid_z0)):
    T_PJIVSM[i] = (grid_z0[i] - grid_z1[i])/1.8 + (grid_z1[i] - grid_z2[i])/2.3  +(grid_z2[i] - grid_z3[i])/3.0
T_PJIVSM =T_PJIVSM*2 # multiply by two to get the two-way travel time


#%% Construct the P-wave two-way travel times estimates for the 3 Lines
if Line_nb=='1':
    T_P = np.concatenate(( np.linspace(-0.5, 1.75, int(8*1/nnb) ), np.ones(int(15*1/nnb))*1.75, np.linspace(1.75, 2.5, int(10*1/nnb) ), np.ones(int(5*1/nnb))*2.5, np.linspace(2.5,1.9, int(5*1/nnb)), np.linspace(1.9,1.7, int(30*1/nnb )) , np.linspace(1.7,.8, int(70*1/nnb ))   ))
    
elif Line_nb=='2':
    T_P = np.concatenate((np.ones(int(5*1/nnb))*1.7, np.linspace(1.7, 1.15, int(14*1/nnb )), np.ones(int(20*1/nnb))*1.15, np.linspace(1.15, 2.5, int(33*1/nnb)), np.ones(int(50*1/nnb))*2.5) )
elif Line_nb=='3':
    T_P = np.concatenate(( np.ones(int(10*1/nnb))*1.3 ,np.linspace(1.3, 2.9, int(10*1/nnb)) ,np.linspace( 2.9,2.2, int(10*1/nnb)), np.ones(int(30*1/nnb))*2.2, np.linspace(2.2, .9, int(33*1/nnb)), np.linspace(.9, .4 , int(30*1/nnb )) ))


distTs = np.linspace(0, len(T_P)*nnb, int(len(T_P) ) ) - minkm # Distance axis for plots


#%%
#  Get Earthquake ACFs P-wave two-way travel time
travel_nbEQ = 1 # number of two-way travel times
v13EQ = -.5 # Search negative peaks within pm 0.65 s (can be changed)
v23EQ = .5 # (can be changed)

mxpk3EQ = np.empty(len(dat_fin2EQ)) # Initiate an array of NaNs that will contain the travel time for each bin
mxpk3EQ[:] = np.NaN
mxpk3allEQ = np.zeros((len(dat_fin2EQ),5))-1 # Array to keep in memory all the negative peaks within the period range of interest
for i in np.arange(len(dat_fin2EQ)):
    if np.mean(dat_fin2EQ[i]) != 0 and T_P[i]!= 'nan':
         peaks3EQ =[]
         peaks3EQ, _ = find_peaks(dat_fin2EQ[i]*-1, height=0,prominence=(0.001, 100)) # Get peaks, multiply by -1 as we want the negative peaks
         temppkEQ =[]
         amp = []
         for ii in np.arange(len(peaks3EQ)):
                 if t[peaks3EQ[ii]] >T_P[i] +v13EQ  and t[peaks3EQ[ii]]< T_P[i] + v23EQ:
                     temppkEQ.append(t[peaks3EQ[ii]])
                     amp.append(np.abs(dat_fin2EQ[i][peaks3EQ[ii]]))
                     
         if len(temppkEQ)==1:
            mxpk3EQ[i] = temppkEQ[0]
         elif len(temppkEQ)>1: # If multiple peaks, get the largest one
             temppkEQ[np.argmax(amp)  ]

mxpk3EQ = mxpk3EQ/travel_nbEQ  # Divide by the number of travel times (Not necessary for EQ ACFs)
mxpk3allEQ = mxpk3allEQ/travel_nbEQ # Divide by the number of travel times (Not necessary for EQ ACFs)


#%% Compute Noise ACF P-wave two-way travel time from 2p3 phases
travel_nb3 = 3 # Number of two-way travel times
v13 = -1.5 # Search negative peaks within pm 2.5 s (can be changed)
v23 = 1.5  #(can be changed)

mxpk3 = np.empty(len(dat_fin2)) # Initiate an array of NaNs that will contain the travel time for each bin
mxpk3[:] = np.NaN

mxpk3all = np.zeros((len(dat_fin2),5))-1 # Array to keep in memory all the negative peaks within the period range of interest

for i in np.arange(len(dat_fin2)):
    if np.mean(dat_fin2[i]) != 0 and T_P[i]!= 'nan':
        peaks3 =[]
        peaks3, _ = find_peaks(dat_fin2[i]*-1, height = 0,prominence=(0.001, 100)) # Get peaks, multiply by -1 as we want the negative peaks
        temppk = []
        amp =[]
        for ii in np.arange(len(peaks3)):
            if t[peaks3[ii]] > T_P[i]*travel_nb3+v13 and t[peaks3[ii]] < T_P[i]*travel_nb3+v23:
                temppk.append(t[peaks3[ii] ])
                amp.append(np.abs(dat_fin2[i][peaks3[ii]]))
        if len(temppk)==1: # Store the negative peak values if any
            mxpk3[i] = temppk[0]
        elif len(temppk)>1: # If multiple peaks, get the one that is the closest to the JIVSM
            mxpk3[i] = temppk[np.argmax(amp) ]
                
mxpk3 = mxpk3/travel_nb3  # Divide by the number of travel times (e.g., 3) to get the two-way travel time
mxpk3all = mxpk3all/travel_nb3  # Divide by the number of travel times (e.g., 3) to get the two-way travel time


#%% Plot the data

# Color of the JIVSM lines
Cestimate = 'peru'
CestimateEQ = 'navy'

col2p = 'dodgerblue'
col2p2 = 'mediumturquoise'
col2p3 = 'orange'



matplotlib.rcParams.update({'font.size': 14})
figo = plt.figure(figsize =(8,14))

# Plot the JIVSM along the Line
axs1 = plt.subplot(4,1,1)
plt.plot(dist_t, grid_z0, 'k') # Surface
plt.plot(dist_t, grid_z1, 'grey') # Layer 1
plt.plot(dist_t, grid_z2, 'grey') # Layer 2
plt.plot(dist_t, grid_z3, 'k' , linewidth = 4) #  Bedrock
plt.plot(dist4, toposta + .175, 'v', markersize = 12, mec = 'k') # Plot station location
plt.grid()
plt.title('Line ' + Line_nb + ': JIVSM P-wave velocity profile' )
plt.ylabel('Depth (km)')
plt.xlim(dist_t[0], dist_t[-1])
plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)
axs1.xaxis.set_minor_locator(AutoMinorLocator())

# Some more parameters that are line specific  
if  Line_nb == '1':
    plt.text(58, -0.32, '1.8')
    plt.text(58, -.9, '2.3')
    plt.text(58, -1.5, '3.0')
    plt.text(47, -3, '$\mathrm{V_p}$ = 5.5 km/s')
    plt.text(-18, .8, '(a)', fontsize = 18)
    axs1.annotate("", xy=(-1.5, -2.9), xytext=(10.5, -2.9), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(2.5, -3.5, 'W')
    axs1.annotate("", xy=(119, -2.9), xytext=(107, -2.9), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(111, -3.5, 'E')
    plt.ylim(-4.01, .7)    
elif  Line_nb == '2':
    plt.text(48, -0.35, '1.8')
    plt.text(48, -.85, '2.3')
    plt.text(48, -1.55, '3.0')
    plt.text(40, -2.8, '$\mathrm{V_p}$ = 5.5 km/s')
    plt.text(-16, .8, '(e)', fontsize = 18)
    axs1.annotate("", xy = (-4, -3.4), xytext=(4, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(-2, -3.9, 'NW')
    axs1.annotate("", xy = (94, -3.4), xytext=(86, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text( 87.7, -3.9, 'SE')
    plt.ylim(-4.01, .7)    
elif  Line_nb == '3':
    plt.text(41, -0.34, '1.8')
    plt.text(41, -1, '2.3')
    plt.text(41, -2.4, '3.0')
    plt.text(41-9, -3.7, '$\mathrm{V_p}$ = 5.5 km/s')
    plt.text(-18, .8, '(i)', fontsize = 18)
    axs1.annotate("", xy=(-2, -3.4), xytext=(10, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(1, -3.9, 'SW')
    axs1.annotate("", xy=(105, -3.4), xytext=(93, -3.4), arrowprops = dict(arrowstyle="simple",facecolor='k') )
    plt.text(97, -3.9, 'NE')
    plt.ylim(-4.01, .7)
    


##%%
# Plot earthquake ACFs
axs3 = plt.subplot(4,1,2)    
plt.plot(dist_t, T_PJIVSM*travel_nbEQ, ':', color = col2p, linewidth =2) # Plot JIVSM two-way travel time
plt.plot(distTs, T_P, color = CestimateEQ)   # Plot estimated P-wave two-way travel time
plt.plot(distTs, T_P+v13EQ, '--', color = CestimateEQ, linewidth = 2)  # Plot estimated P-wave two-way travel time - 0.5 s 
plt.plot(distTs, T_P+v23EQ, '--', color = CestimateEQ, linewidth = 2)  # Plot estimated P-wave two-way travel time + 0.5 s   

for i in np.arange(len(dat_fin2EQ)): # Plot the earthquake ACF waveforms
    factampeq = 4
    if np.mean(dat_fin2EQ[i,:50])==0:
        dat_fin2EQ[i] = np.nan      
    y = dist_t[i]+ dat_fin2EQ[i]*factampeq
    if not np.isnan(mxpk3EQ[i]):
        yminEQ = np.min(dat_fin2EQ[i, int(mxpk3EQ[i]*travel_nbEQ*delta)]*factampeq)
        plt.plot(dist_t[i]+nnb/2+yminEQ ,mxpk3EQ[i]*travel_nbEQ, 'go', mec ='k', markersize=6)
    if not np.isnan(np.sum(y)):
        plt.plot(y ,t,'k', linewidth =.5)
        axs3.fill_betweenx(t , dist_t[i], y , where=(y< dist_t[i]) ,color='gray')
        axs3.fill_betweenx(t , dist_t[i], y , where=(y> dist_t[i]) ,color='k')

if Line_nb == '1':
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth =2)  
     plt.text(124, .4, '2p',  weight="bold", color = col2p, fontsize = 12) 
     plt.text(124, 1.1, 'Est.',  weight="bold", color = CestimateEQ, fontsize = 12) 
     plt.text(124, .8, '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)
     plt.text(-18 ,-.2, '(b)', fontsize = 18) 
elif Line_nb == '2':
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth =2)
     plt.text(96.4, 2.5, '2p',  weight = "bold", color = col2p, fontsize = 12)
     plt.text(96.4, 2.1, 'Est.',  weight = "bold", color = CestimateEQ, fontsize = 12)
     plt.text(96.4, 4.7, '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)
     plt.text(-16 ,-.2, '(f)', fontsize = 18)
elif  Line_nb == '3':
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth = 2)
     plt.text(109.5, .2, '2p',  weight="bold", color = col2p, fontsize = 12)
     plt.text(-18 ,-.2, '(j)', fontsize = 18)
     plt.text(109.5, .6, 'Est.',  weight = "bold", color = CestimateEQ, fontsize = 12) 
     plt.text(109.5, 1., '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)

plt.ylim(5, 0)
plt.xlim(dist_t[0], dist_t[-1])
plt.ylabel('Time (s)')
plt.title('Vertical earthquake ACFs (Filter: ' + str(period1) + '-' + str(period2) + ' s)')
plt.tick_params(which='both', bottom = 1, left = 1, right = 1, top =1)
axs3.xaxis.set_minor_locator(AutoMinorLocator())




# Plot Noise ACFs
axs2 = plt.subplot(4,1,3)
 
plt.plot(dist_t, T_PJIVSM*travel_nb3, ':', color = col2p3, linewidth = 2) # Plot JIVSM three times two-way travel time
plt.plot(distTs, T_P*3, color = Cestimate)  # Plot estimated threee times P-wave two-way travel time
plt.plot(distTs, T_P*3+v13, '--', color = Cestimate, linewidth = 2)  # Plot estimated threee times P-wave two-way travel time - 1.5 s       
plt.plot(distTs, T_P*3+v23, '--', color = Cestimate, linewidth = 2)   # Plot estimated threee times P-wave two-way travel time + 1.5 s        

factval = 3
for i in np.arange(len(dat_fin2)):
    if np.mean(dat_fin2[i,:50])==0:
        dat_fin2[i] = np.nan
    y = dist_t[i] + dat_fin2[i]*factval 
    if not np.isnan(np.sum(y)):
        plt.plot(y, t, 'k', linewidth = .5)
        axs2.fill_betweenx(t, dist_t[i], y, where=(y < dist_t[i]), color='gray')
        axs2.fill_betweenx(t, dist_t[i], y, where=(y > dist_t[i]), color='k')
        disttick = np.zeros((len(mxpk3all[i])))+dist_t[i]     
    if not np.isnan(mxpk3[i]):
        ymin = np.min(dat_fin2[i, int(mxpk3[i]*travel_nb3*delta)]*factval)
        plt.plot(dist_t[i]+nnb/2+ymin, mxpk3[i]*travel_nb3, 'ro', mec ='k', markersize=6)

if  Line_nb== '1':
     plt.plot(dist_t, T_PJIVSM, ':', color = col2p, linewidth =2)
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth =2)
     plt.text(124, .25, '2p',  weight="bold", color = col2p, fontsize = 12)
     plt.text(124, .9, '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)
     plt.text(124, 1.6, '$\mathbf{2p^3}$', weight = "bold", color = col2p3, fontsize = 12)
     plt.text(124, 3.2, 'Est.', weight = "bold", color = Cestimate, fontsize = 12)
     plt.text(-18, -.5, '(c)', fontsize = 18)
elif  Line_nb== '2':
     plt.plot(dist_t, T_PJIVSM, ':', color = col2p, linewidth =2)
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth =2)
     plt.text(96.4, 2.6, '2p',  weight="bold", color = col2p, fontsize = 12)
     plt.text(96.4, 5.1, '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)
     plt.text(96.4, 7.5, '$\mathbf{2p^3}$', weight = "bold", color = col2p3, fontsize = 12)
     plt.text(96.4, 6.7, 'Est.', weight = "bold", color = Cestimate, fontsize = 12)
     plt.text(-16 , -.5, '(g)', fontsize = 18) 
elif  Line_nb== '3':
     plt.plot(dist_t, T_PJIVSM, ':', color = col2p, linewidth =2)
     plt.plot(dist_t, T_PJIVSM*2, ':', color = col2p2, linewidth =2)
     plt.text(109.5, 1.3, '$\mathbf{2p^3}$', weight="bold", color = col2p3, fontsize = 12)
     plt.text(109.5, 0.1, '2p',  weight="bold", color = col2p, fontsize = 12)
     plt.text(109.5, .65, '$\mathbf{2p^2}$',  weight = "bold", color = col2p2, fontsize = 12)
     plt.text(109.5, 2, 'Est.', weight = "bold", color = Cestimate, fontsize = 12)
     plt.text(-19, -.5, '(k)', fontsize = 18)

plt.tick_params(which = 'both', bottom = 1, left = 1, right = 1, top = 1)
plt.ylim(10, 0)
plt.xlim(dist_t[0], dist_t[-1])
plt.ylabel('Time (s)')
plt.title('Vertical noise ACFs (Filter: '+ str(period1)+ '-' + str(period2) + ' s)')
axs2.xaxis.set_minor_locator(AutoMinorLocator())




# Plot P-wave two-way travel time values
axs4 = plt.subplot(4,1,4)
plt.plot(dist_t, T_PJIVSM, color = col2p, linewidth = 3, label = 'JIVSM 2p')
plt.plot(dist_plot+nnb/2,mxpk3, 'ro', label = 'Noise ACFs', mec ='k', markersize = 12, alpha = .75)
plt.plot(dist_plot+nnb/2,mxpk3EQ, 'go', label = 'Earthquake ACFs', mec ='k', markersize = 12, alpha = .75)
plt.plot(dist_t, T_PJIVSM, color = col2p, linewidth = 3)

if  Line_nb == '1':
    plt.text(-18 , -.4, '(d)', fontsize = 18)   
    plt.legend(loc = 4)
elif  Line_nb == '2':
    plt.text(-16 , -.4, '(h)', fontsize = 18)
    plt.legend(loc = 3)
elif Line_nb == '3':
    plt.text(-18 , -.4, '(l)', fontsize = 18)
    plt.legend(loc = 4)

plt.xlim(dist_t[0], dist_t[-1])
plt.ylim(3.5, -0.2)
plt.xlabel('Distance (km)')
plt.title('P-wave two-way travel time')
plt.ylabel('Time (s)')
plt.tick_params(which='both',bottom = 1, left = 1, right = 1, top =1)
minor_ticks = np.arange(0, 4, .5)
axs4.set_yticks(minor_ticks)
axs4.xaxis.set_minor_locator(AutoMinorLocator())
plt.grid()


# Set subplot positions
pos1 = axs1.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0 + 0.09, pos1.width, pos1.height-0.0] 
axs1.set_position(pos2) # set a new position

pos1 = axs3.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0 + 0.0275, pos1.width, pos1.height+0.045] 
axs3.set_position(pos2) # set a new position

pos1 = axs2.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0 + -0.035, pos1.width, pos1.height+0.045] 
axs2.set_position(pos2) # set a new position

pos1 = axs4.get_position() # get the original position 
pos2 = [pos1.x0, pos1.y0 + -0.075, pos1.width, pos1.height+0.02] 
axs4.set_position(pos2) # set a new position


figo.savefig('../Figures/Figure_3_Line_' + Line_nb + '.png', dpi = 100)

