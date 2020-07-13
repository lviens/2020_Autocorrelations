#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:13:47 2020

@author: Loic Viens
Code to reproduce Figure 4 of the paper.
"""

#%%
import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import scipy.io as sio
from matplotlib.colors import ListedColormap
import os
from zipfile import ZipFile


#%% Unzip data in needed

directory = '../Data/Travel_times/'
if not os.path.exists(directory):
    zf = ZipFile('../Data/Travel_times.zip', 'r')
    zf.extractall('../Data/')
    zf.close()

#%%
Vpth = 2.53 # Vp to migrate ACFs to depth (in km/s)


#%% Load Nois ACF travel time
outputf = directory + 'Noise_ACF_travel_times.mat'
datN = sio.loadmat(outputf)
stationN = np.squeeze(datN['station'])  # Station names
latN = np.squeeze(datN['lat']) # Station latitude
lonN = np.squeeze(datN['lon']) # Station longitude
T_PN = np.squeeze(datN['T_P']) # Station JIVSM two-way travel time
JIVSMbasdepth = np.squeeze(datN['depthbas']) # JIVSM basin depth
dataPN = np.squeeze(datN['mxpk']) # Load Noise P-wave two-waytravel time

#%% Load earthquake ACF
input_f = directory + 'EQ_ACF_travel_times.mat'
dat = sio.loadmat(input_f)
station = np.squeeze(dat['station']) # Station names (identical to noise ACFs)
lat = np.squeeze(dat['lat']) # Station latitude (identical to noise ACFs)
lon = np.squeeze(dat['lon']) # Station longitude (identical to noise ACFs)
T_P = np.squeeze(dat['T_P']) # Station JIVSM two-way travel time (identical to noise ACFs)
dataP = np.squeeze(dat['mxpk']) # Load EQ P-wave two-way travel time

#%%
print('Noise: Nb. of stations removed: ' + str(len(np.squeeze(np.where(np.isnan(dataPN)))) ))
print('Earthquake: Nb. of stations removed: ' + str(len(np.squeeze(np.where(np.isnan(dataP)))) ))


#%% Migrate data to depth 
Noisebas = dataPN * Vpth / 2
EQbas = dataP * Vpth / 2

#%% Set the map parameters
north = 36.3
south = 35.1
west = 138.95
east = 141
m = Basemap(projection = 'merc', lat_0 = south, lon_0 = west, llcrnrlat = south, urcrnrlat = north, llcrnrlon = west, urcrnrlon = east, resolution = 'f', epsg = 3832)

#%% Colorbar parameters
# maximum value for the colorbar in Figure 4a-c
vminval = 0
cmaxval = 4.1
cmap_reversed = matplotlib.cm.get_cmap('hot_r')

# maximum value for the colorbar in Figure 4d-e and define new colorbar
cmaxvaldif = 1.2
vminvaldif = -cmaxvaldif
newcolors = [[ 0.,          0.,          0.3,         1. ],
             [ 0.     ,     0.,          0.76117647,  1. ],
             [ 0.33333333,  0.33333333,  1.,          1. ],
             [ 1.,          0.99215686,  0.99215686,  1. ],
             [ 1.,          0.99215686,  0.99215686,  1. ],
             [ 1.,          0.33333333,  0.33333333,  1. ],
             [ 0.82941176,  0.   ,       0.,          1. ],
             [ 0.5,         0.   ,       0.,          1. ]]
Seismic2 = ListedColormap(newcolors, name = 'Seismic2')

#%% Plot Figure 4
fig1 = plt.figure(figsize=(11,14))

# Subplot (a): JIVSM bedrock depth
ax2 = fig1.add_subplot(321)
m.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000) 
x, y = m(lon, lat)
m.scatter(x, y, c = JIVSMbasdepth, marker = 'o', vmin = vminval, vmax = cmaxval, s = 50, cmap = cmap_reversed)
plt.title('JIVSM bedrock depth', fontsize = 13)
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle = 'fancy') # Plot scale
parallels = np.arange(30,40,.25)
m.drawparallels(parallels, labels = [True,False,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(135,142,.5)
m.drawmeridians(meridians, labels = [False,False,False,True], linewidth = .2, fontsize = 12)
plt.annotate('(a)', xy = (-.03, 1.03), xycoords = 'axes fraction', fontsize = 14)

# Basin depth colorbar
cbaxes = fig1.add_axes([0.82, 0.8, 0.04, 0.1]) 
cb = plt.colorbar( cax = cbaxes,orientation = 'vertical', ticks = np.arange(0, 6, 1) )
cb.ax.invert_yaxis()
cb.set_label('km', fontsize = 12)
cb.ax.set_title('Basin depth', fontsize = 12, pad = 10)
cb.ax.tick_params(labelsize = 12)



# Subplot (b):  Plot Noise ACF basin depth
ax3 = fig1.add_subplot(323)
m.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)
x, y = m(lon, lat)
m.scatter(x, y, c = Noisebas, marker = 'o', vmin = vminval, vmax = cmaxval, s = 50, cmap = cmap_reversed)
plt.title('Noise ACF bedrock depth', fontsize = 13)
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle = 'fancy') # Plot scale
parallels =  np.arange(30.00, 40.00, .25) 
m.drawparallels(parallels, labels = [False,True,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(135, 142,.5)
m.drawmeridians(meridians, labels = [False,False,False,True], linewidth = .2, fontsize = 12)
plt.annotate('(b)', xy = (-.03, 1.05), xycoords = 'axes fraction', fontsize = 14)



# Subplot (c): Plot Earthquake ACF basin depth
ax4 = fig1.add_subplot(324)
m.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)
x, y = m(lon, lat)
m.scatter(x, y, c = EQbas, marker = 'o', vmin = vminval, vmax = cmaxval, s = 50, cmap = cmap_reversed)
plt.title('Earthquake ACF bedrock depth', fontsize = 13 )
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle = 'fancy') # Plot scale
parallels = np.arange(30.00, 40.00, .25)
m.drawparallels(parallels, labels = [False,False,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(135, 142,.5)
m.drawmeridians(meridians, labels = [False,False,False,True], linewidth = .2, fontsize = 12)
plt.annotate('(c)', xy = (-.03, 1.05), xycoords = 'axes fraction', fontsize = 14)



# Subplot (d): Plot JIVSM minus noise ACF bedrock depth difference
ax5 = fig1.add_subplot(325)
m.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)
x, y = m(lon, lat)
m.scatter(x, y, c = JIVSMbasdepth - Noisebas, marker ='o', vmin = vminvaldif, vmax = cmaxvaldif, s = 50, cmap = Seismic2, alpha = 1,edgecolors = 'k', linewidths = .1) 
plt.title('JIVSM minus noise ACF bedrock depth \n $\mu = $' + str(np.round(np.nanmean(JIVSMbasdepth-Noisebas), 3)) + ' km, $\sigma = $' + str(np.round(np.nanstd(JIVSMbasdepth-Noisebas),3)) + ' km', fontsize = 13)
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle = 'fancy') # Plot scale
parallels = np.arange(30.00, 40.00, .25)
m.drawparallels(parallels, labels = [False,True,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(135, 142,.5)
m.drawmeridians(meridians, labels = [False,False,False,True], linewidth = .2, fontsize = 12)
plt.annotate('(d)', xy = (-.03, 1.05), xycoords = 'axes fraction', fontsize = 14)



# Subplot (e): Plot JIVSM minus EQ ACF bedrock depth difference
ax6 = fig1.add_subplot(326)
m.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)
x, y = m(lon, lat)
m.scatter(x,y, c=(JIVSMbasdepth-EQbas),marker='o', vmin = vminvaldif, vmax = cmaxvaldif, s = 50,cmap=Seismic2, alpha = 1,edgecolors= 'k',linewidths = .1) 
plt.title( 'JIVSM minus earthquake ACF bedrock depth \n $\mu = $' + str(np.round(np.nanmean(JIVSMbasdepth-EQbas),3)) + ' km, $\sigma = $' + str(np.round(np.nanstd(JIVSMbasdepth-EQbas),3)) + ' km', fontsize = 13 )
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle = 'fancy') # Plot scale
parallels = np.arange(30.00, 40.00, .25)
m.drawparallels(parallels, labels = [False,False,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(135, 142,.5)
m.drawmeridians(meridians, labels = [False,False,False,True], linewidth = .2, fontsize = 12)
plt.annotate('(e)', xy = (-.03, 1.05), xycoords = 'axes fraction', fontsize=14)

# Plot difference colorbar
cbaxes = fig1.add_axes([0.33, 0.04, 0.34, 0.02])
cb2 = plt.colorbar(cax = cbaxes, orientation = 'horizontal', ticks = np.linspace(vminvaldif, cmaxvaldif, 9), drawedges = True) 
cb2.set_label('km', fontsize = 12)
cb2.ax.set_title('Depth difference', fontsize = 12)
cb2.ax.tick_params(labelsize = 12)



#%% Set positions of all the subplots
ax2.set_position([0.27, .63, .44, .44])
ax3.set_position([0.02, .325, .44, .44])
ax4.set_position([0.53, .325, .44, .44])
ax5.set_position([0.02, .0075, .44, .44])
ax6.set_position([0.53, .0075, .44, .44])

fig1.savefig('../Figures/Fig_4.png', dpi = 100)
