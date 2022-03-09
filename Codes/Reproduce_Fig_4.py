#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 11:13:47 2021

@author: Loic Viens
Code to reproduce Figure 4 of Viens et al. (2022, GJI).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This code requires the basemap package and requires Matplotlib 3.1 (or below)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
Vpth =  2.537*1000 # Vp to migrate ACFs to depth (in m/s)

#%% Load Noise and Earthquake ACF travel times
input_f = directory + 'Noise_EQ_ACFs_Travel_times.mat'
dat = sio.loadmat(input_f)
stationN = np.squeeze(dat['station'])  # Station names
lat = np.squeeze(dat['lat']) # Station latitude
lon = np.squeeze(dat['lon']) # Station longitude
depthbas =  np.squeeze(dat['depthbas'])*1000 # JIVSM basin depth
dataPN =  np.squeeze(dat['mxpkN']) # Load Noise P-wave two-waytravel time
dataP = np.squeeze(dat['mxpk']) # Load EQ P-wave two-way travel time
model = 'JIVSM'
#%%
print('Noise: # of stations removed: ' + str(len(np.squeeze(np.where(np.isnan(dataPN)))) ))
print('Earthquake: # of stations removed: ' + str(len(np.squeeze(np.where(np.isnan(dataP)))) ))

#%% Set the map parameters
north = 36.3
south = 35.1
west = 138.95
east = 141
m = Basemap(projection = 'merc', lat_0 = south, lon_0 = west, llcrnrlat = south, urcrnrlat = north, llcrnrlon = west, urcrnrlon = east, resolution = 'f', epsg = 3832)

#%% Colorbar parameters
# maximum value for the colorbar in Figure 4a-c
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
cmap_reversed = matplotlib.cm.get_cmap('hot_r')
vminval = 0
cmaxval = 4100

#%% Set the map parameters
north = 36.3
south = 35.1
west = 138.95
east = 141
m = Basemap(projection='merc', lat_0 = south, lon_0 = west, llcrnrlat=south,urcrnrlat=north, llcrnrlon=west,urcrnrlon=east, resolution='f', epsg=3832)



#%% Plot FIGURE 4
fig1 = plt.figure(figsize=(11,14))

ax2 = fig1.add_subplot(321)
m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000) #ESRI_Imagery_World_2D
x, y = m(lon, lat)

lon_ticks_proj, _=m(np.arange(139,141.5,.5), np.zeros(len(np.arange(139,141.5,.5))))
_, lat_ticks_proj=m(np.zeros(len(np.arange(35.25,36.5,.25))), np.arange(35.25,36.5,.25))


m.scatter(x,y, c=depthbas ,marker='o', vmin = vminval, vmax = cmaxval, s = 50,cmap=cmap_reversed)
plt.title( model   ,fontsize =13)

m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle='fancy')

parallels = np.arange(30,40,.25)
m.drawparallels(parallels,labels=[True,False,False,False], linewidth=.2,fontsize =12)
meridians = np.arange(135,142,.5)
m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=.2,fontsize =12)
plt.annotate('(a)', xy=(-.03, 1.03), xycoords='axes fraction', fontsize=14)

cbaxes = fig1.add_axes([0.82, 0.8, 0.04, 0.1]) 
cb = plt.colorbar( cax = cbaxes,orientation = 'vertical', ticks = np.arange(0,5000, 1000) ) #[vminval2, vminval2+thir , cmaxval2- thir ,cmaxval2 ]
cb.ax.invert_yaxis()

cb.set_label('Meters',fontsize =12)
cb.ax.set_title('Bedrock depth',fontsize =12, pad=10)
cb.ax.tick_params(labelsize=12)
ax2.set_xticks(lon_ticks_proj)
ax2.set_yticks(lat_ticks_proj)
ax2.tick_params(axis='both',which='major')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticklabels([])
ax2.yaxis.set_ticklabels([])



ax3 = fig1.add_subplot(324)
m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000) #ESRI_Imagery_World_2D
x, y = m(lon, lat)
m.scatter(x,y, c=(dataPN)*Vpth/2,marker='o', vmin = vminval, vmax = cmaxval, s = 50,cmap=cmap_reversed) #, alpha = 1,edgecolors= 'k',linewidths = .1

plt.title( 'Noise ACFs',fontsize =13 )
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle='fancy')

parallels =  np.arange(30.00,40.00,.25) 
#m.drawparallelsL(parallels,labels=[False,True,False,False], linewidth=.2,fontsize =12,centerpos='center')
m.drawparallels(parallels,labels=[False,False,False,False], linewidth=.2,fontsize =12)

meridians = np.arange(135,142,.5)
m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=.2,fontsize =12)
plt.annotate('(c)', xy=(-.03, 1.05), xycoords='axes fraction', fontsize=14)
ax3.set_xticks(lon_ticks_proj)
ax3.set_yticks(lat_ticks_proj)
ax3.tick_params(axis='both',which='major')
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticklabels([])
ax3.yaxis.set_ticklabels([])


ax4 = fig1.add_subplot(323)
m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000) #ESRI_Imagery_World_2D

x, y = m(lon, lat)
m.scatter(x,y, c=((dataP)*Vpth/2),marker='o', vmin = vminval, vmax = cmaxval, s = 50,cmap=cmap_reversed) #, alpha = 1,edgecolors= 'k',linewidths = .1
plt.title('Earthquake ACFs',fontsize =13 )
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle='fancy')
parallels = np.arange(30.00,40.00,.25)
m.drawparallels(parallels,labels=[False,True,False,False], linewidth=.2,fontsize =12)
meridians = np.arange(135,142,.5)
m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=.2,fontsize =12)
plt.annotate('(b)', xy=(-.03, 1.05), xycoords='axes fraction', fontsize=14)

ax4.set_xticks(lon_ticks_proj)
ax4.set_yticks(lat_ticks_proj)
ax4.tick_params(axis='both',which='major')
ax4.xaxis.set_ticks_position('both')
ax4.yaxis.set_ticks_position('both')
ax4.xaxis.set_ticklabels([])
ax4.yaxis.set_ticklabels([])




colrnb = 8
vminvaldif = -1000
cmaxvaldif =-vminvaldif
thir = (cmaxvaldif - vminvaldif)/colrnb

ax5 = fig1.add_subplot(326)
m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000) #ESRI_Imagery_World_2D

x, y = m(lon, lat)
m.scatter(x,y, c=(depthbas-((dataPN)*Vpth/2)),marker='o', vmin = vminvaldif, vmax = cmaxvaldif, s = 50,cmap=Seismic2, alpha = 1,edgecolors= 'k',linewidths = .1) #depthbas depthbasJSHIS #discrete_cmap(colrnb, 'seismic')
plt.title( model +' minus noise ACF bedrock depth \n $\mu = $' + str(int(np.round(np.nanmean(depthbas-((dataPN)*Vpth/2)),0)))  +' m, $\sigma = $' + str(int(np.round(np.nanstd(depthbas-((dataPN)*Vpth/2)),0)))+ ' m', fontsize =13 )

m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle='fancy')
parallels = np.arange(30.00,40.00,.25)
#m.drawparallelsL(parallels,labels=[False,True,False,False], linewidth=.2,fontsize =12, centerpos='center')
m.drawparallels(parallels,labels=[False,False,False,False], linewidth=.2,fontsize =12)

meridians = np.arange(135,142,.5)
m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=.2,fontsize =12)
plt.annotate('(e)', xy=(-.03, 1.05), xycoords='axes fraction', fontsize=14)

ax5.set_xticks(lon_ticks_proj)
ax5.set_yticks(lat_ticks_proj)
ax5.tick_params(axis='both',which='major')
ax5.xaxis.set_ticks_position('both')
ax5.yaxis.set_ticks_position('both')
ax5.xaxis.set_ticklabels([])
ax5.yaxis.set_ticklabels([])


ax6 = fig1.add_subplot(325)
m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000) #ESRI_Imagery_World_2D

x, y = m(lon, lat)
m.scatter(x,y, c=(depthbas-((dataP)*Vpth/2)),marker='o', vmin = vminvaldif, vmax = cmaxvaldif, s = 50,cmap=Seismic2, alpha = 1,edgecolors= 'k',linewidths = .1) #depthbas depthbasJSHIS discrete_cmap(colrnb, 'seismic')

plt.title( model +' minus earthquake ACF bedrock depth \n $\mu = $' + str(int(np.round(np.nanmean(depthbas-((dataP)*Vpth/2)),0)))  +' m, $\sigma = $' + str(int(np.round(np.nanstd(depthbas-((dataP)*Vpth/2)),0))) + ' m',fontsize =13 )
m.drawmapscale(140.75, 35.35, 139, 36, 25, barstyle='fancy')

parallels = np.arange(30.00,40.00,.25)
m.drawparallels(parallels,labels=[False,True,False,False], linewidth=.2,fontsize =12)
meridians = np.arange(135,142,.5)
m.drawmeridians(meridians,labels=[False,False,False,True], linewidth=.2,fontsize =12)
plt.annotate('(d)', xy=(-.03, 1.05), xycoords='axes fraction', fontsize=14)

cbaxes = fig1.add_axes([0.29, 0.04, 0.4, 0.02])
cb2 = plt.colorbar( cax = cbaxes,orientation = 'horizontal' ,ticks = np.linspace(vminvaldif,cmaxvaldif ,colrnb+1),drawedges =True ) #ticks =[-1.5 , -.9, -.3  , .3, .9 , 1.5]
cb2.set_label(r'$\leftarrow$' + 'Thicker than JIVSM                Thinner than JIVSM ' +r'$\rightarrow$', fontsize =12)
cb2.ax.set_title('Bedrock depth difference (Meters)',fontsize =12)
cb2.ax.tick_params(labelsize=12)


ax6.set_xticks(lon_ticks_proj)
ax6.set_yticks(lat_ticks_proj)
ax6.tick_params(axis='both',which='major')
ax6.xaxis.set_ticks_position('both')
ax6.yaxis.set_ticks_position('both')
ax6.xaxis.set_ticklabels([])
ax6.yaxis.set_ticklabels([])


#@%%
ax2.set_position([0.28, .63, .44, .44])
ax4.set_position([0.02, .325, .44, .44])
ax3.set_position([0.535, .325, .44, .44])
ax6.set_position([0.02, .0075, .44, .44])
ax5.set_position([0.535, .0075, .44, .44])

fig1.savefig('../Figures/Fig_3.png', dpi = 100)
