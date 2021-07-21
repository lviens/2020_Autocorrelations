#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 10:23:26 2021

@author: Loic Viens
Codes to reproduce the maps of Figure 1 (Viens et al., Submitted to GRL)
The Figure is saved in the "Figures" folder.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This code requires the basemap package and requires Matplotlib 3.1 (or below)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import os
from zipfile import ZipFile


#%% Used to reshape the file with the Kanto Basin latitude, longitude, and depth info
def getIndexPositions(listOfElements, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    indexPosList = []
    indexPos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            indexPos = listOfElements.index(element, indexPos)
            # Add the index position in list
            indexPosList.append(indexPos)
            indexPos += 1
        except ValueError:
            break
    return indexPosList

#%%
dir_ini = './'
sys.path.append(dir_ini)
#%% Directory with station info, Kanto basin depth, earthquake lists,...
directory = '../Data/Other_data/'
if not os.path.exists(directory):
    zf = ZipFile('../Data/Other_data.zip', 'r')
    zf.extractall('../Data/')
    zf.close()

#%% Load station names, latitude, and longitude
crs = open(directory + "/MeSOnet_channel_file", "r")
nm_sta = []
latstaini = []
lonstaini =[]
for columns in ( raw.strip().split() for raw in crs ):  
     if (columns[4]) =='U':
         nm_sta.append(columns[3])
         latstaini.append(columns[13])
         lonstaini.append(columns[14])
crs.close()

#%% Load JIVSM and reshape the file
crs = open(directory + "/Kanto_basin_shape.txt", "r")
latg = []
long =[]
lay3 =[]
for columns in ( raw.strip().split() for raw in crs ):  
    long.append(columns[0])
    latg.append(columns[1])
    lay3.append(columns[15])
crs.close()

# Reshape the data
lay3 = np.array(lay3,dtype=float)
indo = getIndexPositions(long, long[0])
for index, item in enumerate(long):
    long[index] = float(item)
for index, item in enumerate(latg):
    latg[index] = float(item)
    
long1 = np.reshape(long, (int(len(long)/indo[1]),indo[1] ) ).T
latg1 = np.reshape(latg, (int(len(latg)/indo[1]),indo[1] ) ).T
lay31 = np.reshape(lay3, (int(len(lay3)/indo[1]),indo[1] ) ).T/-1000

#%% Load earthquake lists
# All events between May 2017 to April 2020
crs = open(directory + "events_all.txt", "r")
latgeqall = []
longeqall =[]
for columns in ( raw.strip().split(',') for raw in crs ):  
    longeqall.append(np.float(columns[2]) )
    latgeqall.append(np.float(columns[1]) )
crs.close()

# Selected earthquakes
crs = open(directory + "events_used.txt", "r")
latgeq = []
longeq = []
for columns in (raw.strip().split(',') for raw in crs):  
    longeq.append(np.float(columns[2]) )
    latgeq.append(np.float(columns[1]) )
crs.close()

#%% Load Plate boundaries (Pacific and Philippine sea Plate)
crs1 = open(directory +'ps', 'r')
lt = []
ln = []
for columns in ( raw.strip().split() for raw in crs1 ):  
    lt.append(np.float(columns[0]))
    ln.append(np.float(columns[1]))
crs1.close()

crs2 = open(directory + 'pa', 'r')
ln2 = []
lt2 = []
for columns in (raw.strip().split() for raw in crs2):  
    lt2.append(np.float(columns[0]))
    ln2.append(np.float(columns[1]))
crs2.close()
    
#%%  Stations used in this study and get their locations      
     
all_sta = ['E.DSCM', 'E.AM1M', 'E.ANSM', 'E.AZMM', 'E.BKKM', 'E.BTSM', 'E.CF3M', 'E.DAIM', 'E.DGRM', 'E.DICM', 'E.ENZM', 'E.F3KM', 'E.FJSM', 'E.FKCM', 'E.FKMM', 'E.FKOM', 'E.FMNM', 'E.FRAM', 'E.FSMM', 'E.FTKM', 'E.FTPM', 'E.FUNM', 'E.GHGM', 'E.GKSM', 'E.GNZM', 'E.GSJM', 'E.HGCM', 'E.HKBM', 'E.HNOM', 'E.HNPM', 'E.HNSM', 'E.HONM', 'E.HRGM', 'E.HRMM', 'E.HSBM', 'E.HSDM', 'E.HTTM', 'E.HYDM', 'E.HYHM', 'E.HYSM', 'E.ICEM', 'E.IDUM', 'E.IIDM', 'E.IKBM', 'E.IKCM', 'E.IMIM', 'E.IN3M', 'E.INAM', 'E.INOM', 'E.IWNM', 'E.JDJM', 'E.JKPM', 'E.JNUM', 'E.JUWM', 'E.JYHM', 'E.KBAM', 'E.KBRM', 'E.KBTM', 'E.KCBM', 'E.KD9M', 'E.KDKM', 'E.KGKM', 'E.KH2M', 'E.KHBM', 'E.KHDM', 'E.KIMM', 'E.KK2M', 'E.KKHM', 'E.KKSM', 'E.KM5M', 'E.KMGM', 'E.KMHM', 'E.KMKM', 'E.KMNM', 'E.KMRM', 'E.KMTM', 'E.KNDM', 'E.KNJM', 'E.KNMM', 'E.KNZM', 'E.KOHM', 'E.KOYM', 'E.KRCM', 'E.KRPM', 'E.KSGM', 'E.KSHM', 'E.KSNM', 'E.KSOM', 'E.KSRM', 'E.KTOM', 'E.KUDM', 'E.KUYM', 'E.KW8M', 'E.KWHM', 'E.KWWM', 'E.KYNM', 'E.KYTM', 'E.KZ2M', 'E.KZMM', 'E.KZTM', 'E.MBSM', 'E.MBYM', 'E.MD1M', 'E.MDHM', 'E.MGRM', 'E.MHKM', 'E.MKBM', 'E.MKJM', 'E.MKSM', 'E.MNAM', 'E.MNHM', 'E.MNKM', 'E.MNMM', 'E.MNYM', 'E.MOKM', 'E.MRTM', 'E.MSKM', 'E.MSOM', 'E.MSRM', 'E.MYHM', 'E.MYMM', 'E.MZMM', 'E.MZPM', 'E.MZUM', 'E.NARM', 'E.NBKM', 'E.NDOM', 'E.NGOM', 'E.NGSM', 'E.NISM', 'E.NKDM', 'E.NKGM', 'E.NKMM', 'E.NKNM', 'E.NKSM', 'E.NMDM', 'E.NNBM', 'E.NNGM', 'E.NNSM', 'E.NNTM', 'E.NOBM', 'E.NS7M', 'E.NSHM', 'E.NSJM', 'E.NSKM', 'E.NSMM', 'E.NSOM', 'E.NSUM', 'E.NTNM', 'E.OA5M', 'E.OANM', 'E.OBRM', 'E.OBTM', 'E.OG6M', 'E.OHSM', 'E.OJCM', 'E.OKCM', 'E.OKDM', 'E.OMKM', 'E.OMNM', 'E.OMSM', 'E.ONMM', 'E.OOYM', 'E.OYMM', 'E.OYOM', 'E.OYTM', 'E.RKGM', 'E.RMSM', 'E.RYGM', 'E.RYNM', 'E.SBAM', 'E.SDMM', 'E.SEBM', 'E.SECM', 'E.SFHM', 'E.SGOM', 'E.SGWM', 'E.SICM', 'E.SJSM', 'E.SKCM', 'E.SKEM', 'E.SKHM', 'E.SKPM', 'E.SKYM', 'E.SMGM', 'E.SMNM', 'E.SNHM', 'E.SNJM', 'E.SNSM', 'E.SONM', 'E.SR2M', 'E.SRCM', 'E.SRGM', 'E.SRIM', 'E.SRKM', 'E.SRSM', 'E.SRTM', 'E.SSDM', 'E.SSHM', 'E.SSMM', 'E.SSNM', 'E.SSPM', 'E.SSSM', 'E.ST2M', 'E.STHM', 'E.STKM', 'E.SW2M', 'E.SYOM', 'E.SYPM', 'E.SYSM', 'E.TAKM', 'E.TBKM', 'E.TGNM', 'E.TK2M', 'E.TKKM', 'E.TKMM', 'E.TKNM', 'E.TKSM', 'E.TKWM', 'E.TKZM', 'E.TMHM', 'E.TMOM', 'E.TNKM', 'E.TNRM', 'E.TOCM', 'E.TOKM', 'E.TREM', 'E.TRSM', 'E.TSCM', 'E.TSRM', 'E.TTNM', 'E.TTOM', 'E.TWDM', 'E.TYHM', 'E.TYNM', 'E.TYPM', 'E.TYUM', 'E.U1CM', 'E.UHRM', 'E.UNHM', 'E.UNMM', 'E.URMM', 'E.USCM', 'E.UTKM', 'E.YGHM', 'E.YKBM', 'E.YKDM', 'E.YKKM', 'E.YKSM', 'E.YMKM', 'E.YMMM', 'E.YMPM', 'E.YNDM', 'E.YNMM', 'E.YROM', 'E.YRPM', 'E.YSKM', 'E.YSOM', 'E.YSPM', 'E.YSSM', 'E.YT2M', 'E.YTBM', 'E.YUKM', 'E.YYIM', 'E.KYDM', 'E.AYHM', 'E.KSCM' ,'E.ABHM', 'E.MRJM' ,'E.THCM', 'E.FMHM', 'E.OHYM', 'E.SBCM' ,'E.TACM' ,'E.HSMM' ,'E.SIBM', 'OK.AOCM', 'OK.AONM', 'OK.ARMM', 'OK.HRDM', 'OK.KRHM', 'OK.KTGM', 'OK.NHMM', 'OK.NKYM', 'OK.NRAM', 'OK.TKCM']

lonsta_all =[]
latsta_all=[]
for  i in np.arange(len(all_sta)):
    indr = nm_sta.index(all_sta[i])
    latsta_all.append(float(latstaini[indr]))
    lonsta_all.append(float(lonstaini[indr]))
    
#%% Set up map parameters
north = 36.25
south = 35.2
west = 138.975
east = 141
m = Basemap(projection = 'merc', lat_0 = south, lon_0 = west, llcrnrlat = south, urcrnrlat = north, llcrnrlon = west, urcrnrlon = east, resolution = 'f', epsg = 3832)

#%% Set up inset map parameters
north2 = 46
south2 = 30
west2 = 127
east2 = 147
m2 = Basemap(projection='merc', lat_0 = south2, lon_0 = west2, llcrnrlat = south2,urcrnrlat = north2, llcrnrlon =  west2, urcrnrlon = east2, resolution = 'i', epsg = 3832)

#%% Plot Figure 1a
fig1 = plt.figure(figsize = (8,9))
ax0 = fig1.add_subplot(311)

m.arcgisimage(service='World_Shaded_Relief', xpixels = 2000)

x, y = m(lonsta_all, latsta_all)# Plot all stations
m.scatter(x,y,marker='o',s = 40, edgecolors ='k', alpha=1, facecolors='none',linewidths =1)
m.drawmapscale(139.15, 35.85, 139, 36, 25, barstyle='fancy') # Plot scale

# Loop over the 3 lines to plot the stations
for lnb in range(1,4):
    if lnb ==1:
        station = ['E.THCM','E.DICM','E.TYHM','E.TBKM','E.SSHM','E.NKMM','E.TK2M','E.HSDM','E.SRTM','E.NKNM', 'E.TSCM','E.KUDM', 'E.HRGM', 'E.RKGM' ,'E.MNKM','E.FUNM', 'E.KMRM', 'E.NNTM', 'E.STHM','E.TYNM','E.KSCM','E.MDHM', 'E.HKBM','E.TKNM','E.SDMM', 'E.YNMM', 'E.NSMM','E.HGCM','E.HYHM','E.SRCM','E.DAIM','E.OKCM', 'E.FMHM','E.FMNM', 'E.NISM', 'E.FKCM', 'E.OKDM', 'E.KSRM', 'E.IIDM', 'E.SYPM', 'E.YSSM'] 
        az_mean = 100
        mxdist = 122
        col = 'green'
        postext = [139.25, 35.93]
    elif lnb ==2:
        station = ['E.SICM','E.MSKM','E.SSNM','E.YSOM','E.MSRM','E.YGHM','E.SSDM','E.INOM','E.YMMM', 'E.SYPM','E.KMRM', 'E.SR2M', 'E.TGNM', 'E.ABHM','E.NDOM','E.KKHM','E.TNKM', 'E.KUYM','E.MRJM','E.MYHM','E.SGOM', 'E.NNGM', 'E.NKGM','E.YSPM','E.NSUM', 'E.MRTM', 'E.SSMM', 'E.NOBM', 'E.SRSM', 'E.FKOM', 'E.MHKM', 'E.HYDM' ]
        az_mean = 142
        mxdist = 94
        col = 'red'
        postext = [139.7, 36.17]
    elif lnb ==3:
        station = ['E.SYOM','E.FJSM','E.DGRM','E.KMHM','E.SMGM','E.HTTM','E.KSOM','OK.HRDM','E.HNOM','E.BKKM','E.MZPM', 'E.JDJM', 'E.SNHM', 'E.UNMM','E.KMKM','E.MKSM','E.GKSM','E.KDKM','E.ENZM','E.TKMM', 'E.SBAM','E.GNZM','E.YKKM','E.MKJM', 'E.KCBM','E.MZMM', 'E.TKNM','E.MBSM', 'E.YKSM','E.KGKM', 'E.NGSM','E.KKSM', 'E.KWHM', 'E.TNKM','E.KUYM', 'E.IN3M','E.INAM', 'E.KBRM','E.GSJM', 'E.YTBM', 'E.HNPM','E.YKBM','E.TSRM', 'E.TKZM','E.RMSM','E.RYGM', 'E.RYNM']
        az_mean = 34.2
        mxdist = 107
        col = 'blue'
        postext = [140.05, 36.13]
        
    
    # Get lat and lon of the stations along the line
    stalat = []
    stalon = []
    for i in np.arange(len(station)):
        indr = nm_sta.index(station[i])
        stalat.append(float(latstaini[indr]))
        stalon.append(float(lonstaini[indr]))
      
    # Plot stations along each line
    xst, yst = m(stalon, stalat)
    m.scatter(xst, yst, marker='o', s = 40, edgecolors ='k', alpha=.5, linewidths = 1, facecolors = col, zorder = 101) 
    
    # Plot line number
    xp, yp = m(postext[0], postext[1])
    plt.text(xp, yp, 'Line ' + str(lnb),FontSize= 12, fontweight = 'bold')

#%% Set ticks
parallels = np.arange(33, 38, .25)
m.drawparallels(parallels, labels=[False,True,False,False], linewidth = .15, fontsize = 12)
meridians = np.arange(138, 142, .5)
m.drawmeridians(meridians, labels=[False,False,False,True], linewidth = .15, fontsize = 12)

# Plot contour of the Kanto Basin
xx, yy = m(long1, latg1)
m.contour( xx, yy, lay31, np.arange(0.5, 5,.25), linewidths = .8, cmap = plt.cm.hot_r, zorder = .1)
plt.annotate('(a)', xy = (-.065, .96), xycoords = 'axes fraction', fontsize = 14)

# White rectangle and Colorbar
inix = 140.688
x1,y1 = m(inix, 35.25)
x2,y2 = m(inix, 35.595)
x3,y3 = m(inix+.28, 35.595)
x4,y4 = m(inix+.28, 35.25)
poly = Polygon([(x1,y1), (x2,y2), (x3,y3), (x4,y4)], facecolor = 'white', edgecolor = 'none', alpha = .8)
ax0.add_patch(poly)

img = plt.imshow(np.array([[0.5, 5]]), cmap="hot_r")
img.set_visible(False)
cbaxes = fig1.add_axes([0.77, 0.565, 0.02, 0.09]) 
cb = plt.colorbar( cax = cbaxes, orientation = 'vertical', ticks = [1 ,2,3,4,5]) 
cb.ax.set_title('   Bedrock\n   depth', fontsize = 12)
cb.set_label('km', fontsize = 12)
cb.ax.set_yticklabels([1 ,2,3,4,5], fontsize = 12)


#%% Plot the inset map
ax1 = fig1.add_subplot(312)
m2.drawcoastlines()

parallels = np.arange(30, 50, 5)
m2.drawparallels(parallels, labels = [False,False,False,False], linewidth = .2, fontsize = 12)
meridians = np.arange(125, 150, 5)
m2.drawmeridians(meridians, labels = [False,False,False,False], linewidth = .2, fontsize = 12)

# Plot red rectangle around the Kanto Basin
w,o = m2((west ,west),(south, north))
m2.plot(w, o, 'r')
w,o = m2((east ,east),(south, north))
m2.plot(w, o, 'r')
w,o = m2((east ,west),(north, north))
m2.plot(w, o, 'r')
w,o = m2((east ,west),(south, south))
m2.plot(w, o, 'r')

# Plot Plate boundaries
xp, yp = m2(ln, lt)
m2.plot(xp, yp, 'grey' )
xp2, yp2 = m2(ln2, lt2)
m2.plot(xp2, yp2, 'grey' )
xp, yp = m2(135.4, 30.4)
plt.text(xp, yp, 'Philippine\nSea Plate',FontSize = 8)
xp, yp = m2(142.7, 33.6)
plt.text(xp, yp, 'Pacific\n Plate',FontSize = 8)
xp, yp = m2(135, 39)
plt.text(xp, yp, 'Japan', FontSize = 9)


#%% Subplot Figure 1b
ax2 = fig1.add_subplot(313)
# set up orthographic map projection with perspective of satellite looking down at 35.75N, 140 (low resolution coastlines).
center_lat = 35.75
center_lon = 140
width = 22000000
map = Basemap(projection='aeqd', lat_0 = center_lat,lon_0 = center_lon, resolution='l', width=width,height=width)
map.shadedrelief(scale = 0.1)
map.drawcoastlines(linewidth = 0.25)
plt.annotate('(b)', xy = (-0.12, .945), xycoords = 'axes fraction', fontsize=14)

# Plot earthquake locations
xp, yp = map(longeqall, latgeqall)
map.scatter(xp, yp, marker = 'o', s = 10,edgecolors = 'k', alpha=.5, linewidths = .5, facecolors = 'grey', zorder = .1)
xp, yp = map(longeq, latgeq)
map.scatter(xp, yp, marker = 'o', s = 15, edgecolors = 'k', alpha=1, linewidths = .5, facecolors = 'red', zorder = 1)

# plot triangle on the Kanto Basin
xp, yp = map(center_lon, center_lat)
map.scatter(xp, yp, marker = '^', s = 40,edgecolors = 'k', alpha = 1, linewidths = 1, facecolors = 'blue', zorder = 1)

# Plot 30 and 95 degree circles around the Kanto Basin
x,y = map(center_lon, center_lat)
x2,y2 = map(center_lon,center_lat - 30) 
circle1 = plt.Circle((x, y), y2-y, color='k', fill=False, zorder = 101)
plt.gca().add_patch(circle1)

x,y = map(center_lon, center_lat)
x2,y2 = map(center_lon, center_lat - 95) 
circle2 = plt.Circle((x, y), y2-y, color='k', fill=False, zorder = 101)
plt.gca().add_patch(circle2)

# Plot text 35 and 90 degree
xp, yp = map(176, 27)
plt.text(xp, yp, '$30^o$', fontsize=8 )
xp, yp = map(-138.5, 3.6)
plt.text(xp, yp, '$95^o$', fontsize=8 )



#%% Set position of axes, plot, and save
ax0.set_position([0.0775, .41, .78, .7])
ax1.set_position([0.735, .7765, .22, .22])
ax2.set_position([0.25, 0.01, .5, .5])
plt.show()

fig1.savefig('../Figures/Fig_1.png', dpi=300)
