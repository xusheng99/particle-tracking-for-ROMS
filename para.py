import os
import shutil
import sys
import numpy as np
from netCDF4 import Dataset
import math
import random
import time
from multiprocessing import Process,cpu_count
import functools
import platform
import socket
import matplotlib.pyplot as plt
from scipy.interpolate.interpolate import interp1d


######################################  README  ########################################
# this is a 2 dimensional particle tracking model based on ROMS output files
# if a history file is given to drive the model, the surface sigma layer data will be extracted
# this model release a batch of particles all at once, then record their tracks to several files:
# including the *.dat file, one file for one timestep, and a *.npy file record all timesteps.

######################################  PARAS   ########################################

# workdir = 'c:/Users/XUSHENG/Desktop/particle_tracking_results/fast/test/' # path that stores particle tracking results
workdir = '/mnt/i/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow/'
# workdir = 'i:/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow/'

# workdir = './'
os.chdir(workdir) # working directory

## project name, write whatever you want
projectname = 'for bohai community detection' 

## specify revelant files
# domain = 'i:/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow//Grd.nc' # ROMS Grid file name
# datadir = 'i:/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow/' # 
domain = '/mnt/i/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow//Grd.nc' # ROMS Grid file name
datadir = '/mnt/i/application/bh/0910test3/Christmas_new_nowind_TSfixed_nobryflow/' # 

## restart?
L_RST = False # if model restart from a restartfile, True / False
hotfile = '000188.csv' # the restart file

## ROMS output file parameters
his_or_qck = 'quicksave' # drive the particle tracking model using history file(vertical layer) or quicksave file(surface-only)?
if his_or_qck == 'history':
    prefix = 'his_' # the prefix of history files; eg. 'his_0001.nc','his_0002.nc',...,'his_0365.nc' -> prefix = 'his_'
else:
    prefix = 'qck_' # the prefix of quicksave files;
numdigits = 4 # the length of the numbers in a history or quicksave file; eg. 'his_0234.nc' -> numdigits = 4
his_writestep = 24 # how many timesteps are there in a normal history or quicksave file? do remember the timestep numbers in the first outputfile contains the extra intital time step;
                   # those infomations can be accessed by "ncdump -h qkc_0002.nc" on linux
ndefhis = 86400 # how many seconds are there in a normal history or quicksave file? 
DT = 3600 # the time gap of two timesteps in history or quicksave file, unit: seconds


## integrate parameters
timearrow = 'forward' # time flow forwards or backwards
# timearrow = 'backward'
T_cycle = 0 # reset all particle's location every T_cycle time, if T_cycle equals 0, the reset behaviour is disabled.
dt = 600  # integration time step, seconds
release_time = 86400*30  # release time, seconds since roms initialization
total_time = 86400*30 # total tracking time, seconds
write_step = 6 # take a snapshot of all particles to a *.dat file for every write_step*dt seconds, 

## interpolation method: horizontal
# 'IDW': Inverse-Distance-Weighted interpolation, much more faster if the amount of particles is less than 1e5
# 'BATCH_LINEAR': Bi-Linear interpolation, T(n) is about log(n), faster for massive particles
# 'ANA': calculate the flow field in analytical expression instead of reading it from nc files
h_interp_method = 'IDW' # 'IDW' / 'BATCH_LINEAR' / 'ANA'

## land-ocean mask
# 0: recommended, the land-ocean mask is created based on the "mask_rho" in the ROMS grid file
# 2: accurate but extreme slow, based on a closed polygon(coastline.bln)
# 3: analytical method
oceanmask_method = 0 # 0/1/2/3

## the coastline polygonfile, change only when oceanmask_method == 2
coast = '/mnt/i/application/0DATASRC/bh_coastline/bh_google.bln'

## the analytical method to determine land-ocean mask, change only when oceanmask_mask == 3
lon1 = 105.096875
lon2 = 105.103125
lat1 = 0.080625
lat2 = 0.1

## particel behaviour
diffusion = 'off' # if enable particle diffussion(random walk)?
Dh = 0.3 # horizontal diffusion coefficient




## particle release region control
minX = 117 # min longitude
maxX = 127 # max longitude
minY = 35.5 # min latitude
maxY = 41 # max latitude
N_edge = 50 # linear density, how many particles on a one degree segment?










######################################################################################
######################################################################################
######################################################################################

if not L_RST:
    x = np.linspace(minX,maxX,int(N_edge*(maxX-minX)+1))
    y = np.linspace(minY,maxY,int(N_edge*(maxY-minY)+1))   
    [lats,lons] = np.meshgrid(y,x)
    lon = np.squeeze(np.reshape(lons,[(lons.shape[0])*(lons.shape[1]),1]))
    lat = np.squeeze(np.reshape(lats,[(lats.shape[0])*(lats.shape[1]),1]))
    time_ = release_time
    positions = [lon,lat,time_]
else:
    x = np.linspace(minX,maxX,int(N_edge*(maxX-minX)+1))
    y = np.linspace(minY,maxY,int(N_edge*(maxY-minY)+1))   
    [lats,lons] = np.meshgrid(y,x)
    aa = np.loadtxt(hotfile,delimiter=',')
    lon = np.squeeze(np.reshape(np.array(aa[:,0]),[(lons.shape[0])*(lons.shape[1]),1]))
    lat = np.squeeze(np.reshape(np.array(aa[:,1]),[(lons.shape[0])*(lons.shape[1]),1]))
    time_ = aa[0,2]
    positions = [lon,lat,time_]
positions_0 = positions.copy()

if timearrow == 'backward':
    dt = -1 * dt

if T_cycle != 0:
    write_step = T_cycle//dt


## grid 
content2 = Dataset(domain)
mask_rho = content2.variables['mask_rho'][:].data
lon_rho = content2.variables['lon_rho'][:].data
lat_rho = content2.variables['lat_rho'][:].data
content2.close()
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) # grid size
N = int(N_edge*(maxX-minX)+1)*int(N_edge*(maxY-minY)+1) # particle num

from multiprocessing import cpu_count
cpu_count = cpu_count() # read cpu count
timestep_numbers = int(total_time/dt) # how many timesteps to run? 

## device name
host_name = socket.gethostname()
