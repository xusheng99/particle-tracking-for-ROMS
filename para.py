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




######################################  参数   ########################################

nprc = 4 # 想要用来运行的核数, 不得超过本机核数
workdir = 'c:/Users/XUSHENG/Desktop/particle_tracking_results/3dtest/' # 存放粒子追踪结果的文件夹
# workdir = '/mnt/c/Users/XUSHENG/Desktop/particle_tracking_results/fast/'

os.chdir(workdir) # 工作目录

## 项目名称, 仅仅是一个记号
projectname = 'SCS 3D TEST1' 

## 指定各种文件

# domain = 'D:/application/bh/roughs/rough_clm_sst_sss_test/Grd.nc' # ROMS Grid文件
# datadir = 'D:/application/bh/roughs/rough_clm_sst_sss_test/' # 存放ROMS输出流场数据的目录
# domain = '/mnt/d/application/idealised/test6-boundaryflow-highres/Grd.nc' # ROMS Grid文件
# datadir = '/mnt/d/application/idealised/test6-boundaryflow-highres/' # 存放ROMS输出流场数据的目录
domain = 'D:/application/scs_2/DATA/Grd.nc' # ROMS Grid文件
datadir = 'D:/application/scs_2/output/' # 存放ROMS输出流场数据的目录

## 关于水动力模型输出数据
prefix = 'his_' # history文件的前缀
numdigits = 4 # history或quicksave文件的数字编号的位数,含零
his_writestep = 4 # 非初始的history/quicksave文件中的时刻数目, 单位个; 注意, 初始的his/qck_0001.nc比后续的正常的文件多一个时刻. 
ndefhis = 86400 # ROMS每隔多少秒另起一个文件来记录history/quicksave数据
DT = 21600  # ROMS文件中两个相邻时刻的时间差，单位秒

## 粒子运动行为控制
dt = 10800  # 粒子运动步长, 单位秒, 最好远小于DT(ROMS输出时间间隔)
release_time = 86400*367  # seconds since initialization
total_time = 86400*30*3 # 总的模拟时间, 单位秒, 86400秒为一天
write_step = 2 # 每计算多少个步长写入一次数据
release_dep = -30 # 释放时刻粒子深度，负值，单位meters
flowtype = 'uv' # 填'rho'/'uv'

diffusion = 'off' # 是否开启random walk? 'on'/'off'
random_speed_sigma = 0.06 # 随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）

## 判定海陆掩码的岸线文件
# coast = '/mnt/c/Users/XUSHENG/Desktop/roms/track/bh_topo/bh.bln' # 岸线文件
coast = 'D:/application/0DATASRC/bh_coastline/scs.bln'

## 理想实验中海岬的范围
# lon1 = 105.096875
# lon2 = 105.103125
# lat1 = 0.080625
# lat2 = 0.1



## 初始释放粒子范围和数量控制 
# 理想海岬实验-full domain
# minX = 105 # 球面坐标
# maxX = 105.2
# minY = 0
# maxY = 0.1
# N_edge = 2000 # 线密度: 个/度
# 此处N—_edge指每度的宽度上有多少个粒子，对于这种特别小的实验，可以适当调大

## 理想海岬实验-partial domain near headland
# minX = 105.04 # 球面坐标
# maxX = 105.16
# minY = 0.05
# maxY = 0.1
# N_edge = 2000 # 线密度: 个/度

## 测试三维追踪
minX = 112 # 球面坐标
maxX = 115
minY = 15
maxY = 18
N_edge = 40 # 线密度: 个/度








######################################################################################
######################################################################################
######################################################################################

x = np.linspace(minX,maxX,int(N_edge*(maxX-minX)+1))
y = np.linspace(minY,maxY,int(N_edge*(maxY-minY)+1))   
[lats,lons] = np.meshgrid(y,x)
lon = np.reshape(lons,[(lons.shape[0])*(lons.shape[1]),1])
lat = np.reshape(lats,[(lats.shape[0])*(lats.shape[1]),1])
dep = release_dep * np.ones(lon.shape)
time_ = release_time
positions = [lon,lat,dep,time_]

## 网格维度信息
content2 = Dataset(domain)
mask_rho = content2.variables['mask_rho'][:].data
lon_rho = content2.variables['lon_rho'][:].data
lat_rho = content2.variables['lat_rho'][:].data
content2.close()
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) # grid size
N = int(N_edge*(maxX-minX)+1)*int(N_edge*(maxY-minY)+1) # 质点数目
lonstart = lonarray[0]
latstart = latarray[0]
xinterval = (len_xi - 1)/(lonarray[-1] - lonstart)
yinterval = (len_eta - 1)/(latarray[-1] - latstart)
from multiprocessing import cpu_count
cpu_count = cpu_count() # 读取当前计算机核数
timestep_numbers = int(total_time/dt) # 总共运动多少个时间步长数目, 后面要用这个变量

## 设备名
host_name = socket.gethostname()