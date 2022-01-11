import numpy as np
from netCDF4 import Dataset
from para import domain
from interpolation import geth
# coast = 'D:/application/0DATASRC/bh_coastline/bh_google.bln'
# domain = 'd:/application/bh/0910test3/1-0910_bh_CLM/Grd.nc' # Grid文件



content = Dataset(domain)
__mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
content.close()
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) 

# # 对于特殊形状的理想实验，使用更快的方法判断海陆
# def inocean(pnt):
#     [lon,lat] = pnt[0],pnt[1]
#     if lonmin>=lon or lonmax<=lon:
#         return 3 # 左右开边界
#     if latmin>=lat or latmax<=lat or (lon1<lon<lon2 and lat1<lat<lat2):
#         return 2 # 上下固壁边界, 海岬算陆地
#     return 1 

import matplotlib.path as mplPath
from para import coast
a = np.loadtxt(coast,delimiter=',')
## 找出分割符所在行
b = a==1
## 读取其下标
c = np.where(b == True)[0]
## 找到每个多边形的边数
dis = np.array([c[i+1]-c[i] for i in range(len(c)-1)]) - 1
d = np.zeros([len(c)-1,max(dis),2]) # 多边形编号, 节点编号, 经纬
## 提取每个多边形的节点, 存储在不同的行, 不够的位数补零
for i in range(d.shape[0]):
    d[i,0:dis[i],:] = a[c[i]+1:c[i+1],:]
## 创建多个多边形, 用于补齐的零是不计入的
polys = []
for i in range(d.shape[0]): # 创建多个多边形
    crd = d[i,0:dis[i]-1,:]
    poly = mplPath.Path(crd)
    polys.append(poly)


def inpoly(pnt,poly): # 给定点pnt和多边形poly, 返回是否在多边形内部
    # poly = mplPath.Path(poly)
    r = 0
    return (poly.contains_point(pnt,radius=r) or poly.contains_point(pnt,radius=-r))

# 使用多边形判定的方法，比较慢 全岛屿类
# def inocean(pnt):
#     [lon,lat] = pnt[0],pnt[1]
#     if lonmin>=lon or lonmax<=lon or latmin>=lat or latmax<=lat:
#         return 3 # 超出岸线范围，无法判定
#     for i in range(len(polys)):
#         a = inpoly(pnt,polys[i])
#         if a:
#             return 2 # 陆点
#     return 1 # 海点

# # 使用多边形判定的方法，比较慢 第一个为海类
def inocean(pnt):
    [lon,lat] = pnt[0],pnt[1]
    if lonmin>=lon or lonmax<=lon or latmin>=lat or latmax<=lat:
        return 3 # 超出岸线范围，无法判定
    aa = inpoly(pnt,polys[0])
    if not aa:
        return 2
    for i in range(1,len(polys)):
        a = inpoly(pnt,polys[i])
        if a:
            return 2 # 陆点
    return 1 # 海点

from interpolation import lonlat2xy,dep2sigma
def inocean_bottom(lon,lat,dep):
    [x,y] = lonlat2xy(lon,lat)
    sigmas = dep2sigma(x,y,dep)[0] + dep2sigma(x,y,dep)[1]
    flag2 = np.ones(sigmas.shape)
    for i in range(len(sigmas)):
        if sigmas[i]< 0 or dep[i] > 0:
            flag2[i] = 2
    return flag2