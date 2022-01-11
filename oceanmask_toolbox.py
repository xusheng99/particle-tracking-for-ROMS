from locale import Error
import numpy as np
from netCDF4 import Dataset
from numpy.ma.core import concatenate
from para import domain,lon1,lon2,lat1,lat2,oceanmask_method
# coast = 'D:/application/0DATASRC/bh_coastline/bh_google.bln'
# domain = 'd:/application/bh/0910test3/1-0910_bh_CLM/Grd.nc' # Grid文件

def ncread(file,key):
    content = Dataset(file)
    a = content.variables[key][:].data
    return a

content = Dataset(domain)
mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
content.close()
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) 
londel = (lonmax-lonmin)/(len_xi-1)
latdel = (latmax-latmin)/(len_eta-1)


############## SOLUTION 0: SCIPY GRIDDATA #############
## 使用基于C的scipy插值程序，从mask_rho插值得到海陆掩码 ##
xx = np.reshape(lon_rho,[len_xi*len_eta,1])
yy = np.reshape(lat_rho,[len_xi*len_eta,1])
ll = np.concatenate((xx,yy),1)
vv = np.reshape(mask_rho,[len_xi*len_eta,1])
from scipy.interpolate import griddata
def inocean0(pnts):
    [lon,lat] = pnts[0],pnts[1]
    themask = np.squeeze(griddata(ll,vv,(lon,lat),method='nearest'))
    return themask

############## SOLUTION 1: ROMS GRID MASK #############
def cor2gridindex(lon,lat):
    xi_index = np.round((np.float(lon)-lonmin)/londel)
    eta_index = np.round((np.float(lat)-latmin)/latdel)
    return [int(xi_index),int(eta_index)]

def inocean1(pnt):
    [lon,lat] = pnt[0],pnt[1]
    if lonmin>=lon or lonmax<=lon or latmin>=lat or latmax<=lat:
        return 2 # 超出岸线范围，无法判定
    [xi_index,eta_index] = cor2gridindex(lon,lat)
    return int(mask_rho[eta_index,xi_index])

############### SOLUTION 2: land and island polygon,slow but accurate #########
if oceanmask_method == 2:
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

# 使用多边形判定的方法，比较慢，但准确 [海洋闭多边形,bh_google类]
def inocean2(pnt):
    [lon,lat] = pnt[0],pnt[1]
    if lonmin>=lon or lonmax<=lon or latmin>=lat or latmax<=lat:
        return 2 # 超出岸线范围，无法判定
    for i in range(len(polys)):
        a = inpoly(pnt,polys[i])
        if a:
            return 0 # 陆点
    return 1 # 海点

'''    
# 使用多边形判定的方法，比较慢，但准确 [大陆闭多边形,scs.bln类]
def inocean(pnt):
    [lon,lat] = pnt[0],pnt[1]
    if lonmin>=lon or lonmax<=lon or latmin>=lat or latmax<=lat:
        return 2 # 超出岸线范围，无法判定
    pp = polys[0]
    if not inpoly(pnt,pp):
        return 0
    for i in range(1,len(polys)):
        a = inpoly(pnt,polys[i])
        if a:
            return 0 # 陆点
    return 1 # 海点
'''  


############## SOLUTION 3: simple geometric method, limited 
# 对于特殊形状的理想实验，使用更快的方法判断海陆
def inocean3(pnt):
    pass
    [lon,lat] = pnt[0],pnt[1]
    if lonmin>=lon or lonmax<=lon:
        return 2 # 左右开边界
    if latmin>=lat or latmax<=lat or (lon1<lon<lon2 and lat1<lat<lat2):
        return 0 # 上下固壁边界, 海岬算陆地
    return 1 





# switch oceanmask method
if oceanmask_method == 0:
    inocean = inocean0
elif oceanmask_method == 1:
    inocean = inocean1
elif oceanmask_method == 2:
    inocean = inocean2
elif oceanmask_method == 3:
    inocean = inocean3
else:
    raise Error('INVALID OCEANMASK METHOD PARAMETER! CHECK para.py')
