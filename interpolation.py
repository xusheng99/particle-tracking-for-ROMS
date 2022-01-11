from scipy.interpolate.interpolate import interp2d
from para import *
from scipy import interpolate

import os
os.chdir(datadir)

def ncread(file,key):
    content = Dataset(file)
    a = content.variables[key][:].data
    content.close()
    return a

# 某时刻:是什么文件的第几个时间断面
def whichfile(time):
        monika = DT*his_writestep
        if time == 0:
            return prefix + '1'.zfill(numdigits) + '.nc'
        else:
            paul = time%monika
            if paul == 0:
                return prefix + str(int(time//monika+1)).zfill(numdigits) + '.nc',0
                # 实际应为 time//monika 而无需加一，此改动为了后续方便
            else:
                return prefix + str(int(time//monika+1)).zfill(numdigits) + '.nc',paul/DT


s_rho = ncread((datadir+whichfile(release_time)[0]),'s_rho')
Cs_r = ncread((datadir+whichfile(release_time)[0]),'Cs_r')
hc = ncread((datadir+whichfile(release_time)[0]),'hc')
h = ncread(datadir+whichfile(release_time)[0],'h')
Nlev = len(s_rho)

# 三维空间坐标数组构建
lon_rho3 = np.zeros([lon_rho.shape[0],lon_rho.shape[1],Nlev])
for i in range(Nlev):
    lon_rho3[:,:,i] = lon_rho
lat_rho3 = np.zeros([lat_rho.shape[0],lat_rho.shape[1],Nlev])
for i in range(Nlev):
    lat_rho3[:,:,i] = lat_rho
lonlon = np.squeeze(np.reshape(lon_rho3,[lon_rho3.shape[0]*lon_rho3.shape[1]*Nlev,1]))
latlat = np.squeeze(np.reshape(lat_rho3,[lat_rho3.shape[0]*lat_rho3.shape[1]*Nlev,1]))

if np.max(h) < hc:
    trans = s_rho
else:
    trans = Cs_r

# z_rho = np.zeros([lon_rho.shape[0],lon_rho.shape[1],5])
# for i in range(h.shape[0]):
#     for j in range(h.shape[1]):
#         z_rho[i,j] = h[i,j]*trans
# zz = np.squeeze(np.reshape(z_rho,[z_rho.shape[0]*z_rho.shape[1]*5,1]))
ss = np.zeros([lon_rho.shape[0],lon_rho.shape[1],Nlev])
for i in range(Nlev):
    ss[:,:,i] =  i
ss = np.squeeze(np.reshape(ss,[ss.shape[0]*ss.shape[1]*Nlev,1]))

# 水平网格索引
def lonlat2xy(lon,lat):
    # 将经纬度坐标 转化为 格点坐标
    x = (lon - lonstart)*xinterval
    y = (lat - latstart)*yinterval
    return x,y

# 深度索引
geth = interpolate.interp2d(np.arange(0,len_eta,1),np.arange(0,len_xi,1),h.T)
def dep2sigma(x,y,dep):
    # x = 84
    # y = 120
    # dep = -962.668879
    def aa(x,y,dep):
        hh = geth(y,x)
        array = hh*trans # 可以加入对zeta的考虑
        Nlev = lon_rho3.shape[-1]
        if dep < -hh[0]: 
            return -1,0
        if dep <= array[0]:
            return 0,0
        if dep >= array[-1]:
            return Nlev-2,1
        if dep in array:
            return list(array).index(dep),0
        else:
            for i in range(Nlev-1):
                if array[i] < dep <array[i+1]:
                    return i,(dep-array[i])/(array[i+1]-array[i])    
    rsigma = np.zeros(x.shape)
    rsigmastep = np.zeros(y.shape)
    for i in range(len(rsigma)):
        rsigma[i],rsigmastep[i] = aa(x[i],y[i],dep[i])
    return rsigma,rsigmastep


    
# 时间索引    
def time2index(time_):    
    filename,times = whichfile(time_)
    index1 = math.floor(times)-1
    index2 = math.floor(times)
    if index1<0:
        file1 = whichfile(time_-DT)[0]
        rindex1 = his_writestep-1
        file2 = filename
        rindex2 = index2
    else:
        file1 = filename
        file2 = filename
        rindex1 = index1
        rindex2 = index2


# v = v2 + (v3-v2)/*r2


def rho2u(var):
    [eta,xi] = var.shape
    var_u = 0.5*(var[:,0:xi-1]+var[:,1:])
    return var_u

def rho2v(var):
    [eta,xi] = var.shape
    var_v = 0.5*(var[0:eta-1,:]+var[1:,:])
    return var_v

def u2rho(var):
    [Mp,L] = var.shape
    Lp = L+1
    Lm = L-1
    var_rho = np.zeros([Mp,Lp])
    var_rho[:,1:L] = 0.5*(var[:,0:Lm] + var[:,1:L])
    var_rho[:,0] = var_rho[:,1]
    var_rho[:,Lp-1] = var_rho[:,L-1]
    return var_rho

def v2rho(var):
    [M,Lp] = var.shape
    Mp = M + 1
    Mm = M - 1
    var_rho = np.zeros([Mp,Lp])
    var_rho[1:M,:] = 0.5*(var[0:Mm,:] + var[1:M,:])
    var_rho[0,:] = var_rho[1,:]
    var_rho[Mp-1,:] = var_rho[M-1,:]
    return var_rho

# 构造测试数据集，正式运行时务必要注释掉
# lon = positions[0]
# lat = positions[1]
# time_ = positions[3]+7200
# dep = np.zeros([len(lon),1])
# for i in range(len(lon)):
#     dep[i] = np.random.uniform(-3,-geth(lat[i],lon[i])*0.95)


def V(lon,lat,dep,time_):
    # 定位文件和时刻
    filename,times = whichfile(time_)
    index1 = math.floor(times)-1
    index2 = math.floor(times)
    # 确定流场key
    if flowtype == 'rho':
        ukey = 'u_eastward'
        vkey = 'v_northward'
    else:
        ukey = 'u'
        vkey = 'v'
    # 读取前后两个“夹层”时刻的流场
    if index1<0:
        content1 = Dataset(whichfile(time_-DT)[0])
        u_rho1 = content1.variables[ukey][his_writestep-1,:,:].data
        v_rho1 = content1.variables[vkey][his_writestep-1,:,:].data
        w_rho1 = content1.variables['w'][his_writestep-1,:,:].data
        content1.close()
        content2 = Dataset(filename)
        u_rho2 = content2.variables[ukey][index2,:,:].data
        v_rho2 = content2.variables[vkey][index2,:,:].data
        w_rho2 = content2.variables['w'][index2,:,:].data
        content2.close()
    else:
        content1 = Dataset(filename)
        u_rho1 = content1.variables[ukey][index1,:,:,:].data
        v_rho1 = content1.variables[vkey][index1,:,:,:].data
        w_rho1 = content1.variables['w'][index1,:,:].data
        u_rho2 = content1.variables[ukey][index2,:,:,:].data
        v_rho2 = content1.variables[vkey][index2,:,:,:].data
        w_rho2 = content1.variables['w'][index2,:,:].data
        content1.close()
    # 时间线性内插
    u_rho = u_rho1 + (u_rho2 - u_rho1) / (index2 - index1) * (times - 1 - index1)
    v_rho = v_rho1 + (v_rho2 - v_rho1) / (index2 - index1) * (times - 1 - index1)
    w_rho = w_rho1 + (w_rho2 - w_rho1) / (index2 - index1) * (times - 1 - index1)
    [u_rhoo,v_rhoo] = np.zeros([Nlev,len_eta,len_xi]), np.zeros([Nlev,len_eta,len_xi])
    # 剔除陆地掩码点
    u_rho[abs(u_rho) > 999] = 0
    v_rho[abs(v_rho) > 999] = 0
    w_rho[abs(w_rho) > 999] = 0
    # 水平流速 u,v网格转rho网格
    if u_rho.shape != v_rho.shape:
        for i in range(u_rho.shape[0]):
            u_rhoo[i,:,:] = u2rho(u_rho[i,:,:])
            v_rhoo[i,:,:] = v2rho(v_rho[i,:,:])
        u_rho,v_rho = u_rhoo,v_rhoo
    # 垂向流速 w网格转rho网格
    w_rho1 = np.zeros([w_rho.shape[0]-1,w_rho.shape[1],w_rho.shape[2]])
    for i in range(w_rho1.shape[0]):
        w_rho1[i,:,:] = w_rho[i,:,:] + w_rho[i+1,:,:]
    w_rho = 0.5*w_rho1
    # lon,lat,dep -> x,y,sigma   POSITION TO INDEX
    [x,y] = lonlat2xy(lon,lat)
    sigmas = dep2sigma(x,y,dep)[0] + dep2sigma(x,y,dep)[1]
    xi = (sigmas,lat,lon)
    # INTERPOLATION: bi-linear, fill zero for out bounds
    uu = interpolate.interpn((list(np.arange(0,Nlev,1)),latarray,lonarray),u_rho,xi,\
        bounds_error=False,fill_value=0)
    vv = interpolate.interpn((list(np.arange(0,Nlev,1)),latarray,lonarray),v_rho,xi,\
        bounds_error=False,fill_value=0)
    ww = interpolate.interpn((list(np.arange(0,Nlev,1)),latarray,lonarray),w_rho,xi,\
        bounds_error=False,fill_value=0)
    # 单位变换
    uu = uu*(360/(6371000*2*3.1415926*np.cos(lat*(math.pi/180))))
    vv = vv*(360/(6371000*2*3.1415926))
    return [uu,vv,ww]

# dd = V(lon,lat,dep,time_)

# plt.hist(np.reshape(dd[0],[22311,1]))
# plt.hist(np.reshape(dd[1],[22311,1]))
# plt.hist(np.reshape(dd[2],[22311,1]))

# plt.scatter(lon,lat,3,dd[0])
# plt.colorbar()


# plt.scatter(lon,lat,3,dd[1])
# plt.colorbar()

# plt.scatter(lon,lat,3,dd[2])
# plt.colorbar()


















