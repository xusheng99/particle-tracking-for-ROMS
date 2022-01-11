from para import *
from scipy import interpolate
from locale import Error

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

lonlon = np.squeeze(np.reshape(lon_rho,[lon_rho.shape[0]*lon_rho.shape[1],1]))
latlat = np.squeeze(np.reshape(lat_rho,[lat_rho.shape[0]*lat_rho.shape[1],1]))

def IDW_V(lons,lats,time_):
    # 定位文件和时刻 读取
    filename,times = whichfile(time_)
    index1 = math.floor(times)-1
    index2 = math.floor(times)
    if index1<0:
        content1 = Dataset(whichfile(time_-DT)[0])
        u_rho1 = content1.variables['u_sur_eastward'][his_writestep-1,:,:].data
        v_rho1 = content1.variables['v_sur_northward'][his_writestep-1,:,:].data
        content1.close()
        content2 = Dataset(filename)
        u_rho2 = content2.variables['u_sur_eastward'][index2,:,:].data
        v_rho2 = content2.variables['v_sur_northward'][index2,:,:].data
        content2.close()
    else:
        content1 = Dataset(filename)
        u_rho1 = content1.variables['u_sur_eastward'][index1,:,:].data
        v_rho1 = content1.variables['v_sur_northward'][index1,:,:].data
        u_rho2 = content1.variables['u_sur_eastward'][index2,:,:].data
        v_rho2 = content1.variables['v_sur_northward'][index2,:,:].data
        content1.close()

    # 时间线性内插
    u_rho = u_rho1 + (u_rho2 - u_rho1) / (index2 - index1) * (times - 1 - index1)
    v_rho = v_rho1 + (v_rho2 - v_rho1) / (index2 - index1) * (times - 1 - index1)
    u_rho[abs(u_rho)>100] = 0
    v_rho[abs(v_rho)>100] = 0

    # 创建存储空数组
    uu = np.zeros([len(lons)])
    vv = np.zeros([len(lons)])

    for ii in range(len(lons)):
        lon = lons[ii]
        lat = lats[ii]
        # 将经纬度坐标 转化为 格点坐标
        lonstart = lonarray[0]
        latstart = latarray[0]
        xinterval = (len_xi - 1)/(lonarray[-1] - lonstart)
        yinterval = (len_eta - 1)/(latarray[-1] - latstart)
        x = (lon - lonstart)*xinterval
        y = (lat - latstart)*yinterval

        x1 = int(x//1)  #整数部分  
        x2 = x1 + 1
        y1 = int(y//1)
        y2 = y1 + 1

        u = np.zeros([4])
        v = np.zeros([4])

        u[0] = u_rho[y1,x1]
        v[0] = v_rho[y1,x1]
        u[1] = u_rho[y1,x2]
        v[1] = v_rho[y1,x2]
        u[2] = u_rho[y2,x2]
        v[2] = v_rho[y2,x2]
        u[3] = u_rho[y2,x1]
        v[3] = v_rho[y2,x1]

        def distance(m1,n1,m2,n2):  # 欧几里得距离
            return math.sqrt((m1-m2)*(m1-m2)+(n1-n2)*(n1-n2)) 

        d = np.zeros([4])
        d[0] = distance(x,y,x1,y1)
        d[1] = distance(x,y,x2,y1)
        d[2] = distance(x,y,x2,y2)
        d[3] = distance(x,y,x1,y2)

        for point in range(len(d)):  # 正巧卡在了整数编号的网格点上，手动处理分母为零的情况
            if abs(d[point]) < 1e-15:
                return [1*u[point]*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180)))),1*v[point]*(360/(6371000*2*3.14159))]
            
        def func_d(distance_list):  # 每点的权重因子, 平方反比关系
            fd = np.zeros([4])
            for point in range(len(distance_list)):
                fd[point] = (1/(distance_list[point]*distance_list[point]))
            return fd
        fd_list = func_d(d)

        def weight(fd_list):  # 每点的权
            weight = np.zeros([4])
            for point in range(len(fd_list)):
                weight[point] = ((fd_list[point])/(sum(fd_list)))
            return weight
        weight = weight(fd_list)
        
        # 加权平均, 即 权向量u/v 和 速率向量weight 点乘
        uu[ii] = sum(u*weight)
        vv[ii] = sum(v*weight)

    # 扩散项
    if diffusion == 'on':  
        uu = uu + np.sqrt((2*Dh)/dt)*np.random.uniform(-1,1,size=[N,1])
        vv = vv + np.sqrt((2*Dh)/dt)*np.random.uniform(-1,1,size=[N,1])

    # 单位变换
    uu = uu*(360/(6371000*2*3.1415926*np.cos(lat*(math.pi/180))))
    vv = vv*(360/(6371000*2*3.1415926))
    return [uu,vv]


def BATCH_BILINAER_V(lon,lat,time_):

    # 定位文件和时刻
    filename,times = whichfile(time_)
    index1 = math.floor(times)-1
    index2 = math.floor(times)
    if index1<0:
        content1 = Dataset(whichfile(time_-DT)[0])
        u_rho1 = content1.variables['u_sur_eastward'][his_writestep-1,:,:].data
        v_rho1 = content1.variables['v_sur_northward'][his_writestep-1,:,:].data
        content1.close()
        content2 = Dataset(filename)
        u_rho2 = content2.variables['u_sur_eastward'][index2,:,:].data
        v_rho2 = content2.variables['v_sur_northward'][index2,:,:].data
        content2.close()
    else:
        content1 = Dataset(filename)
        u_rho1 = content1.variables['u_sur_eastward'][index1,:,:].data
        v_rho1 = content1.variables['v_sur_northward'][index1,:,:].data
        u_rho2 = content1.variables['u_sur_eastward'][index2,:,:].data
        v_rho2 = content1.variables['v_sur_northward'][index2,:,:].data
        content1.close()

    # 时间线性内插
    u_rho = u_rho1 + (u_rho2 - u_rho1) / (index2 - index1) * (times - 1 - index1)
    v_rho = v_rho1 + (v_rho2 - v_rho1) / (index2 - index1) * (times - 1 - index1)
    u_rho[abs(u_rho)>100] = 0
    v_rho[abs(v_rho)>100] = 0
    
    # 空间插值，双线性插值
    u_rho = np.squeeze(np.reshape(u_rho,[u_rho.shape[0]*u_rho.shape[1],1]))
    uu = interpolate.griddata((lonlon,latlat),u_rho,(lon,lat),fill_value=0)
    v_rho = np.squeeze(np.reshape(v_rho,[v_rho.shape[0]*v_rho.shape[1],1]))
    vv = interpolate.griddata((lonlon,latlat),v_rho,(lon,lat),fill_value=0)

    # 扩散项
    if diffusion == 'on':  
        uu = uu + np.sqrt((2*Dh)/dt)*np.random.uniform(-1,1,size=[N,1])
        vv = vv + np.sqrt((2*Dh)/dt)*np.random.uniform(-1,1,size=[N,1])

    # 单位变换
    uu = uu*(360/(6371000*2*3.1415926*np.cos(lat*(math.pi/180))))
    vv = vv*(360/(6371000*2*3.1415926))
    return [uu,vv]

def ANA_V(lon,lat,time_):
    
if h_interp_method == 'IDW':
    V = IDW_V
elif h_interp_method == 'BATCH_LINEAR':
    V = BATCH_BILINAER_V
elif h_interp_method == 'ANA':
    V = ANA_V
else:
    raise Error('invalid para input: h_interp_method!')