from para import *
from interpolation import V
from oceanmask_toolbox import inocean,inocean_bottom

# 四阶龙格-库塔算法 给定速度计算粒子位置
def gonext(position):
    # position = positions
    t1 = time.perf_counter()
    lon_0 = position[0]
    lat_0 = position[1]
    dep_0 = position[2]
    time_0 = position[3]
    [u1,v1,w1] = V(lon_0,lat_0,dep_0,time_0)
    [u2,v2,w2] = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, dep_0+dt*0.5*w1, time_0+dt*0.5 )
    [u3,v3,w3] = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, dep_0+dt*0.5*w2, time_0+dt*0.5 )
    [u4,v4,w4] = V( lon_0+dt*u3, lat_0+dt*v3, dep_0+dt*w3, time_0+dt)
    u,v,w = (1/6)*(u1+2*u2+2*u3+u4), (1/6)*(v1+2*v2+2*v3+v4), (1/6)*(w1+2*w2+2*w3+w4)
    lon,lat,dep = lon_0 + dt*u, lat_0 + dt*v, dep_0 + dt*w
    t2 = time.perf_counter()
    print('    ','interpolation:',round(t2-t1,2),'s')
    # 掩码计算
    flag = np.zeros([len(lon),1]) # 位置状态flag：海/陆/出边界 对应 1/2/3
    flag2 = inocean_bottom(lon,lat,dep)
    for i in range(len(lon)):
        flag[i] = inocean([lon[i],lat[i]])
        if  flag[i] == 2: # 如果上到了陆地上, 移动到上一步所在位置
            lon[i] = lon_0[i] 
            lat[i] = lat_0[i]
        elif flag[i] == 3: # 如果离开了开边界, 设置为Inf值, 停止计算
            lon[i] = 99999
            lat[i] = 99999
        if flag2[i] == 2: # 如果撞到了海底，回到上一个时刻所在的深度
            dep[i] = dep_0[i]
    t3 = time.perf_counter()
    print('    ','land-ocean masking:',round(t3-t2,2),'s')    
    return [lon,lat,dep,time_0+dt]

# positions = gonext(positions)
# [a,b,c] = V(positions[0],positions[1],positions[2],positions[3])
# positions[0]

# 某时刻粒子位置快照，保存到文本文件方便查看
def snapshot(pos,filename):
    lon = pos[0]
    lat = pos[1]
    dep = pos[2]
    t = pos[3]
    with open(filename,'w',encoding='utf-8') as f:
        for i in range(len(lon)):
            f.write(str(float(lon[i]))+' , ')
            f.write(str(float(lat[i]))+' , ')
            f.write(str(float(dep[i]))+' , ')
            f.write(str(t)+'\n')
    f.close()

# 可视化粒子当前位置分布
def view():
    plt.figure(dpi=300)
    a = positions[0]
    b = positions[1]
    a[a == 99999] = np.NaN
    b[b == 99999] = np.NaN
    plt.scatter(a,b,1)



if __name__ == '__main__':
    print('PROJECT NAME: ',projectname)
    print('\n粒子总数: {} '.format(N))
    print('共计追踪时间: {} hours'.format(round((total_time/3600),3)))
    print('共计追踪时间: {} days'.format(round((total_time/86400),3)))
    print('粒子单步运动步长:{} seconds'.format(dt))

    print(os.getcwd())
    # 清空工作目录
    files = os.listdir(workdir)
    for each in files:
        if '.csv' in each:
            os.remove(workdir+each)

    # 开始追踪
    timenum = int(total_time/dt) 
    savenum = int(total_time/(dt*write_step))
    result = np.empty([len(lon),4,savenum])
    count = 0
    positions = positions
    while positions[3] <= (release_time + total_time):
        if count % write_step == 0:
            result[:,0,count//write_step] = np.squeeze(positions[0])
            result[:,1,count//write_step] = np.squeeze(positions[1])
            result[:,2,count//write_step] = np.squeeze(positions[2])
            result[:,3,count//write_step] = np.squeeze(positions[3]*np.ones([len(lon),1]))
            print('\nDAY',str(round(positions[3]/86400,4)),': SNAPSHOT TAKEN!')
            snapshot(positions,str(count//write_step).zfill(6)+'.csv')
            shutil.move(str(count//write_step).zfill(6)+'.csv',workdir)
        print(count,' of ',timenum)
        positions = gonext(positions)
        count += 1
    np.save('tracking_output.npy',result)

# view()


'''

lonori = positions[0]
latori = positions[1]

lonini = positions[0]
latini = positions[1]

V(positions[0],positions[1],positions[2],positions[3])


# 绘制水平分布
plt.figure(dpi=300)
a = positions[0].copy()
b = positions[1].copy()
a[a == 99999] = np.NaN
b[b == 99999] = np.NaN
plt.scatter(a,b,5)

positions = gonext(positions)

# 绘制深度
plt.figure(dpi=300)
a = positions[0].copy()
b = positions[1].copy()
dd = positions[2].copy()
a[a == 99999] = np.NaN
b[b == 99999] = np.NaN
plt.scatter(a,b,5,dd)
plt.title('depth distribution')
plt.colorbar()

# 绘制upwelling
[a,b,c] = V(positions[0],positions[1],positions[2],positions[3])
plt.figure(dpi=300)
a = positions[0].copy()
b = positions[1].copy()
a[a == 99999] = np.NaN
b[b == 99999] = np.NaN
plt.scatter(a,b,5,c)
plt.title('vertical flow')
plt.colorbar()

'''


# plt.figure(dpi=300)
# plt.hist((a-lonini)/(0.2/160))