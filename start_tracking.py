from para import *
from interpolation import V
from oceanmask_toolbox import inocean
import warnings 
warnings.filterwarnings('ignore')

# RK4
def gonext(position):
    position = positions
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    t1 = time.perf_counter()

    ## land-ocean mask 1
    flag = inocean([lon_0,lat_0])
    t2 = time.perf_counter()
    print('    ',round(t2-t1,2))

    ## particle movement
    [u1,v1] = V(lon_0,lat_0,time_0)
    [u2,v2] = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, time_0+dt*0.5 )
    [u3,v3] = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, time_0+dt*0.5 )
    [u4,v4] = V( lon_0+dt*u3, lat_0+dt*v3, time_0+dt)
    u,v = (1/6)*(u1+2*u2+2*u3+u4), (1/6)*(v1+2*v2+2*v3+v4)
    u = u*flag
    v = v*flag
    lon,lat = lon_0 + dt*u, lat_0 + dt*v
    t3 = time.perf_counter()
    print('    ', round(t3-t2,2))
    
    ## land-ocean mask 2
    flag = inocean([lon_0,lat_0])
    for i in range(len(lon)):
        if  flag[i] == 0: # if it went ashore, back
            lon[i] = lon_0[i] 
            lat[i] = lat_0[i]
        elif flag[i] == 2: # if it left the open boundary, kill it
            lon[i] = 99999
            lat[i] = 99999
            # lon[i] = lon_0[i]
            # lat[i] = lat_0[i]
    t4 = time.perf_counter()
    print('    ',round(t4-t3,2))
    return [lon,lat,time_0+dt]



# take a snap shot for all particles to a dat file 
def snapshot(pos,filename):
    lon = pos[0]
    lat = pos[1]
    t = pos[2]
    with open(filename,'w',encoding='utf-8') as f:
        for i in range(len(lon)):
            f.write(str(float(lon[i]))+' , ')
            f.write(str(float(lat[i]))+' , ')
            f.write(str(t)+'\n')
    f.close()

# visualize the particle location
def view():
    plt.figure(dpi=300)
    a = positions[0]
    b = positions[1]
    a[a == 99999] = np.NaN
    b[b == 99999] = np.NaN
    plt.scatter(a,b,1)



if __name__ == '__main__':
    T1 = time.perf_counter()
    print('PROJECT NAME: ',projectname)
    print('\nparticle count: {} '.format(N)) 
    print('tracking time: {} hours'.format(round((total_time/3600),3)))
    print('tracking time: {} days'.format(round((total_time/86400),3)))
    print('time per step:{} seconds'.format(dt))

    timenum = int(total_time/abs(dt)) 
    savenum = int(total_time/(abs(dt)*write_step))+1
    result = np.empty([len(lon),3,savenum])
    count = int((positions[2] - release_time)//abs(dt))
    positions = positions
    if timearrow == 'forward':
        while positions[2] <= (release_time + total_time):
            if count % write_step == 0:
                result[:,0,count//write_step] = np.squeeze(positions[0])
                result[:,1,count//write_step] = np.squeeze(positions[1])
                result[:,2,count//write_step] = np.squeeze(positions[2]*np.ones([len(lon),1]))
                print('\nDAY',str(round(positions[2]/86400,4)),': SNAPSHOT TAKEN!')
                snapshot(positions,str(count//write_step).zfill(6)+'.csv')
                # shutil.move(str(count//write_step).zfill(6)+'.csv',workdir)
            print(count,' of ',timenum)
            if T_cycle != 0:
                if positions[2] - release_time == T_cycle: # 
                    positions[0] = positions_0[0]
                    positions[1] = positions_0[1]
                    print('#################')
                    print('all particles are reset')
                    print('#################\n')
            positions = gonext(positions)
            count += 1
    elif timearrow == 'backward':
        while positions[2] >= (release_time - total_time):
                if count % write_step == 0:
                    result[:,0,count//write_step] = np.squeeze(positions[0])
                    result[:,1,count//write_step] = np.squeeze(positions[1])
                    result[:,2,count//write_step] = np.squeeze(positions[2]*np.ones([len(lon),1]))
                    print('\nDAY',str(round(positions[2]/86400,4)),': SNAPSHOT TAKEN!')
                    snapshot(positions,'backward_'+str(count//write_step).zfill(6)+'.csv')
                    # shutil.move('backward_'+str(count//write_step).zfill(6)+'.csv',workdir)
                print(count,' of ',timenum)
                positions = gonext(positions)
                count += 1
    T2 = time.perf_counter()
    np.save('tracking_output.npy',result)
    T3 = time.perf_counter()

    print('TOTAL TIME CONSUME: {} seconds'.format(round(T3-T1),4))
    print('WRITING OUTPUTFILE: {} seconds'.format(round(T3-T2),4))
    # shutil.move('tracking_output.npy',workdir)






# view()


'''

for i in range(1000):
    positions = gonext(positions)
    print(i,'!!!!!')
    positions[2] = positions_0[2]


lonori = positions[0]
latori = positions[1]

lonini = positions[0]
latini = positions[1]


dt = -2400
positions = gonext(positions)   

plt.figure(dpi=300)
a = positions[0]
b = positions[1]
# a[a == 99999] = np.NaN
# b[b == 99999] = np.NaN
plt.scatter(a,b,0.1)
print(positions[-1])


dt = 240
for i in range(10):
    positions = gonext(positions)
plt.figure(dpi=300)
a = positions[0]
b = positions[1]
# a[a == 99999] = np.NaN
# b[b == 99999] = np.NaN
plt.scatter(a,b,0.00001)


'''


# plt.figure(dpi=300)
# plt.hist((a-lonini)/(0.2/160))