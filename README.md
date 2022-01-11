# particle-tracking-for-ROMS
particle tracking model, works with the ROMS output file(qck.nc, his.nc)

## description
the 3-dimension version of the model, particle moves vertically.

## features and numerical schemes
* 4th order RK method is applied to the particle movement progress.
* support time backward tracking.
* diffusion is implemented by random walk algorithm.


## usage
1. edit the parameter is para.py, the meaning are explained in the comment
2. cd to the path where "start_tracking.py" is in, run:
~~~
python ./start_tracking_3d.py
~~~
if on linux platform:
~~~
python3 ./start_tracking_3d.py
~~~
