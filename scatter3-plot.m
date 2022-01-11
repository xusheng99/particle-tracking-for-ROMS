h = ncread('D:\application\scs_2\DATA\Grd.nc','h');
lon = ncread('D:\application\scs_2\DATA\Grd.nc','lon_rho');
lat = ncread('D:\application\scs_2\DATA\Grd.nc','lat_rho');

aa = dir('*.csv')

for i = 1:2:121;
figure()

subplot(1,2,1);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
view(-10,15)

subplot(1,2,2);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
axis equal
view(0,90)

showmax
savefig2png([aa(i).name,'.png'])
end
%%
figure()
i = 1
subplot(1,2,1);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
view(-10,15)
subplot(1,2,2);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
axis equal
view(0,90)
scatter3(115,15,0,0.00001,'filled')
showmax;title('20170101')

savefig2png(['20170101.png'])




%%

figure()
i = 60
subplot(1,2,1);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
view(-10,15)
subplot(1,2,2);
mesh(lon',lat',-log10(h'));shading interp;colorbar;
zlabel('log(h)')
hold on
% mesh(lon',lat',zeros(size(h')));
a = load(aa(i).name);
scatter3(a(:,1),a(:,2),-log10(-a(:,3)),4,'filled');
axis equal
view(0,90)
scatter3(115,15,0,0.00001,'filled')
showmax;title('20170115')

savefig2png(['20170115.png'])

%%
ncst = load('all_china_sea_i.mat');
scatter(ncst.ncst(:,1),ncst.ncst(:,2),1,'filled');
hold on
a(a(:,3)>-200) = NaN;

scatter(a(:,1),a(:,2),2,a(:,3),'filled');colorbar;
axis([105,125,5,25]);
title('depth 20170115');
axis equal


