latdtu18=[53.0000000000:0.0166666600:66.9999944000];
londtu18=[8.0000000000:0.0166666600:30.9999908000]';
fid=fopen('mss18.grd','r');
t=textscan(fid, '%f');
mssdtuu18=reshape(t{1,1},[1381,841]);


[lon_geo, lat_geo, z_geo]=grdread2('NKG2015.grd');
ndtu18=griddata(lon_geo, lat_geo, z_geo,londtu18,latdtu18);


[C,h] =contourf(londtu18,latdtu18,(mssdtuu18'-ndtu18));
set(h,'LineColor','none')
hold on
plot(coast(:,2),coast(:,1),'.k')

