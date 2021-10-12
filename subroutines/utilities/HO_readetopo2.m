function [zt,lon,lat] = HO_readetopo2(etopo_path)
%%
data_lon = ncread(etopo_path,'x');
data_lon = double(data_lon + 180.0);
data_lat = ncread(etopo_path,'y');

nx = length(data_lon); ny = length(data_lat);

z = ncread(etopo_path,'z');
data_z = zeros(ny,nx);

for i = 1:nx
    for j = 1:ny

        data_z(j,i) = z(i,j);
        
    end
end
data_z = [data_z(:,5401:end) data_z(:,1:5400)];
%%
lonL = 100; lonR = 300; latB = -15; latT = 65;
%lonL = 1; lonR = 360; latB = -65; latT = 65;
lenx = (lonR - lonL + 1)*1; leny = (latT - latB + 1)*1;
%
lon = linspace(lonL,lonR,lenx); lat = linspace(latB,latT,leny);
ilon = zeros(1,lenx); ilat = zeros(1,leny);
for i = 1:lenx
    [~, Min] = min(abs(data_lon-lon(i)));
    ilon(i) = Min;
end
for j = 1:leny
    [~, Min] = min(abs(data_lat-lat(j)));
    ilat(j) = Min;
end

zt = zeros(lenx,leny);
for i = 1:lenx
    for j = 1:leny
        
        zt(i,j) = data_z(ilat(j),ilon(i));
        
    end
end

end


    
    
                    
                    