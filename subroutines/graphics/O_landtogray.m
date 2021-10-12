function [lon_land_box,lat_land_box]=O_landtogray(alon,alat,fsm,im,jm)

alon_new=zeros(im+1,jm+1);
alat_new=zeros(im+1,jm+1);

alon_new(2:im,2:jm)=(alon(1:im-1,1:jm-1)+alon(2:im,1:jm-1)+alon(2:im,2:jm)+alon(1:im-1,2:jm))/4;
alat_new(2:im,2:jm)=(alat(1:im-1,1:jm-1)+alat(2:im,1:jm-1)+alat(2:im,2:jm)+alat(1:im-1,2:jm))/4;

alon_new(2:im,1)=alon_new(2:im,2);
alat_new(2:im,1)=(alat(1:im-1,1)+alat(2:im,1))-alat_new(2:im,2);

alon_new(2:im,jm+1)=alon_new(2:im,jm);
alat_new(2:im,jm+1)=(alat(1:im-1,jm)+alat(2:im,jm))-alat_new(2:im,jm);

alon_new(1,2:jm)=(alon(1,1:jm-1)+alon(1,2:jm))-alon_new(2,2:jm);
alat_new(1,2:jm)=alat_new(2,2:jm);

alon_new(im+1,2:jm)=(alon(im,1:jm-1)+alon(im,2:jm))-alon_new(im,2:jm);
alat_new(im+1,2:jm)=alat_new(im,2:jm);

alon_new(1,1)      =alon_new(1,2);     alat_new(1,1)      =alat_new(2,1);
alon_new(1,jm+1)   =alon_new(1,jm);    alat_new(1,jm+1)   =alat_new(2,jm+1);
alon_new(im+1,1)   =alon_new(im+1,2);  alat_new(im+1,1)   =alat_new(im,1);
alon_new(im+1,jm+1)=alon_new(im+1,jm); alat_new(im+1,jm+1)=alat_new(im,jm+1);

%-------------------------------------------------------------------------
    allpoints=im*jm;
    land_grid_no=allpoints-sum(sum(fsm));
    lon_land_box=zeros(4,land_grid_no);
    lat_land_box=zeros(4,land_grid_no);
    k_grid=0; 
    lon_min=alon(1,1);lat_min=alat(1,1);
    lon_max=alon(im,jm);lat_max=alat(im,jm);
    for ii=1:im
        for jj=1:jm
if (fsm(ii,jj)==0 && alon(ii,jj)>=lon_min-0.5 && alon(ii,jj)<=lon_max+0.5 && alat(ii,jj)>=lat_min-0.5 && alat(ii,jj)<=lat_max+0.5)
                k_grid=k_grid+1;
                lon_land_box(1:4,k_grid)=[alon_new(ii,jj);alon_new(ii+1,jj);alon_new(ii+1,jj+1);alon_new(ii,jj+1)];
                lat_land_box(1:4,k_grid)=[alat_new(ii,jj);alat_new(ii+1,jj);alat_new(ii+1,jj+1);alat_new(ii,jj+1)];
            end
        end
    end
    fill_h=fill(lon_land_box(:,1:k_grid),lat_land_box(:,1:k_grid),[0.5 0.5 0.5]);
    set(fill_h, 'LineStyle', 'none');
    hold on;
    return
