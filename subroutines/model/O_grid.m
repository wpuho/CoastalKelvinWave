function [dx,dy,are,aru,arv,h,fsm,dum,dvm,cor,dt,xl,yl,zt,lon,lat]=O_grid  ...
                                                (small,g,dtfac,nperx,npery,etopo_path)
%
%Given xy-domain length (xl,yl), specify rectangular model domain
%   grid, topography & Coriolis
%
%Grid and domain array sizes
hc=10.0;                         %channel depth (if constant)
%{
xl=9.e5; yl=6.e5;                %channel x & y lengths 
dxc= 2.e4; im = floor(xl/dxc)+1;    %1-12 & 16-17,19;
dyc= 1.e4; jm = floor(yl/dyc)+1; %1-16 & 18;
dx(1:im,1:jm) = dxc; %---constant dx for this problem
dy(1:im,1:jm) = dyc; %---constant dy
%}
radx = 6378000.0;
d2r = pi/180.0;
%progdir='./';etopo_path = [progdir,'ETOPO2v2c_f4_netCDF/ETOPO2v2c_f4.nc'];
[zt, lon, lat] = HO_readetopo3(etopo_path);
im = length(lon);
jm = length(lat);
%zt = [zt; zt(1,:); zt(2,:)];
dy = zeros(im,jm); dx = zeros(im,jm);
for i = 1:im
    for j = 1:jm
        dy(i,j) = radx*abs(0.3333)*d2r;
    end
end
%dy(:,jm) = dy(:,jm-1);
for i = 1:im
    for j = 1:jm
        dx(i,j) = radx*cos(lat(j)*d2r)*(0.3333*d2r);
    end
end
%dx(im,:) = dx(1,:);
%dx(im-1,:) = dx(im-2,:);

xl = zeros(1,jm); yl = zeros(1,im);
for j = 1:jm
    xl(j) = radx*cos(lat(j)*d2r)*((lon(end)-lon(1))*d2r);
end
for i = 1:im
    yl(i) = radx*abs((lat(end)-lat(1)))*d2r;
end
%---grid areas centered at el & (ua,va) points
imm1 = im-1; jmm1 = jm-1;
are(1:im,1:jm)=dx(1:im,1:jm).*dy(1:im,1:jm);
aru(2:im,1:jm)=.25*(dx(1:imm1,1:jm)+dx(2:im,1:jm)).*(dy(1:imm1,1:jm)+dy(2:im,1:jm)); 
aru(1,1:jm)=aru(2,1:jm);
arv(1:im,2:jm)=.25*(dx(1:im,1:jmm1)+dx(1:im,2:jm)).*(dy(1:im,1:jmm1)+dy(1:im,2:jm)); 
arv(1:im,1)=arv(1:im,2);
%
%Specify channel topography, h (= depth when there is no motion):
hw = hc; he = hc;       %h@west  and east
hs = hc; hn = hc;       %h@south and north
%{
dhx= (he-hw)/xl;        %---depth-increament for x-linear h (below)
dhy= (hn-hs)/yl;        %---depth-increament for y-linear h (below)
%}
dhx = 0; dhy = 0;
h=zeros(im,jm)+small;   %---initialize w/land points everywhere
for i = 1:im;   h(i,:) = hw+dhx.*(i-1)*dx(i,:); end; 
for j = 1:jm;   h(:,j) = hs+dhy.*(j-1)*dy(:,j); end;
%
%h(1,:) = small; h(im,:) = small; %coast west  & east
 %h(:,1) = small; %h(:,jm) = small; %coast south & north
 for i = 1:im
     for j = 1:jm
         
         if (zt(i,j)>=-1000)
             h(i,j) = small;
         else
             h(i,j) = hc;
         end
         
     end
 end
%
%Define land masks for elevation & velocity, using the "h" above: 
%   fsm & dum =1 at water points, =0 at land/"wall" points
%   fsm is for elevation point; dum is for ua-point
fsm(1:im,1:jm)=0.0; dum=fsm; dvm=fsm;         %---initialize
for j=1:jm; for i=1:im; if(h(i,j)>small)fsm(i,j)=1.0; end;end;end;
%
%Take care of periodic boundary to be compatible with the above mask:
if (nperx == 1 && im > 3);  %x-periodic: 2 west cells are land if 1 is
for j=1:jm; if (fsm(1,j) == 0.0); fsm(2,j)=0.0; end; end; 
end;
if (npery == 1 && jm > 3);  %y-periodic: 2 south cells are land if 1 is
for i=1:im; if (fsm(i,1) == 0.0); fsm(i,2)=0.0; end; end; 
end;
%
dum(2:im,:) = fsm(1:imm1,:).*fsm(2:im,:); dum(1,:)=dum(2,:);
dvm(:,2:jm) = fsm(:,1:jmm1).*fsm(:,2:jm); dvm(:,1)=dvm(:,2);
%
h(:,:)=fsm(:,:).*h(:,:)+small; %Redo "h" again to make sure lands are set
%{
%Coriolis parameter:
lat0 =43.5; %25.0;              %---Model domain center latitude
                                %   lat0=43.5 gives f = 1.e-4
cor0 = 4.*pi*sind(lat0)/86400.; %---Coriolis @domain center
                                %   maybe modified below
cor0 = 10./100000.;             %set to be c/100km, so Rossby a = 100km,
                                % where "c" is from O_modparam
                                % So, lat0 (above) is ignored, but it is:
                                % lat0=asind(86400*cor0/(4*pi)); =43.4364N
cor0 = 1.*cor0;                 %---Make earth spin faster
cor=zeros(im,jm)+cor0;          %---Coriolis
%}
cor = zeros(im,jm);
omega = (360./86400.);
for i = 1:im
    for j = 1:jm
        
        cor(i,j) = 2*(omega*d2r)*sin(lat(j)*d2r);
        
    end
end
%Calculate time step "dt" satisfying CFL:
%   CFL condition is: sqrt(gh)*dt/dx <= 1/2 --> 2d version:
%       dtcfl = 1./(2.*sqrt(hmax*g)*sqrt(1./dxmin^2+1./dymin^2));
hmax = O_maxval(h); dxmin = O_minval(dx); dymin = O_minval(dy);
dtcfl = dtfac*(1./(2.*sqrt(hmax*g)*sqrt(1./dxmin^2+1./dymin^2)));
dt = floor(dtcfl); %dt = floor(0.5*dxmin/sqrt(hmax*g)); %-->1d CFL
%
return;

%%
%{
%Given xy-domain length (xl,yl), specify rectangular model domain
%   grid, topography & Coriolis
%
%Grid and domain array sizes
hc=10.0;                         %channel depth (if constant)
xl=9.e5; yl=6.e5;                %channel x & y lengths 
dxc= 1.e4; im = floor(xl/dxc)+1;  %01,02,05,06
%dxc= 2.5e3; im = floor(xl/dxc)+1; %03,04
%ForcedKW:dxc= 1.25e3; im = floor(xl/dxc)+1; %15%23
dyc= 1.e4; jm = floor(yl/dyc)+1; %01,02,05,06
%dyc= 2.5e3; jm = floor(yl/dyc)+1; %04
dx(1:im,1:jm) = dxc; %---constant dx for this problem
dy(1:im,1:jm) = dyc; %---constant dy
%---grid areas centered at el & (ua,va) points
imm1 = im-1; jmm1 = jm-1;
are(1:im,1:jm)=dx(1:im,1:jm).*dy(1:im,1:jm);
aru(2:im,1:jm)=.25*(dx(1:imm1,1:jm)+dx(2:im,1:jm)).*(dy(1:imm1,1:jm)+dy(2:im,1:jm)); 
aru(1,1:jm)=aru(2,1:jm);
arv(1:im,2:jm)=.25*(dx(1:im,1:jmm1)+dx(1:im,2:jm)).*(dy(1:im,1:jmm1)+dy(1:im,2:jm)); 
arv(1:im,1)=arv(1:im,2);
%
%Specify channel topography, h (= depth when there is no motion):
hw = hc; he = hc;       %h@west  and east
hs = hc; hn = hc;       %h@south and north
dhx= (he-hw)/xl;        %---depth-increament for x-linear h (below)
dhy= (hn-hs)/yl;        %---depth-increament for y-linear h (below)
h=zeros(im,jm)+small;   %---initialize w/land points everywhere
for i = 1:im;   h(i,:) = hw+dhx.*(i-1)*dx(i,:); end; 
for j = 1:jmm1; h(:,j) = hs+dhy.*(j-1)*dy(:,j); end;
%
%{
%progname=['KelvinWaveCoast06'];
%Put h(x) cosine bump centered at 2*im/3 (x=600km), wavelength=600km;
wln=6.e5/dxc; a=6.0;
x=zeros(im,jm); y=zeros(im,jm); x1=0.0:1.0:im-1; y1=0.0:1.0:jm-1; 
for j=1:jm; x(:,j)=x1(:); end; for i=1:im; y(i,:)=y1(:); end;
ic=round(2*im/3); xc=ic-1; ys=y1(2); yn=y1(jm-1); 
xmax=0.25*wln; %(K*W*02); %cos pulse S&N 
%cosx=zeros(im,jm); cosx=cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax); %06-11
%12,13,14,15,17,18,19,20,21,22,23: Gaussian, next line:
ic=round(im/2); xc=ic-1; cosx=zeros(im,jm); cosx=exp(-((x-xc)/(0.25*wln)).^2); 
%16: Gaussian, next line:
%ic=round(im/2); xc=ic-1; cosx=zeros(im,jm); cosx=exp(-((x-xc)/(0.5*wln)).^2); 
%h=hc-0.8*hc*cosx.*(exp((ys-y)/a)); %KelvinWaveCoast06
h=hc-0.8*hc*cosx;                   %KelvinWaveCoast11,12,13,14,15,16,17,18,19,20,21,22,23
%}
%
%h(1,:) = small; h(im,:) = small; %coast west  & east
 h(:,1) = small; %ForcedKW:h(:,jm) = small; %coast south & north %01,02,03,04,07
%h(:,1) = small; h(:,jm) = small; %coast south & north %05,06,08
%
%Define land masks for elevation & velocity, using the "h" above: 
%   fsm & dum =1 at water points, =0 at land/"wall" points
%   fsm is for elevation point; dum is for ua-point
fsm(1:im,1:jm)=0.0; dum=fsm; dvm=fsm;         %---initialize
for j=1:jm; for i=1:im; if(h(i,j)>small)fsm(i,j)=1.0; end;end;end;
%
%Take care of periodic boundary to be compatible with the above mask:
if (nperx == 1 && im > 3);  %x-periodic: 2 west cells are land if 1 is
for j=1:jm; if (fsm(1,j) == 0.0); fsm(2,j)=0.0; end; end; 
end;
if (npery == 1 && jm > 3);  %y-periodic: 2 south cells are land if 1 is
for i=1:im; if (fsm(i,1) == 0.0); fsm(i,2)=0.0; end; end; 
end;
%
dum(2:im,:) = fsm(1:imm1,:).*fsm(2:im,:); dum(1,:)=dum(2,:);
dvm(:,2:jm) = fsm(:,1:jmm1).*fsm(:,2:jm); dvm(:,1)=dvm(:,2);
%
h(:,:)=fsm(:,:).*h(:,:)+small; %Redo "h" again to make sure lands are set
%
%Coriolis parameter:
lat0 =43.5; %25.0;              %---Model domain center latitude
                                %lat0=43.5 gives f = 1.e-4
cor0 = 4.*pi*sind(lat0)/86400.; %---Coriolis @domain center
cor0 = 1.15740740741/60000.;    %set to be c/60km, so Rossby a = 60km,
                                % where "c" is from O_modparam
                                % So, lat0 (above) is ignored, but it is:
                                % lat0=asind(86400*cor0/(4*pi)); =7.6215N
%cor0 = 10./100000.;             %09:lat0 = 43.4364307146557
cor0 = 1.*cor0;                 %---Make earth spin faster
cor=zeros(im,jm)+cor0;          %---Coriolis
%
%Calculate time step "dt" satisfying CFL:
%   CFL condition is: sqrt(gh)*dt/dx <= 1/2 --> 2d version:
%       dtcfl = 1./(2.*sqrt(hmax*g)*sqrt(1./dxmin^2+1./dymin^2));
hmax = O_maxval(h); dxmin = O_minval(dx); dymin = O_minval(dy);
dtcfl = dtfac*(1./(2.*sqrt(hmax*g)*sqrt(1./dxmin^2+1./dymin^2)));
dt = floor(dtcfl); %dt = floor(0.5*dxmin/sqrt(hmax*g)); %-->1d CFL
%
return;
%}
