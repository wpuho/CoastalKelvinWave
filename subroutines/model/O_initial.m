function [elf,uaf,vaf,elb,eln,uab,uan,vab,van,d,bndwe,bndsn]=O_initial(h)
%
[im,jm]=size(h);
%
elb(1:im,1:jm)=0.0; eln=elb; elf=elb; 
uab(1:im,1:jm)=0.0; uan=uab; uaf=uab; 
vab(1:im,1:jm)=0.0; van=vab; vaf=vab; 
d(1:im,1:jm)=h(1:im,1:jm)+eln(1:im,1:jm);
%{
%KelvinWavesCoast:
%The followings assume dxc=dyc=10km, Rossby rad a=60km & hc=10m in O_grid.m
%   So a (=a/dy) = 6 is used below;  Also c=1.15740740741 m/s
a = 60/10; hc=10; c=1.15740740741;
%a = 60/2.5;  %04: dy=2.5km
el0=0.0; ua0=(el0*c/hc); wln=300/10; %wln=wavelength; el0=amplitude in m
%wln=300/2.5; %03 This has dx=2.5km instead of 10km
%ForcedKW:wln=300/1.25; %This has dx=1.25km instead of 10km
%wln=300/0.625; %This has dx=0.625km instead of 10km
x=zeros(im,jm); y=zeros(im,jm); x1=0.0:1.0:im-1; y1=0.0:1.0:jm-1; 
for j=1:jm; x(:,j)=x1(:); end; for i=1:im; y(i,:)=y1(:); end;
ic=round(im/2); ic=round(im/4);
xc=ic-1; ys=y1(2); yn=y1(jm-1); 
%KW @ S & N. Just remove the 2nd exp(...) for S only, i.e. set ncoast=0
ncoast=0.0; %for S only
%ncoast=1.0; %for S & N
xmax=0.25*wln; %(K*W*02); %cos pulse S&N 
%xmax=1.e18; %progname=['KelvinWaveCoast01'] & K*W*03; %cosine wave S&N
cosx=zeros(im,jm); cosx=cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax); %cosx=0*(abs(x-xc)>=xmax) + cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax);
elb=el0*cosx.*(exp((ys-y)/a)-ncoast*exp((y-yn)/a)); %elb=el0*cos(2*pi*(xmxc)/wln).*(exp((ys-y)/a)-exp((y-yn)/a));
%cosx=zeros(im,jm); cosx=cos(2*pi*(x-xc-0.5)/wln).*(abs(x-xc-0.5)<xmax);
uab=ua0*cosx.*(exp((ys-y)/a)+ncoast*exp((y-yn)/a)); %uab=ua0*cos(2*pi*(xmxc)/wln).*(exp((ys-y)/a)+exp((y-yn)/a));
uab(1:im-1,:)=0.5*(uab(1:im-1,:)+uab(2:im,:));
eln=elb; uan=uab; d=h+eln;
%
%---figure(1); contour(eln,'ShowText','on');
%}
%%
%Specify boundary elev/vels (here = initial cond), to be used in O_bc.m:
bndwe=zeros(4,jm); bndsn=zeros(im,4); %1:2 for elev, 3:4 for velo
bndwe(1,:)=eln(1,:); bndwe(2,:)=eln(im,:);
%
%%return;

