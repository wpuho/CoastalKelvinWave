function [ela,wus,wvs,uwnd,vwnd]=O_surfforce(dx,dy,n,dt)
%
%   ela = inverse barometer (m)
%   wus = minus u-windstress (m/s)^2 @u point, LHS of x-momentum eqn;
%   wvs = minus v-windstress (m/s)^2 @v point, LHS of y-momentum eqn;
%   uwnd = + u-wind (m/s) @el point;
%   vwnd = + v-wind (m/s) @el point;

[im,jm]=size(dx); ela=zeros(im,jm); wus=zeros(im,jm); wvs=zeros(im,jm);
                  uwnd=zeros(im,jm); vwnd=zeros(im,jm);
                  uwnd(:,31:61) = 5.0;
%{
%ForcedKW:Comment out next 4 lines originally from StormSurgeLocal
%ForcedKW:%x-wind = 5m/s ramping over 1 day: 
%ForcedKW:uwnd(:,:)=5.0;
%ForcedKW:%rampt=(n*dt/86400.); if (rampt<=1.0); uwnd(:,:)=rampt*5.0; end;
%ForcedKW:rampt=(n*dt/86400.); if (rampt<0.0); uwnd(:,:)=rampt*5.0; end;
%
%Put wind nwindday days after t=0 to allow initial eddy to form:
nmoveday=0;     %02-08:
%nmoveday=3*86400/dt;  %09: wind moves after one day
nwindday=0; %ForcedKWCoast01:nwindday=999*86400/dt;
if (n>nwindday);
%
ic=round(im/4); %wln=300/10; %ForcedKWCoast02:wln & ic same as in *initial*
wln=3.e5; %has to be in meters
utc=0.0;
if (n>nmoveday);
g=0.133959190672; hc=10.; 
%g=10; hc=10.; %09:
utc=sqrt(g*hc); %ForcedKWCoast02:g & hc same as in *modparam* & *grid*
end;
%utc=0; %08:
%
x=zeros(im,jm); x1=0.0:1.0:im-1; %y=zeros(im,jm); y1=0.0:1.0:jm-1;
for j=1:jm; x(:,j)=x1(:).*dx(:,j); end; %for i=1:im; y(i,:)=y1(:).*dy(i,:); end;
%xc=(ic-1)+(n-nwindday)*dt*utc; %ys=y1(2); yn=y1(jm-1); 
xc=(ic-1)*dx+(n-nwindday)*dt*utc;
xmax=0.25*wln; %(K*W*02); %cos pulse S&N 
uwnd=1.581*cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax); %06,07:stress=5.e-6
%uwnd=4.*4.64717194000958*cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax); 
        %09:Increase wind speed to get the same response for g=10:
        %   1.581*(10/0.133959190672)^(1/4)
%uwnd=1.581*2*cos(2*pi*(x-xc)/wln).*(abs(x-xc)<xmax); %08:stress=2.e-5
%uwnd=1.5*(abs(x-xc)<xmax); %02,03,04,05:stress=4.5e-6
%}
%{
%Inverse Barometer & Wind: 
vtc=+0.0;
utc=+10.0; rtc=1.5e5;     %TC translation vel & radius
pc=950.; p0=1010.;        %TC center and ambient pressures
rtcsqi=1./rtc^2;
for j=1:jm; for i=1:im; 
        xc=(5)*dx(i,j)+(n-nwindday)*dt*utc; yc=(round(jm/2)- 1)*dy(i,j);
%***    xc=(im-3)*dx(i,j)+(n-nwindday)*dt*utc; yc=(round(jm/2)- 1)*dy(i,j);
        x=(i-1)*dx(i,j); y=(j-1)*dy(i,j); %assume constant dx & dy
%***    expon=exp(-((x-xc)^2 + (y-yc)^2)*rtcsqi); 
%***    ela(i,j)=p0+(pc-p0)*expon;                %pressure, Gaussian;
%***    uwnd(i,j)= (pc-p0)*rtcsqi*2.*(y-yc)*expon; %uwnd=-dp/dy
%***    vwnd(i,j)=-(pc-p0)*rtcsqi*2.*(x-xc)*expon; %vwnd=+dp/dx
if (abs(x-xc)<=2) uwnd(i,j)= 10.0; end;
end; end;
%ela=zeros(im,jm); %-(1.e-4)*ela;   %inverse barometer, where 1.e-4=1/(rhos*g)
%ela=-(1.e-4)*ela;   %inverse barometer, where 1.e-4=1/(rhos*g)
ela=-(0.0746)*ela;   %inverse barometer, where 0.0746=(1/(rhos*g))*(g/g')
                     % i.e. use g' is more correct
%Adjust wind speed so Max= AtkinsonHoliday*MWR1977; 0.514444*6.7 kt to m/s
%***vmxadjust=(3.4468*(1010-pc)^0.644)/(O_maxval(sqrt(uwnd.^2+vwnd.^2))+1.e-12);
%***uwnd=uwnd*vmxadjust; vwnd=vwnd*vmxadjust;
%
%}
%Windstress:ForcedKW:Keep yinoey option from StormSurgeLocal
yinoey=0; %=1 if Yin & Oey (2006) formula is used
if (yinoey == 0);
%Use simple formula, --> 0.5e-4 kinematic stress (m/2)^2 for 5m/s wind
for j=1:jm; for i=1:im;
 spd = sqrt( uwnd(i,j)^2 + vwnd(i,j)^2 );
 wus(i,j) = -2.e-6 * spd * uwnd(i,j);
 wvs(i,j) = -2.e-6 * spd * vwnd(i,j);
 end; end;
else;
%Use Yin & Oey:
for j=1:jm; for i=1:im;
 spd = sqrt( uwnd(i,j)^2 + vwnd(i,j)^2 );
 if (spd <= 11.)      cda = 0.0012;
 elseif (spd <= 19.)  cda = 0.00049 + 0.000065*spd;
 elseif (spd <= 100.) cda = 0.001364 + 0.0000234*spd - 2.31579e-7*spd^2;
 else                 cda = 0.00138821;
 end;
%1.19675654e-3 = rhoa/rhow = 1.225/1023.6
wus(i,j) = -1.19675654e-3 * cda * spd * uwnd(i,j);
wvs(i,j) = -1.19675654e-3 * cda * spd * vwnd(i,j);
end; end;
end;
%
%Transfer to u,v points:
%--uwnd(2:im,:)=0.5*(uwnd(1:im-1,:)+uwnd(2:im,:)); %u-point
%--vwnd(:,2:jm)=0.5*(vwnd(:,1:jm-1)+vwnd(:,2:jm)); %v-point
wus(2:im,:)=0.5*(wus(1:im-1,:)+wus(2:im,:)); %u-point
wvs(:,2:jm)=0.5*(wvs(:,1:jm-1)+wvs(:,2:jm)); %v-point
%
%end;
%
return

