%Changes from ATOP linus version (2019/May/25), search "%lyoChanged"
% default is nps=-2, i.e. png instead of ps; png does not work on windows
function [plotfile]=O_xycontour04(x,y,z,     ...
	   xrange,yrange,zrange,xi,yi,fsm,nps,pltitle,plfile,pvel,fig,wh)
%
% A more advanced version of leoxycontour03.m w/pvel
%
%---------------------Plot xy-contours of "z"-------------------------
%
  % x,y,z must be defined on input:
% x = 1d array (= row vector of length 'im', i.e. w/im columns);
% y = 1d array (= row vector of length 'jm', i.e. w/jm columns);
% OR x & y are 2d (im,jm) = positions of "z(im,jm)" %lyo20200430:
% z(im,jm) = 2d array with im rows & jm columns;
%
% Other variables have below defaults if undefined (NaN or empty string):
% xrange(1:2) = min & max of x; (1:2) MUST be undefined (or defined);
% yrange(1:2) = min & max of y; (1:2) MUST be undefined (or defined);
% eg if xrange(1) is (un)defined, then xrange(2) must also be (un)defined
% xi,yi = 10 equal x,y axes tickmarks (intervals);
% fsm(im,jm) = land mask to shade land gray, deafult is no shading
% nps =-1 ps, =-2 png; default -1; =0 none; =1,2,.. "animation", -append
% pltitle = plot title, default = 'O_xycontour04';
% plfile = plot output filename, default = 'O_xycontour04'(.*)
% pvel = vel vector input as complex array = u + 1i*v
% fig = must be given for "animation" if nps > 0; otherwise error
% wh = h (ocean depth) is used for 3D surf plot only
%
%---------------------User's functions required:--------------------- 
% 1. O_landtogray.m;
% 2. O_maxval.m;
% 3. O_minval.m;
%--------------------------------------------------------------------
%
%Set x,yrange, pltitle etc, if necessary:
[im,jm]=size(z); wrk=zeros(size(z));
if (isnan(xrange));
xrange(1:2)=0.0; xrange(1)=O_minval(x); xrange(2)=O_maxval(x); end;
if (isnan(yrange));
yrange(1:2)=0.0; yrange(1)=O_minval(y); yrange(2)=O_maxval(y); end;
if (isnan(zrange));
%zrange(1:2)=0.0; zrange(1)=round(O_minval(z)); zrange(2)=round(O_maxval(z)); end;
zrange(1:2)=0.0; zrange(1)=(O_minval(z)); zrange(2)=(O_maxval(z)); end;
if (isnan(xi)) xi=(xrange(2)-xrange(1))/10.; end;
if (isnan(yi)) yi=(yrange(2)-yrange(1))/10.; end;
if (length(pltitle) == 0); pltitle=['O_xycontour04']; end;
if (isnan(nps)); nps=-2; end; %lyoChanged
if (length(plfile) == 0 ); plfile=['./plots/' 'O_xycontour04']; end;
%
%Contour interval:
zmax=zrange(2); zmin=zrange(1);
nc=8; ci(1:nc+1)=0.0; dz=(zmax-zmin)/nc;
for i=1:nc+1; ci(i)=zmin+(i-1)*dz; end;
%Get rid of value near zero:
cin=zeros(size(ci))+99999;  
for i=1:nc+1;if(abs(ci(i))/abs(zmax)>0.05);cin(i)=ci(i);end;end;
kc=0;for i=1:nc+1;if(cin(i)~=99999);kc=kc+1;cit(kc)=cin(i);end;end;%
if(exist('cit'));ci=cit;end;
cmap=O_lbmap(nc,'jrcmap');
%
%lyo20200430:beg---
if (numel(x) == numel(z)); 
    xx=x; yy=y;
else %Create matrices w/rows = x & cols = y, to match z(im,jm):
    [yy xx] = meshgrid(y,x);
end;
%lyo20200430:end---
%
fonts=11; %16; %if(nps > 0);fonts=24;end;%larger on screeen but doesn't work

%pause;
%S = ['leoxycontour03, fsm(1,1) = ' num2str(fsm(1,1)) '******************']; disp(S)
%xrange(1),xrange(2),yrange(1),yrange(2),zrange(1),zrange(2)
%pause;
%???if (fsm(1,1)==9);
if (exist('wh')) ;
wrk=z; if (jm == 3) wrk(:,3)=wrk(:,2); end;
surf(xx,yy,wrk,'FaceAlpha',0.75);  
%Above is default: view([azimuth elevation]); azimuth=-37.5, elevation=30
%--surf(xx,yy,z,'FaceAlpha',0.75);
shading interp; caxis([zrange(1) zrange(2)]); %Nhi:caxis([-3.0 3.0]); %
%daspect([1.2 1 0.01]); %Nhi:colormap(jet); c = colorbar; c.Label.String = 'eta (m)';
daspect([1 1 0.01]);    %Tsunami1D: to enlarge "z-surface" nicer
%daspect([1 1 0.1]);    %Tide1D:
%view(18,-30);  %
%---view(18,20); %view([azimuth elevation]);
view(144,20);

axis([xrange(1) xrange(2) yrange(1) yrange(2) zrange(1) zrange(2)]);
set(gca,'ztick',[zrange(1):2.*dz:zrange(2)]);

hold on;
wrk=wh; if (jm == 3) wrk(:,3)=wrk(:,2); end;
surf(xx,yy,wrk,'FaceAlpha',0.95); shading interp;
colormap(cmap);
if (exist('pvel')) ;
uscalei=1./real(pvel(1,1)); vscalei=1./imag(pvel(1,1));
pvel(1,1)=pvel(1,2);
xxs=NaN(im,jm);yys=NaN(im,jm);zzs=NaN(im,jm);uus=NaN(im,jm);vvs=NaN(im,jm);
iskip=2; jskip=2; if (im>180) iskip=5; end; if (jm>180) jskip=5; end;
%for i=1:2:im; for j=1:2:jm; %skip every 2 grids
for i=1:iskip:im; for j=1:jskip:jm; %skip every 5 grids in x for %18
xxs(i,j)=xx(i,j);yys(i,j)=yy(i,j);zzs(i,j)=z(i,j);
uus(i,j)=real(pvel(i,j));vvs(i,j)=imag(pvel(i,j));
end; end;
quiver3(xxs,yys,zzs,uus,vvs,zeros(im,jm),2.,'-k');
%    plotvec=O_psliceuv(x,y,pvel,3,2,'black');
%   plotvec=O_psliceuv(x,y,pvel,1,.2,'black');
pvel=real(pvel)*uscalei+1i*imag(pvel)*vscalei;
spdmax=O_maxval(abs(pvel));
%For "view(18,20)" in L86:
%---text(xrange(1)+.05*(xrange(2)-xrange(1)), ...
%--- yrange(2)-.05*(yrange(2)-yrange(1)), ...
%--- zrange(2)-.10*(zrange(2)-zrange(1)), ...
%--- ['MaxSpeed = ' num2str(spdmax,'%10.2e\n') ' m/s'],'Color','black', ...
%--- 'FontSize',12);end;
%For "view(144,20)" in L87:
%zrange(2)=2.5:%text(xrange(2)*1.1, yrange(2)-.05*(yrange(2)-yrange(1)), 2.5*zrange(2), ...
%zrange(2)=2:
text(xrange(2)*1.1, yrange(2)-.05*(yrange(2)-yrange(1)), 2.8*zrange(2), ...
 ['MaxSpeed = ' num2str(spdmax,'%10.2e\n') ' m/s'],'Color','black', ...
 'FontSize',12);end;
%
else;
contourf(xx,yy,z,[ci],'LineColor','none');
%--hold on; contour(xx,yy,z,[2.*ci],'LineColor','k','ShowText','on'); %[4.*ci]
%--hold on; contour(xx,yy,z,[2.*ci],'LineColor','k'); %[4.*ci]
hold on; [c,h]=contour(xx,yy,z,[2.*ci],'LineColor','k'); 
h.LevelList=round(h.LevelList,2);  %rounds levels to 3rd decimal place
clabel(c,h);
%contourf(xx,yy,z,[ci],'LineColor','k');
%fonts=16; if(nps > 0);fonts=16;end;%larger on screeen but doesn't work
%set(gca,'fontsize',fonts)  %x & y axis legends
%title(pltitle, 'fontsize',fonts);
axis([xrange(1) xrange(2) yrange(1) yrange(2)]);
if (abs(xrange(2)-xrange(1)) == abs(yrange(2)-yrange(1))); axis equal; 
    fonts=1.0*fonts; end;  %fonts=0.7*fonts; end;
set(gca,'xtick',[xrange(1):xi:xrange(2)]);
set(gca,'ytick',[yrange(1):yi:yrange(2)]);
colormap(cmap);
%
%Plot velocity vectors if 'pvel' exists:
if (exist('pvel')) ;
uscalei=1./real(pvel(1,1)); vscalei=1./imag(pvel(1,1));
pvel(1,1)=pvel(1,2);
    plotvec=O_psliceuv(xx,yy,pvel,3,2,'black');
%   plotvec=O_psliceuv(x,y,pvel,1,.2,'black');
pvel=real(pvel)*uscalei+1i*imag(pvel)*vscalei;
spdmax=O_maxval(abs(pvel));
text(xrange(1)+.05*(xrange(2)-xrange(1)), ...
 yrange(2)-.05*(yrange(2)-yrange(1)), ...
 ['MaxSpeed = ' num2str(spdmax,'%10.2e\n') ' m/s'],'Color','black', ...
 'FontSize',12);
end;
%
end;
%
set(gca,'fontsize',fonts)  %x & y axis legends
if (exist('fig')) ;
title(pltitle, 'fontsize',fonts);
else;
%--titpos=1.07*yrange(2);
%--text(0.5*(xrange(1)+xrange(2)),titpos,pltitle, ...
%--    'HorizontalAlignment','center','FontSize',fonts);
annotation('textbox',[.001 .05 .96 .96],'String',pltitle,'EdgeColor','none');
end;
colorbar('southoutside'); caxis([zmin zmax]); 
hold on;
%
if (isnan(fsm));
%do nothing:
else;
if (fsm(1,1) ~= 9);  %(fsm(1,1)==0 | fsm(1,1)==1);
[lon_land_box,lat_land_box]=O_landtogray(xx,yy,fsm,im,jm);
fill_h=fill(lon_land_box,lat_land_box,[0.5 0.5 0.5]); %0.5=gray for land
set(fill_h, 'linestyle','none');
end;
end;
%
if (nps <= 0);
%
%nps <= 0: Invidual, e.g. ps (png or no) file-----------------------
if (nps == 0 ); plotfile = [ 'No Plot File' ]; return; end;
if (nps == -1); plotfile = [plfile '.ps' ];
print('-dpsc', '-loose', plotfile);
end;
if (nps == -2); plotfile = [plfile '.png'];
print('-dpng', '-loose', plotfile);
end;
else;       %---if (nps <= 0);
%
%nps > 0: Animation file--------------------------------------------
%GIF containing a series of images all combined into one file, see:
%https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab
frame = getframe(fig); lm = frame2im(frame); [imind,cm] = rgb2ind(lm,256); 
if nps==1;  %---1st time in this function
imwrite(imind,cm,plfile,'gif', 'Loopcount',inf); 
else;
imwrite(imind,cm,plfile,'gif','WriteMode','append'); 
end;        %---if nps==1;
plotfile = [plfile];
end;        %---if (nps <= 0);
%
hold off;
%
return;
%
