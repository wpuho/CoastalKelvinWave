function [plotfile]=O_xycontour01(x,y,z)
%
% Code works if called from MAIN e.g.:
%    [plotfile] = O_xycontour01(xpl,tpl,elout);
%
%---------------------Plot xy-contours of "z"-------------------------
%
  % x,y,z must be defined on input:
% x = 1d array (= row vector of length 'im', i.e. w/im columns);
% y = 1d array (= row vector of length 'jm', i.e. w/jm columns);
% z(im,jm) = 2d array with im rows & jm columns;
%
% Other variables take on the below defaults if undefined:
% xrange(1:2) = min & max of x; (1:2) MUST be undefined (or defined);
% yrange(1:2) = min & max of y; (1:2) MUST be undefined (or defined);
% eg if xrange(1) is (un)defined, then xrange(2) must also be (un)defined
% xi,yi = 10 equal x,y axes tickmarks (intervals);
% pltitle = plot title, default = 'O_xycontour01';
% fsm(im,jm) = land mask to shade land gray, deafult is no shading
% nps = 1 ps, set = 2 for an additional png; default is 1;
% plfile = plot output filename, default = 'O_xycontour01'(.*)
%
%---------------------User's functions required:--------------------- 
% 1. O_landtogray.m;
% 2. O_maxval.m;
% 3. O_minval.m;
%--------------------------------------------------------------------
%
plotdir=['./plots/']; mkdir(plotdir);
%Set x,yrange, pltitle etc, if necessary:
[im,jm]=size(z); wrk=zeros(size(z));
if (exist('xrange') == 0); 
xrange(1:2)=0.0; xrange(1)=O_minval(x); xrange(2)=O_maxval(x); end;
if (exist('yrange') == 0); 
yrange(1:2)=0.0; yrange(1)=O_minval(y); yrange(2)=O_maxval(y); end;
if (exist('xi') == 0) xi=(xrange(2)-xrange(1))/10.; end;
if (exist('yi') == 0) yi=(yrange(2)-yrange(1))/10.; end;
if (exist('pltitle') == 0); pltitle=['O_xycontour01']; end;
if (exist('nps') == 0); nps=2; end;
if (exist('plfile') == 0); plfile=['O_xycontour01']; end;
%
%Contour interval:
%???zmax=round(O_maxval(z)); zmin=round(O_minval(z));
zmax=(O_maxval(z)); zmin=(O_minval(z));
nc=8; ci(1:nc+1)=0.0; dz=(zmax-zmin)/nc;
for i=1:nc+1; ci(i)=zmin+(i-1)*dz; end;
%
%Create matrices w/rows = x & cols = y, to match z(im,jm):
[yy xx] = meshgrid(y,x);
%???[xx yy] = meshgrid(x,y);
%
contourf(xx,yy,z,[ci],'LineColor','none');
colorbar; caxis([zmin zmax]);
set(gca,'fontsize',16)  %x & y axis legends
title(pltitle, 'fontsize',16);
axis([xrange(1) xrange(2) yrange(1) yrange(2)]);
set(gca,'xtick',[0:xi:xrange(2)]);
set(gca,'ytick',[0:yi:yrange(2)]);
hold on;
%
if (exist('fsm') == 1);
[lon_land_box,lat_land_box]=O_landtogray(xx,yy,fsm,im,jm);
fill_h=fill(lon_land_box,lat_land_box,[0.5 0.5 0.5]); %0.5=gray for land
set(fill_h, 'linestyle','none');
end;
%
if (nps == 1);
plotfile = [plfile '.ps']; plotfile=[plotdir plotfile];
print('-dpsc', '-loose', plotfile);
end;
if (nps == 2);
plotfile = [plfile '.png']; plotfile=[plotdir plotfile];
print('-dpng', '-loose', plotfile);
%%if (nps == 2); print('-dpng', '-loose', [plfile '.png']); end;
%if (nps == 2); saveas(gcf,[plfile '.png']); end;
end;
%
hold off;
%
return;
%
