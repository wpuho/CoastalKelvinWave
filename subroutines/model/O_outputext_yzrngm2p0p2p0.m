function [nfig]=O_outputext(progname,elout,uaout,vaout,outtim, ...
                            cor0,dx,dy,fsm,h,ana,fileid)
%% Plot & Print Outputs for "external" or "barotropic" Variables:
%
% Search for "How to". Some turned off using "%{ ... %}".
% Plot outputs in subfolder "./plots"
%
plotdir=['./plots/']; mkdir(plotdir); 
%
[im,jm,nm]=size(elout); imm1=im-1; jmm1=jm-1;
%Plot parameters & constants:
nfig = 0;       %initialize figure# displayed on screen: figure(nfig+1) etc
%Tsunami2D:---beg
zscale=3.0;      %4.0;     %5.; elevation min/max scale for plots
%Tsunami2D:---end
dayi=1./86400.;
jc=round(jm/2); xl=sum(dx(2:im,jc)); 
ic=round(im/2); yl=sum(dy(ic,2:jm));

%Move velocity to elevation point:
uaout(1:imm1,:,:) = 0.5*(uaout(1:imm1,:,:)+uaout(2:im,:,:)); 
uaout(im    ,:,:) = uaout(imm1,:,:);
vaout(:,1:jmm1,:) = 0.5*(vaout(:,1:jmm1,:)+vaout(:,2:jm,:)); 
vaout(:,jm    ,:) = vaout(:,jmm1,:);
%
%Define plot axes: x, y in km @elevation points, & time 'tpl' in days:
xpl(1:im)=0.0; ypl(1:jm)=0.0; tpl(1:nm)=0.0; 
%xpl(2)=0.5*dx(2);for i=3:im;xpl(i)=xpl(i-1)+0.5*(dx(i-1,jc)+dx(i,jc));end; %x=0@ua(2,:)
for i=2:im;xpl(i)=xpl(i-1)+0.5*(dx(i-1,jc)+dx(i,jc));end; %x=0@el(1,:) i=1 el-point
%?-for i=2:im;xpl(i)=xpl(i-1)+dx(i-1,jc);end; %x=0@ua(1,:) i=1 ua-point
xpl(:)=xpl(:)*1.e-3;
%ypl(2)=0.5*dy(2);for j=3:jm;ypl(j)=ypl(j-1)+0.5*(dy(ic,j-1)+dy(ic,j));end; %y=0@va(:,2)
for j=2:jm;ypl(j)=ypl(j-1)+0.5*(dy(ic,j-1)+dy(ic,j));end; %y=0@el(:,1) j=1 el-point
%?-for j=2:jm;ypl(j)=ypl(j-1)+dy(ic,j-1);end; %y=0@va(:,1) j=1 va-point
ypl(:)=ypl(:)*1.e-3;
for n=1:nm; tpl(n)=(n-1)*outtim*dayi; end;
%
%-------------------------------------------------------------------------
%How to use O_xycontour01: simple contour, @time near tpl(nm/2) day
%{
nfig=nfig+1; figure(nfig);
pelout=zeros(im,jm);n=round(nm/2);pelout(:,:)=elout(:,:,n);pelout(fsm==0)=NaN;
[plotfile]=O_xycontour01(xpl,ypl,pelout); %-> O_xycontour01.png
%}
%-------------------------------------------------------------------------
%How to use O_xycontour04 to make xt-contours:
nfig=nfig+1; figure(nfig);
pelout=zeros(im,nm);fsmxt=zeros(im,nm);j=2; %round(jm/2);
for i=1:im;for n=1:nm;fsmxt(i,n)=fsm(i,j);pelout(i,n)=elout(i,j,n);end;end;
pelout(fsmxt==0)=NaN;
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl);xrng(2)=O_maxval(xpl);
yrng(1)=O_minval(tpl);yrng(2)=O_maxval(tpl);
zrng(1)=O_minval(pelout);zrng(2)=O_maxval(pelout);
%zrng(1)=-2.;zrng(2)=+2.;
wrk(:,:)=fsmxt(:,:);
plfil=[progname 'xt-' 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(1),'%10.3f\n') ...
'-' num2str(tpl(nm),'%10.3f\n')]; pltit = plfil;
plfil=[plotdir plfil];
[plotfile] = O_xycontour04(xpl,tpl,pelout, ...
	    xrng,yrng,zrng,xi,yi,wrk,nps,pltit,plfil);%,''); 
            %No pvelxy (vectors), No fig (animation), No wh (surf)
if exist(ana,'file') == 2;
[ael,aua]=O_anatidechannel(elout,xpl,tpl,dx,h);
nfig=nfig+1; figure(nfig);
plfil=[progname 'ana_xt-' 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(1),'%10.3f\n') ...
'-' num2str(tpl(nm),'%10.3f\n')]; pltit = plfil;
plfil=[plotdir plfil];
[plotfile] = O_xycontour04(xpl,tpl,ael, ...
	    xrng,yrng,zrng,xi,yi,wrk,nps,pltit,plfil);%,''); 

end;
%
%-------------------------------------------------------------------------
%How to use O_xycontour04 w/o velocity vectors:
%{
nfig=nfig+1; figure(nfig);
pelout=zeros(im,jm);n=round(nm/2);pelout(:,:)=elout(:,:,n);pelout(fsm==0)=NaN;
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl);xrng(2)=O_maxval(xpl);
yrng(1)=O_minval(ypl);yrng(2)=O_maxval(ypl);
wrk(:,:)=fsm(:,:);
plfil=[progname 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n),'%10.3f\n')]; pltit = plfil;
plfil=[plotdir plfil];
[plotfile] = O_xycontour04(xpl,ypl,pelout, ...
	    xrng,yrng,zrng,xi,yi,wrk,nps,pltit,'');%,''); %-> O_xycontour04.png
            %No pvelxy (vectors), No fig (animation), No wh (surf)
%}
%-------------------------------------------------------------------------
%How to use O_xycontour04 w/vectors to make a series of *png (or *ps) files 
% which then maybe used to make animation gif file w/this linux command:
% "convert -crop 0x0 -delay 50 -loop 20 WhateverFile.ps MaybeOtherFile.gif", or
% "convert -delay 120 -loop 0 WhateverFile.gif MaybeOtherFile.gif"
%{%
%The follwings make animation by *ps etc, but replaced by direct gif next
nfig=nfig+1; figure(nfig);
pelout=zeros(im,jm);
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl);xrng(2)=O_maxval(xpl);
yrng(1)=O_minval(ypl);yrng(2)=O_maxval(ypl);
zrng(1)=-zscale;zrng(2)=+zscale/3.;
zrng(1)=-1.6;zrng(2)=1.6;
xpl2=zeros(im,jm); ypl2=zeros(im,jm);
for j=1:jm; xpl2(:,j)=xpl(:); end; for i=1:im; ypl2(i,:)=ypl(:); end;
wrk(:,:)=fsm(:,:);
nps = -2; %png files;  
%nps = -1; %ps files;
velsc=2.*max(O_maxval(sqrt(uaout.^2+vaout.^2))); %=0.5;
uscale=(xl./O_maxval(dx))/velsc; %1 grid cell is 0.5 m/s;
vscale=(yl./O_maxval(dy))/velsc; %1 grid cell is 0.5 m/s;
n1=1 ;ni= 4;n2=n1+ni;   %n1=nm ;ni= 8;n2=nm; 
for n=n1:ni:n2;
pelout(:,:)=elout(:,:,n); pelout(fsm==0)=NaN;
plfil=[progname 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n),'%10.3f\n')]; pltit = plfil;
plfil=[plotdir plfil];
%velsc=2.*max(O_maxval(sqrt(uaout.^2+vaout.^2))); %=0.5;
%uscale=(xl./O_maxval(dx))/velsc; %1 grid cell is 0.5 m/s;
%vscale=(yl./O_maxval(dy))/velsc; %1 grid cell is 0.5 m/s;
pvelxy(:,:)=uscale*uaout(:,:,n)+vscale*1i*vaout(:,:,n); %the w-component = 0;
pvelxy(1,1)=uscale+1i*vscale;
[plotfile] = O_xycontour04(xpl2,ypl2,pelout, ...
	    xrng,yrng,zrng,xi,yi,wrk,nps,pltit,plfil,pvelxy);%,'');
            %Yes pvelxy (vectors), No fig (animation), No wh (surf)
end;     %---for n=n1:...
%}%
%
%-------------------------------------------------------------------------
%How to use O_xycontour04 to directly make perspective surf animation gif:
%   The animation *gif can be split into invidual gif files by uploading 
%   it into this online program: https://ezgif.com/split
%
nfig=nfig+1; figure(nfig);
pelout=zeros(im,jm);
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl);xrng(2)=O_maxval(xpl);
yrng(1)=O_minval(ypl);yrng(2)=O_maxval(ypl);
zrng(1)=-1.6;zrng(2)=1.6; %zrng(1)=-zscale;zrng(2)=+zscale/3.;
zrng(1)=-2.0;zrng(2)=2.0;
xpl2=zeros(im,jm); ypl2=zeros(im,jm);
for j=1:jm; xpl2(:,j)=xpl(:); end; for i=1:im; ypl2(i,:)=ypl(:); end;
wrk(:,:)=fsm(:,:);
%
%???wrk(1,1)=9; %---Set (1,1) to 9 for surface plot animation
wh=h; wh(fsm==0)=NaN;  wh=-0.75*zscale*wh./O_maxval(wh);
%pause;
%S=['Channel2d_02, wrk(1,1) = ' num2str(wrk(1,1)) '******************']; disp(S)
%pause;
n1=1;ni=1;n2=nm;    %n1=1;ni=8;n2=nm; 
fig = figure(nfig);
axis tight manual % this ensures that getframe() returns a consistent size
plfil=[progname 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n1),'%10.3f\n') ...
'to' num2str(tpl(n2),'%10.3f\n') 'ni' num2str(ni)];
%                              
plfil=[plotdir plfil];
animfile = [plfil '.gif'];
velsc=2.*max(O_maxval(sqrt(uaout.^2+vaout.^2))); %=0.5;
uscale=(xl./O_maxval(dx))/velsc; %1 grid cell is 0.5 m/s;
vscale=(yl./O_maxval(dy))/velsc; %1 grid cell is 0.5 m/s;
%vscale=(yl./O_maxval(dy))/(0.2*velsc); %1 grid cell is 0.5 m/s;
for nps=n1:ni:n2; n=nps-n1+1; 
pelout(:,:)=elout(:,:,nps); pelout(fsm==0)=NaN;
plfil=[progname 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(nps),'%10.3f\n')]; %pltit=plfilt;
pltit = ['t=' num2str((n-1)*outtim*dayi) ' d,' ...
    ' Max=' num2str(O_maxval(pelout),'%8.1f\n') ', Min=' num2str(O_minval(pelout),'%8.1f\n')]; 
%--[plotfile] = O_xycontour04(xpl2,ypl2,pelout, ...
%--	    xrng,yrng,zrng,xi,yi,wrk,n,pltit,animfile,'',fig,wh);
            %No pvelxy (vectors), Yes fig (animation), Yes wh (surf)
%velsc=2.*max(O_maxval(sqrt(uaout.^2+vaout.^2))); %=0.5;
%uscale=(xl./O_maxval(dx))/velsc; %1 grid cell is 0.5 m/s;
%vscale=(yl./O_maxval(dy))/velsc; %1 grid cell is 0.5 m/s;
pvelxy(:,:)=uscale*uaout(:,:,nps)+vscale*1i*vaout(:,:,nps); %the w-component = 0;
pvelxy(1,1)=uscale+1i*vscale;
[plotfile] = O_xycontour04(xpl2,ypl2,pelout, ...
	    xrng,yrng,zrng,xi,yi,wrk,n,pltit,animfile,pvelxy,fig,wh);
            %Yes pvelxy (vectors), Yes fig (animation), Yes wh (surf)
end;
%
%-------------------------------------------------------------------------
%How to make elevation vs. x plots at different times:
%{
nfig=nfig+1; figure(nfig);
hpl(:,:)=h(:,:);
n2=nm; n1=n2-4; ni=1; %plot last 12hours, every ni*3hours
yrng(1)=-zscale*0.5; yrng(2)=zscale*0.5; hpl(1,:) = -yrng(2); xplmin=0.9*xpl(2);
for nps=n1:ni:n2; 
pelout(:,:)=elout(:,:,round(nps)); pelout(fsm==0)=NaN;
plot(xpl(:),pelout(:,2),'b','LineWidth',2);
hold on; text(0.1*(xplmin+xpl(1)),pelout(2,2),num2str(tpl(nps)));
%hold on; plot(xpl(2:im),-hpl(2:im,2),'k','LineWidth',2);
hold on; plot([xplmin xpl(2:im)],-hpl(1:im,2),'k','LineWidth',2);
hold on; plot([xpl(1) xplmin],[-hpl(1,2) -hpl(1,2)],'k','LineWidth',2);
msl=zeros(im); %(0.0:1.0:im-1)*0.0;
hold on; plot(xpl(1:im),msl(1:im),'k--','LineWidth',1);
axis([xrng(1) xrng(2) yrng(1) yrng(2)]);
end
%pltit = ['EL Last Period=' num2str(per/86400,'%10.1f\n') ' day'];
pltit = [progname 'EL_Last_Period'];
title(pltit, 'fontsize',16);
pltit=[plotdir pltit];
saveas(gcf,[pltit '.png']); %print('-dpng', '-loose', [pltit '.png']);
%saveas(gcf,[pltit '.ps']); %print('-dpsc', '-loose', [pltit '.ps']);
%}
%-------------------------------------------------------------------------
%Suppress plots after 3D surf animation (i.e. xz & yz section plots):BEG
%{
%How to plot xz-section contours of ua w/vectors, @j=jm/2 & t=tpl(nm/2)day
%Plot to follow time-dependent free surface using sigma-coordinate
nfig=nfig+1; figure(nfig);
%Define (xpl2,zpl2)=(x,z)-coordinate of the u-point at time n=tpl(nm/2):
j=round(jm/2); nfix=nm; %round(nm/2);
j=2;
kb=5; sige=-(0.0:1.0:kb-1)/(kb-1); %Use sigma-coord == POM's "z"
sigu(1:kb-1)=0.5*(sige(1:kb-1)+sige(2:kb)); %sigma at u-point == POM "zz"
sigu=sige; %define at edge of sigma cell to fill up color -> nicer plot
xpl2=zeros(im,kb); zpl2=zeros(im,kb); pua2=zeros(im,kb); fsm2=zeros(im,kb);
for k=1:kb; fsm2(:,k)=fsm(:,j); end; for k=1:kb; xpl2(:,k)=xpl(:); end;
%Time-dependent zpl2 & pua2:
for n=nfix:nfix; %1:nm;
for k=1:kb; zpl2(:,k)=(h(:,j)+elout(:,j,n))*sigu(k)+elout(:,j,n); end;
%for k=1:kb-1; zpl2u(:,k)=(h(:,j)+elout(:,j,n))*sigu(k)+elout(:,j,n); end;
for k=1:kb; pua2(:,k)=uaout(:,j,n); end; pua2(fsm2==0)=NaN;
%Contour pua2 @(xpl2,zpl2): contourf(xpl2,zpl2,pua2);
xrng=NaN(1,2);yrng=xrng;zrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pua2));nps=xi;pltit='';plfil='';
xrng(1)=O_minval(xpl2);xrng(2)=O_maxval(xpl2);%O_maxval(xpl);%nps=-1;%ps
%yrng(1)=round(O_minval(zpl2)*1.2); yrng(2)=zscale;  %Tide1D:
yrng(1)=-zscale; yrng(2)=zscale/3.;                     %Tsunami1D:
yrng(1)=-2.0; yrng(2)=2.0;

zrng(1)=-O_maxval(abs(uaout));zrng(2)=-zrng(1);% +O_maxval(abs(uaout));
wrk=fsm2; %--- array to put gray shading at land points
%---pvel=pua2+1i*0.0; %the imaginary part "w"-component = 0;
pvel=zeros(im,kb)+1i*0.0;
plfil=[progname 'xz-' 'xl' num2str(xl*1.e-3) 'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n),'%10.3f\n')]; pltit = plfil;
plfil=[plotdir plfil];
[plotfile]=O_xycontour04(xpl2,zpl2,pua2,xrng,yrng,zrng,xi,yi,wrk,nps, ...
    pltit,plfil);%,pvel);%plfil);
            %No pvel: vectors, No fig: animation, No wh: surf
%    pltit,plfil,pvel);%plfil);
            %Yes pvel: vectors, No fig: animation, No wh: surf
%   pltit,'');%plfil); %this works, 'pvel' does not exist in *contour04
end;
%
%-------------------------------------------------------------------------
%Make xz-animation following the free surface & superimpose u-contours
%
%How to make animation xz-section contours of ua w/vectors, ...
%   The animation *gif can be split into individual gif files by uploading 
%   it into this online program: https://ezgif.com/split
%
nfig=nfig+1; figure(nfig);
j=round(jm/2); %nfix=round(nm/2);
j=2;
kb=5; sige=-(0.0:1.0:kb-1)/(kb-1); %Use sigma-coord == POM's "z"
sigu(1:kb-1)=0.5*(sige(1:kb-1)+sige(2:kb)); %sigma at u-point == POM "zz"
sigu=sige; %define at edge of sigma cell to fill up color -> nicer plot
xpl2=zeros(im,kb); zpl2=zeros(im,kb); pua2=zeros(im,kb); fsm2=zeros(im,kb);
for k=1:kb; fsm2(:,k)=fsm(:,j); end; for k=1:kb; xpl2(:,k)=xpl(:); end;
%{%
pelout=zeros(im,jm); %Tsunami1D: uncomment this to correct Tide1D version
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl2);xrng(2)=O_maxval(xpl2);
%--yrng(1)=O_minval(ypl);yrng(2)=O_maxval(ypl);
%--zrng(1)=-zscale;zrng(2)=+zscale;
zrng(1)=-O_maxval(abs(uaout));
zrng(2)=-zrng(1);% +O_maxval(abs(uaout));
wrk=fsm2;
%3dsurface:wrk(1,1)=9; %---Set (1,1) to 9 for surface plot animation
%3dsurface:wh=h; wh(fsm==0)=NaN;  wh=-0.75*zscale*wh./O_maxval(wh);
%pause;
%S=['Channel2d_02, wrk(1,1) = ' num2str(wrk(1,1)) '******************']; disp(S)
%pause;
n1=65;ni= 1;n2=nm;%70; 
n1=1;ni= 1;n2=nm;%70; 
fig = figure(nfig);
axis tight manual % this ensures that getframe() returns a consistent size
plfil=[progname '-uxz-' ...
                   'xl' num2str(xl*1.e-3) ...
                   'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n1),'%10.3f\n') 'to' num2str(tpl(n2),'%10.3f\n')];
%                              
plfil=[plotdir plfil];
animfile = [plfil '.gif'];
for nps=n1:ni:n2; n=nps-n1+1; 
for k=1:kb; zpl2(:,k)=(h(:,j)+elout(:,j,nps))*sigu(k)+elout(:,j,nps); end;
%for k=1:kb-1; zpl2u(:,k)=(h(:,j)+elout(:,j,nps))*sigu(k)+elout(:,j,n); end;
for k=1:kb; pua2(:,k)=uaout(:,j,nps); end; pua2(fsm2==0)=NaN;
%Tsunami1D: uncomment next line to correct Tide1D version
pelout(:,:)=elout(:,:,nps); pelout(fsm==0)=NaN;
%{
plfilt=[progname '-uxz-' ...
                    'xl' num2str(xl*1.e-3) ...
                    'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(nps),'%10.3f\n')]; %pltit=plfilt;
%}
pltit = ['u, t=' num2str((n-1)*outtim*dayi) ' d,' ...
    ' Max=' num2str(O_maxval(pelout),'%8.1f\n') ', Min=' num2str(O_minval(pelout),'%8.1f\n')]; 
yrng(1)=round(O_minval(zpl2)*1.2); yrng(1)=-zscale; %Tsunami1D:
yrng(2)=zscale/3.;
yrng(1)=-2.0; yrng(2)=2.0;


%--zrng(1)=-O_maxval(abs(uaout));zrng(2)=-zrng(1);% +O_maxval(abs(uaout));
%---pvel=pua2+1i*0.0; %the w-component = 0;
pvel=zeros(im,kb)+1i*0.0;
[plotfile] = O_xycontour04(xpl2,zpl2,pua2, ...
	    xrng,yrng,zrng,xi,yi,wrk,n,pltit,animfile,pvel,fig);%,wh);
            %Yes pvel: vectors, Yes fig: animation, No wh: surf
end;
%
%-------------------------------------------------------------------------
%Make xz-animation following the free surface & superimpose v-contours
%
nfig=nfig+1; figure(nfig);
j=round(jm/2); %nfix=round(nm/2);
j=2;
kb=5; sige=-(0.0:1.0:kb-1)/(kb-1); %Use sigma-coord == POM's "z"
sigu=sige; %define at edge of sigma cell to fill up color -> nicer plot
xpl2=zeros(im,kb); zpl2=zeros(im,kb); pua2=zeros(im,kb); fsm2=zeros(im,kb);
for k=1:kb; fsm2(:,k)=fsm(:,j); end; for k=1:kb; xpl2(:,k)=xpl(:); end;
%{%
pelout=zeros(im,jm); %Tsunami1D: uncomment this to correct Tide1D version
xrng=NaN(1,2);yrng=xrng;xi=NaN;yi=xi;wrk=NaN(size(pelout));nps=xi;pltit='';plfil='';
zrng=xrng;
xrng(1)=O_minval(xpl2);xrng(2)=O_maxval(xpl2);
%--yrng(1)=O_minval(ypl);yrng(2)=O_maxval(ypl);
zrng(1)=-0.2; %-O_maxval(abs(vaout));  %<---FIX this for Eddy2D:
zrng(2)=-zrng(1);% +O_maxval(abs(uaout));
%--velsca=max(abs(O_minval(vaout)),O_maxval(vaout));
%--zrng(1)=-1.1*velsca;zrng(2)=+1.1*velsca;
wrk=fsm2;
%3dsurface:wrk(1,1)=9; %---Set (1,1) to 9 for surface plot animation
%3dsurface:wh=h; wh(fsm==0)=NaN;  wh=-0.75*zscale*wh./O_maxval(wh);
%pause;
%S=['Channel2d_02, wrk(1,1) = ' num2str(wrk(1,1)) '******************']; disp(S)
%pause;
n1=1;ni= 1;n2=nm;%70; 
fig = figure(nfig);
axis tight manual % this ensures that getframe() returns a consistent size
plfil=[progname '-vxz-' ...
                   'xl' num2str(xl*1.e-3) ...
                   'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(n1),'%10.3f\n') 'to' num2str(tpl(n2),'%10.3f\n')];
%                              
plfil=[plotdir plfil];
animfile = [plfil '.gif'];
for nps=n1:ni:n2; n=nps-n1+1; 
for k=1:kb; zpl2(:,k)=(h(:,j)+elout(:,j,nps))*sigu(k)+elout(:,j,nps); end;
for k=1:kb; pua2(:,k)=vaout(:,j,nps); end; pua2(fsm2==0)=NaN;
pelout(:,:)=elout(:,:,nps); pelout(fsm==0)=NaN;
%{
plfilt=[progname '-vxz-' ...
                    'xl' num2str(xl*1.e-3) ...
                    'yl' num2str(yl*1.e-3) ...
'cor' num2str(cor0,'%10.1e\n') 'day' num2str(tpl(nps),'%10.3f\n')]; %pltit=plfilt;
%}
pltit = ['v, t=' num2str((n-1)*outtim*dayi) ' d,' ...
    ' Max=' num2str(O_maxval(pelout),'%8.1f\n') ', Min=' num2str(O_minval(pelout),'%8.1f\n')]; 
yrng(1)=round(O_minval(zpl2)*1.2); yrng(1)=-zscale; %Tsunami1D:
yrng(2)=zscale/3.;
yrng(1)=-2.0; yrng(2)=2.0;


%---pvel=pua2+1i*0.0; %the w-component = 0;
pvel=zeros(im,kb)+1i*0.0;
[plotfile] = O_xycontour04(xpl2,zpl2,pua2, ...
	    xrng,yrng,zrng,xi,yi,wrk,n,pltit,animfile,pvel,fig);%,wh);
            %Yes pvel: vectors, Yes fig: animation, No wh: surf
end;
%
%}
%Suppress plots after 3D surf animation (i.e. xz & yz section plots):END
%-------------------------------------------------------------------------
fprintf(fileid,['--------Printouts from O_outputext---------\n']);
fprintf(fileid,['cor0 =' num2str(cor0) '\n']);
fprintf(fileid,['MaxEl=' num2str(O_maxval(elout)) ', MinEl=' num2str(O_minval(elout)) '\n']);
fprintf(fileid,['MaxUa=' num2str(O_maxval(uaout)) ', MinUa=' num2str(O_minval(uaout)) '\n']);
fprintf(fileid,['MaxVa=' num2str(O_maxval(vaout)) ', MinVa=' num2str(O_minval(vaout)) '\n']);
%
return
%


