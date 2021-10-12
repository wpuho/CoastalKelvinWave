function [ael,aua]=O_anatidechannel(elout,xpl,tpl,dx,h)
%
%% Analytical solution for tide in 1D channel
%
elevamp=0.2; ome=1.454441043328608e-04; %amp & freq from O_bc.m
g=10.;                                  %gravity from O_modparam.m
%
[im,jm,nm]=size(elout); ic=round(im/2); jc=round(jm/2);
ael=zeros(im,nm); aua=zeros(im,nm);
%
%analytic x-point for elev: (Note: input xpl in km & tpl in days)
ax=xpl*1.e3; %first define ax to be same as xpl, both at el-point
ax(2:im)=ax(2:im)-0.5*dx(1,jc); %adjust origin to left-wall u@i=2 ax(1)=0
al = sum(dx(2:im))-0.5*dx(im);   %channel length
%Adjust al so amplitude near x=0 matches numr soln - necessary becasue 
%   there is ambiguity of the channel length on staggered grid
al = al*1.0194; % Slightly longer (=96.8km instead of 95km) for dx=10km
%al = al*1.012; % if dx=5km
%
%phase speed, period, wavelength & wavenumber:
ac = sqrt(g*h(ic,jc)); ap = 2.*pi/ome; awvl = ap*ac; ak = 2.*pi/awvl; 
%
%Solutions, equations (1.2.44a & b):
aelamp = elevamp/cos(ak*al); auaamp = ac*aelamp/h(ic,jc); omed=ome*86400.
for n=1:nm; for i=1:im;
        ael(i,n) = aelamp*cos(ak*ax(i))*sin(omed*tpl(n));
        aua(i,n) = auaamp*sin(ak*ax(i))*(1-cos(omed*tpl(n)));
    end; end;
return
%

