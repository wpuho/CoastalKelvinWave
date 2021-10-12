function [ff]=O_advt2d(fb,fbclim,sw,dx,dy,u,v,art,aru,arv,  ...
                        fsm,dum,dvm,aam,wriv,dt2)
%%
%   Solves d(D)/dt + d(uD)/dx + d(vD)/dy = Diffusion_of_(D)
%   using the Smolarkiewicz scheme to minimize grid oscillations
%   This MATLAB funciton is adapted after Oey [2005] West&Dry code:
%   c:../atop/mpiPOM/wrk_lyo_mpiPOM_pom_wetdry/proc/ ...
%           wadproc01_mode2_wper=5_drysea_cor=c_nowind_advt2d.f
%           See notes in wadproc01*...............advt2d-COPY.f
%   Note: nitera(e.g. =3) is directly specified below;
%           nitera=1 is upwind; nitera >=3 all give v.similar results;
%         sw (=smol in O_modparam.m) should be ~=1; 
%           0 < sw < 1 gives smoother solutions when nitera > 1.
%
%%
[im,jm]=size(fb); imm1=im-1; jmm1=jm-1;
ff=zeros(im,jm);
%Local work arrays:
xflux=zeros(im,jm);yflux=zeros(im,jm);fbmem=zeros(im,jm);
xmassflux=zeros(im,jm);ymassflux=zeros(im,jm);
%Horizontal mass fluxes:
xmassflux(2:im,:)=0.5*(dy(1:imm1,:)+dy(2:im,:)).*u(2:im,:);
ymassflux(:,2:jm)=0.5*(dx(:,1:jmm1)+dx(:,2:jm)).*v(:,2:jm);
fbmem(:,:)=fb(:,:);
%Start Smolarkiewicz scheme:
for itera=1:3; %3;          %itera=1:nitera, where nitera=3;
%Upwind advection scheme:
xflux(2:im,2:jm)=0.5e0*                                                 ...
 ((xmassflux(2:im,2:jm)+abs(xmassflux(2:im,2:jm))).*fbmem(1:imm1,2:jm)+ ...
  (xmassflux(2:im,2:jm)-abs(xmassflux(2:im,2:jm))).*fbmem(2:im,2:jm));
yflux(2:im,2:jm)=0.5e0*                                                 ...
 ((ymassflux(2:im,2:jm)+abs(ymassflux(2:im,2:jm))).*fbmem(2:im,1:jmm1)+ ...
  (ymassflux(2:im,2:jm)-abs(ymassflux(2:im,2:jm))).*fbmem(2:im,2:jm));
%Add net advective fluxes and step forward in time:
ff(1:imm1,1:jmm1)=xflux(2:im,1:jmm1)-xflux(1:imm1,1:jmm1)+ ...
                  yflux(1:imm1,2:jm)-yflux(1:imm1,1:jmm1);
ff(:,:)=fsm(:,:).*(fbmem(:,:).*art(:,:)-dt2*ff(:,:))./(art(:,:));
%
%Calculate antidiffusion velocity - given ff, update xmassflux & ymassflux:
O_smoladif;
%
fbmem(:,:)=ff(:,:);
end;
%End of Smolarkiewicz scheme loop
%
%Add horizontal diffusive fluxes:
fbmem(:,:)=fb(:,:)-fbclim(:,:); %For fb=D, remove topography h=fbclim, so:
                                % diffusion based on Laplacian(elb) only;
                                % diffusion then "nudges" ff(=D) back to h
xmassflux=zeros(im,jm);ymassflux=zeros(im,jm);
xmassflux(2:im,2:jm)=0.5e0*(aam(2:im,2:jm)+aam(1:imm1,2:jm));
ymassflux(2:im,2:jm)=0.5e0*(aam(2:im,2:jm)+aam(2:im,1:jmm1));
xflux(2:im,2:jm)=-xmassflux(2:im,2:jm)                                  ...
     .*(fbmem(2:im,2:jm)-fbmem(1:imm1,2:jm)).*dum(2:im,2:jm)            ...
     .*(dy(2:im,2:jm)+dy(1:imm1,2:jm))./(dx(2:im,2:jm)+dx(1:imm1,2:jm));
yflux(2:im,2:jm)=-ymassflux(2:im,2:jm)                                  ...
     .*(fbmem(2:im,2:jm)-fbmem(2:im,1:jmm1)).*dvm(2:im,2:jm)            ...
     .*(dx(2:im,2:jm)+dx(2:im,1:jmm1))./(dy(2:im,2:jm)+dy(2:im,1:jmm1));
%Add net horizontal fluxes and step forward in time:
ff(1:imm1,1:jmm1)=ff(1:imm1,1:jmm1)-dt2*(                  ...
          ( xflux(2:im,1:jmm1)-xflux(1:imm1,1:jmm1)         ...
           +yflux(1:imm1,2:jm)-yflux(1:imm1,1:jmm1) )       ...
           ./art(1:imm1,1:jmm1)  +  wriv(1:imm1,1:jmm1) );
%
%%
function O_smoladif
%
%   Inputs:  ff, xmassflux & ymassflux
%   Outputs: xmassflux & ymassflux
%
%Calculates the antidiffusive velocity used to reduce the numerical 
%   diffusion associated with the upstream differencing scheme
%
%Local variables:
value_min=1.e-9; epsilon=1.0e-14;
%Recalculate mass fluxes with antidiffusion velocity:
for j=2:jmm1; for i=2:im;
    if (ff(i,j)<value_min | ff(i-1,j)<value_min); xmassflux(i,j)=0.e0;
    else;
    udx=abs(xmassflux(i,j));
    u2dt=dt2*xmassflux(i,j)*xmassflux(i,j)/(aru(i,j));
    mol=(ff(i,j)-ff(i-1,j))/(ff(i-1,j)+ff(i,j)+epsilon);
    xmassflux(i,j)=(udx-u2dt)*mol*sw;
    abs_1=abs(udx); abs_2=abs(u2dt);
    if(abs_1<abs_2);xmassflux(i,j)=0.e0;end;
    end;
    end;      end;
%
for j=2:jm;   for i=2,imm1;
    if(ff(i,j)<value_min | ff(i,j-1)<value_min); ymassflux(i,j)=0.e0;
    else;
    vdy=abs(ymassflux(i,j));
    v2dt=dt2*ymassflux(i,j)*ymassflux(i,j)/(arv(i,j));
    mol=(ff(i,j)-ff(i,j-1))/(ff(i,j-1)+ff(i,j)+epsilon);
    ymassflux(i,j)=(vdy-v2dt)*mol*sw;
    abs_1=abs(vdy); abs_2=abs(v2dt);
    if(abs_1<abs_2);ymassflux(i,j)=0.e0;end;
    end;
    end;      end;
end;
%
%%
end

