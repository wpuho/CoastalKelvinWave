function [aam]=O_horvis(smagc,dx,dy,u,v,fsm)
%
%Calculate horizontal viscosity/diffusivity coefficient aam (m^2/s)
%%
[im,jm]=size(u); imm1=im-1; jmm1=jm-1; imm2=im-2; jmm2=jm-2;
aam=zeros(im,jm);
%
aam(2:imm1,2:jmm1)=smagc*dx(2:imm1,2:jmm1).*dy(2:imm1,2:jmm1)         ...
  .*sqrt( ((u(3:im,2:jmm1)-u(2:imm1,2:jmm1))./dx(2:imm1,2:jmm1)).^2   ...
  +((v(2:imm1,3:jm)-v(2:imm1,2:jmm1))./dy(2:imm1,2:jmm1)).^2          ...
  +.5*( .25*(u(2:imm1,3:jm)  +u(3:im,3:jm)                            ...
            -u(2:imm1,1:jmm2)-u(3:im,1:jmm2))./dy(2:imm1,2:jmm1)      ...
       +.25*(v(3:im,2:jmm1)  +v(3:im,3:jm)                            ...
            -v(1:imm2,2:jmm1)-v(1:imm2,3:jm))./dx(2:imm1,2:jmm1)).^2);
aam(1,:)=aam(2,:); aam(im,:)=aam(imm1,:);
aam(:,1)=aam(:,2); aam(:,jm)=aam(:,jmm1);
aam(:,:)=aam(:,:).*fsm(:,:);
%
return;



