function [elf,uaf,vaf]=O_bc(nflag,elf,uaf,vaf,fsm,dum,dvm,n,dt,g,h, ...
                            eln,nperx,npery,bndwe,bndsn)
%
[im,jm]=size(elf); imm1=im-1; jmm1=jm-1; imm2=im-2; jmm2=jm-2;
elbw(1,:)=bndwe(1,:); elbe(1,:)=bndwe(2,:); 
uabw(1,:)=bndwe(3,:); uabe(1,:)=bndwe(4,:);
elbs(:,1)=bndsn(:,1); elbn(:,1)=bndsn(:,2); 
vabs(:,1)=bndsn(:,3); vabn(:,1)=bndsn(:,4);
%--uabw=zeros(1,jm); elbw=zeros(1,jm)+0.5; uabe=zeros(1,jm); elbe=zeros(1,jm)-0.5;
%--vabs=zeros(im,1); elbs=zeros(im,1); vabn=zeros(im,1); elbn=zeros(im,1);
%
%Tide1D: elev amplitude & tidal frequency at im:
%--elevamp=0.2; ome=1.454441043328608e-04; %per=12.*3600.; ome=2.*pi/per;
%
if (nflag == 1);        %Elevation boundary condition
if (nperx == 1);        %x-periodic
elf(im,:)=elf(3,:); elf(1,:)=elf(imm2,:); elf(2,:)=elf(imm1,:);
else;                   %regular x-BCs
%--elf(im,:)=elevamp*sin(ome*n*dt);  %Tide1D: 
elf(1,:)=elf(2,:); elf(im,:)=elf(imm1,:);
end;
if (npery == 1);        %y-periodic
elf(:,jm)=elf(:,3); elf(:,1)=elf(:,jmm2); elf(:,2)=elf(:,jmm1);
else;                   %regular y-BCs
elf(:,1)=elf(:,2); elf(:,jm)=elf(:,jmm1);
end;
%
elf(:,:)=elf(:,:).*fsm(:,:);
end;
%
if (nflag == 2); %uaf & vaf boundary condition
if (nperx == 1);        %x-periodic
uaf(im,:)=uaf(3,:); uaf(1,:)=uaf(imm2,:); uaf(2,:)=uaf(imm1,:);
vaf(im,:)=vaf(3,:); vaf(1,:)=vaf(imm2,:); vaf(2,:)=vaf(imm1,:);
else;                   %regular x-BCs
%west & east radiation if open:
uaf(2,:)=uabw(1,:)-sqrt(g./h(2,:)).*(eln(2,:)-elbw(1,:)); 
%--vaf(1,:)=0.0; uaf(1,:)=uaf(2,:); %---Just fill this west-end value
vaf(1,:)=vaf(2,:); uaf(1,:)=uaf(2,:); 
uaf(im,:)=uabe(1,:)+sqrt(g./h(imm1,:)).*(eln(imm1,:)-elbe(1,:)); 
vaf(im,:)=vaf(imm1,:); %0.0;        %---Just fill this east-end value
end;
if (npery == 1);        %y-periodic
uaf(:,jm)=uaf(:,3); uaf(:,1)=uaf(:,jmm2); uaf(:,2)=uaf(:,jmm1);
vaf(:,jm)=vaf(:,3); vaf(:,1)=vaf(:,jmm2); vaf(:,2)=vaf(:,jmm1);
else;                   %regular y-BCs
%south & north radiation if open:
vaf(:,2)=vabs(:,1)-sqrt(g./h(:,2)).*(eln(:,2)-elbs(:,1)); 
uaf(:,1)=0.0; vaf(:,1)=vaf(:,2); %---Just fill this south-end value
vaf(:,jm)=vabn(:,1)+sqrt(g./h(:,jmm1)).*(eln(:,jmm1)-elbn(:,1)); 
uaf(:,jm)=0.0;                   %---Just fill this north-end value
end;
%
uaf(:,:)=uaf(:,:).*dum(:,:);     %---Note this is u-wall boundary condition
vaf(:,:)=vaf(:,:).*dvm(:,:);     %---Note this is v-wall boundary condition
end;
%
return;


