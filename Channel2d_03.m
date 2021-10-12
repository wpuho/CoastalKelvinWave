%% ------- Building Computer Models Of Ocean Atmos -------
%               Two-Dimensional Channel Codes
%
%------------ Report bugs and/or improvements to: ------------
%               Leo Oey lyo@alumni.princeton.edu 
%               (20200506; 20190518)
%
%---------------User's functions required:-------------------- 
% 1. Subfolder "./subroutines/model/"       contains model functions 
% 2. Subfolder "./subroutines/graphics/"    contains plot  functions 
% 3. Subfolder "./subroutines/utilities/"   contains other functions 
%
%------------------------------------------------------------- 
% Model sea-level elevation (el) and depth-averaged current (ua,va) in a 
%   2d channel x=(0,xl), y=(0,yl) of dimension xl*yl
%
%   In the Lecture Note: 
%       http://oeylectures.pbworks.com/w/file/136377039/
%       BuildingComputerModelsOfOceanAtmos.pdf
%   el = Greek-eta; ua = u; and va = v.
% Governing equations (equation 1.3.5 + Coriolis + friction):
%   d(el)/dt = -d(d*ua)/dx -d(d*va)/dy;
%   D(d*ua)/Dt = -(gd)d(el-ela)/dx + f*va + WindX - FrictionX + ViscousX
%   D(d*va)/Dt = -(gd)d(el-ela)/dy - f*ua + WindY - FrictionY + ViscousY
%   where d/dt, d/dxd/dy are partial, but D/Dt = Material Derivative,
%   el=elevation, (ua,va)=depth-averaged (x,y)-velocity, 
%   ela=-pas/(rhos*g) = the surface elevation of an inverse barometer
%                       pas & rhos = surface atmos pressure & water density
%   d=h+el & h=channel depth
%
%------------------------------------------------------------- 
% Numerical Method:
%   Leap-Frog time stepping & central-space difference on staggered grid
%
%--------------------Model Plot Outputs----------------------- 
% Are in subfolder "./plots"
%
%------------------------------------------------------------- 
%% Program begins: 
% Clear/close everything to start fresh:
clear all; close all; clc; fclose all;
% Tunr off annoying warnings:
warning off;
% Set program dir: directory where this program is
progdir='./';
% Set leading filename of output plot & text files
progname=['ForcedKWCoast07'];     %<<< 
%-------------------------------------------|
txtdir=['./txtout/']; mkdir(txtdir);        %
fileid=fopen([txtdir progname '.txt'],'w'); %
%-------------------------------------------|
matdir=['./matout/']; mkdir(matdir);        %
matdatfile=[ matdir progname ];             %
%-------------------------------------------|
% Set paths:
addpath(strcat(progdir,'subroutines/model/'));
addpath(strcat(progdir,'subroutines/graphics/'));
addpath(strcat(progdir,'subroutines/utilities/'));
addpath(strcat(progdir,'subroutines/anasoln/'));
etopo_path = [progdir,'ETOPO2v2c_f4_netCDF/ETOPO2v2c_f4.nc'];
%
%% Do model calculation only if matdatfile does NOT exist:
%if exist([matdatfile '.mat']) ~= 2;
if exist([matdatfile 'test.mat']) ~= 2; %Override: always do model calc.
%
%% Specify model params & constants (MKS units, unless otherwise stated):
%tottim = model total time for this run
%outtim = model output time interval
[smoth,smol,dtfac,fnln,small,g,pi,tottim,outtim,                ...
                                nperx,npery,smagc,nadv]=O_modparam;
%
%% Specify model domain (xl,yl), grid, topography, Coriolis & time step:
[dx,dy,are,aru,arv,h,fsm,dum,dvm,cor,dt,xl,yl,zt,lon,lat]=O_grid(small,g,dtfac, ...
                                                        nperx,npery,etopo_path);
[im,jm]=size(dx); imm1=im-1; jmm1=jm-1; cor0=cor(round(im/2),round(jm/2));
%Adjust "dt" so that outtim below = interger multiple of dt, = nout
dt   = outtim/ceil(outtim/dt);          %--- <= dt based on CFL
nmod   = tottim/dt;                     %---Total# model time steps
%
%% Output Parameters & Arrays:
nout   = outtim/dt;                     %---Output every nout time steps
nm     = tottim/(nout*dt)+1;            %---Total# outputs
ncount = 1;
%Output arrays:
elout(1:im,1:jm,1:nm) = 0.0; %---output arrays; (:,:,1) for t=0 fields
uaout(1:im,1:jm,1:nm) = 0.0; %---1d OK; 2d/3d maybe too memory-intensive
vaout(1:im,1:jm,1:nm) = 0.0; %---1d OK; 2d/3d maybe too memory-intensive
%
%% Specify Model Arrays & Initial Conditions:
%Integration arrays, (b,n,f) = (backward,now,future) time levels:
[elf,uaf,vaf,elb,eln,uab,uan,vab,van,d,bndwe,bndsn]=O_initial(h);
%
%% Solve, Time-Stepping:
elout(:,:,ncount) = eln(:,:).*fsm(:,:); 
uaout(:,:,ncount) = uan(:,:).*dum(:,:);
vaout(:,:,ncount) = van(:,:).*dvm(:,:);
for n = 1:nmod ;
%
%[elf,uaf,vaf,elb,eln,uab,uan,vab,van,d]=O_solveext      ...
[elf,uaf,vaf,elb,eln,uab,uan,vab,van,d,ela,uwnd,vwnd]=O_solveext      ...
            (elb,eln,uab,uan,vab,van,d,                               ...
             h,dx,dy,dt,are,aru,arv,cor,fsm,dum,dvm,                  ...
             n,g,fnln,smoth,smol,nperx,npery,smagc,                   ...
             nadv,bndwe,bndsn);
if (mod(n,nout) == 0); ncount=ncount+1; 
elout(:,:,ncount) = eln(:,:); 
uaout(:,:,ncount) = uan(:,:);
vaout(:,:,ncount) = van(:,:);
elaout(:,:,ncount) = -(1.e4)*ela(:,:); %sea-level pressure (hPa)
uwndout(:,:,ncount) = uwnd(:,:);       %x-wind @el-point (m/s)
vwndout(:,:,ncount) = vwnd(:,:);       %y-wind @el-point (m/s)
end;
%
end;
%{
%% Plots & Prints:
ana = strcat(progdir,'subroutines/anasoln/O_anatsunami.m'); %Tsunami1D:
[nfig]=O_outputext(progname,elout,uaout,vaout,outtim,cor0,dx,dy,fsm,h,  ...
                    ana,fileid);
%
%% Save matfile:
save( matdatfile,'elout', 'uaout',  'vaout',                    ...
                 'elaout','uwndout','vwndout','outtim','cor0',  ...
                 'dx','dy','fsm','h','ana','fileid');
[matdatfile '.mat saved']
%
%% End Program:
fprintf(fileid,['--------Printouts from MAIN---------\n']);
fprintf(fileid,['smoth =' num2str(smoth) '\n']);
fprintf(fileid,['smol =' num2str(smol) '\n']);
fprintf(fileid,['dtfac =' num2str(dtfac) '\n']);
fprintf(fileid,['fnln =' num2str(fnln) '\n']);
fprintf(fileid,['small =' num2str(small) '\n']);
fprintf(fileid,['g =' num2str(g) '\n']);
fprintf(fileid,['tottim =' num2str(tottim) '\n']);
fprintf(fileid,['outtim =' num2str(outtim) '\n']);
fprintf(fileid,['nperx =' num2str(nperx) '\n']);
fprintf(fileid,['npery =' num2str(npery) '\n']);
fprintf(fileid,['smagc =' num2str(smagc) '\n']);
fprintf(fileid,['nadv =' num2str(nadv) '\n']);
fprintf(fileid,['-------------------------------------------\n']);
fprintf(fileid,['Model' progname ' ends @day=' num2str(tottim/86400.) ';\n']);
fprintf(fileid,['Last fig#=' num2str(nfig) '\n']);
fprintf(fileid,['-------------------------------------------\n']);
fclose(fileid);
%
else; 
%% matdatfile does exist, so load & plot it only (i.e. w/o model calc):
load(matdatfile);
['Finished loading ' matdatfile '.mat; now plot:']
ana = strcat(progdir,'subroutines/anasoln/O_anatsunami.m'); %Tsunami1D:
[nfig]=O_outputext(progname,elout,uaout,vaout,outtim,cor0,dx,dy,fsm,h,  ...
                    ana,fileid);
%}
HO_plot2(elout,uaout,vaout,fsm,outtim,lon,lat)
end;
%%
return
%