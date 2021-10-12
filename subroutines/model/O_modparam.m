function [smoth,smol,dtfac,fnln,small,g,pi,tottim,outtim,       ...
                             nperx,npery,smagc,nadv]=O_modparam
%% Model Params & Constants (MKS units, unless otherwise stated):
%
smoth=0.1;                       %temporal filter constant Leap-Frog Scheme
smol=0.0;                        %=1.0 for Smolarkiewicz, otherwise 0.
nperx=0; npery=0;                %=1 if x- or y-periodic, otherwise =0.
smagc=0.0;                       %Smagorinsky cons; <0 for aam=-smagc m^2/s
nadv=1;                          %n-time interval for calc. advect+diffuse 
dtfac = 1.0;                     %=1 use dtcfl, <1 for more stability, e.g.
                                 % if smol.ne.0, dt<dtcfl is often needed;
                                 % but w/alpha=0.225 (O_solveext.m) then 
                                 % dtcfl=1 is stable for smol.ne.0.
linear=1; fnln=(1.-real(linear));%linear=1 linear, =0 nonlin  %19%22
small=1.e-6;                     %program constant
%g=0.133959190672; %=10.;         %gravity, chosen=c^2/H, H=10m in O_grid.m,
                                 % The "c" chosen so KW travels 500km in 5d
                                 % c*5*86400=500km, c=1.15740740741,
                                 % g = 1.33959190672/10 = 0.133959190672
                                 
g=10.;                           %09:Change from above reduced gravity, 
                                 % Also change tottim & outtim below, AND:
                                 % O_surfforce: uwnd=4.64717194000958*..
                                 % O_grid: cor0 = 10./100000.; a=100km
pi=acos(-1.);
tottim = 1.0*10.000*86400.;       %---total model integration time
outtim = 1.0* 0.125*86400.;       %---output interval
%tottim = 1.0*4.000*86400.;        %09
%outtim = 1.0* 0.025*86400.;       %09
%
return;
%