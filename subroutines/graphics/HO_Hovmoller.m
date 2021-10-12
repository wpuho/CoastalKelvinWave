function HO_Hovmoller(lon,lat,elout,uaout,vaout,outim)

tper = size(elout,3); ylen = size(elout,2); xlen = size(elout,1);
lonL = 26; lonR = 26; latB = 31; latT = 41;
%lonL = 151; lonR = 176; latB = 36;
xl = lonR - lonL + 1;
yl = latT - latB + 1;

Havmoller1 = zeros(tper,yl); Havmoller2 = zeros(tper,yl);
Havmoller = zeros(tper,yl);
time_interval = outim;
ytick = string(tper);
%{
for t = 1:tper
    Havmoller(t,:) = elout(lonL:lonR,latB,t);
    if (t==1)
        ytick(t) = [' 0 hr'];
    else
        ytick(t) = [num2str(outim/3600),' hr'];
        outim = outim + time_interval;
    end
end
%}
day = 0;
for t = 1:tper
    Havmoller1(t,:) = uaout(lonL,latB:latT,t);
    Havmoller2(t,:) = vaout(lonL,latB:latT,t);
    %Havmoller(t,:) = elout(lonL,latB:latT,t);
    if (t==1)
        %titlee = ['Time = 0 day 0 hr'];
         ytick(t) = [' 0 day 0 hr'];
    else
        hr = mod(outim/3600,24);
        if (hr==0)
            day = day + 5;
        end
        %titlee = ['Time = ',num2str(day),' day ',num2str(hr),' hr'];
        ytick(t) = [num2str(day),' day ',num2str(hr),' hr'];
        outim = outim + time_interval;
    end
end
    
load('mycmap.mat','mymap');

Havmoller = sqrt(Havmoller1.^2 + Havmoller2.^2);
imagesc(Havmoller);
set(gca,'Ydir','normal');

cmin = -0.5; cmax = 0.5;
colormap(mymap)
caxis([cmin cmax])
cbh = colorbar; % Create Colorbar
cbh.Ticks = linspace(cmin, cmax, 9); % Create 9 ticks from zero to 1
cbh.TickLabels = num2cell(cmin:((cmax-cmin)-0)/8:cmax); % Replace the labels of these 9 ticks with the numbers 1 to 8
cbh.LineWidth = 1.;
%set(get(cbh,'title'),'string','elevation','FontSize',15);
%
set(gca,'ytick',[1:20:tper])
set(gca,'yticklabel',ytick(1,1:(tper-1)/10:tper),'FontSize',12)    
%
set(gca,'xtick',[1:2:yl])
set(gca,'xticklabel',num2cell([lat(latB):2:lon(latT)]),'FontSize',16)
%}
title('Elevation')
%{
set(gca,'ytick',[1:20:tper])
set(gca,'yticklabel',ytick(1,1:(tper-1)/6:tper),'FontSize',12)    
%
set(gca,'xtick',[1:3:xl])
set(gca,'xticklabel',num2cell([lonL:3:lonR]),'FontSize',16)

title('Elevation')
%}
   
