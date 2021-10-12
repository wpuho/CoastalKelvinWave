function HO_plot2(elout,uaout,vaout,fsm,outim,lon,lat)
%%
%{
mycolors = [0 0 1; 1 1 1; 1 0 0];
num = 4;
mycolors(num,1:3) = 1;
for ii = 1:num
    mycolors(ii,1) = 0;
    mycolors(ii,2) = 1/num*(ii);
    mycolors(ii,3) = 1; 
    mycolors(ii+num,1) = 1;
    mycolors(ii+num,2) = 0;
    mycolors(ii+num,3) = 1/num*(num-ii);
end
%}
%%
lonL = min(lon); lonR = max(lon); latB = min(lat); latT = max(lat);
lenx = length(lon); leny = length(lat);
cmin = -10; cmax = 10;
lonn = meshgrid(lon,lat); latt = meshgrid(lat,lon);
%%
dim = size(fsm);
xl = dim(1); yl = dim(2);
%
for i = 1:xl
    for j = 1:yl
        
        if (fsm(i,j) == 0)
            fsm(i,j) = missing;
            %elout(i,j,5) = missing;
            %vaout(i,j,5) = missing;
            %uaout(i,j,5) = missing;
        end
        
    end
end
%}
%%
load('mycmap.mat','mymap');

res = 2;
nn = size(elout(:,:,:));
nn = nn(3);
time_interval = outim;
day = 0;
for n = nn:nn

    Vmax = max(max(sqrt(uaout(:,:,n).^2+vaout(:,:,n).^2)));
    %figure('Visible','on');
    figure('units','normalized','outerposition',[0 0 1 1],'Visible','on')
    
    imAlpha=ones(size(elout(:,:,n)'));
    imAlpha(isnan(fsm(:,:)'))=0;
    imagesc(elout(:,:,n)','AlphaData',imAlpha);
    %
    colormap(mymap);
    %
    caxis([cmin cmax])
    cbh = colorbar; % Create Colorbar
    cbh.Ticks = linspace(cmin, cmax, 9); % Create 9 ticks from zero to 1
    cbh.TickLabels = num2cell(cmin:((cmax-cmin)-0)/8:cmax); % Replace the labels of these 9 ticks with the numbers 1 to 8
    cbh.LineWidth = 1.;
    cbh.Location = 'southoutside';
    %}
    set(gca,'color',[190/255 190/255 190/255]);
    set(gca,'Ydir','normal');
    %
    hold on
    [~, Min] = min(abs(lat-0));
    eqx = [1 xl]; %eqy = [floor(yl/2)+1 floor(yl/2)+1];
    eqy = [Min Min];
    line(eqx,eqy,'Color','r','LineWidth',1);
    hold off
    %}
    %{
    hAx1=gca;
    hAx2=axes('Position',hAx1.Position,'XAxisLocation','top','YAxisLocation','right','color','none');

    set(hAx2,'xtick',[]); set(hAx2,'ytick',[])
    set(hAx2,'xticklabel',[]); set(hAx2,'yticklabel',[])
    hold(hAx2,'on')
    h = quiver(110.0,67.0,5,0,0,'k','Autoscale','off'); text(115.5,67.0,'5 m/s'); text(127.0,67.0,['Maximum Current Speed = ',num2str(Vmax)]);
    set(h,'MaxHeadSize',50);
    quiver(lonn(1:res:end-1,1:res:end-1),latt(1:res:end-1,1:res:end-1)',uaout(1:res:end-1,1:res:end-1,n)',vaout(1:res:end-1,1:res:end-1,n)','k','Autoscale','off');
    %quiver(uaout(1:res:end-1,1:res:end-1,n)',vaout(1:res:end-1,1:res:end-1,n)','k','Autoscale','off');
    hold(hAx2,'off')
    
    %}
    %
    hold on
    quiver(uaout(:,:,n)',vaout(:,:,n)','k','Autoscale','off');
    h = quiver(10.0,75.0,5,0,0,'k','Autoscale','off'); text(15.5,75.0,'5 m/s'); text(27.0,75.0,['Maximum Current Speed = ',num2str(Vmax)]);
    set(h,'MaxHeadSize',50);
    [~, Min] = min(abs(lat-0));
    eqx = [1 xl]; %eqy = [floor(yl/2)+1 floor(yl/2)+1];
    eqy = [Min Min];
    line(eqx,eqy,'Color','r','LineWidth',1);
    hold off
    %}
    %
        hAx1=gca;
    set(hAx1,'ytick',[1:(yl-1)/8:leny])
    set(hAx1,'yticklabel',num2cell([latB:(latT-latB)/8:latT]),'FontSize',16)    
    set(hAx1,'xtick',[1:(lenx-1)/8:lenx])
    set(hAx1,'xticklabel',num2cell([lonL:(lonR-lonL)/8:lonR]),'FontSize',16)
    %}
    %{
    if (n==1)
        titlee = ['Time = 0 hr'];
    else
        titlee = ['Time = ',num2str(outim/3600),' hr'];
        outim = outim + time_interval;
    end
    title(titlee,'FontSize',24)
    %}
    if (n==1)
        titlee = ['Time = 0 day 0 hr'];
    else
        hr = mod(outim/3600,24);
        if (hr==0)
            day = day + 1;
        end
        titlee = ['Time = ',num2str(day),' day ',num2str(hr),' hr'];
        outim = outim + time_interval;
    end
    title(titlee,'FontSize',24)
    
    %frames = getframe(gcf);
    
    %set(gcf,'unit','centimeters','position',[5 3 25 17])
    
    % Make GIF file
    %{   
    image=frame2im(frames);
    [im,cm]=rgb2ind(image,256);
    if n==1
        imwrite(im,cm,'Kelvinwave_test2.gif','gif','LoopCount',Inf,'DelayTime',0.125);
    else
        imwrite(im,cm,'Kelvinwave_test2.gif','gif','WriteMode','append','DelayTime',0.125);
    end
    %}

end

return