function HO_plot(elout,uaout,vaout,fsm,outim)

nn = size(elout(:,:,:));
nn = nn(3);
time_interval = outim;

for n = 1:nn
    
    figure('Visible','on');
    hold on
    
    contourf(elout(:,:,n)');
    colorbar; caxis([-10 10]);
    quiver(uaout(:,:,n)',vaout(:,:,n)','k');
    contour(fsm','k')
    
    if (n==1)
        titlee = ['Time = 0 hr'];
    else
        titlee = ['Time = ',num2str(outim/3600),' hr'];
        outim = outim + time_interval;
    end
    title(titlee)
    
    frames(n) = getframe(gcf);
    hold off
    
    % Make GIF file
    %{
    image=frame2im(frames(n));
    [im,cm]=rgb2ind(image,256);
    
    if n==1
        imwrite(im,cm,'Kelvinwave_5day.gif','gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(im,cm,'Kelvinwave_5day.gif','gif','WriteMode','append','DelayTime',0.5);
    end
    %}
        
end

return

