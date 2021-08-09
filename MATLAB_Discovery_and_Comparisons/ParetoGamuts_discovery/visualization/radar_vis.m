function [] = radar_vis(buffers_c, rD, rd, rA, labels, extents)
%RADAR_VIS interactive visualization for the radar plot of design variables
%over time
%   
    fig = uifigure();%('Position',[100 100 350 275]);
%     fig = figure;
    ax = uiaxes(fig);
    x = 1:10;
    hplot = plot(ax, x, 0*x);
    sld = uislider(fig,...
               'Position',[100 75 120 3],...
               'ValueChangingFcn',...%@(sld,event) sliderMoving(event,cg));
                @(sld,event) makeplot(sld, event, x, hplot));
end


function sliderMoving(event,cg)
    cg.Value = event.Value;
end


% function myslider
%     x = 1:10;
%     hplot = plot(x,0*x);
%     h = uicontrol('style','slider','units','pixel','position',[20 20 300 20]);
%     addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,x,hplot));
% end
function makeplot(hObject,event,x,hplot)
    n = get(hObject,'Value');
    set(hplot,'ydata',x.^n);
    drawnow;
end