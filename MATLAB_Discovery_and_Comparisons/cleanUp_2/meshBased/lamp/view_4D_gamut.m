% GENERATING THE LAMP GAMUT IMAGE
function [] = view_4D_gamut(l_buffer, l_pts)
    lpoints_rescaled = l_pts;
    lpoints_rescaled(:, l_buffer.rD+l_buffer.rA+1) = lpoints_rescaled(:, l_buffer.rD+l_buffer.rA+1) / 3;
    l_perf1idx = l_buffer.rD + l_buffer.rA + 1;
    names = l_buffer.mFunc.varNames;

    
    figure;
    scatter3(lpoints_rescaled(:, l_perf1idx), lpoints_rescaled(:, l_perf1idx+1), lpoints_rescaled(:, l_perf1idx+2), 40, lpoints_rescaled(:, l_perf1idx-1));
    xlabel(names(l_perf1idx)); 
    ylabel(names(l_perf1idx + 1)); 
    zlabel(names(l_perf1idx + 2)); 
    hcb = colorbar;
    colormap jet;
    colorTitleHandle = get(hcb,'Title');
    titleString = 'Focal Height';
    set(colorTitleHandle ,'String',titleString);
end