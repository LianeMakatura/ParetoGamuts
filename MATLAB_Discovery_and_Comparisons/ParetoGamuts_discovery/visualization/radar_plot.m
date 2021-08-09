function [] = radar_plot(buffers_c, rD, rd, rA, axis_labels, axis_limits)
    numSteps = size(buffers_c, 1);
    
    axes_labels = axis_labels;
    axes_interval = 10;
    axes_precision = 'none';
    axes_limits = axis_limits;
    
    for step=1:numSteps
        desPts = [];
        perfPts = [];
        
        b = buffers_c(step);
        paretoIndices = b.getParetoInd;

        for j = 1:length(paretoIndices)
            entry = b.Buff(paretoIndices(j));
            patch = entry.bestPatch;
            if  patch > 0 
                desPts = [desPts; entry.minD(1:rD)];
                perfPts = [perfPts; entry.minF];
            end
        end
        spider_plot(desPts, axes_labels, axes_interval, axes_precision, axes_limits);
%         spider_plot(perfPts);
        drawnow;
%         keyboard;
    end
    

end

