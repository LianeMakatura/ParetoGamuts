function [] = vanilla_visualization(buffers_c, rD, rd, rA)
    numSteps = size(buffers_c, 1);
    
    % ===== pick indices to visualize
    des1_idx = 1;           % between 1 and rD
    des2_idx = 2;
    app_idx = rD + 1;       % between rD+1 and rD+rA
    perf1_idx = rD+rA+1;     % between rD+rA+1 and rD+rA+rd 
    perf2_idx = rD+rA+2;
    disp([des1_idx, des2_idx, app_idx, perf1_idx, perf2_idx])
    
    hold on;
    allPoints = [];
    for step=1:numSteps
        b = buffers_c(step);
        
        paretoIndices = b.getParetoInd;

        for j = 1:length(paretoIndices)
            entry = b.Buff(paretoIndices(j));
            patch = entry.bestPatch;
            if  patch > 0 
                % Store all the points in a single matrix
                newPoint = [entry.minD entry.minF];
                allPoints = [allPoints; newPoint];
            end
        end
    end
    %% Plot all of the discretized pareto fronts, overlayed on top of the previous 3D figure.
    scatter3(allPoints(:,perf1_idx), allPoints(:,perf2_idx), allPoints(:,app_idx));

    figure; hold on;
    scatter3(allPoints(:,des1_idx), allPoints(:,des2_idx), allPoints(:,app_idx));
   
end
