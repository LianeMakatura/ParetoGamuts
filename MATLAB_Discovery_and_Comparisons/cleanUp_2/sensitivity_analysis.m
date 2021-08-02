close all;
rng(4);
addPaths;
warning('off','all');

%% First run the original algorithm with the application variable unconstrained. 
patchToPoint = discover_global_front(0.5, true);

%% Then run the algorithm at fixed, incremental values of the application variables.
allPoints = [];

for i = 1:10
    % empirically determined that the true front lies in this range
    zFixed = i*0.1;
    [bufferFixedZ, figures, mFunc] = single_z_front(21, zFixed, true);
    paretoIndices = bufferFixedZ.getParetoInd;

    for j = 1:length(paretoIndices)
        entry = bufferFixedZ.Buff(paretoIndices(j));
        patch = entry.bestPatch;
        if  patch > 0 
            % Store all the points in a single matrix
            newPoint = [entry.minD entry.minF];
            allPoints = [allPoints; newPoint];
        end
    end
end

%% Plot all of the discretized pareto fronts, overlayed on top of the previous 3D figure.
scatter3(allPoints(:,4), allPoints(:,5), allPoints(:,3));

figure; hold on;
scatter3(allPoints(:,1), allPoints(:,2), allPoints(:,3));

%% Finally visualize the sensitvity analysis.
% visualize_sensitvity(patchToPoint, allPoints, true, mFunc.rD, mFunc.rd);