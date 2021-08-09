function [ patchToPoint ] = discover_global_front(z, plotGlobalFront)
% Discover the global front for the unconstrained sum of sines problem.
% ----- Params -----
% @z: Float between [0, 1] denoting the starting value of the application
% variable in the unconstrained sum of sines problem.
% @plotGlobalFront: Boolean denoting whether or not to plot the global front.
% ----- Returns -----
% @patchToPoint: Dictionary mapping patch indices (integers) to globally optimal, 
%   continuous patches (as vectors).

% Params
functionID = 23; % Unconstrained function ID
[buffer, ~, ~] = single_z_front(functionID, z, false);

% We are only interested in the Pareto-optimal front
paretoIndices = buffer.getParetoInd;

%% Create two plots from this data:
% a) The global 2D pareto front
% b) The global 3D pareto front with the application variable itself serving as the 
% last objective function

% Construct a dictionary mapping patch IDs to points
patchToPoint = containers.Map;
for i = 1:length(paretoIndices)
    entry = buffer.Buff(paretoIndices(i));
    patch = entry.bestPatch;
    if  patch > 0 
        patchKey = num2str(patch);

        newPoint = [entry.minD entry.minF];
        if isKey(patchToPoint, patchKey)
            patchToPoint(patchKey) = [patchToPoint(patchKey); newPoint];
        else
            patchToPoint(patchKey) = newPoint;
        end
    end
end

if (plotGlobalFront)
    % Plot the 2D front
    figure; hold on;
    k = keys(patchToPoint);
    v = values(patchToPoint);
    for j = 1 : length(k)
        point = v(j);
        point = point{:};
        scatter(point(:,4), point(:,5));
    end

    % Plot the 3D front
    figure; hold on; grid on;
    for j = 1 : length(k)
        point = v(j);
        point = point{:};
        scatter3(point(:,4), point(:,5), point(:,3));
    end
    
    figure; hold on; grid on;
    for j = 1 : length(k)
        point = v(j);
        point = point{:};
        scatter3(point(:,1), point(:,2), point(:,3));
    end
end

end

