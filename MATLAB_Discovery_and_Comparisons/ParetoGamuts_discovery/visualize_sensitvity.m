function visualize_sensitvity(patchToPoint, allPoints, expandFromGlobalFront, rD, rd)
% Visualize the expansion of a patch from a previously discovered point,
% compared to the global front and constrained fronts.
% -----
% @patchToPoint: Dictionary mapping patch indices (integers) to globally optimal, 
%   continuous patches (as vectors). (unconstrained z)
% @allPoints: [n, rD+rd]-sized unsorted vector of points representing a number of 
%   piecewise-continuous Pareto-optimal points, each individually
%   constrained by (different) application variables.
% @expandFromGlobalFront: Boolean denoting whether to expand from one of the
%   globally-optimal points (i.e. patchToPoint) or one of the locally optimal
%   points (allPoints).
% @rD: Integer dimension of design space.
% @rd: Integer dimension of performance space.

%% Parameters and initialization

% Mapping Functions for objective & gradient evaluations
unconstrainedMF = MappingFunction(23, 0.5);
%unconstrainedMF_expansion = MappingFunction(22, 0.5);

globalPatches = values(patchToPoint);

% Arbitrary values
patchIdx = 4;%3 (kursawe); %4 (sumSines);
pointIdxInPatch = 40;%20 (kursawe); %40(sumSines);

locallyOptimalPointIdx = 2000;

patchSize = 30;
stepSize = 0.01;
z = 0.5; % should do nothing

% Fixed values

% Extract matrix from cell array
optimalPatch = globalPatches(patchIdx);
optimalPatch = optimalPatch{:};

if expandFromGlobalFront
    globallyOptimalPoint = optimalPatch(pointIdxInPatch, :);    
else
    globallyOptimalPoint = allPoints(locallyOptimalPointIdx, :);
end

%% Get the exploration direction(s)

expansionParams = struct('rD', rD, 'rd', rd, 'constrain_z', false);

% test with the old expansion methods (Schulz 2018)
if false
    explorationDirection = getExplorationDirection(globallyOptimalPoint(1:rD),unconstrainedMF, expansionParams);      

    explorationDirection = explorationDirection(:,1);
    expandedPoints = [];
    center = globallyOptimalPoint(1:rD);
    for i = -10:10
        dx = i*0.01*explorationDirection;
        expandedPoints = [expandedPoints; center + dx'];
    end

    figure; hold on;
    for j = 1 : length(globalPatches)
        patch = globalPatches(j);
        patch = patch{:};
        scatter3(patch(:,1), patch(:,2), patch(:,3), 'red');
    end

    scatter3(allPoints(:,1), allPoints(:,2), allPoints(:,3), 'green');
    scatter3(globallyOptimalPoint(1), globallyOptimalPoint(2), globallyOptimalPoint(3), 250, 'black', 'filled');
    scatter3(expandedPoints(:,1), expandedPoints(:,2), expandedPoints(:,3), 'yellow');

elseif true
    % test with the new expansion method
    rA = 1;
    symb = sym('x', [1, rD]);
    symb_z = symb(rD);
    symb_x = symb(1:rD-rA);
    
    x_vals = globallyOptimalPoint(1:rD-rA);
    z_vals = globallyOptimalPoint(rD-rA+1:rD);
    funcs = {unconstrainedMF.f1, unconstrainedMF.f2};

    explorationDirection = getExplorationDirections_xz(funcs, ...
                                symb_x, x_vals, symb_z, z_vals, ...
                                rD-rA, rd, rA);



%     explorationDirection = getExplorationDirection(globallyOptimalPoint(1:rD),unconstrainedMF_expansion, expansionParams);
    
    % Throw away the alpha_primes
    explorationDirection = explorationDirection(rd+1:rD-1+rd+rA, :);
%     explorationDirection = orth(explorationDirection')';
    dir1 = explorationDirection(:,1)';
    dir2 = explorationDirection(:,2)';   

    expandedPoints = zeros((patchSize+1)^2, rD);

    % Center point around which to expand
    centerPt = globallyOptimalPoint(1:rD);
    count = 1;

    %% Expand from the center point in the calculated directions (assumed to be 2)

    for ii = -patchSize:patchSize
        dx1 = ii * dir1 * stepSize;
        for jj = -patchSize:patchSize
            dx2 = jj * dir2 * stepSize;
            offsetPt = centerPt + dx1 + dx2;
            
            expandedPoints(count, :) = offsetPt;
            count = count + 1;
        end
    end

    % Borrowed from expandPatch
    validation = unconstrainedMF.validate(expandedPoints);
    validIndices = find(validation);

    evaluatedPoints = zeros(size(expandedPoints,1), rd);

    evaluatedPoints(validIndices,:) = unconstrainedMF.eval(expandedPoints(validIndices,:), z);

    patchPts = [expandedPoints evaluatedPoints];
    patchPts = patchPts(validIndices, :);

    %% Plot the results

    % Plot Performance space
    figure; hold on;
    for j = 1 : length(globalPatches)
        patch = globalPatches(j);
        patch = patch{:};
        scatter3(patch(:,4), patch(:,5), patch(:,3), 'red');
    end

    scatter3(patchPts(:,4), patchPts(:,5), patchPts(:,3), 'blue');
    scatter3(allPoints(:,4), allPoints(:,5), allPoints(:,3), 'green');
    scatter3(globallyOptimalPoint(4), globallyOptimalPoint(5), globallyOptimalPoint(3), 250, 'black', 'filled');

    % Plot Design space
    figure; hold on;
    for j = 1 : length(globalPatches)
        patch = globalPatches(j);
        patch = patch{:};
        scatter3(patch(:,1), patch(:,2), patch(:,3), 'red');
    end
    scatter3(patchPts(:,1), patchPts(:,2), patchPts(:,3), 'blue');
    scatter3(allPoints(:,1), allPoints(:,2), allPoints(:,3), 'green');
    scatter3(globallyOptimalPoint(1), globallyOptimalPoint(2), globallyOptimalPoint(3), 250, 'black', 'filled');

    end
end
