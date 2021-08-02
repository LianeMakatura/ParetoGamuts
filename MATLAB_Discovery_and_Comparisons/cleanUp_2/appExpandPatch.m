function [patchPts, connectivityList, expDirs] = appExpandPatch(centerPt, explorationDirection, ...
                            mFunc, buffer, extents, userParams)
    rD = mFunc.rD;
    rd = mFunc.rd;
    rA = mFunc.rA;
                            
    % Degenerate Case 1: More than d+A-1 directions to explore 
    if size(explorationDirection,2) > rd + rA - 1
        dispstat('More than rd-1+rA directions; choosing subset.', 'keepprev', 'keepthis')
        % Pick d-1 random directions from the selection
        indices = randsample(size(explorationDirection,2),rd+rA-1);
        explorationDirection = explorationDirection(:,indices);
    end

    % Degenerate Case 2: Fewer than d+A-1 directions to explore
    if size(explorationDirection,2) < rd+rA-1
        % Do not expand the point
        dispstat('Degenerate; fewer than rd-1+rA directions.', 'keepprev', 'keepthis')
        patchPts = []; connectivityList = []; expDirs = explorationDirection;
        return
    end
    
    % Orthonormalize the directions
%     explorationDirection = orth(explorationDirection);
    if size(explorationDirection,2) < rd+rA-1
        % Do not expand the point
        dispstat('Orthonormalized directions degenerate; fewer than rd-1+rA directions.', 'keepprev', 'keepthis')
        patchPts = []; connectivityList = []; expDirs = explorationDirection;
        return
    end
    expDirs = explorationDirection; % update it with the final directions
     
    % Only keep patches that haven't already been discovered
    alreadyExplored = buffer.checkIfAlreadyExpanded(centerPt, explorationDirection);
    if alreadyExplored
        dispstat('Patch already found.', 'keepprev', 'keepthis')
        patchPts = []; connectivityList = [];
        return
    end
                        
%     stepSize = userParams.stepSize;  
    stepSize =  extents .* (userParams.stepSize/5); 
    centerDesApp = centerPt(1:rD+rA);
        
    % ====== Expand the design patch ======
    % create the range for each dimension
    numDim = rd + rA - 1;
    ranges = cell(numDim, 1);
    range_len = zeros(numDim, 1);
    for i=1:numDim
        rangei = -extents(i):stepSize(i):extents(i); % fix this to have proper min/max for each dim
        ranges{i} = rangei;
        range_len(i) = length(rangei); % number of points in each dimension
    end
        
    numpoints = prod(range_len);
    stepIdx = cell(1, numDim);
    expandedPoints = zeros(numpoints, rD+rA);
    connPts = zeros(numpoints, numDim); % to store the offsetgrid for connectivity
    for ptidx = 1:numpoints
        offsetPt = centerDesApp; %reset the origin
        [stepIdx{:}] = ind2sub(range_len, ptidx);  % get index into each dimension
        connPts(ptidx, :) = [stepIdx{:}];

        % get design point; index into each dim range, add offsets to center
        for i=1:numDim
            r = ranges{i};
            step = r(stepIdx{i});
            offsetPt = offsetPt + step * explorationDirection(:,i)'; % row vector
        end
        
        expandedPoints(ptidx, :) = offsetPt;
    end
    

    % Filter out any points that violate constraints
    validation = mFunc.validate(expandedPoints);
    validIndices = find(validation);
   
    %evaluate the final patch
    evaluatedPoints = zeros(size(expandedPoints,1), rd);
    evaluatedPoints(validIndices,:) = mFunc.eval(expandedPoints(validIndices,1:rD), ...
                                                expandedPoints(validIndices,rD+1:rD+rA)); 
    patchPts = [expandedPoints evaluatedPoints];
    patchPts = patchPts(validIndices, :);
    
    
    % Build a connectivity list that only includes tris/tets between valid
    % vertices
    connPts = connPts(validIndices, :);
    if isempty(connPts)
        connectivityList = [];
        return
    end
    rankDeficient = (rank(connPts(2:end,:) - connPts(1,:)) < numDim);
        
    if userParams.fillHoles_bufferProjection && (numDim == 2 || numDim == 3) && size(validIndices, 1) >= 3 ...
        && ~rankDeficient
        connectivityList = delaunay(connPts);
    else
        connectivityList = [];
    end
    
end
