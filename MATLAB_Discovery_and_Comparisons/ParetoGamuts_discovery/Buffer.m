classdef Buffer < handle
%% Instance variables
properties
    % Static variables
    rD
    rd
    numCells
    validBufferIndices
    nNeighbours
    delta
    tolerance
    directionDelta
    improvementBreakScore
    distanceBreakScore
    numRepeatedNoImprovement

    % Dynamic variables
    Buff
    emptyCells
    currentDistance
    directionScores
    nPatches
    nNeighboursOpt
    selectedPatches
    finalScores
    finalDist
    converged

    expandedCenters
    expandedDirections
end

%% Functions
methods
    %% Constructor.
    function obj = Buffer(userParams, fixedParams)
        if nargin==0 % needed for Object arrays
            return
        end
        
        obj.expandedCenters = [];
        obj.expandedDirections = [];

        nNeighbours = userParams.neighForScore;
        nNeighboursOpt = userParams.neighForOpt;
        obj.nNeighboursOpt = nNeighboursOpt;
        obj.directionDelta = userParams.directionDelta;
        % Number of discrete cells in buffer
        numCells = userParams.numCells;
        obj.numCells = numCells;
        % Termination metrics
        obj.improvementBreakScore = userParams.breakScoreRelative;
        obj.distanceBreakScore = userParams.breakScoreDistance;
        obj.numRepeatedNoImprovement = userParams.numRepeatedNoImprovement;
        obj.finalScores = zeros(1, userParams.nMaxRuns);
        obj.finalDist  = zeros(1, userParams.nMaxRuns); 
        obj.converged = false;

        % Dimensionality of design/objective space
        obj.rD = fixedParams.rD;
        obj.rd = fixedParams.rd;
        obj.delta = 1/numCells;
        obj.tolerance= userParams.tolerance;

        % info in each Buff cell
        basicInfo = struct('solutions', [], 'minD', zeros(1,obj.rD), 'minF', zeros(1,obj.rd), 'bestPatch', 0,  'neighbours', 0 );

        if (obj.rd == 2)
            obj.validBufferIndices=[1:obj.numCells]';
            obj.Buff = repmat(basicInfo, 1, numCells);
            obj.emptyCells = ones(1, numCells);
            obj.currentDistance = inf*ones(1, numCells);
            obj.directionScores = zeros(1, numCells);

            % assign neighbors to each cell, nNeighbours in each direction
            % except for corner cases
            for i=(1+nNeighbours):(numCells-nNeighbours)
                obj.Buff(i).neighbours = [i-nNeighbours:i+nNeighbours]';
            end

            for i=1:nNeighbours
                obj.Buff(i).neighbours = [1:i+nNeighbours]';
            end  

            for i=(numCells-nNeighbours +1):numCells
                obj.Buff(i).neighbours = [i-nNeighbours:numCells]';
            end 
        else
            [X, Y] = meshgrid(1:numCells, 1:numCells);
            sumXY = X+Y;
            obj.validBufferIndices=find(sumXY <= numCells +1);              
            obj.Buff = repmat(basicInfo, numCells, numCells);
            obj.emptyCells = zeros(numCells, numCells);
            obj.emptyCells(obj.validBufferIndices) = 1;
            obj.currentDistance = zeros(numCells, numCells);
            obj.currentDistance(obj.validBufferIndices) = inf;
            obj.directionScores = zeros(numCells, numCells);

            [X, Y] = meshgrid(-nNeighbours:nNeighbours,-nNeighbours:nNeighbours);
            SUMXY  = X+Y;
            for i=1:obj.numCells
                for j = 1:obj.numCells
                    X_pos = X+i;
                    Y_pos = Y+j;
                    pass1 = find(X+i >0);
                    pass2 = find(Y_pos(pass1) >0);
                    positivePos = pass1(pass2);
                    pass3 = find( SUMXY(positivePos)  <=  numCells -i - j +1);
                    finalPos = positivePos(pass3);
                    obj.Buff(i,j).neighbours = [X_pos(finalPos) Y_pos(finalPos)];
                end
            end                         
        end     
        disp('Instantiated buffer');
    end
        
    %% Update the buffer with a new point on a given patch.
    % @param idx: Integer (potentially 2D) cell index to be updated
    % @param pointD: rD-sized value of the new point in design space
    % @param pointF: rd-sized value of the new point in objective space
    % @param dist: Float value of the distance  from the origin of the new point
    % @param patchID: Integer ID of the patch associated with the point
    % ---
    % Appropriately updates the solution(s) at the given cell index.
    function updateBuffer(obj, idx, pointD, pointF, dist, patchID)
        % Obtain the proper int-valued cell index
        if(length(idx) > 1)
            idx = sub2ind(size(obj.Buff), idx(1), idx(2));
        end
        
        
        solution = struct('D', pointD , 'F', pointF, 'dist', dist, 'patch', patchID);

       % If the new point is within our tolerance threshold
       if (abs(obj.currentDistance(idx) - dist) < obj.tolerance)
           % Do not add a lone point if there exist other solutions already
           if(patchID ~=0) 
                obj.Buff(idx).solutions = [obj.Buff(idx).solutions solution];
           end   
           % If there is only one point and it is within tolerance, remove point
           if (obj.Buff(idx).solutions(1).patch == 0) 
                obj.Buff(idx).solutions = solution;
                if (size(obj.Buff(idx).solutions) > 1)
                    % TODO: Remove this? I don't think this is reachable
                    disp('Something went wrong. Press any key to unpause...')
                    pause;
                end
           end
       else 
            obj.Buff(idx).solutions = solution;
       end

       % Update the optimal solutions for the cell
       if obj.currentDistance(idx) > dist
           obj.currentDistance(idx) = dist;
           obj.Buff(idx).minD = pointD;
           obj.Buff(idx).minF = pointF;
           obj.Buff(idx).bestPatch = patchID;
       end  
       
       % Close out the empty index for the cell (only relevant for the
       % first point/patch added).
       obj.emptyCells(idx) = 0;
    end
        
    %% Add a single point to the buffer.
    % @param point: length (rD + rd) vector representing the point to be
    % added
    function addPoint(obj, point)
        pointD = point(1:obj.rD);
        pointF = point(obj.rD+1:obj.rD+obj.rd);
        [idx, dist] = obj.getCellIndex(pointF);
        % Only add the lone point if it is better than all other solutions.
        if (obj.currentDistance(idx) > (dist))
            updateBuffer(obj, idx, pointD, pointF, dist, 0)
        end
    end
        
    %% Get the cell indices and alpha values for a line of points.
    function [indices, alphas] = getIndicesInLine(obj, projA, projB)
            if(	projA < projB)
                projMin = projA;
                projMax = projB;
            else
                projMin = projB;
                projMax = projA;
            end

            minIdx = ceil((projMin/obj.delta) + 1/2);
            maxIdx = ceil((projMax/obj.delta) + 1/2) -1;
            indices = minIdx:maxIdx;
            alphas = zeros(length(indices), 2); 
            for i = 1: length(indices)
                idx = indices (i);
                projCenter = idx*obj.delta - obj.delta/2;
                alpha = (projB - projCenter)/(projB - projA);
                alphas(i, :) = [alpha, (1 - alpha)];                    
            end
    end
        
    %% Get the cell indices and alpha values for a simplex of points.
    function [indices, alphas] = getIndicesInTriangle(obj, projA, projB, projC)
        indices = [];
        alphas = [];

        x_vals = [projA(1), projB(1), projC(1)];
        y_vals = [projA(2), projB(2), projC(2)];

        minIdx = ceil((min(x_vals)/obj.delta) + sqrt(2)/2);
        maxIdx = ceil((max(x_vals)/obj.delta) + sqrt(2)/2) -1;
        minIdy = ceil((min(y_vals)/obj.delta) + sqrt(2)/2);
        maxIdy = ceil((max(y_vals)/obj.delta) + sqrt(2)/2) -1;

        for i= minIdx:maxIdx
            for j= minIdy:maxIdy
                projCenter = [i*obj.delta - sqrt(2)/2*obj.delta, j*obj.delta - sqrt(2)/2*obj.delta];         
                % get baricentric coordinates
                u = projB - projA;
                v = projC - projA;
                X = linsolve([u' v'], [projCenter - projA]');
                if (min(X) > 0 && sum(X)< 1)
                    indices = [indices; i j];
                    alphas = [alphas; X' (1-sum(X))];
                end
            end
        end
    end
        
    %% Add an expanded patch to the buffer.
    function addPatch(obj, patch, patchID, validation, direction, center)
        obj.expandedDirections = [obj.expandedDirections; direction'];
        obj.expandedCenters = [obj.expandedCenters; center(1:obj.rD)];

        if obj.rd == 2
            [projections, distances] = obj.projectPoint(patch(:,obj.rD+1:obj.rD+obj.rd));

            for p=1:(size(patch,1)-1)
                if (validation(p)*validation(p+1) > 0)
                    aPointF = patch(p,obj.rD+1:obj.rD+obj.rd);
                    bPointF = patch(p +1,obj.rD+1:obj.rD+obj.rd);
                    aPointD = patch(p,1:obj.rD);
                    bPointD = patch(p +1,1:obj.rD);
                    projA = projections(p,:);
                    projB = projections(p +1,:);
                    distA = distances(p);
                    distB = distances(p +1);                        
                    [indices, alphas] = obj.getIndicesInLine(projA, projB);

                    for i = 1: length(indices)
                        idx = indices (i);
                        optF = aPointF*(alphas(i, 1)) + bPointF*(alphas(i, 2));
                        optD = aPointD*(alphas(i, 1)) + bPointD*(alphas(i, 2));
                        optDist = distA*(alphas(i, 1)) + distB*(alphas(i, 2));
                        if (obj.currentDistance(idx) > (optDist - obj.tolerance))
                            updateBuffer(obj, idx, optD, optF, optDist, patchID);
                        end
                    end
                end
            end
        else
            patchPerm = permute(patch,[3 1 2]);
            patchRows =  patchPerm(:, :)';

            [projections, distances] = obj.projectPoint(patchRows(:,obj.rD+1:obj.rD+obj.rd));

            SIZ = [ size(patch,1) size(patch,2)];
            for p1=1:(size(patch,1)-1)
                for p2=1:(size(patch,2)-1)

                    tri1IsValid = validation(p1, p2)*validation(p1+1, p2)*validation(p1, p2+1);
                    tri2IsValid = validation(p1 +1, p2 +1)*validation(p1+1, p2)*validation(p1, p2+1);

                    if  ((tri1IsValid > 0) || (tri2IsValid > 0))

                        pa = (p2 -1)*SIZ(1) + p1;
                        pb = (p2 -1)*SIZ(1) + p1 + 1;
                        pc = (p2)*SIZ(1) + p1;
                        pd = (p2 )*SIZ(1) + p1 + 1;

                        aPointF = patchRows(pa, obj.rD+1:obj.rD+obj.rd)';
                        bPointF =  patchRows(pb, obj.rD+1:obj.rD+obj.rd)';
                        cPointF =  patchRows(pc, obj.rD+1:obj.rD+obj.rd)';
                        dPointF =  patchRows(pd, obj.rD+1:obj.rD+obj.rd)';                   
                        aPointD = patchRows(pa, 1:obj.rD)';
                        bPointD = patchRows(pb, 1:obj.rD)';
                        cPointD = patchRows(pc, 1:obj.rD)';
                        dPointD = patchRows(pd, 1:obj.rD)';

                        [projA] = projections(pa, :);
                        [projB] = projections(pb, :);
                        [projC] = projections(pc, :);
                        [projD] = projections(pd, :);

                        [distA] = distances(pa, 1);
                        [distB] = distances(pb, 1);
                        [distC] = distances(pc, 1);
                        [distD] = distances(pd, 1);
                        if (tri1IsValid > 0)
                            [indices, alphas] = obj.getIndicesInTriangle(projA, projB, projC);
                            for i = 1: size(indices,1)
                                idx = sub2ind(size(obj.Buff), indices(i,1), indices(i,2));
                                optF = aPointF*(alphas(i, 1)) + bPointF*(alphas(i, 2)) + cPointF*(alphas(i, 3));
                                optD = aPointD*(alphas(i, 1)) + bPointD*(alphas(i, 2))+ cPointD*(alphas(i, 3));
                                optDist = distA*(alphas(i, 1)) + distB*(alphas(i, 2))+ distC*(alphas(i, 3));
                                if (obj.currentDistance(idx) > (optDist - obj.tolerance))
                                    updateBuffer(obj, idx, optD', optF', optDist, patchID);
                                end
                            end
                        end
                        if (tri2IsValid > 0)
                            [indices, alphas] = obj.getIndicesInTriangle(projD, projB, projC);
                            for i = 1: size(indices,1)
                                idx = sub2ind(size(obj.Buff), indices(i,1), indices(i,2));
                                optF = dPointF*(alphas(i, 1)) + bPointF*(alphas(i, 2)) + cPointF*(alphas(i, 3));
                                optD = dPointD*(alphas(i, 1)) + bPointD*(alphas(i, 2))+ cPointD*(alphas(i, 3));
                                optDist = distD*(alphas(i, 1)) + distB*(alphas(i, 2))+ distC*(alphas(i, 3));
                                if (obj.currentDistance(idx) > (optDist - obj.tolerance))
                                    updateBuffer(obj, idx, optD', optF', optDist, patchID);
                                end
                            end
                        end  
                    end
                end
            end
        end 
    end
        
    %% Determine whether a point and direction have already been expanded.
    function [haveExpanded] = checkHaveExpanded(obj, center, direction)
        if size(obj.expandedCenters,1) == 0
            haveExpanded = false;
            obj.expandedDirections = [obj.expandedDirections; direction'];
            obj.expandedCenters = [obj.expandedCenters; center];   
            return;
        end

        if obj.rd == 2
           for i = 1:size(obj.expandedCenters, 1)  
                storedDirection = obj.expandedDirections(i,:)';
                dirMatrix = [storedDirection direction];
                [~, S, ~] = svd(dirMatrix);
                singVal = S(2,2);
                if singVal < 0.05 %0.0001 before tweaking for rebuttal *********************
                    storedCenter = obj.expandedCenters(i,:);
                    centerDisplacement = storedCenter - center;
                    dirPointMatrix = [storedDirection centerDisplacement'];

                    [~, S, ~] = svd(dirPointMatrix);
                    singVal2 = S(2,2);
                    if singVal2 < 0.05
                        haveExpanded = true;
                        return;
                    end
                end
            end

        elseif obj.rd == 3   
            for i = 0:size(obj.expandedCenters, 1)-1

                storedDirection = obj.expandedDirections(i*2+1:i*2+2,:)';
                dirMatrix = [storedDirection direction];
                [~, S, ~] = svd(dirMatrix);
                singVal = S(3,3);
                if singVal < 0.0001
                    storedCenter = obj.expandedCenters(i+1,:);
                    centerDisplacement = storedCenter - center;
                    dirPointMatrix = [storedDirection centerDisplacement'];

                    [~, S, ~] = svd(dirPointMatrix);
                    singVal2 = S(3,3);
                    if singVal2 < 0.0001
                        haveExpanded = true;
                        return;
                    end
                end
            end
        end
        haveExpanded = false;
    end
                
    %% Update the scores of the buffer based on discovered solutions.
    function [bufferScore, overallDist] = updateScores(obj)
        sumBufferScore = 0;
        for idx=obj.validBufferIndices'
            if(obj.rd > 2)
                neighInd = sub2ind(size(obj.Buff), obj.Buff(idx).neighbours(:,1), obj.Buff(idx).neighbours(:,2));
            else
                neighInd = obj.Buff(idx).neighbours;
            end
            emptyArray = obj.emptyCells(neighInd);
            popScore = sum(emptyArray)/length(emptyArray);

            nFilledNeighbours = length(emptyArray) -sum(emptyArray);
            discontinuityScore = 0;
            if (nFilledNeighbours > 1)
                    distArray = obj.currentDistance(neighInd);
                    averageDistance = (sum(distArray))/nFilledNeighbours;
                    if (obj.currentDistance(idx) > averageDistance)
                        discontinuityScore = (obj.currentDistance(idx) - averageDistance)/obj.currentDistance(idx);
                    end
            else
                discontinuityScore = 1;
            end
            obj.directionScores(idx) =  popScore + discontinuityScore;
            sumBufferScore = sumBufferScore + obj.directionScores(idx);
        end

        bufferScore =  sumBufferScore/length(obj.validBufferIndices);
        fullCells =  find(obj.emptyCells(obj.validBufferIndices(:)) == 0);
        fullCells = obj.validBufferIndices(fullCells);
        overallDist = sum(obj.currentDistance(fullCells))/length(fullCells);

        disp(['Improvement Score: ' num2str(bufferScore) ' Distance score: '  num2str(overallDist)]);
    end

    %% Determine the optimization targets on the buffer for a list of samples.
    function [targets, dirs] = getRandomDirections(obj, samples)
        dirDelta = obj.directionDelta;
        targets = zeros(size(samples));
        dirs = zeros(size(samples));
        for i=1:size(samples,1)
            idx = obj.getCellIndex(samples(i,:));

            neighFound = false;
            while(~neighFound)
                neighFound = true;
                optInd = idx + randi([-obj.nNeighboursOpt, obj.nNeighboursOpt], size(idx));
                if (any(optInd< 1) || any(optInd > obj.numCells) || sum(optInd) >  obj.numCells +1) 
                    neighFound = false;
                end                    
            end

            goal = 2*getCenter(obj, optInd);
            origin = samples(i,:);
            dir = goal - origin/norm(origin);
            dir = dir/norm(dir);

            amountToPush  = dirDelta*norm(origin);

            targets(i,:) = samples(i,:) + amountToPush*dir;     
            dirs(i, :) = dir;
        end                     
    end
        
    %% Stochastically obtain a list of points to be optimized in the next iteration of the main loop.
    function [randomPoints] = getPointsToImprove(obj, nSamples, mFunc, z)     
        idx = find((obj.emptyCells(obj.validBufferIndices) == 0));
        idx = obj.validBufferIndices(idx); % all valid, non-empty cells
        n = length(idx);
        filtering = ceil(rand(nSamples,1)*n); %pick a random idx # for each sample
        idxActual = idx(filtering);
        randomPoints = zeros(nSamples, obj.rD);
        
        i = 1;
        while i < nSamples+1
            val = idxActual(i);
            sigma = 1./2.^(10*rand(1, 1));
            randPt = obj.Buff(val).minD - sigma*(2*rand(size(obj.Buff(val).minD))-1); 
            randPt(obj.rD) = obj.Buff(val).minD(obj.rD);% don't change the app var
            randomPoints(i,:) = randPt;
            i = i+1;
        end
        randomPoints = min(randomPoints, 1);
        randomPoints = max(randomPoints, 0.00001);     
    end
    
    %% Check if Buffer is empty
    function [empty] = isEmpty(obj)
        filled = find(obj.emptyCells(obj.validBufferIndices)==0);
        if isempty(filled)
            empty = true;
        else
            empty = false;
        end
    end

    %% Project a point in space onto the buffer
    function [proj, dist] = projectPoint(obj, q)
        dist = sqrt(sum(q.^2,2));          
        p = q./dist; 
        m = pi/2;

        if(obj.rd >2)
            t2 = atan(p(:,2)./p(:,1))/m;
            t2 = t2 - floor(t2);
            t2(find(p(:,1) < 0.000001)) = 1;
            t2(find(p(:,2) < 0.000001)) = 0;
            proj = [1 - acos(p(:,3))/m,  acos(p(:,3))/m.*t2 ];
        else
            proj =  1 - atan(p(:,2)./p(:,1))/m;                    
            proj(find(p(:,1) < 0.000001)) = 0;
            proj(find(p(:,2) < 0.000001)) = 1;
        end
    end
        
    %% Get the center point of a cell by index
    function center = getCenter(obj, idx)
        alpha = (idx-1/2)/obj.numCells;

        if (obj.rd > 2)              
            % mirror because of buffer shape
            m_1 = 1 -alpha(1); 
            m_2 = alpha(2);

            % Rescale triangle
            a = (m_1)/2;
            b = m_2 - m_1./2;

            % rotate to first quadrant
            rot_a = (a - b);
            rot_b = (a + b);

            alpha = [rot_a rot_b];                
        end      
        
        pushUpVector = ones(1, obj.rd)/obj.rd;
        center = [alpha 1-sum(alpha)] - pushUpVector;
    end
         
    %% Get the cell index and distance-from-origin of a point.
    function [idx, dist] = getCellIndex(obj, q)
        [proj_q, dist] = obj.projectPoint(q);
        idx = min(floor(proj_q/obj.delta) +1, obj.numCells);
        idx = max(idx, 1);
    end
        
    %% Check whether the convergence criteria have been met and the main loop can terminate.
    function converged = terminate(obj, iteration)
        [obj.finalScores(iteration), obj.finalDist(iteration)]  = obj.updateScores();
        N = obj.numRepeatedNoImprovement;
        improved = false;
        close = false;
        
        if(iteration> N + 1)
            if(abs(obj.finalScores(iteration) - obj.finalScores(iteration - N)) < obj.improvementBreakScore)
                improved = true;
            end
            if(abs(obj.finalDist(iteration) - obj.finalDist(iteration - N)) < obj.distanceBreakScore)
                close = true;
            end            
        end
        converged = improved && close;
    end
    
    function hv = getHyperVolume(obj)
        % Add HyperVolume calculation
        nadirPt = ones(obj.rd, 1);
        idx = obj.getParetoInd();
        pointsF = [];
        for i=1:size(idx, 1)
            if ~obj.emptyCells(idx(i))
                pointsF = [pointsF; obj.Buff(idx(i)).minF];
            end
        end
        hv = lebesgue_measure(pointsF', nadirPt);
    end
    
    %% Get the indices of Pareto-optimal cells in the buffer.    
    function ind = getParetoInd(obj)
        validCells = obj.validBufferIndices(find(obj.emptyCells(obj.validBufferIndices) == 0));

        A = zeros(size(validCells,1), obj.rd);
        i = 1;
        for idx = validCells'
            A(i,:) = obj.Buff(idx).minF;
            i = i+ 1;
        end
        pInd = getParetoIndices(A);
        ind = validCells(find(pInd));
    end

    %% Unlabel all patches in the buffer.
    function unlabelPatches(obj)
        validIndices = obj.validBufferIndices(obj.emptyCells(obj.validBufferIndices) == 0);
        for i = 1:length(validIndices)
            cell = obj.Buff(validIndices(i));
            [~, I] = min([cell.solutions.dist]);
            sol = cell.solutions(I);
            obj.Buff(validIndices(i)).minD = sol.D;
            obj.Buff(validIndices(i)).minF = sol.F;
            obj.Buff(validIndices(i)).bestPatch = sol.patch;
        end
    end
        
    %% Fill the empty cells post-graphcut.
    % Post processing function to re-expand patches around outlier points that 
    % were incorrectly labeled due to graph-cut imprecision.
    function fillEmptyCells(obj, stepSize, mFunc, z)
        if obj.rd == 2
            return;
        end

        validCellIndices = obj.validBufferIndices(find(obj.emptyCells(obj.validBufferIndices) == 0));
        nCells = length(validCellIndices);
        % Iterate through all the valid cells in the buffer
        for i = 1:nCells 
            cellIdx = validCellIndices(i);
            buff = obj.Buff(cellIdx); % Get the cell object
            [ii, jj]= ind2sub(size(obj.Buff),  cellIdx);
            % Speedup alternative to `[X,Y] = meshgrid(-1:1,-1:1)`
            X = [-1 0 1; -1 0 1; -1 0 1];
            Y = [-1 -1 -1; 0 0 0; 1 1 1];

            p = [X(:)+ ii, Y(:)+ jj]; % Get all the neighbors
            p(ceil(size(p,1)/2),:) = []; % Get rid of the middle entry (which is the original point)
            patchesToConsider = zeros(size(p,1),1);
            isOutlier = true; % Dynamic indicator boolean
            minDistFromNeighbor = inf; 

            % For all neighbors
            for j = 1:size(p,1) 
                point = p(j,:);
                if max(point) > obj.numCells || min(point) < 1
                    % Do nothing if the neighbor falls outside the
                    % bounds
                    continue
                end

                neighborPatch = obj.selectedPatches(point(1),point(2));

                distFromNeighbor = norm(buff.minD - obj.Buff(point(1),point(2)).minD);
                if distFromNeighbor < minDistFromNeighbor
                    minDistFromNeighbor = distFromNeighbor;
                end

                if neighborPatch == obj.selectedPatches(ii,jj)
                    % Break if the original point borders another point within the same patch
                    %(i.e. is not an outlier)
                    isOutlier = false;
                end
                patchesToConsider(j) = neighborPatch;
            end

            if obj.selectedPatches(ii,jj) == 0 || isOutlier || minDistFromNeighbor > .125

                options = unique(patchesToConsider);
                sortedOptions = [options, histc(patchesToConsider(:),options)];

                if size(sortedOptions,1) == 1 && sortedOptions(1,1) ~= 0
                    meanD = zeros(1,obj.rD);
                    meanF = zeros(1,obj.rd);
                    for w = 1:size(p,1)
                        neighborBuff = obj.Buff(p(w,1),p(w,2));
                        meanD = meanD + neighborBuff.minD;
                        meanF = meanF + neighborBuff.minF;
                    end
                    obj.Buff(ii,jj).minD = meanD/size(p,1);
                    obj.Buff(ii,jj).minF = meanF/size(p,1);
                    obj.selectedPatches(ii,jj) = sortedOptions(1,1);
                    continue;
                end

                haveFilled = false;
                for k = 1:size(sortedOptions,1)
                    patchToCheck = sortedOptions(k,1);
                    if patchToCheck == 0
                        continue
                    end

                    for h = 1:size(p,1)
                        neigh = p(h,:);
                        if max(neigh) > obj.numCells || min(neigh) < 1
                            % Do nothing if the neighbor falls outside the
                            % bounds
                            continue
                        end                            
                        if obj.selectedPatches(neigh(1),neigh(2)) == patchToCheck
                            centerForExpansion = [obj.Buff(neigh(1),neigh(2)).minD obj.Buff(neigh(1),neigh(2)).minF];
                            break
                        end
                    end
                    % Re-expand the patch
                    directions = obj.expandedDirections(patchToCheck*2-1:patchToCheck*2,:);

                    userParams = struct('stepSize', stepSize/2, 'patchSize', 10);
                    [patch, validation] = expandPatch(directions', centerForExpansion, mFunc, obj, userParams, true, z);  
                    validIndices = find(validation);

                    % Do a bunch of reformating/validation of the patch
                    patchPerm = permute(patch,[3 1 2]);
                    patchRows =  patchPerm(:, :)';
                    patchRowsForProjection = patchRows(:,obj.rD+1:obj.rD+obj.rd);

                    [projections, ~] = obj.projectPoint(patchRowsForProjection);

                    realLocation = [ii jj];

                    % Loop through the points in the new (higher
                    % resolution) patch and attempt to fill the empty
                    % buffer cell.
                    for z = 1:size(validIndices,1)
                        validIndex = validIndices(z);
                        projectedPoint = projections(validIndex,:);
                        projInd = min(floor(projectedPoint/obj.delta)+1,obj.numCells);
                        projInd = max(projInd,1);
                        if isequal(realLocation,projInd)
                            replacementPoint = patchRows(validIndex,:);
                            newD = replacementPoint(1:obj.rD);
                            newF = replacementPoint(obj.rD+1:obj.rD+obj.rd);
                            obj.Buff(cellIdx).minD = newD;
                            obj.Buff(cellIdx).minF = newF;
                            obj.selectedPatches(ii,jj) = patchToCheck;

                            haveFilled = true;
                            break;
                        end
                    end
                    if haveFilled
                        break;
                    end
                end
                % If we cannot find a point, fill the empty cells using
                % the best stored distance.
                if ~haveFilled
                    [~, distIdx] = min([buff.solutions.dist]);
                    sol = buff.solutions(distIdx);
                    obj.Buff(ii,jj).minD = sol.D;
                    obj.Buff(ii,jj).minF = sol.F;
                    obj.selectedPatches(ii,jj) = sol.patch;
                end
            end

        end
    end
        
    %% Use graph-cut to optimally label the buffer patches.
    function [labels, chosenLabels]=  labelPatches(obj, nPatches, nSize)
        addpath graphCut

        obj.selectedPatches = zeros(size(obj.Buff));
        if (nargin <3)
           nSize = 2;
        end 

        validCells = obj.validBufferIndices(find(obj.emptyCells(obj.validBufferIndices) == 0));

        nNodes = length(validCells);
        nLabels = nPatches;
        segclass = zeros(nNodes,1);
        pairwise = sparse(nNodes,nNodes);
        unary = 10*ones(nLabels, nNodes);
        validIndToNodeMap =zeros(max(validCells), 1);

        for i = 1:nNodes
            idx = validCells(i);
            patches = zeros(size(obj.Buff(idx).solutions)); 
            distances = zeros(size(obj.Buff(idx).solutions)); 
            for s=1:length(obj.Buff(idx).solutions)
                patches(s) = obj.Buff(idx).solutions(s).patch; 
                distances(s) = obj.Buff(idx).solutions(s).dist; 
            end
             minError = min(distances);
             maxError = max(distances) - minError;
             if( maxError < 0.000001)
                 unary(patches + 1, i) = 0;
             else           
                unary(patches + 1, i) =  (distances - minError)/ obj.tolerance;
             end
             validIndToNodeMap(idx) = i;
        end


        if (obj.rd == 2)
            for i = 1:nNodes-1
                pairwise(i, i+1) = 1;
                pairwise(i+1, i) = 1;
            end
        else            
            [X, Y] = meshgrid(-nSize:nSize,-nSize:nSize);
            for idx = 1:nNodes
                [i, j]= ind2sub(size(obj.Buff),  validCells(idx));
                p = [X(:)+ i, Y(:)+ j];
                p = max(p, 0);
                aboveMax = find(p  > obj.numCells);
                p(aboveMax) = 0; 
                s = p(:,1).*p(:,2);

                f = find(s);
                neigh = sub2ind(size(obj.Buff), X(f)+ i, Y(f)+ j);
                validN = neigh(ismember(neigh, validCells));
                pairwise(idx, validIndToNodeMap(validN)) = 1;
            end
        end 

        labelcost = 1 - eye(nLabels,nLabels);
        unary = real(unary);

        labels = GCMex(segclass, single(unary), pairwise, single(labelcost),0);

        chosenLabels = unique(labels);

        obj.nPatches = size(chosenLabels,1);
        iPatch = 1;
        for el = chosenLabels'
            cells = find(labels == el);
            for idx= cells'
                cell = validCells(idx);
                obj.Buff(cell).bestPatch = iPatch;
                foundSolutionWithAssignedPath = false;
                for solution = obj.Buff(cell).solutions
                    if (solution.patch == el)  
                        obj.selectedPatches(cell) = el;
                        if(foundSolutionWithAssignedPath)
                            if(obj.currentDistance(cell) > solution.dist)
                                obj.Buff(cell).minD = solution.D;
                                obj.Buff(cell).minF = solution.F;
                                obj.currentDistance(cell) = solution.dist;
                            end
                        else
                            obj.Buff(cell).minD = solution.D;
                            obj.Buff(cell).minF = solution.F;
                            obj.currentDistance(cell) = solution.dist;
                            foundSolutionWithAssignedPath = true;
                        end
                    end
                end
                if(~foundSolutionWithAssignedPath)
                   obj.Buff(cell).bestPatch = 0; 
                end
            end
            iPatch = iPatch+1;
        end
    end
end 

end
