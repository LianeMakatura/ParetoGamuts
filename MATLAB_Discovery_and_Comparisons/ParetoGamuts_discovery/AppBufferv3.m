classdef AppBufferv3 < handle
properties
    % === static variables
    % problem specific (fixedparams)
    mFunc
    rD
    rd
    rA
    varNames
    
    % optimization specific (userparams)
    numCells 
    directionDelta
    numRepeatedNoImprovement
    improvementTol
    fillHoles_gamutExtraction
    gamutFillTol
    nMaxRuns
    
    % === dynamic variables
    patchArray
    cylCellDims
    cylPerfArray
    cellCenters
    
    % termination criteria
    HV_history
    converged
end

methods
    % Constructor
    function obj = AppBufferv3(mFunc, userParams, fixedParams)
        if nargin==0 % needed for Object arrays
            return
        end
        
        % dimensionality of design / perf / application space
        obj.mFunc = mFunc;
        obj.rD = fixedParams.rD;
        obj.rd = fixedParams.rd; 
        obj.rA = fixedParams.rA;
        obj.varNames = fixedParams.varNames;
        
        obj.numCells = userParams.numCells;
        obj.directionDelta = userParams.directionDelta;

        obj.patchArray = {};

        numRadialPerfCells = obj.numCells^(obj.rd - 1);
        numContextCells = obj.numCells^(obj.rA);
        obj.cylCellDims = repmat({obj.numCells}, obj.rd -1 + obj.rA, 1);
        
        obj.cylPerfArray = cell([numRadialPerfCells*numContextCells, 1]); % cylindrical performance array; radial for d-1 perf dimensions, cartesian for A application dims

        obj.cellCenters = obj.getCellCenters();
        
        % termination metrics
        obj.converged = false;
        obj.HV_history = zeros([1, numContextCells]); % initial HV is 0, no samples
        obj.nMaxRuns = userParams.nMaxRuns;
        obj.numRepeatedNoImprovement = userParams.numRepeatedNoImprovement;
        obj.improvementTol = userParams.improvementTol;
        obj.fillHoles_gamutExtraction = userParams.fillHoles_gamutExtraction;
        obj.gamutFillTol = userParams.gamutExtractCellPercent; % how much perturbation around the center, as a fraction of cell size
    end
    
    %% add patch to the buffer and interval structures
    function [] = addPatch(obj, patchID, centerPt, explorationDirections, ... 
                            extents, validPatchPts, connectivityList, fillHoles)
        desPts = validPatchPts(:, 1:obj.rD);
        appPts = validPatchPts(:, obj.rD+1:obj.rD+obj.rA);
        perfPts = validPatchPts(:, obj.rD+obj.rA+1:end);

        % convert to (normalized by pi/2) spherical coords (alpha space) and get span on first
        % rd-1 elements
        numCylDims = length(obj.cylCellDims);
        numPts = size(validPatchPts, 1);
        
        % put the points in the buffer
        if fillHoles && (numCylDims == 2 || numCylDims == 3) && ~isempty(connectivityList)
            % fill in all the holes in the patch 
            radii = zeros(numPts, 1);
            cylPatchPts = zeros(numPts, numCylDims);
            for ptID=1:numPts
                [phi, r] = cartesian2polar(perfPts(ptID, :));
                phi = phi / (pi / 2);                       % normalize all components to [0,1]

                % collect the points for post processing only
                cylPatchPts(ptID, :) = [phi, appPts(ptID, :)];
                radii(ptID, :) = r;
            end
            
             % fill in all cells within the patch expanse, and save the
             % augmented patch points 
            [validPatchPts] = obj.placePatchInCylBuffer(patchID, cylPatchPts, radii, validPatchPts, connectivityList);
            
        else
            % only put in the specific points we found (don't yet know how
            % to fill in the holes in higher dimensions)
            for ptID=1:numPts
                [phi, r] = cartesian2polar(perfPts(ptID, :));
                phi = phi / (pi / 2);                       % normalize all components to [0,1]

                % place point into radialPerfBuffer
                obj.placePtInCylBuffer(patchID, ptID, phi, r, appPts(ptID, :));
            end
        end
        
        % patch array keeps all patch info and points in cartesian coords
        newPatch = struct( ...
            'patchID', patchID, ...
            'centerPt', centerPt, ...
            'expDirs', explorationDirections, ...
            'dirExtents', extents, ...
            'patchPts', validPatchPts ...
        );
        obj.patchArray{patchID} = newPatch; 
    end
    
    
    
    % @param cylPts - nx(rd-1+A) array of cylindrical points [phi, app]
    % @param radii - nx1 array of radii
    % @return cartestianPatchPts - mx(rD+rA+rd) array of points in the patch
    function [cartesianPatchPts] = placePatchInCylBuffer(obj, patchIdx, cylPatchPts, radii, cartesianPatchPtsin, patchConnectivity)
        keepOrigPts = true;
        
        % get triangulation/tetrahedralization connectivity
%         T = delaunayTriangulation(cylPatchPts);
        T = triangulation(patchConnectivity, cylPatchPts);
        
        % test all cell centers for their containing cell & barycentric 
        numCylDims = length(obj.cylCellDims);
        if numCylDims == 2
            % faster version
            [triIdx,baryCoords] = pointLocationQuadTree(T,obj.cellCenters);
        else
            [triIdx, baryCoords] = pointLocation(T, obj.cellCenters); % get tri/tet id for each querypt
        end
        cellIDsLin = (1:length(triIdx))';
%         validIdcs = ~isnan(triIdx);
        validIdcs = boolean(~isnan(triIdx) .* ~isnan(baryCoords(:,1)));
        triIdx = triIdx(validIdcs, :); 
        
        if isempty(triIdx)
            cartesianPatchPts = cartesianPatchPtsin; % keep the original (probably a line in radial space)
            return;
        end
        
        baryCoords = baryCoords(validIdcs, :);
        cellIDsLin = cellIDsLin(validIdcs, :);
        
        % for each point, toss it in to buffer and the patchpts list
        if keepOrigPts
            ptIdx = size(cylPatchPts, 1) + 1;
        else
            ptIdx = 1;
        end
        newCartesianPoints = zeros(length(triIdx), obj.rD+obj.rA+obj.rd);
        for i=1:length(triIdx)
            cellIDlin = cellIDsLin(i);
            p=triIdx(i);
            if isnan(p) %outside the patch
                continue;
            end

            % interpolate the Barycentric coordinates in all spaces
            b = baryCoords(i, :);
            triIDs = T.ConnectivityList(p, :);
            
            tDesPts = cartesianPatchPtsin(triIDs, 1:obj.rD);
            newDesPt = (tDesPts' * b')'; % row with interpolated design point (assume the same coords)
            
            tAppPts = cartesianPatchPtsin(triIDs, obj.rD+1:obj.rD+obj.rA);
            newAppPt = (tAppPts' * b')';
            
%             tPerfPts = cartesianPatchPtsin(triIDs, obj.rD+obj.rA+1:end);
%             newCartPerfPt = (tPerfPts'*b')'; % row with interpolated pt
%             newRadius = radii(triIDs, :)' * b';

            newCartPerfPt = obj.mFunc.eval(newDesPt, newAppPt);
            [~, newRadius] = cartesian2polar(newCartPerfPt);

            newCartesianPoints(i, :) = [newDesPt, newAppPt, newCartPerfPt];
            obj.updateCylBuffer(cellIDlin, patchIdx, ptIdx, newRadius);
            ptIdx = ptIdx + 1;
        end
        if keepOrigPts
            cartesianPatchPts = [cartesianPatchPtsin; newCartesianPoints];
        else
            cartesianPatchPts = newCartesianPoints;
        end
    end
    
    
    
    function [] = placePtInCylBuffer(obj, patchIdx, ptIdx, ptPhi, ptRadius, ptApp)
        cellIDlin = obj.getLinearCellIdx([ptPhi, ptApp], obj.cylCellDims);
        obj.updateCylBuffer(cellIDlin, patchIdx, ptIdx, ptRadius);
    end
    
    function [] = updateCylBuffer(obj, cellIDlin, patchIdx, ptIdx, ptRadius)
        numToKeep = 50;
        % place a reference to this point id in the appropriate buffer cell 
        newEntry = [patchIdx, ptIdx, ptRadius];
        obj.cylPerfArray{cellIDlin} = [obj.cylPerfArray{cellIDlin}; newEntry];
        
        if size(obj.cylPerfArray{cellIDlin}, 1)> numToKeep
            radiusColIdx = 3;                       % each row of pts is [patchID, ptID(in patchArray.patchPts), radius]
            pts = sortrows(obj.cylPerfArray{cellIDlin}, radiusColIdx);      % sorted in fixedcontext cell from best to worst
            obj.cylPerfArray{cellIDlin} = pts(1:numToKeep, :);
        end
    end
    
    
    % Identify the cell (within the interval structure) that a certain point 
    % belongs in
    function [cellID, clamped] = getSubscriptedCellIdx(obj, vals, dimCell)
        rem = floor(vals .* cell2mat(dimCell)');
        if vals < 1         % account for the upper edge of the very last cell (would otherwise go to the next cell up)
            cellID = rem + 1; 
        else
            cellID = rem;
        end
        
        % make sure it's in the valid range, warn if not
        clamped = false;
        if any(cellID > obj.numCells)
%             fprintf("CellID out of bounds. Clamping to max cell.\n");
            cellID = min(cellID, obj.numCells);
            clamped = true;
        end
        if any(cellID < 1)
%             fprintf("CellID out of bounds. Clamping to min cell.\n");
            cellID = max(cellID, 1);
            clamped = true;
        end
    end
    
    
    function [cellIDlin, clamped] = getLinearCellIdx(obj, vals, dimCell)
        [cellIDsub, clamped] = obj.getSubscriptedCellIdx(vals, dimCell);
        cellIDlin = obj.getLinIDfromSubID(cellIDsub, dimCell);
    end
    
    function [cellIDlin] = getLinIDfromSubID(obj, cellIDsub, dimCell)
        if length(cellIDsub) ~= length(dimCell)
            disp('Mismatch in indexing function!');
            keyboard;
        end
        if length(dimCell) > 1
            sz = cell2mat(dimCell);
            cellIDsub = num2cell(cellIDsub);        % have to input each dim as separate arg to sub2ind, only possible with cell array
            cellIDlin = sub2ind(sz, cellIDsub{:});
        else
            cellIDlin = cellIDsub;
        end
    end
    
    function [cellIDsub] = getSubIDfromLinID(obj, cellIDlin, dimCell)
        if length(dimCell) > 1
            sz = cell2mat(dimCell);
            cellIDsub = ind2sub(sz, cellIDlin);
        else
            cellIDsub = cellIDlin;
        end
    end
    
    function [] = addPoint(obj, point) % decided not to do this
    end
    
    function [randomPoints] = getPointsToImprove(obj, numSamples)
        % for now, just pick some random point in the buffer
        % and draw from a distribution around it
        
        numPatches = length(obj.patchArray);
        randomPoints = zeros(numSamples, obj.rD+obj.rA);
        
        i=1;
        while i < numSamples+1
            patchID = ceil(rand(1,1) * numPatches);
            
            % if we didn't want to store the points, just walk along the 
            % expansion plane to random point and draw from distribution
            % around that
            pts = obj.patchArray{patchID}.patchPts;
            ptID = ceil(rand(1,1) * size(pts, 1));
            pt = pts(ptID, :);
            pt = pt(1:obj.rD+obj.rA); % get only design and app var vals
            
            sigma = 1./2.^(10*rand(1, 1));
%             randPt = pt - sigma*(2*rand(size(pt))-1); 
            randPt = pt - sigma*(2*rand(size(pt))-1); 
            
            if obj.mFunc.hasLinearConstraints %check if it violates constraints
                [numViolated, ~] = obj.mFunc.p.checkLinearConstraints(randPt(1, 1:obj.rD), randPt(obj.rD+1:obj.rD+obj.rA));
                if numViolated == 0 %valid sample
                    randomPoints(i, :) = randPt;
                    i = i+1;
                end
            else
                randomPoints(i,:) = randPt; %always accept
                i = i+1;
            end
        end
        
        randomPoints = min(randomPoints, 1);
        randomPoints = max(randomPoints, 0.00001);   
    end
    
    
    % get target points and directions for scalarization within single
    % contextual configuration (app values not affected)
    % @param samples - sample to optimize in PERFORMANCE space (only rd params; no app var values)
    function [targets, dirs] = getRandomDirections(obj, samples)
        dirDelta = obj.directionDelta;
        numSamples = size(samples, 1);
        targets = zeros(size(samples));
        dirs = zeros(size(samples));
        
        for i=1:numSamples
            % identify a goal point 
            origin = samples(i,:);
            distToSample = norm(origin);
            
            % TODO: add randomness by taking direction from neighboring
            % area 
            goal = obj.getGoalPt(origin);
            dir = goal - origin/distToSample;
            dir = dir/norm(dir);

            amountToPush  = dirDelta*distToSample;

            targets(i,:) = samples(i,:) + amountToPush*dir;     
            dirs(i, :) = dir;
        end
    end
    
    % internal function only (used by getRandomDirections)
    function [goal] = getGoalPt(obj, sample)
        % assume pt given in Cartesian coords
        % get direction toward front in this area (see Fig 5 of Schulz 2018)
        pt = sample / sum(sample); % sums to 1, gives ~"alpha" / perf tradeoff
        goal = pt - ( ones(size(sample))/obj.rd ); % offset so the (0.5)^rd vector gives 0^rd
        goal = 2 * goal; 
    end
    
    function [centers] = getCellCenters(obj)
        N = length(obj.cylCellDims);
        nCells = size(obj.cylPerfArray, 1);
        centers = zeros([nCells, N]);
        
        cylDims = [obj.cylCellDims{:}];
        
        cellCenterCoords = cell(1, N);
        for ptidx = 1:nCells
            % get the subscripted indices for the fixed context
           [cellCenterCoords{:}] = ind2sub(cylDims, ptidx);
           centerPt = ([cellCenterCoords{:}] - 0.5) ./ cylDims;
           centers(ptidx, :) = centerPt;
        end
    end
    
    
    function alreadyExplored = checkIfAlreadyExpanded(obj, centerPt, expDirs)
        alreadyExplored = false;
        singValTol = 0.0001;
        
        for patchID = 1:length(obj.patchArray)
            existingDirs = obj.patchArray{patchID}.expDirs; % each col an [x', z'] direction
            
            % dirs of length rd+rA, so this is max possible rank; if this 
            % singular value is non-zero, dirs linearly independent. 
            maximalRank = obj.rd+obj.rA; 
            
            dirMatrix = [existingDirs expDirs];
            [~, S, ~] = svd(dirMatrix);
            singVal = S(maximalRank, maximalRank); 
            if singVal < singValTol
                % lin dep dirs --> planes parallel; check if identical by
                % checking whether direction between existing centerpt to
                % newcenterpt is contained in old expansion plane.
                
                existingCenter = obj.patchArray{patchID}.centerPt(1:obj.rD+obj.rA);
                center = centerPt(1:obj.rD+obj.rA);
                centerDisplacement = existingCenter - center;
                dirPointMatrix = [existingDirs centerDisplacement'];

                [~, S, ~] = svd(dirPointMatrix);
                singVal2 = S(maximalRank, maximalRank);
                if singVal2 < singValTol
                    alreadyExplored = true;
                    return;
                end
            end
        end
    end
    
    
    function [stop, converged, percentTargStability] = terminate(obj, iterNum, useHV)      
        if useHV
            approximate = false;
        else
            approximate = true;
        end
        % compute Gamut hypervolume and iterations without improvement
        HVs = obj.getGamutHypervolume(approximate);
        obj.HV_history = [obj.HV_history; HVs];
        numStableIters = obj.countItersNoImprovement(obj.improvementTol);
        
        % check if we've converged or hit limits
        percentTargStability = numStableIters / obj.numRepeatedNoImprovement; % used to decide when to randomly sample again
        if numStableIters >= obj.numRepeatedNoImprovement
            converged = true;
            stop = true;
        elseif iterNum > obj.nMaxRuns
            stop = true;
            converged = false;
        else
            stop = false; 
            converged = false;
        end
        obj.converged = converged;
    end
    
    function empty = isEmpty(obj)
        if isempty(obj.patchArray)
            empty = true;
        else
            empty = false;
        end
    end
    
    %% ==================== TERMINATION CRITERIA ===============================
    function [numItersNoImprovement] = countItersNoImprovement(obj, tol)
        numItersNoImprovement = 0;
        numMeasurements = size(obj.HV_history, 1);
        if numMeasurements <2 %not enough data points yet
            return
        end
        
        % define the improvement function over the gamut
        aggFunc = @(curr, prev)((curr - prev)*(curr - prev)'); %sum of squared differences for individual FC hypervolumes
        
        % test improvement over previous iterations; stops as soon as
        % there's progress
        currHV = obj.HV_history(numMeasurements, :);
        for prevIdx=(numMeasurements-1):-1:1
            improvement = aggFunc(currHV, obj.HV_history(prevIdx, :));
            if improvement > tol
                break; %sufficiently substantial progress 
            else
                numItersNoImprovement = numItersNoImprovement+1;
            end
        end
        
    end
    
    
    % compute hypervolume of gamut
    function [fcHyperVols] = getGamutHypervolume(obj, approximate)
%         contextDims = obj.cylCellDims{obj.rd:obj.rd+obj.rA}; %pick out the contextual dimensions
        contextDims = [];
        for i = obj.rd-1+1:obj.rd-1+obj.rA
            contextDims = [contextDims, obj.cylCellDims{i}];
        end
%         contextDims = obj.cylCellDims{obj.rd-1+1:obj.rd-1+obj.rA}; %pick out the contextual dimensions
        numContexts = prod(contextDims);
        
        % get the hv for each slice
        fcHyperVols = zeros(numContexts, 1);
        contextCell = cell([1, obj.rA]);
        for fcIDlinear=1:numContexts
            % get the subscripted indices for the contextual values
           [contextCell{:}] = ind2sub(contextDims, fcIDlinear);
           fcHyperVols(fcIDlinear) = obj.getFixedContextHypervolume(contextCell, approximate);
        end
        fcHyperVols = fcHyperVols';
    end
    
    function hv = getFixedContextHypervolume(obj, contextCellSub, approximate)
        nadirPt = ones(obj.rd, 1); % upper bound for all objectives
        
        FC_cells = obj.getFixedContextCellIDs(contextCellSub);
        
        % gather the non-dominated (Pareto optimal) points in this fixed context
        F = [];
        for i=1:length(FC_cells)
            pts = obj.cylPerfArray{FC_cells(i)}; 
            if isempty(pts)
                continue;
            end
            
            radiusColIdx = 3;                       % each row of pts is [patchID, ptID(in patchArray.patchPts), radius]
            [m, minidx] = min(pts(:, radiusColIdx));           % find index of point with minimum radius in this cell 
            minPtRef = pts(minidx, :);
            
            % get the actual [Cartesian] point from patch Arrays
            patchID = minPtRef(1); 
            ptID =minPtRef(2);
            p = obj.patchArray{patchID}.patchPts(ptID, :);    
            
            F = [F; p];
        end
        
        if isempty(F)
            hv = 0;
            return;
        end
        
        % must check for dominated points in each context
        % pass performance points in, to get non-dominate points
        nonDominated = getParetoIndices(F(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
        paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices
        F = F(paretoIndices, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd);
        
        if approximate
            hv = approximate_hypervolume_ms(F', nadirPt);
        else
            hv = lebesgue_measure(F', nadirPt); % transpose F because need each point's performance in a column
        end
    end
    
    
    
    %% ==================== IDENTIFY PARETO OPTIMALITY =========================    
    function [finalPts, finalLabels, filledPatchPts] = getParetoGamut(obj)
        numPtsToKeep = 1;
        candPts = []; candLabels = []; candContexts = [];
        
        % loop over a particular cell in the cylindrical performance space
        for linearCellID = 1:length(obj.cylPerfArray)
            % sort the points by radius, so pts(1) is closest to origin in cell, and points(end) is furthest
            % radius is only wrt it's fixed-context radial array
            pts = obj.cylPerfArray{linearCellID}; 
            if isempty(pts)
                continue;
            end
            
            radiusColIdx = 3;                       % each row of pts is [patchID, ptID(in patchArray.patchPts), radius]
            pts = sortrows(pts, radiusColIdx);      % sorted in fixedcontext cell from best to worst
            
            % extract best N points in this cell of the cylindrical perf
            % (alpha, app) buffer, and get the actual [Cartesian] point from patch Arrays
            N = min(size(pts, 1), numPtsToKeep);
            for i=1:N
                patchID = pts(i,1); 
                ptID =pts(i, 2);
                p = obj.patchArray{patchID}.patchPts(ptID, :);
                
                % figure out the appropriate contextual cell for this point
                if obj.rA > 0
                    appVals = p(obj.rD+1:obj.rD+obj.rA);
                    appCellID = obj.getLinearCellIdx(appVals, obj.cylCellDims(obj.rd-1+1:end));
                else
                    appCellID = 1;
                end
                
                candPts = [candPts; p];
                candLabels = [candLabels; patchID];
                candContexts = [candContexts; appCellID];
            end
        end
        
        % must check for dominated points in each context
        finalPts = []; finalLabels = [];
        if obj.rA == 0
            numContexts = 1;
        else
            context_dims = [];
            for dim = obj.rd-1+1:length(obj.cylCellDims)
                context_dims = [context_dims, obj.cylCellDims{dim}];
            end
            numContexts = prod(context_dims);
        end
        for i=1:numContexts
            singlepts_idx = find(candContexts == i);
            singlePts = candPts(singlepts_idx, :);
            singleLabels = candLabels(singlepts_idx, :);
            if isempty(singlePts)
                continue;
            end
            
            % pass performance points in, to get non-dominate points
            nonDominated = getParetoIndices(singlePts(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
            paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices

            finalPts = [finalPts; singlePts(paretoIndices, :)];
            finalLabels = [finalLabels; singleLabels(paretoIndices, :)];
        end
        
        if obj.fillHoles_gamutExtraction %&& obj.rd -1 + obj.rA== 2
            lastParetoCheck = true; % should be true to prevent adding in dominated points

            uniqueLabels = unique(finalLabels);
            numPatchIDs = length(uniqueLabels);

            mappedLabels = zeros(size(finalLabels));
            for i=1:numPatchIDs
                ids = find(finalLabels == uniqueLabels(i));
                mappedLabels(ids) = i;
            end
            
            pts = finalPts;
            % get all samples within the desired boundary of the
            % points, and test for pareto optimality within all of
            % these patch points (in case there are independent pieces of a contiguous patch piece eg ZDT3)
            filledPatchPts = cell(numPatchIDs, 1);
            for i=1:numPatchIDs
                origPatchNum = uniqueLabels(i); %get the original patch ID of this point
                patchIds = find(mappedLabels == i);
                dirs = obj.patchArray{origPatchNum}.expDirs'; % subspace directions, transposed so each dir on row
                augDesPts = pts(patchIds, 1:obj.rD+obj.rA);
                if size(augDesPts, 1) < 3
                    continue;
                end
                projPatchPts = proj_patch2subspace(augDesPts, dirs); %project to expansion subspace to get 2 or 3d space for triangulation
                allPts = obj.patchArray{origPatchNum}.patchPts;
                projAllPts = proj_patch2subspace(allPts(:, 1:obj.rD+obj.rA), dirs);

                % find the boundary elements of the patch / volume that
                % was deemed optimal
                boundaryIdx = boundary(projPatchPts, 1);

                % test all points for whether they're in this boundary
                polyx = projPatchPts(boundaryIdx, 1);
                polyy = projPatchPts(boundaryIdx, 2);
                queryx = projAllPts(:, 1);
                queryy = projAllPts(:, 2);
                [inPoly] = inpolygon(queryx,queryy,polyx,polyy); %will only work for 2d; gives logical array
                allPts = allPts(inPoly, :);

                % get Pareto optimal points out of this set
                finalFullPatch = [];
                appDims = cell(1, obj.rA);
                [appDims{:}] = obj.cylCellDims{obj.rd-1+1:obj.rd-1+obj.rA};
                appDims = [appDims{:}];
                numContexts = prod(appDims);
                cellSize = 1./appDims;
                cellIdx = cell(1, obj.rA); 
                for c=1:numContexts % must check pareto optimality by context
                    % get center of context, and find all points within
                    % cell
                    if length(appDims) == 1
                        center = (c - 0.5) * cellSize;
                        % find all points near that context
                        FCidx = find( (abs(allPts(:, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*obj.gamutFillTol ) < 0 );
                        FCpoints = allPts(FCidx, :);
                    else
                       [cellIdx{:}] = ind2sub(appDims, c);
                       center = ([cellIdx{:}] - 0.5) .* cellSize;

                       FCpoints = [];
                       for pInd=1:size(allPts, 1)
                           keep = max(abs(allPts(pInd, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*obj.gamutFillTol ) < 0;
                           if keep
                               FCpoints = [FCpoints; allPts(pInd, :)];
                           end
                       end
                    end

                    % compute pareto dominance in fixed context
                    if ~isempty(FCpoints)
                        nonDominated = getParetoIndices(FCpoints(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
                        paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices
                        finalFullPatch = [finalFullPatch; FCpoints(paretoIndices, :)];
                    end
                end
                filledPatchPts{i} = finalFullPatch;
            end
            
            % update the finalPts and finalLabels
            finalPts = [];
            finalLabels = [];
            for i=1:numPatchIDs
                finalPts = [finalPts; filledPatchPts{i}];
                numPts = size(filledPatchPts{i}, 1);
                finalLabels = [finalLabels; uniqueLabels(i)*ones([numPts, 1])];
            end

            if lastParetoCheck
                % ============================================================
                % check pareto dominance one more time (patches may be
                % internally consistent but not with others)
                finalParetoPoints = [];
                finalParetoLabels = [];
                appDims = cell(1, obj.rA);
                [appDims{:}] = obj.cylCellDims{obj.rd-1+1:obj.rd-1+obj.rA};
                appDims = [appDims{:}];
                numContexts = prod(appDims);
                cellSize = 1./appDims;
                cellIdx = cell(1, obj.rA); 

                for c=1:numContexts % must check pareto optimality by context
                    % get center of context, and find all points within
                    % cell
                    if length(appDims) == 1
                        center = (c - 0.5) * cellSize;
                        % find all points near that context
                        FINALidx = find( (abs(finalPts(:, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*obj.gamutFillTol ) < 0 );
                        FINALpoints = finalPts(FINALidx, :);
                        FINALlabels = finalLabels(FINALidx, :);
                    else
                       [cellIdx{:}] = ind2sub(appDims, c);
                       center = ([cellIdx{:}] - 0.5) .* cellSize;

                       FINALpoints = [];
                       FINALlabels = [];
                       for pInd=1:size(finalPts, 1)
                           keep = max(abs(finalPts(pInd, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*obj.gamutFillTol ) < 0;
                           if keep
                               FINALpoints = [FINALpoints; finalPts(pInd, :)];
                               FINALlabels = [FINALlabels; finalLabels(pInd, :)];
                           end
                       end
                    end

                    % compute pareto dominance in fixed context
                    if ~isempty(FINALpoints)
                        nonDominated = getParetoIndices(FINALpoints(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
                        paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices
                        finalParetoPoints = [finalParetoPoints; FINALpoints(paretoIndices, :)];
                        finalParetoLabels = [finalParetoLabels; FINALlabels(paretoIndices, :)];
                    end
                end
                finalPts = finalParetoPoints;
                finalLabels = finalParetoLabels;
                % ======================================
            end

        else
            filledPatchPts = {};
        end
        
    end

    % context cell is already a subscripted index
    function [linearFCids] = getFixedContextCellIDs(obj, contextCell)
        % return all cell ids that are associated with this fixed context
        % loop over all
        radialDims=cell(obj.rd-1, 1);
        [radialDims{:}] = obj.cylCellDims{1:obj.rd-1}; %pick out the alpha dimensions 
        radialDims = [radialDims{:}];
        numFCpoints = prod(radialDims);
        
        linearFCids = zeros([1, numFCpoints]);
        radialCoords = cell(1, obj.rd-1);
        for ptidx = 1:numFCpoints
            % get the subscripted indices for the fixed context
           [radialCoords{:}] = ind2sub(radialDims, ptidx);
           
           % get the linear index within the full cylindrical perf array
           cylDims = cell2mat(obj.cylCellDims);
           linearFCids(ptidx) = sub2ind(cylDims, radialCoords{:}, contextCell{:});
        end
    end

    
    
    %% ================== VISUALIZATION ==========================
    function [design, perf] = visualizeAppBuff(obj, idxToVis, groundTruthOnly, providedPts, providedLabels)
        plotPoints = true;
        plotSurfaces = false;
        plotTriangular = false; % !!!! can only do this if rd-1+rA == 2 or 3! 
        fillHoles = false; % !!!! can only do this if rd-1+rA == 2 or 3! 
        
        colorByContext = true; % if false, color by patch
        monochrome = false;
        contextColorRange = [-0.3, 1.1];
        
        numPatches = length(obj.patchArray);
        filledPatchPts = {};
        
        if nargin < 2 % only obj 
            d1=1; d2=2; a1=1; p1=1; p2=2; % default values
            groundTruthOnly = false;
        else
            d1=idxToVis.designIndex1;
            d2=idxToVis.designIndex2; 
            a1=idxToVis.appIndex;
            p1=idxToVis.perfIndex1;
            p2=idxToVis.perfIndex2;
        end
        
        if d1 < 1 || d1 > obj.rD || d2 < 1 || d2 > obj.rD
            fprintf('Design indices must be between 1 and %d', obj.rD);
            return
        elseif p1 < 1 || p1 > obj.rd || p2 < 1 || p2 > obj.rd
            fprintf('Performance indices must be between 1 and %d', obj.rd);
            return
        elseif obj.rA > 0 && (a1 < 1 || a1 > obj.rA)
            fprintf('Contextual indices must be between 1 and %d', obj.rA);
            return
        end
        
        % ===== pick indices to visualize
        des1_idx = d1;
        des2_idx = d2;
        app_idx = obj.rD + a1;
        perf1_idx = obj.rD+obj.rA+p1;
        perf2_idx = obj.rD+obj.rA+p2;
        
            
        % ==== set up figures
        perf = figure('Name', 'Performance Space'); 
        hold on;    
        % always keep the same proportions
        xlim([0,1]); xlabel(obj.varNames(perf1_idx))
        ylim([0,1]); ylabel(obj.varNames(perf2_idx))
        if obj.rA > 0
            zlim([0,1]); zlabel(obj.varNames(app_idx))
        end
        axis equal;
        set(gca,'FontSize', 14);

        design = figure('Name', 'Design Space'); 
        hold on;
        % always keep the same proportions
        xlim([0,1]); xlabel(obj.varNames(des1_idx))
        ylim([0,1]); ylabel(obj.varNames(des2_idx))
        if obj.rA > 0
            zlim([0,1]); zlabel(obj.varNames(app_idx))
        end
        axis equal; 
        set(gca,'FontSize', 14);

        if exist('providedPts', 'var') && ~groundTruthOnly
            uniqueLabels = unique(providedLabels);
            numPatchIDs = length(uniqueLabels);

            mappedLabels = zeros(size(providedLabels));
            for i=1:numPatchIDs
                ids = find(providedLabels == uniqueLabels(i));
                mappedLabels(ids) = i;
            end
            
            pts = providedPts;
            color = hsv(numPatchIDs);
            markerColors = color(mappedLabels, :); % assign correct color to each point
            
            markerSize = 25;
            
            if plotTriangular
                % get the connectivity for each patch
                patchConnectivity = cell(numPatches, 1);
                for i=1:numPatchIDs
                    origPatchNum = uniqueLabels(i); %get the original patch ID of this point
                    patchIds = find(mappedLabels == i);
                    dirs = obj.patchArray{origPatchNum}.expDirs'; % subspace directions, transposed so each dir on row
                    augDesPts = pts(patchIds, 1:obj.rD+obj.rA);
                    if size(augDesPts, 1) < 3
                        continue;
                    end
                    projPatchPts = proj_patch2subspace(augDesPts, dirs); %project to expansion subspace to get 2 or 3d space for triangulation

                    % find the boundary elements of the patch / volume 
                    boundaryIdx = boundary(projPatchPts, 1);
                    T = delaunay(projPatchPts);
                    
                    %remove all-edge triangles
                    keep = logical(ones(size(T, 1), 1));
                    for triID=1:size(T, 1)
                        tri = T(triID, :);
                        if (ismember(tri(1), boundaryIdx) && ismember(tri(2), boundaryIdx) && ismember(tri(3), boundaryIdx))
                            keep(triID) = 0;
                        end
                    end
                    T = T(keep, :);
                    patchConnectivity{i} = T;
                end
            end
            
            if fillHoles
                % get all samples within the desired boundary of the
                % points, and test for pareto optimality within all of
                % these patch points (in case there are independent pieces of a contiguous patch piece eg ZDT3)
                filledPatchPts = cell(numPatches, 1);
                for i=1:numPatchIDs
                    origPatchNum = uniqueLabels(i); %get the original patch ID of this point
                    patchIds = find(mappedLabels == i);
                    dirs = obj.patchArray{origPatchNum}.expDirs'; % subspace directions, transposed so each dir on row
                    augDesPts = pts(patchIds, 1:obj.rD+obj.rA);
                    if size(augDesPts, 1) < 3
                        continue;
                    end
                    projPatchPts = proj_patch2subspace(augDesPts, dirs); %project to expansion subspace to get 2 or 3d space for triangulation
                    allPts = obj.patchArray{origPatchNum}.patchPts;
                    projAllPts = proj_patch2subspace(allPts(:, 1:obj.rD+obj.rA), dirs);
                    
                    % find the boundary elements of the patch / volume that
                    % was deemed optimal
                    boundaryIdx = boundary(projPatchPts, 1);
                    
                    % test all points for whether they're in this boundary
                    polyx = projPatchPts(boundaryIdx, 1);
                    polyy = projPatchPts(boundaryIdx, 2);
                    queryx = projAllPts(:, 1);
                    queryy = projAllPts(:, 2);
                    [inPoly] = inpolygon(queryx,queryy,polyx,polyy); %will only work for 2d; gives logical array
                    allPts = allPts(inPoly, :);
                    
                    % get Pareto optimal points out of this set
                    finalFullPatch = [];
                    appDims = obj.cylCellDims{obj.rd:obj.rd+obj.rA};
                    numContexts = prod(appDims);
                    cellSize = 1./appDims;
                    for c=1:numContexts % must check pareto optimality by context
                        % get center of context, and find all points within
                        % cell
                        if length(appDims) == 1
                            center = (c - 0.5) * cellSize;
                            % find all points near that context
                            FCidx = find( (abs(allPts(:, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*0.5 ) < 0 );
                            FCpoints = allPts(FCidx, :);
                        else
                           [cellIdx{:}] = ind2sub(appDims, c);
                           center = (cellIdx - 0.5) * cellSize;
                           
                           FCpoints = [];
                           for pInd=1:size(allPts, 1)
                               keep = max(abs(allPts(pInd, obj.rD+1:obj.rD+obj.rA) - center) - cellSize*0.5 ) < 0;
                               if keep
                                   FCpoints = [FCpoints; allPts(pInd, :)];
                               end
                           end
                        end
                        
                        % compute pareto dominance in fixed context
                        if ~isempty(FCpoints)
                            nonDominated = getParetoIndices(FCpoints(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
                            paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices
                            finalFullPatch = [finalFullPatch; FCpoints(paretoIndices, :)];
                        end
                    end
                    filledPatchPts{i} = finalFullPatch;
                end
            end
            
            
            figure(perf); % create or get performance space fig
            if obj.rA > 0
                % plot the observed data
                if plotPoints
                    if colorByContext
                        plottingColors = pts(:, app_idx); % color by context value
                        caxis(contextColorRange);
                        if monochrome
                            colormap gray;
                        else
                            colormap jet;
                        end
                    else
                        plottingColors = markerColors;
                    end
                    scatter3(pts(:, perf1_idx), pts(:,perf2_idx), pts(:,app_idx),...
                        markerSize, plottingColors);
                end
                if fillHoles
                    for i=1:numPatchIDs
                        allpts = filledPatchPts{i};
                        if isempty(allpts)
                            continue;
                        end
                        if colorByContext
                            allptcolors = allpts(:, app_idx); % color by context value
                            caxis(contextColorRange);
                            if monochrome
                                colormap gray;
                            else
                                colormap jet;
                            end
                        else
                            allptcolors = color(i, :); % color by family; pick out the hue
                        end
                        
                        scatter3(allpts(:, perf1_idx), allpts(:,perf2_idx), allpts(:,app_idx), ...
                            markerSize, allptcolors);
                    end
                end
                
                % plot the interpolated surfaces for each
                if plotSurfaces
                    for i=1:numPatchIDs
                        patchIds = find(mappedLabels == i);
                        patchPts = pts(patchIds, [perf1_idx, perf2_idx, app_idx]);

                        [XX, YY, ZZ] = obj.interpolatePatch(patchPts);
                        if ~isempty(ZZ)
                            c = color(i, :);
                            cgrid = zeros(size(ZZ,1), size(ZZ,2), 3);
                            cgrid(:,:,1) = c(1); cgrid(:,:,2) = c(2); cgrid(:,:,3) = c(3);
                            surf(XX,YY,ZZ, cgrid);
                        end
                    end
                end
                if plotTriangular
                    for i=1:numPatchIDs
                        patchIds = find(mappedLabels == i);
                        patchPts = pts(patchIds, [perf1_idx, perf2_idx, app_idx]);
                        if size(patchPts, 1) < 3
                            continue;
                        end
                        
                        T = patchConnectivity{i};
                        t = trisurf(T,patchPts(:,1), patchPts(:,2), patchPts(:,3));
                    end
                end
                
            else
                scatter(pts(:, perf1_idx), pts(:,perf2_idx),...
                    markerSize, markerColors);
                % plot the interpolated curves for each
%                 for i=1:numPatchIDs
%                     patchIds = find(mappedLabels == i);
%                     patchPts = pts(patchIds, [perf1_idx, perf2_idx]);
%                     
%                     plot(patchPts(:, 1), patchPts(:, 2))
%                 end
            end

            figure(design); % create or get design space fig
            if obj.rA > 0
                if plotPoints
                    if colorByContext
                        plottingColors = pts(:, app_idx); % color by context value
                        caxis(contextColorRange);
                        if monochrome
                            colormap gray;
                        else
                            colormap jet;
                        end
                    else
                        plottingColors = markerColors;
                    end
                    scatter3(pts(:, des1_idx), pts(:,des2_idx), pts(:,app_idx),...
                        markerSize, plottingColors);
                end
                if fillHoles
                    for i=1:numPatchIDs
                        allpts = filledPatchPts{i};
                        if isempty(allpts)
                            continue;
                        end
                        if colorByContext
                            allptcolors = allpts(:, app_idx); % color by context value
                            caxis(contextColorRange);
                            if monochrome
                                colormap gray;
                            else
                                colormap jet;
                            end
                        else
                            allptcolors = color(i, :); % color by family; pick out the hue
                        end
                        scatter3(allpts(:, des1_idx), allpts(:,des2_idx), allpts(:,app_idx), ...
                            markerSize, allptcolors);
                    end
                end
                
                % plot the interpolated surfaces for each
                if plotSurfaces
                    for i=1:numPatchIDs
                        patchIds = find(mappedLabels == i);
                        patchPts = pts(patchIds, [des1_idx, des2_idx, app_idx]);

                        [XX, YY, ZZ] = obj.interpolatePatch(patchPts);
                        if ~isempty(ZZ)
                            c = color(i, :);
                            cgrid = zeros(size(ZZ,1), size(ZZ,2), 3);
                            cgrid(:,:,1) = c(1); cgrid(:,:,2) = c(2); cgrid(:,:,3) = c(3);
                            surf(XX,YY,ZZ, cgrid);
                        end
                    end
                end
                
            else
                scatter(pts(:, des1_idx), pts(:,des2_idx),...
                    markerSize, markerColors);
                
                % plot the interpolated curves for each
%                 for i=1:numPatchIDs
%                     patchIds = find(mappedLabels == i);
%                     patchPts = pts(patchIds, [des1_idx, des2_idx]);
%                     
%                     plot(patchPts(:, 1), patchPts(:, 2))
%                 end
            end
                
        elseif ~exist('providedPts', 'var') && ~groundTruthOnly
            color=hsv(numPatches);
            for i=1:numPatches
                patch = obj.patchArray{i};
                pts = patch.patchPts;
                cpt = patch.centerPt;

                % plot the patch points (design, perf)
                figure(perf); % create or get performance space fig
                if obj.rA > 0
                    scatter3(pts(:, perf1_idx), pts(:,perf2_idx), pts(:,app_idx),...
                        'MarkerFaceColor', color(i,:), 'MarkerEdgeColor', color(i, :));
                    scatter3(cpt(perf1_idx), cpt(perf2_idx), cpt(app_idx),...
                        250, 'black', 'filled');
                else
                    scatter(pts(:, perf1_idx), pts(:,perf2_idx),...
                        'MarkerFaceColor', color(i,:), 'MarkerEdgeColor', color(i, :));
                    scatter(cpt(perf1_idx), cpt(perf2_idx),...
                        250, 'black', 'filled');
                end

                figure(design);
                if obj.rA > 0
                    scatter3(pts(:, des1_idx), pts(:,des2_idx), pts(:,app_idx),...
                        'MarkerFaceColor', color(i,:), 'MarkerEdgeColor', color(i, :));
                    scatter3(cpt(des1_idx), cpt(des2_idx), cpt(app_idx),...
                        250, 'black', 'filled');
                else
                    scatter(pts(:, des1_idx), pts(:,des2_idx),...
                        'MarkerFaceColor', color(i,:), 'MarkerEdgeColor', color(i, :));
                    scatter(cpt(des1_idx), cpt(des2_idx),...
                        250, 'black', 'filled');
                end
            end
        end
        
        if obj.mFunc.hasGroundTruth && groundTruthOnly
            figure(perf);
            [pts_global, pts_gamut] = obj.mFunc.getGroundTruth();
            disp('Plotting ground truth')
            
            % plot global ground truth (lower envelope of gamut)
            numFrontPieces = length(pts_global);
            for i=1:numFrontPieces
                piece = pts_global{i};
                if obj.rd == 2
                    plot(piece(:, 1),piece(:, 2), 'black','LineWidth',4);
                end
            end
            
            % plot gamut ground truth
            numFrontPieces = length(pts_gamut);
            C = [];
            if groundTruthOnly && exist('providedPts', 'var')
                queryPts = [providedPts(:, perf1_idx), providedPts(:,perf2_idx), providedPts(:,app_idx)];
            end
            for i=1:numFrontPieces
                piece = pts_gamut{i};
                if obj.rd == 2 && obj.rA == 1
                    % compute the distance of scatter plot to ground truth
                    % (if plotting ground truth only)
                    if groundTruthOnly && exist('providedPts', 'var')
                        if isstruct(piece)
                            triStruct = true;
                            s = piece;
                        else
                            triStruct = false;
                            s = surf(piece{1},piece{2},piece{3});
                            s.EdgeColor = 'none';
                        end
                        C = [C, colorByDistance(s, queryPts, triStruct)];
                        markerSize = 40;
%                         set(gca,'ColorScale','log')
                    end
                end
            end
            
            if groundTruthOnly && exist('providedPts', 'var')
                % compute the closest distances over all surface pieces (each column result of different surface patch)
                if size(C, 2) > 1
                   C =  min(C,[],2);
                end

                % color the scatter plot by found points.
                figure;
                scatter3(queryPts(:, 1), queryPts(:, 2), queryPts(:, 3), markerSize, C, 'filled');
                colormap(jet)
                colorbar
                caxis([0 max(C)]);
                xlabel('f1');
                ylabel('f2');
                zlabel('z');
                xlim([0,1]); ylim([0,1]); zlim([0,1]);
                set(gca, 'FontSize', 18);                
            end
            
        end
        
        figure(perf); hold off;
        figure(design); hold off;
    end
    
    
    
    
    
    
    
    % https://www.mathworks.com/help/gads/plot-3-d-pareto-set.html
    % output of this can be fed to surf(XX,YY,ZZ,'LineStyle','none')
    % @param n x 3 matrix of data, where column 1,2,3 are X,Y,Z coords 
    function [XX, YY, ZZ] = interpolatePatch(obj, data)
        eps = 0.005;
        
        minx = min(data(:, 1)); maxx = max(data(:, 1));
        miny = min(data(:, 2)); maxy = max(data(:, 2));
        minz = min(data(:, 3)); maxz = max(data(:, 3));
        
        if abs(minx - maxx) > eps && abs(miny - maxy) > eps %proper function, expand as normal
            sgr = linspace(minx,maxx);
            ygr = linspace(miny,maxy);
            [XX,YY] = meshgrid(sgr,ygr);
            [XX, YY] = obj.maskPatch(data(:, 1), data(:, 2), XX, YY);
            
            F = scatteredInterpolant(data(:,1),data(:,2),data(:,3),'linear','none');
            ZZ = F(XX,YY);
            
        elseif abs(minx - maxx) < eps && abs(miny - maxy) > eps && abs(minz - maxz) > eps
            % vertical plane with x fixed
            sgr = linspace(miny,maxy);
            ygr = linspace(minz,maxz);
            [YY,ZZ] = meshgrid(sgr,ygr);
            [YY, ZZ] = obj.maskPatch(data(:, 2), data(:, 3), YY, ZZ);
            
            XX = minx * ones(size(YY));
            
        elseif abs(minx - maxx) > eps && abs(miny - maxy) < eps && abs(minz - maxz) > eps
            % vertical plane with y fixed
            sgr = linspace(minx,maxx);
            ygr = linspace(minz,maxz);
            [XX,ZZ] = meshgrid(sgr,ygr);
            [XX, ZZ] = obj.maskPatch(data(:, 1), data(:, 3), XX, ZZ);
            
            YY = miny * ones(size(XX));
            
        elseif abs(minx - maxx) < eps && abs(miny - maxy) < eps && abs(minz - maxz) > eps
            % vertical line where only z varies
            ZZ = linspace(minz,maxz);
            ZZ = [ZZ', ZZ'+0.00001];
            YY = miny * ones(size(ZZ));
            XX = minx * ones(size(ZZ));
        elseif abs(minx - maxx) < eps && abs(miny - maxy) > eps && abs(minz - maxz) < eps
            % line where only y varies
            YY = linspace(miny,maxy);
            YY = [YY', YY'+0.00001];
            ZZ = minz * ones(size(YY));
            XX = minx * ones(size(YY));
        elseif abs(minx - maxx) > eps && abs(miny - maxy) < eps && abs(minz - maxz) < eps
            % line where only x varies
            XX = linspace(minx,maxx);
            XX = [XX', XX'+0.00001];
            YY = miny * ones(size(XX));
            ZZ = minz * ones(size(XX));
        elseif abs(minx - maxx) < eps && abs(miny - maxy) < eps && abs(minz - maxz) < eps
            % single point
            XX = []; YY = []; ZZ = [];
        else
            disp('interpolation case reached but not implemented!')
            fprintf('x fixed: %d; y fixed: %d; z fixed: %d;', abs(minx - maxx) < eps, abs(miny - maxy) < eps, abs(minz - maxz) < eps);
        end
        
    end

    
    function [XX, YY] = maskPatch(obj, xdata, ydata, xgrid, ygrid)
        if numel(xdata) < 3 || numel(ydata) < 3 % need at least this many for convex hull
            XX = xgrid;
            YY = ygrid;
            return;
        end
        % implement the meshgrid mask according to convex hull of data points
        % in patch (so it doesn't extend blindly)
        queryx = xgrid(:);
        queryy = ygrid(:);
    
        xlin = xdata(:); ylin = ydata(:);
        [k,~] = convhull(xlin, ylin);
        polyx = xlin(k);
        polyy = ylin(k);

        [inPoly] = inpolygon(queryx,queryy,polyx,polyy); %logical array
        
        % mask out parts of the mesh grid that exceed the patch
        queryx(~inPoly) = NaN;
        queryy(~inPoly) = NaN;
        
        XX = reshape(queryx, size(xgrid));
        YY = reshape(queryy, size(ygrid));
    end
    
end % end of methods
   

end


        