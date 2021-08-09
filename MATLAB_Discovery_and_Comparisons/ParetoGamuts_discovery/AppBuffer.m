classdef AppBuffer < handle
properties
    % static variables
    % problem specific (fixedparams)
    mFunc
    rD
    rd
    rA
    varNames
    
    % optimization specific (userparams)
    numCells 
    deltaCell
    directionDelta
    nMaxRuns
    
    % dynamic variables
    patchArray
    desInts
    perfInts
    appInts
    alphaInts
    
    converged
end

methods
    % Constructor
    function obj = AppBuffer(mFunc, userParams, fixedParams)
        if nargin==0 % needed for Object arrays
            return
        end
        
        % dimensionality of design / perf / application space
        obj.mFunc = mFunc;
        obj.rD = fixedParams.rD;
        obj.rd = fixedParams.rd; 
        obj.rA = fixedParams.rA;
        obj.varNames = fixedParams.varNames;
        
        % Interval structures, to track span of each patch (accelerated searching)
        obj.numCells = userParams.numCells;
        obj.deltaCell = 1.0 / userParams.numCells;
        obj.directionDelta = userParams.directionDelta;
        
        cell = struct('rangeMin', 0, 'rangeMax', 0, 'patchIDs', []);
        cellArray = repmat(cell, 1, obj.numCells);
        obj.desInts = repmat(cellArray, obj.rD, 1); % one cell Array for each design dim
        obj.perfInts = repmat(cellArray, obj.rd, 1);
        obj.appInts = repmat(cellArray, obj.rA, 1);
        obj.alphaInts = repmat(cellArray, obj.rd-1, 1);
        
        obj.patchArray = {};
        
        % termination metrics
        obj.converged = false;
        obj.nMaxRuns = userParams.nMaxRuns;
    end
    
    %% add patch to the buffer and interval structures
    function [] = addPatch(obj, patchID, centerPt, explorationDirections, ... 
                            extents, validPatchPts, connectivityList)
                        
        desPts = validPatchPts(:, 1:obj.rD);
        appPts = validPatchPts(:, obj.rD+1:obj.rD+obj.rA);
        perfPts = validPatchPts(:, obj.rD+obj.rA+1:end);
        
        dSpan = obj.getRangeSpanned(obj.rD, desPts);
        pSpan = obj.getRangeSpanned(obj.rd, perfPts);
        aSpan = obj.getRangeSpanned(obj.rA, appPts);
        
        % convert to (normalized by pi/2) spherical coords (alpha space) and get span on first
        % rd-1 elements
        numPts = size(validPatchPts, 1);
        sphericalPts = zeros(numPts, obj.rd + obj.rA); % 1:rd-1 phi (directions), rd is radius, rd+1:rd+rA is app
        for i=1:numPts
            [phi, r] = cartesian2polar(perfPts(i, :));
            phi = phi / (pi / 2);
            sphericalPts(i, :) = [phi, r, appPts(i, :)];
        end
        alphaSpan = obj.getRangeSpanned(obj.rd-1, sphericalPts(:, 1:obj.rd-1));
        
        % patch array keeps all patch info and points in cartesian coords
        newPatch = struct( ...
            'patchID', patchID, ...
            'centerPt', centerPt, ...
            'expDirs', explorationDirections, ...
            'dirExtents', extents, ...
            'perfSpan', pSpan, ...
            'appSpan', aSpan, ...
            'desSpan', dSpan, ...
            'alphaSpan', alphaSpan, ...
            'patchPts', validPatchPts, ...
            'sphericalPatchPts', sphericalPts ...
        );
        
        obj.patchArray{patchID} = newPatch;  
        
        % place patchID in relevant cells of the interval structure
        obj.placeInIntervalStructs(obj.rD, dSpan, 'desInts', patchID);
        obj.placeInIntervalStructs(obj.rd, pSpan, 'perfInts', patchID);
        obj.placeInIntervalStructs(obj.rA, aSpan, 'appInts', patchID);
        obj.placeInIntervalStructs(obj.rd-1, alphaSpan, 'alphaInts', patchID);
    end
    
    %% Determine the range of values spanned by the patch in question
    % @param numDims - number of axes along which to measure min/max (rD, rd, rA)
    % @params pts - NxnumDims matrix of points in patch. 
    %               Assumes each row is a pt, each col is a particular axis
    function [rangeStruct] = getRangeSpanned(obj, numDims, pts)
        minVal = zeros(1, numDims);
        maxVal = zeros(1, numDims);
        for i=1:numDims
            vals = pts(:, i); %all values for the ith performance metric
            minVal(i) = min(vals);
            maxVal(i) = max(vals);
        end
        rangeStruct.minVals = minVal;
        rangeStruct.maxVals = maxVal;
    end
    
    %% Place patch ID into all cells that it crosses
    % @param numDims - number of axes along which to place the points (rD, rd, rA)
    % @param xSpan - rangeStruct with the min/max extent of the patch in each dim
    % @param xInts - string with the Ints struct to update (eg, 'desInts')
    % @param patchID - ID of patch to place
    function [] = placeInIntervalStructs(obj, numDims, xSpan, xInts, patchID)
        for dim=1:numDims
            minc = obj.getIntCell_singleAxis(xSpan.minVals(dim));
            maxc = obj.getIntCell_singleAxis(xSpan.maxVals(dim));
            minc = max(1, minc);
            maxc = min(obj.numCells, maxc);
            % append to all relevant interval cells
            for c=minc:maxc
                cell = obj.(xInts)(dim, c); 
                patchIDList = [cell.patchIDs; patchID];
                obj.(xInts)(dim, c).patchIDs = patchIDList;
            end
        end
    end
    
    % Identify the cell (within the interval structure) that a certain point 
    % belongs in
    function [cellID] = getIntCell_singleAxis(obj, val)
        rem = floor(val / obj.deltaCell); 
        cellID = rem + 1; 
        
        % make sure it's in the valid range, warn if not
        if cellID > obj.numCells
            fprintf("CellID %d out of bounds. Clamping to max cell.", cellID);
            cellID = obj.numCells;
        elseif cellID < 1
            fprintf("CellID %d out of bounds. Clamping to min cell.", cellID);
            cellID = 1;
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
    
    
    function converged = terminate(obj, iterNum)
        if iterNum > obj.nMaxRuns
            converged = true;
        else
            converged = false;
        end
    end
    
    function empty = isEmpty(obj)
        if isempty(obj.patchArray)
            empty = true;
        else
            empty = false;
        end
    end
    
    
    
    % ==================== IDENTIFY PARETO OPTIMALITY =========================    
    
    %% Identify all patches that intersect a given cell in some space
    % @param numDims - number of axes to intersect for desired cell (eg: rD, rd, rA)
    % @param xInts - string with the Ints struct to search (eg: 'desInts')
    % @param xVals - 1 x numDims array, with the point in x space ; if [],
    %                   use cellIDs instead
    % @param cellIDs - 1 x numDims array, with cell # in x space; only used
    %                   if xVals==[]
    function patches = getRelevantPatches_SingleSpace(obj, numDims, xInts, xVals, cellIDs)
        if numDims == 0 % happens if no application variable, return all patches
            patches = 1:length(obj.patchArray);
            return
        end
        
        for dim=1:numDims
            if length(xVals) > 0 % use point in space
                val = xVals(dim);
                cellID = obj.getIntCell_singleAxis(val);
            else
                cellID = cellIDs(dim);
            end
            cell = obj.(xInts)(dim, cellID);
            p = cell.patchIDs;
            
            if dim==1 
                patches = p; % accept all patches in initial dimension
            else
                patches = intersect(patches,p); % find patches in common
            end
        end
    end

    % Retrieve the pareto front at a particular context
    % @param appVals - 1 x rA array, with the point in application space ; if [],
    %                   use cellIDs instead
    % @param appCells - 1 x rA array, with cell # in application space; only used
    %                   if appVals==[]
    % @param samplesPerDim - number of samples for each dimension in alpha
    %                   (hyperspherical performance) space
    %
    % @return pts - single best point for each of the alpha cells (in cartesian coords, des+app+perf)
    function [pts, patchIDs] = getSingleContextFront(obj, appVals, appCells, samplesPerDim)        
        % identify patches that go through this application value/cell
        if length(appVals) > 0
            appPatches = obj.getRelevantPatches_SingleSpace(obj.rA, 'appInts', appVals, []);
        else
            appPatches = obj.getRelevantPatches_SingleSpace(obj.rA, 'appInts', [], appCells);
        end

        % for every cell in alpha space (rd-dim, context fixed)
        numSamples = samplesPerDim^(obj.rd-1); 
        spaceSize = ones(1, obj.rd-1) * samplesPerDim;
        alphaCellIDs = cell(1, obj.rd-1);
        
        pts = zeros(numSamples, obj.rD+obj.rA+obj.rd);
        patchIDs = zeros(numSamples, 1);
        for sampleID = 1:numSamples
            % get rd values for this cell
            [alphaCellIDs{:}] = ind2sub(spaceSize, sampleID);
            alphaCells = cell2mat(alphaCellIDs);
            
            % get candidate patches, by intersecting app and perf spaces
            alphaPatches = obj.getRelevantPatches_SingleSpace(obj.rd-1, 'alphaInts', [], alphaCells);
            patches = intersect(alphaPatches, appPatches);
            
            % select best point to represent this element
            [pt, patchID] = obj.getBestPointForCell(alphaCells, appVals, appCells, patches);
            pts(sampleID, :) = pt;
            patchIDs(sampleID) = patchID;
        end
        
    end
    
    % Get the best point in a particular cell
    % @param perfCellIDs - 1 x (rd-1) array, with the cell IDs in alpha space
    % @param appVals - 1 x rA array with values in application space; if
    %                   empty, use the appCell values instead
    % @param appCells - 1 x rA array, with cell # in application space; only used
    %                   if appVals==[]
    % @param patchIDs - list of patchIDs that intersect the desired cell
    % @return pt - best point for cell, 1 x (rD+rA+rd) in cartesian coords
    function [bestPt, bestPatchID] = getBestPointForCell(obj, perfCellIDs, appVals, appCellIDs, patchIDs)
        if length(appVals) > 0 % use point in space
            for dim=1:length(appVals)
                val = appVals(dim);
                appCellIDs(dim) = obj.getIntCell_singleAxis(val);
            end
        end
        
        % get representative point for this cell
        numPatches = length(patchIDs);
        bestDist = inf;
        bestPt = ones(1, obj.rD + obj.rA + obj.rd);
        bestPatchID = 0;
        for id=1:numPatches
            [minPatchPt, dist] = obj.getMinPatchPointWithinCell(perfCellIDs, appCellIDs, patchIDs(id));
            if dist < bestDist
                bestDist = dist; 
                bestPt = minPatchPt;
                bestPatchID = patchIDs(id);
            end
        end
    end
    
    % Get the points of a patch that intersect with a particular cell
    % @param perfCellIDs - 1 x (rd-1) array, with the cell IDs in alpha space
    % @param appCells - 1 x rA array, with cell # in application space
    % @param patchID - ID of patch to intersect with cell
    % @return minpt - 1 x (rD+rA+rd) best (min radius) patch point in cell; in
    %               cartesian coords
    % @return dist - distance of minpt to the origin
    function [minpt, dist] = getMinPatchPointWithinCell(obj, perfCellIDs, appCellIDs, patchID)
        patchPts = obj.patchArray{patchID}.sphericalPatchPts;
        numpts = size(patchPts, 1);
        perfPts = patchPts(:, 1:obj.rd-1);
        radPts = patchPts(:, obj.rd); 
        appPts = patchPts(:, obj.rd+1:obj.rd+obj.rA);
        
        tol = (1 / obj.numCells);% / 2;
        
        perfCellCenter = obj.deltaCell * (perfCellIDs - 0.5); 
        appCellCenter = obj.deltaCell * (appCellIDs - 0.5); 

        ptIndices = [];
        for i=1:numpts
            dist_from_center = abs(perfPts(i,:) - perfCellCenter);
            perfPtValid = all(dist_from_center<tol); % within bounds on all dims
            
            dist_from_center = abs(appPts(i,:) - appCellCenter);
            appPtValid = all(dist_from_center<tol); % within bounds on all dims
            
            if perfPtValid && appPtValid
                ptIndices= [ptIndices, i];
            end
        end
        
        % get index of min point (minimum radius), out of the valid indices
        [rad, idx] = min(radPts(ptIndices)); %idx is location of min in valid pts
        orig_idx = ptIndices(idx);              % get idx of min in original pts
        
        % get min pt from cartesian points.
        minpt = obj.patchArray{patchID}.patchPts(orig_idx, :);
        
        % for benefit of Javascript visualization -- rewrite the application value
        %minpt(obj.rD+1:obj.rD+obj.rA) = appCellCenter;
        
        dist = rad;
    end
    
    function [pts, patchIDs] = getParetoGamut(obj)
        % for now, we'll visualize the same number of cells as we store for
        % each interval; may want to revisit that. 
        samplesPerAppDim = obj.numCells;
        samplesPerPerfDim = obj.numCells;
        
        % for every cell in app space (rA-dim)
        numSamples = samplesPerAppDim^(obj.rA);
        spaceSize = ones(1, obj.rA) * samplesPerAppDim;
        appCellIDs = cell(1, obj.rA);

        pts = [];
        patchIDs = [];
        for sampleID = 1:numSamples
            % get rd values for this cell
            [appCellIDs{:}] = ind2sub(spaceSize, sampleID);
            
            % gives single best point for each of the alpha cells (in cartesian coords, des+app+perf)
            [singlepts, singlePatchIDs] = obj.getSingleContextFront([], cell2mat(appCellIDs), samplesPerPerfDim);
            
            % pass performance points in, to get non-dominate points
            nonDominated = getParetoIndices(singlepts(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
            paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices

            
            pts = [pts; singlepts(paretoIndices, :)];
            patchIDs = [patchIDs; singlePatchIDs(paretoIndices, :)];
        end
        
    end
    
    
    
    
    
    % ================== VISUALIZATION ==========================
    function [design, perf] = visualizeAppBuff(obj, idxToVis, providedPts, providedLabels)
        numPatches = length(obj.patchArray);
        
        if nargin < 2 % only obj 
            d1=1; d2=2; a1=1; p1=1; p2=2; % default values
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
        perf = figure('Name', 'Performance Space'); hold on;    
        % always keep the same proportions
        xlim([0,1]); xlabel(obj.varNames(perf1_idx))
        ylim([0,1]); ylabel(obj.varNames(perf2_idx))
        if obj.rA > 0
            zlim([0,1]); zlabel(obj.varNames(app_idx))
        end
        axis equal;

        design = figure('Name', 'Design Space'); hold on;
        % always keep the same proportions
        xlim([0,1]); xlabel(obj.varNames(des1_idx))
        ylim([0,1]); ylabel(obj.varNames(des2_idx))
        if obj.rA > 0
            zlim([0,1]); zlabel(obj.varNames(app_idx))
        end
        axis equal; 

        if exist('providedPts', 'var')
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
            
            figure(perf); % create or get performance space fig
            if obj.rA > 0
                scatter3(pts(:, perf1_idx), pts(:,perf2_idx), pts(:,app_idx),...
                    markerSize, markerColors);
            else
                scatter(pts(:, perf1_idx), pts(:,perf2_idx),...
                    markerSize, markerColors);
            end

            figure(design); % create or get design space fig
            if obj.rA > 0
                scatter3(pts(:, des1_idx), pts(:,des2_idx), pts(:,app_idx),...
                    markerSize, markerColors);
            else
                scatter(pts(:, des1_idx), pts(:,des2_idx),...
                    markerSize, markerColors);
            end
                
        else
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
        
        % ===== visualize ground truth ======
        if obj.mFunc.hasGroundTruth
            figure(perf);
            [pts_global, pts_gamut] = obj.mFunc.getGroundTruth();
            
            % plot global ground truth (lower envelope of gamut)
            numFrontPieces = length(pts_global);
            for i=1:numFrontPieces
                piece = pts_global{i};
                if obj.rd == 2
                    disp('Plotting ground truth')
                    plot(piece(:, 1),piece(:, 2), 'black','LineWidth',4);
                end
            end
            
            % plot gamut ground truth
            numFrontPieces = length(pts_gamut);
            for i=1:numFrontPieces
                piece = pts_gamut{i};
                if obj.rd == 2 && obj.rA == 1
                    disp('Plotting ground truth')
                    surf(piece{1},piece{2},piece{3});
                end
            end
        end
        
        figure(perf); hold off;
        figure(design); hold off;
    end
end % end of methods
    
end


        