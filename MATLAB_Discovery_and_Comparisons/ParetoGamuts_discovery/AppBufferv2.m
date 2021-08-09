classdef AppBufferv2 < handle
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
    directionDelta
    nMaxRuns
    
    % dynamic variables
    patchArray
    radialCellDims
    radialPerfArray
    
    appCellDims
    
    converged
end

methods
    % Constructor
    function obj = AppBufferv2(mFunc, userParams, fixedParams)
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
        obj.directionDelta = userParams.directionDelta;

        obj.patchArray = {};

        numRadialPerfCells = obj.numCells^(obj.rd - 1);
        obj.radialCellDims = repmat({obj.numCells}, obj.rd -1, 1);
        obj.radialPerfArray = cell([numRadialPerfCells, 1]);
        
        obj.appCellDims = repmat({obj.numCells}, obj.rA, 1);
        
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

        % convert to (normalized by pi/2) spherical coords (alpha space) and get span on first
        % rd-1 elements
        numPts = size(validPatchPts, 1);
        sphericalPts = zeros(numPts, obj.rd + obj.rA); % 1:rd-1 phi (directions), rd is radius, rd+1:rd+rA is app
        for ptID=1:numPts
            [phi, r] = cartesian2polar(perfPts(ptID, :));
            phi = phi / (pi / 2);                       % normalize all components to [0,1]
            
            % store it -- might not actually be necessary
            sphericalPts(ptID, :) = [phi, r, appPts(ptID, :)];
            
            % place point into radialPerfBuffer
            obj.placePtInRadialBuffer(patchID, ptID, phi, r);
        end
        
        % patch array keeps all patch info and points in cartesian coords
        newPatch = struct( ...
            'patchID', patchID, ...
            'centerPt', centerPt, ...
            'expDirs', explorationDirections, ...
            'dirExtents', extents, ...
            'patchPts', validPatchPts, ...
            'sphericalPatchPts', sphericalPts ...
        );
        
        obj.patchArray{patchID} = newPatch; 
        
    end
    
    function [] = placePtInRadialBuffer(obj, patchIdx, ptIdx, ptPhi, ptRadius)
        cellIDlin = obj.getLinearCellIdx(ptPhi, obj.radialCellDims);
        
        % place a reference to this point id in the appropriate buffer cell 
        newEntry = [patchIdx, ptIdx, ptRadius];
        
        obj.radialPerfArray{cellIDlin} = [obj.radialPerfArray{cellIDlin}; newEntry];
    end
    
    
    % Identify the cell (within the interval structure) that a certain point 
    % belongs in
    function [cellID] = getSubscriptedCellIdx(obj, vals)
        rem = floor(vals * obj.numCells);
        if vals < 1         % account for the upper edge of the very last cell (would otherwise go to the next cell up)
            cellID = rem + 1; 
        else
            cellID = rem;
        end
        
        % make sure it's in the valid range, warn if not
        if any(cellID > obj.numCells)
%             fprintf("CellID out of bounds. Clamping to max cell.\n");
            cellID = min(cellID, obj.numCells);
        end
        if any(cellID < 1)
%             fprintf("CellID out of bounds. Clamping to min cell.\n");
            cellID = max(cellID, 1);
        end
    end
    
    
    function [cellIDlin] = getLinearCellIdx(obj, vals, dimCell)
        cellIDsub = obj.getSubscriptedCellIdx(vals);
        
        if length(cellIDsub) > 1
            sz = cell2mat(dimCell);
            cellIDsub = num2cell(cellIDsub);        % have to input each dim as separate arg to sub2ind, only possible with cell array
            cellIDlin = sub2ind(sz, cellIDsub{:});
        else
            cellIDlin = cellIDsub;
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
    
    
    
    
    
    %% ==================== IDENTIFY PARETO OPTIMALITY =========================    
    function [finalPts, finalLabels] = getParetoGamut(obj)
        numPtsToKeep = 1;
        candPts = []; candLabels = []; candContexts = [];
        
        % loop over a particular cell in the spherical performance space
        for linearCellID = 1:length(obj.radialPerfArray)    
            % sort the points by radius, so pts(1) is closest to origin in cell, and points(end) is furthest
            pts = obj.radialPerfArray{linearCellID}; 
            
            if isempty(pts)
                continue;
            end
            
            radiusColIdx = 3;                       % each row is [patchID, ptID(in patchArray.patchPts), radius]
            pts = sortrows(pts, radiusColIdx);
            
            contextCells = cell([obj.numCells^(obj.rA), 1]);
            numContexts = length(contextCells);
            
            % loop over all points in this cell of the spherical
            % performance (alpha) buffer
            for i=1:size(pts, 1)
                % get the actual point from patch Arrays
                patchID = pts(i,1); ptID =pts(i, 2);
                patch = obj.patchArray{patchID};
                p = patch.patchPts(ptID, :);
                
                % figure out the appropriate contextual cell for this point
                if obj.rA > 0
                    appVals = p(1, obj.rD+1:obj.rD+obj.rA);
                    appCellID = obj.getLinearCellIdx(appVals, obj.appCellDims);
                else
                    appCellID = 1;
                end
                
                % if the context has n points, this is the (n+1)th best
                % point we've found for this performance cell; save it as
                % one of the best if (n+1) < numToKeep (will be whittled to 1 in graph cut)
                if size(contextCells{appCellID}, 1) < numPtsToKeep
                    contextCells{appCellID} = [contextCells{appCellID}; p, patchID, appCellID];
                end
            end
            
            % parse the context cells into the final return buffer
            for i=1:numContexts
                points = contextCells{i};
                numCols = size(points, 2);
                if ~isempty(points)
                 
                    p = points(:, 1:numCols-2); % design, application, and performance point
                    l = points(:, numCols-1);     % patch label for this point
                    c = points(:, numCols);         % context ID for this point

                    candPts = [candPts; p];
                    candLabels = [candLabels; l];
                    candContexts = [candContexts; c];
                end
            end
        end
        
        % must check for dominated points in each context
%         % sort by the given context
%         [candContexts, idx] = sortrows(candContexts, 1);
%         candLabels = candLabels(idx, :);
%         candPts = candPts(idx, :);
        finalPts = []; finalLabels = [];
        for i=1:numContexts
            singlepts_idx = find(candContexts == i);
            singlePts = candPts(singlepts_idx, :);
            singleLabels = candLabels(singlepts_idx, :);
            
            % pass performance points in, to get non-dominate points
            nonDominated = getParetoIndices(singlePts(:, obj.rD+obj.rA+1:obj.rD+obj.rA+obj.rd)); %gives logical array
            paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices

            finalPts = [finalPts; singlePts(paretoIndices, :)];
            finalLabels = [finalLabels; singleLabels(paretoIndices, :)];
        end
    end

    
    
    
    
    
    
    %% ================== VISUALIZATION ==========================
    function [design, perf] = visualizeAppBuff(obj, idxToVis, groundTruthOnly, providedPts, providedLabels)
        plotPoints = true;
        plotSurfaces = false;
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
        set(gca,'FontSize', 14);

        design = figure('Name', 'Design Space'); hold on;
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
            
            figure(perf); % create or get performance space fig
            if obj.rA > 0
                % plot the observed data
                if plotPoints
                    scatter3(pts(:, perf1_idx), pts(:,perf2_idx), pts(:,app_idx),...
                        markerSize, markerColors);
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
            else
                scatter(pts(:, perf1_idx), pts(:,perf2_idx),...
                    markerSize, markerColors);
                % plot the interpolated curves for each
                for i=1:numPatchIDs
                    patchIds = find(mappedLabels == i);
                    patchPts = pts(patchIds, [perf1_idx, perf2_idx]);
                    
                    plot(patchPts(:, 1), patchPts(:, 2))
                end
            end

            figure(design); % create or get design space fig
            if obj.rA > 0
                if plotPoints
                    scatter3(pts(:, des1_idx), pts(:,des2_idx), pts(:,app_idx),...
                        markerSize, markerColors);
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
                for i=1:numPatchIDs
                    patchIds = find(mappedLabels == i);
                    patchPts = pts(patchIds, [des1_idx, des2_idx]);
                    
                    plot(patchPts(:, 1), patchPts(:, 2))
                end
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
        
        if obj.mFunc.hasGroundTruth
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
            queryPts = [providedPts(:, perf1_idx), providedPts(:,perf2_idx), providedPts(:,app_idx)];
            for i=1:numFrontPieces
                piece = pts_gamut{i};
                if obj.rd == 2 && obj.rA == 1
                    s = surf(piece{1},piece{2},piece{3});
                    
                    % compute the distance of scatter plot to ground truth
                    % (if plotting ground truth only)
                    if groundTruthOnly 
                        C = [C, colorByDistance(s, queryPts)];
                        markerSize = 40;
%                         set(gca,'ColorScale','log')
                    end
                end
            end
            
            if groundTruthOnly
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


        