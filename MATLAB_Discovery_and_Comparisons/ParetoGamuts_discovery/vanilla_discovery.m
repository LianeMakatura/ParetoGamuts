% we want to take the typical mapping function with constrained z

function [buffers_c] = vanilla_discovery(constrained_functionID, unconstrained_functionID, ...
                                        z_min, z_max, numSteps, interContextSharing, contextSharingRadius) 
    all_linear_patch_pts = {};
    linear_idx = 0;
                                    
    %% Construct the problem statement as a set of mapping functions 
    % which deals with evaluations of the appropriate objective functions
    mFunc_uc = MappingFunction(unconstrained_functionID); % might need dummy zval
    if constrained_functionID == 25 || constrained_functionID == 27
        disp('found rA')
        properZsep = true; % indicate that theres proper separation between x, z
        rA = mFunc_uc.rA;
    else
        properZsep = false;
        rA = 1;
    end
    
    stepSize = (z_max - z_min) / numSteps;
    mFuncs_c(numSteps+1, 1) = MappingFunction();         % create object array for mFuncs
    for i=0:numSteps
        z = z_min + i*stepSize;
        mFunc_c = MappingFunction(constrained_functionID, z);
        mFuncs_c(i+1) = mFunc_c;
    end
    
    if interContextSharing % create func cell array for expansion dirs xz
        funcsCell = {mFunc_uc.f1, mFunc_uc.f2};
        rD = mFuncs_c(1).rD;
        if properZsep
            xSymb = sym('x', [1 rD-1]);
            zSymb = sym('z', [1 rA]);
        else
            xSymb = sym('x', [1 rD]);
            zSymb = xSymb(1, rD-rA+1:rD); % get last element of the xSymb array
            xSymb = xSymb(1, 1:rD-rA); % shrink x symb vec appropriately
        end
        sharingOffset = linspace(-contextSharingRadius, contextSharingRadius, 2*contextSharingRadius+1);
    end

    %% Initialize parameters fixed by the problem statement
    % & user adjustable params
    
    % unconstrained problem
    imageTitle = mFunc_uc.getImageTitle();
    [rD,rd] = mFunc_uc.getDimensions();
    numCells = mFunc_uc.getNumCells(); 
    [minBound, maxBound] = mFunc_uc.getBounds();
    fixedParams_uc = struct('rD', rD, 'rd', rd, 'imageTitle', imageTitle, 'minBound', minBound, ...
                         'maxBound', maxBound, 'constrain_z', false);
 
    userParams_uc = struct('tolerance', 0.01, 'directionDelta', 0.3, 'numCells', numCells, 'patchSize', 0.5*numCells,...
        'stepSize', 1/numCells, 'numInitialSamples', 10, 'numNewSamples', 10, ...
        'numRepeatedNoImprovement', 3, 'breakScoreRelative', 0.1/numCells, 'breakScoreDistance', 0.0001, ...
        'neighForScore', ceil(0.05*numCells),'neighForOpt', ceil(0.2*numCells), 'nMaxRuns', 15);

    % constrained problem -- I believe the same set of params can be used
    % for all constrained instances; may have to change
    
    imageTitle = mFunc_c.getImageTitle();
    [rD,rd] = mFunc_c.getDimensions();
    numCells = mFunc_c.getNumCells(); 
    [minBound, maxBound] = mFunc_c.getBounds();
    fixedParams_c = struct('rD', rD, 'rd', rd, 'imageTitle', imageTitle, 'minBound', minBound, ...
                         'maxBound', maxBound, 'constrain_z', true);

    userParams_c = struct('tolerance', 0.01, 'directionDelta', 0.3, 'numCells', numCells, 'patchSize', 0.5*numCells,...
                        'stepSize', 1/numCells, 'numInitialSamples', 10, 'numNewSamples', 10, ...
                        'numRepeatedNoImprovement', 3, 'breakScoreRelative', 0.1/numCells, 'breakScoreDistance', 0.0001, ...
                        'neighForScore', ceil(0.05*numCells),'neighForOpt', ceil(0.2*numCells), 'nMaxRuns', 15);
    disp([rd, rD, rA])
%    if ~properZsep
    rD = rD - rA;
 %   end
    
    %% Initialize buffer data structure
    buffer_uc = Buffer(userParams_uc, fixedParams_uc);
    
    buffers_c(numSteps+1, 1) = Buffer(); % create object array for buffers
    for i=0:numSteps
        b = Buffer(userParams_c, fixedParams_c);
        buffers_c(i+1) = b;
    end
    
    %% Initialize figures class for display
    %figures = showFigures(rD, rd, buffer, functionID);

    %% Initialize a random population of samples
    if interContextSharing % only use a fraction of the samples in each slice
        initSamplesPerZ = ceil(userParams_c.numInitialSamples / (numSteps+1));
        newSamplesPerZ = ceil(userParams_c.numNewSamples / (numSteps+1));
    else % use all samples in each slice
        initSamplesPerZ = userParams_c.numInitialSamples;
        newSamplesPerZ = userParams_c.numNewSamples;
    end


    %% Start main loop
    tic;
    disp('Starting main loop');
    patchIDs = zeros(numSteps+1);
    
    num_zero_projz = 0;
    num_proj_z = 0;
    for run=1:userParams_c.nMaxRuns
        disp('----------------------------------------------------------------------');
        disp(['Iteration ' num2str(run)]);
        
        % get samples for first slice, either randomly or focusing on
        % under-explored areas
        % CHANGE THIS BACK (seeing if constrained z the issue)
        curr_z_val = z_min + stepSize;
        curr_z_idx = 2;
        if buffers_c(curr_z_idx).isEmpty()
            disp('buff is empty')
            randomSamples = mFuncs_c(curr_z_idx).sample(initSamplesPerZ, curr_z_val);
        else
            [randomPoints] = buffers_c(curr_z_idx).getPointsToImprove(newSamplesPerZ);
            randomSamples = [randomPoints mFuncs_c(curr_z_idx).eval(randomPoints, curr_z_val)];
        end
        
        % Loop over all "slices" of the z value
        for step=2:numSteps %% +1 (AVOIDING DIVISION BY 0 FOR LBRACKET)
            fprintf('curr z idx=%d, value=%0.3f\n', curr_z_idx, curr_z_val);
            %% Step 1: Add each random point to the buffer first in case we get lucky
            
            % NOTE: This is for experimentation with 'z' as a design variable.
            %[numRandomSamples, ~] = size(randomSamples);
    %         randomSamples(:,rD) = ones(numRandomSamples,1)*z;

            % this slice already converged; keep going
            if buffers_c(curr_z_idx).converged
                if curr_z_idx < numSteps+1
                    curr_z_val = curr_z_val + stepSize;
                    curr_z_idx = curr_z_idx + 1;
                    if buffers_c(curr_z_idx).isEmpty()
                        randomSamples = mFuncs_c(curr_z_idx).sample(initSamplesPerZ, curr_z_val);
                    else
                        [randomPoints] = buffers_c(curr_z_idx).getPointsToImprove(newSamplesPerZ);
                        randomSamples = [randomPoints mFuncs_c(curr_z_idx).eval(randomPoints, curr_z_val)];
                    end
                end
                continue;
            end
            
            for i = 1:size(randomSamples,1)
                buffers_c(curr_z_idx).addPoint(randomSamples(i,:));
            end


            %% Step 2: Locally optimize current samples
            % Get target directions for our random samples
            [targets, dirs] = buffers_c(curr_z_idx).getRandomDirections(randomSamples(:,rD+rA+1:rD+rA+rd));

            % Attempt to push each sample at its corresponding target point
            [optimizedSamples, optSampAlphas, optSampBetas] = ...
                optimize(randomSamples, targets, dirs, mFuncs_c(curr_z_idx), fixedParams_c, curr_z_val);

            disp('Finished local optimizations');

            %% Step 3: Perform local exploration on optimized points
            for v = 1:size(optimizedSamples,1)
                if interContextSharing
                    disp('SHARING!-----');
                    centerPoint = optimizedSamples(v,:);
                    centerAlpha = optSampAlphas(v, :)';
                    centerBeta = optSampBetas(v, :)';
%                     test_uc = mFunc_uc.eval(centerPoint(1:rD), curr_z_val); %should be same as orig center point 3,4,5
%                     test_c = mFuncs_c(curr_z_idx).eval(centerPoint(1:rD), curr_z_val); %should be same as orig center point 3,4,5
                        
                    % Get the direction(s) around which to locally explore each sample
                    xVals = centerPoint(1:rD);
                    zVals = centerPoint(rD+1:rD+rA);
                    [xz_explorationDirection_orig, numActiveC] = getExplorationDirections_xz(funcsCell, ...
                                                xSymb, xVals, zSymb, zVals, ...
                                                centerAlpha, centerBeta, ...
                                                rD, rd, rA);

                    % Separate x'z' from alpha' (and numActiveC beta', if constrained)
                    alphaPrimes = xz_explorationDirection_orig(1:rd, :);
                    betaPrimes = xz_explorationDirection_orig(rd+1:rd+numActiveC, :);
                    xz_explorationDirection = xz_explorationDirection_orig(rd +numActiveC+1:rd+numActiveC+rD+rA, :);
                             
                    % find trustworthy expansion in design performance space
                    [des_extents, perf_extents, principal_dirs, principal_extents] = getExtents(...
                                        centerPoint, centerAlpha, centerBeta, ...
                                        alphaPrimes, betaPrimes, ...
                                        xz_explorationDirection, mFunc_uc, rD, rd, rA);
                        
                    % VISUALIZATION: initialize arrays to collect predicted/optimized pts
                    principal_lines = [expandLine(centerPoint, principal_dirs(:, 1), principal_extents(1), mFunc_uc, rD, rd, rA);
                                        expandLine(centerPoint, principal_dirs(:, 2), principal_extents(2), mFunc_uc, rD, rd, rA)];
                    
                    predictedSharedPts = [];
                    correctedPts = [];
                    % ==== expand patch by given directions
                    patchPts = patchFromDirections(centerPoint, xz_explorationDirection, userParams_uc, ...
                                            mFunc_uc, rD, rA, rd, des_extents, false);
                    linearPatchPts_center = patchFromDirections(centerPoint, xz_explorationDirection, userParams_uc, ...
                                            mFunc_uc, rD, rA, rd, perf_extents, true);
                    linear_idx = linear_idx + 1;
                    all_linear_patch_pts{linear_idx} = linearPatchPts_center;
                    
                    % Push point to other fronts in sharing radius
                    for idx=1:2*contextSharingRadius + 1
                        offset = sharingOffset(idx);
                        nearby_z_idx = curr_z_idx + offset;
                        if nearby_z_idx < 1 ...
                                || nearby_z_idx > numSteps + 1 ...
                                || offset == 0
                            continue; % out of bounds or self (offset=0)
                        end
                        fprintf('Sharing from idx %d to nearby idx %d', curr_z_idx, nearby_z_idx)
                        % push to nearby front 
                        % ================ NEW APPROACH TO PUSHING ========
                        targ_dir = [zeros(1, rD),1]; % move toward app var (Assumes single app var)
                        
                        % project onto exploration subspace
                        proj_dir = proj_vec2subspace(targ_dir', xz_explorationDirection');
                        proj_dir = proj_dir'; % want row vector
                        
                        % follow this direction until z=c
                        nearby_z_val = curr_z_val + offset * stepSize;
                        centerPoint_z = centerPoint(rD+rA); %% ASSUMES 1 Z VALUE
                        proj_dir_z = proj_dir(rD+rA);
                        if proj_dir_z==0
                            num_zero_projz = num_zero_projz + 1;
                            keyboard;
                            continue;
                        end
                        num_proj_z = num_proj_z + 1;
                        t = (nearby_z_val - centerPoint_z) / proj_dir_z;
                        new_des_pt = centerPoint(1:rD+rA) + t.*proj_dir;
                        
                        new_perf_pt = mFuncs_c(nearby_z_idx).eval(new_des_pt, nearby_z_val); %eval new design point
                        
                        predictedSharedPts = [predictedSharedPts; [new_des_pt, new_perf_pt]];
                        
                        % =============END NEW APPROACH TO PUSHING ========
                        
                        % ---- optimize in that front -----
                        % Get target directions for our random samples
                        [targets, dirs] = buffers_c(nearby_z_idx).getRandomDirections(new_perf_pt);
                        % Attempt to push each sample at its corresponding target point
                        [optimizedPt, optPtAlpha, optPtBeta] = ...
                            optimize([new_des_pt, new_perf_pt], targets, dirs, mFuncs_c(nearby_z_idx), fixedParams_c, nearby_z_val);

                        % test without additional optimization correction
%                         optimizedPt = [new_des_pt, new_perf_pt];
                        correctedPts = [correctedPts; optimizedPt];


                        % ---- expand in neighboring front ---- 
                        % NOT SURE ABOUT THE +rA
                        explorationDirection = getExplorationDirection(optimizedPt(1:rD+rA),mFuncs_c(nearby_z_idx), fixedParams_c);

                        % Perform the exploration
                        [patch, validation] = expandPatch(explorationDirection, optimizedPt, ...
                                                mFuncs_c(nearby_z_idx), buffers_c(nearby_z_idx), ...
                                                userParams_c, false, nearby_z_val);

                        % Add the patch to the buffer if the former is nonempty
                        if (size(patch,1) > 0)
                            patchIDs(nearby_z_idx) = patchIDs(nearby_z_idx) +1;
                            buffers_c(nearby_z_idx).addPatch(patch, patchIDs(nearby_z_idx), ...
                                                validation, explorationDirection, optimizedPt);
%                             if nearby_z_idx < numSteps + 1
%                                 linearPatchPts = linearizedPerf(mFunc_uc, optimizedPt, patch, rD, rd, rA);
%                                 linear_idx = linear_idx + 1;
%                                 all_linear_patch_pts{linear_idx} = linearPatchPts;
%                             end
                        else
                            disp('patch non existent')
                        end
                        
                    end
                    
                    % ----- Expand point in curr front (as before) ------
                    explorationDirection = getExplorationDirection(centerPoint(1:rD+rA),mFuncs_c(curr_z_idx), fixedParams_c);

                    % Perform the exploration
                    [patch, validation] = expandPatch(explorationDirection, centerPoint, ...
                                            mFuncs_c(curr_z_idx), buffers_c(curr_z_idx), ...
                                            userParams_c, false, curr_z_val);

                    % Add the patch to the buffer if the former is nonempty
                    if (size(patch,1) > 0)
                        patchIDs(curr_z_idx) = patchIDs(curr_z_idx) +1;
                        buffers_c(curr_z_idx).addPatch(patch, patchIDs(curr_z_idx), ...
                                            validation, explorationDirection, centerPoint);
                                        
%                         linearPatchPts = linearizedPerf(mFunc_uc, centerPoint, patch, rD, rd, rA);
%                         linear_idx = linear_idx + 1;
%                         all_linear_patch_pts{linear_idx} = linearPatchPts; 
                    end
                    
                    
                    
                    % display the results every so often
                    if mod(step + v, 2) == 0 %0 (3 never happens)
                        intermediateVis(buffers_c, centerPoint, patchPts, predictedSharedPts,...
                                    correctedPts, z_min, z_max, rD, rd, rA, linearPatchPts_center,...
                                    principal_lines);
                    end
                    fprintf('\n--- percentage of zero proj_dir_z is %0.6f', num_zero_projz / num_proj_z)

                % not sharing, do only the in-context expansion
                else
                    centerPoint = optimizedSamples(v,:);
                    % Get the direction(s) around which to locally explore each sample
                    explorationDirection = getExplorationDirection(centerPoint(1:rD+rA),mFuncs_c(curr_z_idx), fixedParams_c);

                    % Perform the exploration
                    [patch, validation] = expandPatch(explorationDirection, centerPoint, ...
                                            mFuncs_c(curr_z_idx), buffers_c(curr_z_idx), ...
                                            userParams_c, false, curr_z_val);

                    % Add the patch to the buffer if the former is nonempty
                    if (size(patch,1) > 0)
                        patchIDs(curr_z_idx) = patchIDs(curr_z_idx) +1;
                        buffers_c(curr_z_idx).addPatch(patch, patchIDs(curr_z_idx), ...
                                            validation, explorationDirection, centerPoint);
                    end
                end
            end

            disp('Finished expanding patches');

            %% Step 4: Check stopping criteria
            if buffers_c(curr_z_idx).terminate(run)
                fprintf('Single-context buffer for z=%0.3f converged!\n', curr_z_val);
                buffers_c(curr_z_idx).converged = true;
            end

            %% Step 5: Update z and get new samples for next iteration
            if curr_z_idx < numSteps+1
                curr_z_val = curr_z_val + stepSize;
                curr_z_idx = curr_z_idx + 1;
                if buffers_c(curr_z_idx).isEmpty()
                    randomSamples = mFuncs_c(curr_z_idx).sample(initSamplesPerZ, curr_z_val);
                else
                    [randomPoints] = buffers_c(curr_z_idx).getPointsToImprove(newSamplesPerZ);
                    randomSamples = [randomPoints mFuncs_c(curr_z_idx).eval(randomPoints, curr_z_val)];
                end
            end
        end
        % finished all z-values
        % check convergence for whole problem
        gamut_converged = true;
        for i=1:numSteps+1
            if ~buffers_c(i).converged
                gamut_converged = false;
                fprintf('slice %d not converged yet\n', i)
                break;
            end
        end
        
        % if whole gamut converged, stop optimizing
        if gamut_converged
            disp('Full gamut converged!');
            break;
        end
    end
    
    % plot all patches
    figure; hold on;
    for i=1:linear_idx
        pts = all_linear_patch_pts{i};
        scatter3(pts(:, 4), pts(:, 5), pts(:, 3));
    end
    hold off;
    
    %% Post processing
%     buffer.labelPatches(patchID + 1);
%     buffer.fillEmptyCells(userParams.stepSize, mFunc);
%     disp('Finished labeling patches');
    toc;
    fprintf('^ timing for main loop only, sharing = %d', interContextSharing)
end



function [] = intermediateVis(buffers, centerPt, patchPts, predictedSharedPts, ...
                                correctedPts, zmin, zmax, rD, rd, rA, linearPatchPts,...
                                principal_lines)
    perf = figure('Name', 'Performance Space'); hold on;    
    % always keep the same proportions
    xlim([0,1])
    ylim([0,1])
    zlim([zmin, zmax])
    zlabel('application var')
    axis equal;
    
    design = figure('Name', 'Design Space'); hold on;
    % always keep the same proportions
    xlim([0,1])
    ylim([0,1])
    zlim([zmin, zmax])
    zlabel('application var')
    axis equal; 
    
    % ===== pick indices to visualize
    des1_idx = 1;           % between 1 and rD
    des2_idx = 2;
    app_idx = rD + 1;       % between rD+1 and rD+rA
    perf1_idx = rD+rA+1;     % between rD+rA+1 and rD+rA+rd 
    perf2_idx = rD+rA+2;
    disp([des1_idx, des2_idx, app_idx, perf1_idx, perf2_idx])
                           
    
    % ==== put all (pareto optimal) buffer points in single array
    numBuffs = length(buffers);
    allParetoPts = [];
    
    for buffIdx=1:numBuffs
        b = buffers(buffIdx);
        paretoIndices = b.getParetoInd;

        for j = 1:length(paretoIndices)
            entry = b.Buff(paretoIndices(j));
            patch = entry.bestPatch;
            if  patch > 0 
                % Store all the points in a single matrix
                newPoint = [entry.minD entry.minF];
                allParetoPts = [allParetoPts; newPoint];
            end
        end
    end
   
    
    % ==== draw all pareto optimal points
    figure(perf) % create or get performance space fig
    scatter3(allParetoPts(:,perf1_idx), allParetoPts(:,perf2_idx), allParetoPts(:,app_idx),...
            'green', 'filled'); % perf space
    scatter3(centerPt(perf1_idx), centerPt(perf2_idx), centerPt(app_idx), 250, ...
            'black', 'filled');

    
    figure(design); % create or get design space fig
    scatter3(allParetoPts(:,des1_idx), allParetoPts(:,des2_idx), allParetoPts(:,app_idx),...
            'green', 'filled'); % design space
    
    % ==== draw center point in distinct manner
    scatter3(centerPt(des1_idx), centerPt(des2_idx), centerPt(app_idx), 250, ...
        'black', 'filled');
    
    % ==== display predicted shared & corrected points
    figure(perf) % create or get performance space fig
    scatter3(patchPts(:,perf1_idx), patchPts(:,perf2_idx), patchPts(:,app_idx),...
            'magenta', 'filled');
    scatter3(predictedSharedPts(:,perf1_idx), predictedSharedPts(:,perf2_idx), predictedSharedPts(:,app_idx), 100,...
            'red', 'filled');
    scatter3(correctedPts(:,perf1_idx), correctedPts(:,perf2_idx), correctedPts(:,app_idx), 100,...
            'blue', 'filled');
    scatter3(principal_lines(:,perf1_idx), principal_lines(:,perf2_idx), principal_lines(:,app_idx), ...
            'black', 'filled');
        
    if exist('linearPatchPts', 'var')
        scatter3(linearPatchPts(:,perf1_idx), linearPatchPts(:,perf2_idx), linearPatchPts(:,app_idx),...
            'cyan', 'filled');
    end


    figure(design); % create or get design space fig
    scatter3(patchPts(:,des1_idx), patchPts(:,des2_idx), patchPts(:,app_idx), ...
            'magenta', 'filled'); % design space
    scatter3(predictedSharedPts(:,des1_idx), predictedSharedPts(:,des2_idx), predictedSharedPts(:,app_idx), 100,...
            'red', 'filled'); % design space
    scatter3(correctedPts(:,des1_idx), correctedPts(:,des2_idx), correctedPts(:,app_idx), 100,...
            'blue', 'filled'); % design space
    scatter3(principal_lines(:,des1_idx), principal_lines(:,des2_idx), principal_lines(:,app_idx), ...
            'black', 'filled');
    if exist('linearPatchPts', 'var')
        scatter3(linearPatchPts(:,des1_idx), linearPatchPts(:,des2_idx), linearPatchPts(:,app_idx),...
            'cyan', 'filled');
    end
    
    drawnow;
    keyboard;
    figure(perf);
    close;
    figure(design);
    close;
end

function linepts = expandLine(centerPt, dir, extent, mFunc_uc, rD, rd, rA)
    stepSize=0.01;
    steps = -extent:stepSize:extent;
    
    centerDes = centerPt(1:rD+rA);
    linepts = zeros(length(steps), rD+rA);
    for i=1:length(steps)
        linepts(i, :) = centerDes + steps(i)*dir';
    end
    
        % Borrowed from expandPatch
    validation = mFunc_uc.validate(linepts);
    validIndices = find(validation);

    evaluatedPoints = zeros(size(linepts,1), rd);
    if mFunc_uc.id == 26 %% id where x and z are separated
        evaluatedPoints(validIndices,:) = mFunc_uc.eval(linepts(validIndices,1:rD), ...
                                                    linepts(validIndices,rD+1:rD+rA)); 
    else
        % the z parameter (0 here) doesn't mean anything; not used, just
        % required
        evaluatedPoints(validIndices,:) = mFunc_uc.eval(linepts(validIndices,:), 0);
    end
    
    linepts = [linepts, evaluatedPoints];
    linepts = linepts(validIndices, :);
end

function patchPts = patchFromDirections(centerPt, explorationDirection, userParams_uc, ...
                            mFunc_uc, rD, rA, rd, extents, linearized)
%     patchSize = 50; %userParams_uc.patchSize; %30; %
    stepSize = 0.01; %userParams_uc.stepSize; %0.01;%
    dir1 = explorationDirection(:,1)';
    dir2 = explorationDirection(:,2)';  

    range_dir1 = -extents(1):stepSize:extents(1);
    range_dir2 = -extents(2):stepSize:extents(2);
    
    expandedPoints = zeros(length(range_dir1)*length(range_dir2), rD+rA);
    count = 1;
    
    centerDes = centerPt(1:rD+rA);

    %% Expand from the center point in the calculated directions (assumed to be 2)
    for dx1 = range_dir1
        %dx1 = ii * dir1 * stepSize;
        for dx2 = range_dir2
            %dx2 = jj * dir2 * stepSize;
            offsetPt = centerDes + dx1*dir1 + dx2*dir2;
            
            expandedPoints(count, :) = offsetPt;
            count = count + 1;
        end
    end

    % Borrowed from expandPatch
    validation = mFunc_uc.validate(expandedPoints);
    validIndices = find(validation);
    
    if linearized
        patchPts = evalLinearized(mFunc_uc, centerPt, expandedPoints, rD, rd, rA, extents);
    else
        evaluatedPoints = zeros(size(expandedPoints,1), rd);
        if mFunc_uc.id == 26 %% id where x and z are separated
            evaluatedPoints(validIndices,:) = mFunc_uc.eval(expandedPoints(validIndices,1:rD), ...
                                                        expandedPoints(validIndices,rD+1:rD+rA)); 
        else
            % the z parameter (0 here) doesn't mean anything; not used, just required
            evaluatedPoints(validIndices,:) = mFunc_uc.eval(expandedPoints(validIndices,:), 0);
        end
        patchPts = [expandedPoints evaluatedPoints];
    end

    patchPts = patchPts(validIndices, :);
end


function [linearPatchPts] = evalLinearized(mFunc_uc, centerPt, patchPts, rD, rd, rA, extents)
    DF = jacobian(mFunc_uc.f);
    numPts = size(patchPts, 1);
    linearPatchPts = zeros(numPts, rD+rA+rd);
    
    % evaluate the jacobian on the center point
    DFevald = evalAtPt(DF, centerPt(1:rD), centerPt(rD+1:rD+rA));
    
    % compute other points (inefficient, for visualization only)
    x0 = centerPt(1:rD+rA);    % design and app vars
    F0 = centerPt(rD+rA+1:end);           % F(x,z)
    for i=1:numPts
        x = patchPts(i, 1:rD+rA);  % get design and application vars
        f_x = F0' + DFevald * (x - x0)';
        linearPatchPts(i, :) = [x, f_x'];
    end
    % return, then plot in the intermediate vis to see how well we do 
end

