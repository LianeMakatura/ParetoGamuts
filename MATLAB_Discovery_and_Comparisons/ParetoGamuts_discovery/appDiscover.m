function [stats, buffer, des_figh, perf_figh, pts, labels, filledPtsArray, idxToVis] = appDiscover(problemID, visualize)
    stats = struct();
    totalTime_start = tic;
    addPaths;
    dispstat('','init');
    
    if nargin < 2
        visualize = true;
    end
    
    % path to current directory (ParetoGamuts_discovery/), 
    % IMPT: make sure to include final /
    basewd = "<working directory including final />";
    if strcmp(basewd, "<working directory including final />")
        disp("FATAL (appDiscover.m): placeholder working directory detected. Please update.")
    end
    
    % for meshBased problems only
    prepJSvis = false;
    mesh_provider = "jscad";

    %% We construct the problem statement as a mapping function which deals with evaluations of the appropriate objective functions
    problemSetupTime_start = tic;
    mFunc = AppMappingFunction(problemID, basewd);

    %% Initialize parameters fixed by the problem statement
    imageTitle = mFunc.getImageTitle();
    [rD,rd,rA] = mFunc.getDimensions();
    numCells = mFunc.getNumCells(); 
    [minBound, maxBound] = mFunc.getBounds();
    fixedParams = struct('rD', rD, ...
                         'rA', rA, ...
                         'rd', rd, ...
                         'varNames', mFunc.varNames, ...
                         'imageTitle', imageTitle, ...
                         'minBound', minBound, ...
                         'maxBound', maxBound);

    %% Initialize user-adjustable parameters
    % num(Initial/New)Samples
    if mFunc.functionID == 12
        stepSize = 0.1;
    elseif rd + rA < 4 
        stepSize = 0.01;
    else
        stepSize = 0.03; % 0.1, reduced for better resolution lamp example but takes a while
    end
    
    if rd + rA == 3
        fillHolesBufferProjection = true;
        fillHolesGamutExtraction = true;
        useExactHV = true;
    else
        fillHolesBufferProjection = false;
        fillHolesGamutExtraction = false;
        useExactHV = false;
    end
    
    userParams = struct('tolerance', 0.01, ...
                        'directionDelta', 0.3, ... 
                        'numCells', numCells, ...
                        'stepSize', stepSize, ...
                        'fd_eps', 1e-3, ... 
                        'fd_type', 'central', ...
                        'numInitialSamples', 10, ...
                        'numNewSamples', 10, ...
                        'numRepeatedNoImprovement', 3, ...
                        'improvementTol', 0.001, ...
                        'fillHoles_bufferProjection', fillHolesBufferProjection, ...
                        'fillHoles_gamutExtraction', fillHolesGamutExtraction,...
                        'gamutExtractCellPercent', 0.8,...
                        'useHV', useExactHV, ...
                        'nMaxRuns', 50);
    stats.time_problemSetup = toc(problemSetupTime_start);


    %% Initialize buffer data structure + random samples
    % samples in (rD + rA) space (and their image is in (rd+rA) space)
    discoveryTime_start = tic;
    buffer = AppBufferv3(mFunc, userParams, fixedParams);
    randomSamples = mFunc.sample(userParams.numInitialSamples);
    
    %% Start main loop
    disp('Starting main loop');
    patchID = 0;
    for run=1:userParams.nMaxRuns
        disp('----------------------------------------------------------------------');
        disp(['Iteration ' num2str(run)]);

        %% Step 1: Locally optimize current samples (within single context)
        samplesPerf = randomSamples(:,rD+rA+1:rD+rA+rd);
        [targets, dirs] = buffer.getRandomDirections(samplesPerf);
        [optimizedSamples, optSampAlphas, optSampBetas] ...
            = appOptimize(randomSamples, targets, dirs, mFunc, fixedParams, userParams);

        disp('Finished local optimizations');

        %% Step 2: Perform local exploration on optimized points
        dispstat('Exploring around optimized points...','keepthis','keepprev','timestamp'); 
        for v = 1:size(optimizedSamples,1)
            % update progress bar - this may remove a bit of the previous
            % output if fprintf or disp are used inside this loop...
            % generally not an issue, but change to dispstat(...,'keepprev', 'keepthis') if needed (either this call, or the disp/fprintf in between) 
%             dispstat(sprintf('Processing patch %d / %d (%0.1f%% complete)',v, size(optimizedSamples,1), v / size(optimizedSamples,1) * 100), 'timestamp');
            
            centerPt = optimizedSamples(v,:);
            centerAlpha = optSampAlphas(v, :)';
            centerBeta = optSampBetas(v, :)';
            
            % Get the direction(s) around which to locally explore each sample
            xVals = centerPt(1:rD);
            zVals = centerPt(rD+1:rD+rA);
            [expDirsOrig, numActiveC] ...
                = appGetExplorationDirections_xz(mFunc, ...
                        xVals, zVals, centerAlpha, centerBeta, ...
                        rD, rd, rA, userParams);
                    
            % Separate x'z' from alpha' (and numActiveC beta', if constrained)
%             alphaPrimes = expDirsOrig(1:rd, :); % not used currently, here for reference
%             betaPrimes = expDirsOrig(rd+1:rd+numActiveC, :); % not used currently, here for reference
            expDirs1 = expDirsOrig(rd +numActiveC+1:rd+numActiveC+rD+rA, :);
            [expDirs, degenerate] = canonizeDirections(expDirs1, rd, rD, rA);
            if degenerate
                disp('Canonized directions are degenerate');
                continue;
            end
                             
            % find amount to expand along each design space direction
            [des_extents] = appGetExtents(centerPt, expDirs, mFunc, rD, rd, rA, userParams);
            
            % Perform the exploration, reassign expDirs (can be modified)
            [validPatchPts, connectivityList, expDirs] = appExpandPatch(centerPt, expDirs, ...
                                            mFunc, buffer, des_extents, userParams);
                        
            % Add the patch to the buffer if the former is nonempty
            if (size(validPatchPts,1) > 3)
                patchID = patchID +1;
                buffer.addPatch(patchID, centerPt, expDirs, ... 
                            des_extents, validPatchPts, connectivityList, ...
                            userParams.fillHoles_bufferProjection); % may need to change how extents are passed
            end
             
        end
        fprintf('Finished expanding patches.\n\n');

        %% Step 3: Check stopping criteria
        [stop, converged, percentTargStability] = buffer.terminate(run, userParams.useHV);
        if stop
            if converged
                disp('Solution converged!')
            else
                disp('Maximum number of iterations or function evaluations reached.')
            end
            break;
        end
        
        %% Step 4: Randomly sample around buffer cells for next iteration
        if percentTargStability > 0.3 %we've had at least half the required number of iters without improvement
            randomPoints = mFunc.sample(userParams.numNewSamples);
        else
            % Could modify
            [randomPoints] = buffer.getPointsToImprove(userParams.numNewSamples);% 11 for bicopter, 13 for camera
        end
        randPtsDes = randomPoints(:, 1:rD);
        randPtsApp = randomPoints(:, rD+1:rD+rA);
        randPtsPerf = mFunc.eval(randPtsDes, randPtsApp);
        randomSamples = [randPtsDes,  randPtsApp, randPtsPerf];
    end
    %% End main loop
    stats.time_discovery = toc(discoveryTime_start);
    fprintf('Time to discover: %s\n', durationString(stats.time_discovery));

    %% Step 5: Post processing + stats update
    idxToVis = struct('designIndex1', 1, ...
                      'designIndex2', 2, ...
                      'appIndex', 1, ...
                      'perfIndex1', 1, ...
                      'perfIndex2', 2);

    disp('Extracting Pareto Gamut')
    gamutExtractionTime_start = tic;
    [pts, labels, filledPtsArray] = buffer.getParetoGamut();
    stats.time_gamutExtraction = toc(gamutExtractionTime_start);
    fprintf('Time to extract gamut: %s\n', durationString(stats.time_gamutExtraction));
    
    stats.time_totalAlg_noVis = toc(totalTime_start);
    
    stats.numFuncEvals =  mFunc.numFuncEvals;
    stats.numJacobianEvals = mFunc.numJacEvals;
    stats.numHessianEvals = mFunc.numHessEvals;
    stats.numScalarizationFuncEvals = mFunc.numScalarizationFuncEvals;
    stats.hypervolumeHistory_fixedContext = buffer.HV_history;
    stats.hypervolumeFinal_fixedContext = buffer.HV_history(end,:);
    stats.hypervolumeFinal_gamutNorm = norm(stats.hypervolumeFinal_fixedContext, 2); %norm of final FC HV vector. Our notion of gamut HV
    stats.numIters = run;
    stats.converged = converged;
    stats.numPtsOnGamut = size(pts, 1);
    stats.numPatches = size(buffer.patchArray);
    
    %% visualization prep
    des_figh = 0;
    perf_figh = 0;
    if visualize
        [des_figh, perf_figh] = buffer.visualizeAppBuff(idxToVis, false, pts, labels);
    end
    
%     if mFunc.hasGroundTruth() && visualize
%         [~, ~] = buffer.visualizeAppBuff(idxToVis, true, pts, labels); % visualize ground truth only
%     end

    if prepJSvis && mFunc.hasMeshes()
        mesh_radius = 10; %should be in each mFunc

        prepForInteractiveVis(mFunc, pts, labels, mesh_provider, mesh_radius);
        disp('Meshes and fronts ready for visualization. Start ParetoExplorer/vis in VSCode, and GoLive to view.')
    end

    stats.time_totalForAppDiscover = toc(totalTime_start);
end




