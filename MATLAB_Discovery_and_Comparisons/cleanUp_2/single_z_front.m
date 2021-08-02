function [buffer, figures, mFunc, hv_fixed_context, n_evals_fixed_context, n_scalar_evals_fixed_context] = single_z_front(functionID, z, constrain_z) 
    %clearvars -except functionID z;
    %% We construct the problem statement as a mapping function which deals with evaluations of the appropriate objective functions
    addPaths();
    mFunc = MappingFunction(functionID, z);

    % ADDED -- sample the function to get a feel for the evaluated range
    % may need this to determine scaling factor
%     massSamples = mFunc.sample(100000, z);
%     [rD,rd] = mFunc.getDimensions();
%     figure; hold on;
%     scatter(massSamples(:,rD+1),massSamples(:,rD+2));
%     figure; hold on;
%     scatter3(massSamples(:,1),massSamples(:,2),massSamples(:,3));
%     pause;
    
    %% Initialize parameters fixed by the problem statement
    imageTitle = mFunc.getImageTitle();
    [rD,rd] = mFunc.getDimensions();
    numCells = mFunc.getNumCells(); 
    [minBound, maxBound] = mFunc.getBounds();
    fixedParams = struct('rD', rD, 'rd', rd, 'imageTitle', imageTitle, 'minBound', minBound, ...
                         'maxBound', maxBound, 'constrain_z', constrain_z);
       
    %% Initialize user-adjustable parameters
    userParams = struct('tolerance', 0.01, ...
                        'directionDelta', 0.3, ...
                        'numCells', numCells, ...
                        'patchSize', 1*numCells,...
                        'stepSize', 1/numCells, ...
                        'numInitialSamples', 10, ...
                        'numNewSamples', 10, ...
                        'numRepeatedNoImprovement', 3, ...
                        'breakScoreRelative', 0.1/numCells, ...
                        'breakScoreDistance', 0.0001, ...
                        'neighForScore', ceil(0.05*numCells),...
                        'neighForOpt', ceil(0.2*numCells), ...
                        'nMaxRuns', 15);
    
    %% Initialize buffer data structure
    buffer = Buffer(userParams, fixedParams);
    
    %% Initialize figures class for display
    figures = showFigures(rD, rd, buffer, functionID);

    %% Initialize a random population of samples
    randomSamples = mFunc.sample(userParams.numInitialSamples, z);
%     massSamples = mFunc.sample(10000, z);
%     figure;
%     scatter(massSamples(:,4), massSamples(:,5));

    %% Start main loop
    disp('Starting main loop');
    
    patchID = 0;
    
    for run=1:userParams.nMaxRuns
        disp('----------------------------------------------------------------------');
        disp(['Iteration ' num2str(run)]);
        
        %% Step 1: Add each random point to the buffer first in case we get lucky
        
        % NOTE: This is for experimentation with 'z' as a design variable.
        [numRandomSamples, ~] = size(randomSamples);
%         randomSamples(:,rD) = ones(numRandomSamples,1)*z;

        for i = 1:size(randomSamples,1)
            buffer.addPoint(randomSamples(i,:));
        end
        

        %% Step 2: Locally optimize current samples
        
        % Get target directions for our random samples
        [targets, dirs] = buffer.getRandomDirections(randomSamples(:,rD+1:rD+rd));
        
        % Attempt to push each sample at its corresponding target point
        optimizedSamples = optimize(randomSamples, targets, dirs,  mFunc, fixedParams, z);
        
        disp('Finished local optimizations');

        %% Step 3: Perform local exploration on optimized points
        for v = 1:size(optimizedSamples,1)
            centerPoint = optimizedSamples(v,:);
            % Get the direction(s) around which to locally explore each sample
            explorationDirection = getExplorationDirection(centerPoint(1:rD),mFunc, fixedParams);
            % Perform the exploration
            [patch, validation] = expandPatch(explorationDirection, centerPoint, mFunc, buffer, userParams, false, z);

            % Add the patch to the buffer if the former is nonempty
            if (size(patch,1) > 0)
                patchID = patchID +1;
                buffer.addPatch(patch, patchID, validation, explorationDirection, centerPoint);
            end
        end
        
        disp('Finished expanding patches');

        %% Step 4: Check stopping criteria
        terminate = buffer.terminate(run);
        if terminate %|| mFunc.numFuncEvals > max_eval_constraints
            disp('Solution converged!')
            break;
        end

        %% Step 5: Randomly sample around buffer cells for next iteration
        [randomPoints] = buffer.getPointsToImprove(userParams.numNewSamples, mFunc, z);
        randomSamples = [randomPoints mFunc.eval(randomPoints, z)];
    end
    hv_fixed_context = buffer.getHyperVolume();
    n_evals_fixed_context = mFunc.numFuncEvals;
    n_scalar_evals_fixed_context = mFunc.numScalarizationFuncEvals;
    %% Post processing
    buffer.labelPatches(patchID + 1);
    buffer.fillEmptyCells(userParams.stepSize, mFunc, z);
    disp('Finished labeling patches');
end