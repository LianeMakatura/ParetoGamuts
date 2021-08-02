close all;
clear all;
rng(4);
addPaths;
warning('off','all');

%% Specify the problem statement
functionID = 1;
%%% Function ID to Problem Key:
% 1: ZDT1
% 2: ZDT2
% 3: ZDT3
% 4: ZDT4
% 5: ZDT6
% 6: DTLZ1
% 7: DTLZ2
% 8: DTLZ3
% 9: Sum of sines 2d
% 10: Sum of sines 3d
% 11: Kursawe
% 12: CAD Lamp %% haven't tested yet
% 13: CAD Brake Hub
% 14: CAD Wrench
%%%

%% We construct the problem statement as a mapping function which deals with evaluations of the appropriate objective functions
mFunc = MappingFunction(functionID);

%% Initialize parameters fixed by the problem statement
imageTitle = mFunc.getImageTitle();
[rD,rd] = mFunc.getDimensions();
numCells = mFunc.getNumCells(); 
[minBound, maxBound] = mFunc.getBounds();
fixedParams = struct('rD', rD, 'rd', rd, 'imageTitle', imageTitle, 'minBound', minBound, ...
                     'maxBound', maxBound);

%% Initialize user-adjustable parameters
userParams = struct('tolerance', 0.01, 'directionDelta', 0.3, 'numCells', numCells, 'patchSize', 0.5*numCells,...
                    'stepSize', 1/numCells, 'numInitialSamples', 10, 'numNewSamples', 10, ...
                    'numRepeatedNoImprovement', 3, 'breakScoreRelative', 0.1/numCells, 'breakScoreDistance', 0.0001, ...
                    'neighForScore', ceil(0.05*numCells),'neighForOpt', ceil(0.2*numCells), 'nMaxRuns', 250);

%% Initialize buffer data structure
buffer = Buffer(userParams, fixedParams);
%% Initialize figures class for display
figures = showFigures(rD, rd, buffer, functionID);

%% Initialize a random population of samples
randomSamples = mFunc.sample(userParams.numInitialSamples);
%% Start main loop
disp('Starting main loop');
patchID = 0;
for run=1:userParams.nMaxRuns
    disp('----------------------------------------------------------------------');
    disp(['Iteration ' num2str(run)]);
    %% Step 1: Add each random point to the buffer first in case we get lucky
    for i = 1:size(randomSamples,1)
        buffer.addPoint(randomSamples(i,:));
    end

    %% Step 2: Locally optimize current samples
    targets = buffer.getRandomDirections(randomSamples(:,rD+1:rD+rd));
    optimizedSamples = optimize(randomSamples, targets, mFunc, fixedParams);

    % Display the samples before & after local optimization
    % figures.showCurrentIteration(randomSamples, optimizedSamples)

    disp('Finished local optimizations');

    %% Step 3: Perform local exploration on optimized points
    for v = 1:size(optimizedSamples,1)
        centerPoint = optimizedSamples(v,:);
        % Get the direction(s) around which to locally explore each sample
        explorationDirection = getExplorationDirection(centerPoint(1:rD),mFunc, fixedParams);
        % Perform the exploration
        [patch, validation] = expandPatch(explorationDirection, centerPoint, mFunc, buffer, userParams, false);

        % Add the patch to the buffer if the former is nonempty
        if (size(patch,1) > 0)
            patchID = patchID +1;
            buffer.addPatch(patch, patchID, validation, explorationDirection, centerPoint);
        end
    end
    disp('Finished expanding patches');

    %% Step 4: Check stopping criteria
    if buffer.terminate(run)
        disp('Solution converged!')
        break;
    end

    %% Step 5: Randomly sample around buffer cells for next iteration
    [randomPoints] = buffer.getPointsToImprove(userParams.numNewSamples);
    randomSamples = [randomPoints mFunc.eval(randomPoints)];

end
%% End main loop

%% Step 6: Post processing
buffer.labelPatches(patchID + 1);
buffer.fillEmptyCells(userParams.stepSize, mFunc);
disp('Finished labeling patches');
figures.showLabeledParetoFront(true);
save([imageTitle '.mat']);
