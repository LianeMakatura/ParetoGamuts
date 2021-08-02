rng(4);
close all;
functionID = 23; % Unconstrained function ID
z = 0.5;
% [buffer,figures, mFunc] = single_z_front(functionID, z, false);

expandedCenters = buffer.expandedCenters;
evaluatedCenters = mFunc.eval(expandedCenters, z);

paretoIndices = buffer.getParetoInd;
paretoPoints = [];
for i = 1:length(paretoIndices)
    idx = paretoIndices(i);
    entry = buffer.Buff(idx);
    patch = entry.bestPatch;
    
    if patch > 0 
        patchKey = num2str(patch);

        newPoint = [entry.minD entry.minF];
        paretoPoints = [paretoPoints; newPoint];
    end
end

randPt = paretoPoints(200,:);
randPt = GP;
expansionParams = struct('rD', 3, 'rd', 2, 'constrain_z', false);
explorationDirection = getExplorationDirection(randPt(1:3), mFunc, expansionParams)
explorationDirection = explorationDirection(:,1);
% explorationDirection = explorationDirection(:,1);
expPts = [];
for j = -50:50
    dx = j*0.01*explorationDirection;
    expPts = [expPts; randPt(1:3) + dx'];
end

evalPts = mFunc.eval(paretoPoints(:,1:3));
% figures.showLabeledParetoSet(true); hold on;
% scatter3(expandedCenters(:,1), expandedCenters(:,2), expandedCenters(:,3), 100, 'black')
figure; hold on;
scatter3(paretoPoints(:,1), paretoPoints(:,2), paretoPoints(:,3), 25, 'green')
scatter3(randPt(1),randPt(2),randPt(3), 100, 'blue');
scatter3(expPts(:,1),expPts(:,2),expPts(:,3), 50, 'red');

expPtsEval = mFunc.eval(expPts);
% figures.showLabeledParetoFront(true);
figure; hold on;
% scatter(expandedCenters(:,1), evaluatedCenters(:,2), 100, 'black')
scatter(paretoPoints(:,4), paretoPoints(:,5), 100, 'green')
% scatter(evalPts(:,1), evalPts(:,2), 'blue');
scatter(expPtsEval(:,1), expPtsEval(:,2), 'blue');
