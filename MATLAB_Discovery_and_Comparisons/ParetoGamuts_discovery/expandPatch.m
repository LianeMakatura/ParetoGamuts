function [ allPoints, validation ] = expandPatch( direction, centerPoint, mFunc, buffer, userParams, alwaysExpand, z)
[rD, rd] = mFunc.getDimensions();
allPoints = [];
validation = [];
stepSize = userParams.stepSize;
patchSize = userParams.patchSize;

%% Degenerate Case 1: More than d-1 directions to explore 
if size(direction,2) >= rd 
    % Pick d-1 random directions from the selection
    indices = randsample(size(direction,2),rd-1);
    direction = direction(:,indices);
end

%% Degenerate Case 2: Fewer than d-1 directions to explore
if size(direction,2) < rd-1
    % Do not expand the point
    disp('degenerate')
    return
end

%% Do not explore if we have already expanded the same patch
if ~alwaysExpand 
    if buffer.checkHaveExpanded(centerPoint(1:rD),direction)
        disp('found already')
        return
    end
end

if rd == 2
    searchSpace = linspace(-patchSize, patchSize, 2*patchSize+1);
    displacements = direction * stepSize * searchSpace;
    originalSample = repmat(centerPoint(1:rD)',1,2*patchSize+1);
    newSamples = (originalSample + displacements)';
else
    stepSize = stepSize * 2;
    searchSpace = linspace(-patchSize, patchSize, 2*patchSize+1);
    displacements = direction(:,1) * stepSize * searchSpace;

    originalSample = repmat(centerPoint(1:rD)',1,2*patchSize+1);
    firstSamples = (originalSample + displacements)';
    newSamples = zeros((patchSize*2+1)^2,rD);
    % Might be a cleaner way of doing this with a mesh grid or something,
    % but basically computes the samples in the first direction, then
    % computes the samples in the other direction
    for l = 1:patchSize*2 + 1
       samp = firstSamples(l,:);
       delta = direction(:,2) * stepSize * searchSpace;
       nv = samp + delta';
       newSamples((l-1)*(patchSize*2+1)+1:l*(patchSize*2+1),:) = nv;
    end    
    
    
end

if rd == 3
    validation = mFunc.validate(newSamples);
    evaluatedPoints = mFunc.eval(newSamples, z);
    allPoints = [newSamples evaluatedPoints];
else
    validation = mFunc.validate(newSamples);

    valPoints = find(validation);

    evaluatedPoints = zeros(size(newSamples,1), rd);

    evaluatedPoints(valPoints,:) = mFunc.eval(newSamples(valPoints,:), z);

    allPoints = [newSamples evaluatedPoints];
    allPoints = allPoints(valPoints, :);
    validation = validation(valPoints, :);
end
%% Reshape the patch to a rectangle for the 3-d case
if rd == 3
    allPoints = reshape(allPoints, patchSize*2+1,patchSize*2+1,rd+rD);
    validation = reshape(validation, patchSize*2+1,patchSize*2+1);
end

end
