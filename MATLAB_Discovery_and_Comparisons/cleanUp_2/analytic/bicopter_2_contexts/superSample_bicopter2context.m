
%% Set the parameters
anchorPtsPath = './analytic/bicopter/bicopter_anchor.mat';
numPtsToTest = 1000000;
N = 16; % number of timesteps (set inside bicopter)


%% Get random inputs of correct dimensions
numVariablePts = N*2; % actuation for each prop, every time step
despts = rand(numPtsToTest, numVariablePts); %normalized between 0, 1
length = rand(numPtsToTest, 1)';
density = rand(numPtsToTest, 1)';
oursIn = despts';


%% Compute results and find min/max for scaling
c = bicopter_2_contexts(pwd);

[distanceF, ~, ~] = c.performanceMetric1(); % already a matlab function
oursDist = distanceF(oursIn, [length;density]);
fprintf("Distance min: %0.6f,  max: %0.6f\n", min(oursDist), max(oursDist));

[energyF, ~, ~] = c.performanceMetric2(); % already a matlab function
oursEnergy = energyF(oursIn, [length; density]);
fprintf("Energy min: %0.6f,  max: %0.6f\n", min(oursEnergy), max(oursEnergy));
