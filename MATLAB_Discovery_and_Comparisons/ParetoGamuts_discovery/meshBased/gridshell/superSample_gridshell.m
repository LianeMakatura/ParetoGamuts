
%% Set the parameters
NX = 27; % total number of points along x, including fixed ends
NY = 27; % total number of points along y, including fixed ends

numPtsToTest = 1000000;


%% Get random inputs of correct dimensions
numVariablePts = (NX-2)*(NY-2); %endpoints fixed
z = rand(numPtsToTest, numVariablePts); %normalized between 0, 1
theta = rand(numPtsToTest, 1);
oursIn = mat2cell(z, numPtsToTest, ones(1, numVariablePts));


%% Compute results and find min/max for scaling
% c = Gridshell(pwd, N, N, [1,2]);
c = Gridshell_nonsymb(pwd, NX, NY, [1,2]);

morningF = matlabFunction(c.performanceMetric1());
oursMorning = morningF(oursIn{:}, theta);
fprintf("Morning min: %0.6f,  max: %0.6f\n", min(oursMorning), max(oursMorning));

eveningF = matlabFunction(c.performanceMetric2());
oursEvening = eveningF(oursIn{:}, theta);
fprintf("Evening min: %0.6f,  max: %0.6f\n", min(oursEvening), max(oursEvening));

smoothF = matlabFunction(c.smoothness());
smoothness = smoothF(oursIn{:});
fprintf("smoothness min: %0.6f,  max: %0.6f\n", min(smoothness), max(smoothness));

powerF = matlabFunction(c.powerOutput(c.morningSun));
power = powerF(oursIn{:}, theta);
fprintf("power min: %0.6f,  max: %0.6f\n", min(power), max(power));
