classdef MappingFunction < handle
%% Instance variables
properties     
    % Visualization
    imageTitle
    % For Gamut's experimentation
    includeZ
    % For C++ 
    useCPP
    cppCalls     
    % Static constants
    minBound
    maxBound
    rD
    rd
    functionID
    % Function evaluation variables
    f
    
    % not in all
    f1
    f2
    rA
    id
    
    fun
    deriv
    hess
    fhSymb
    symFreeAnalytic
    numFuncEvals
    app_var
    numScalarizationFuncEvals
end
%% Functions
methods        
    %% Constructor.
    function obj = MappingFunction(id, z) % z only necessary sometimes
        if nargin==0 % needed for Object arrays
            return
        end
        
        obj.numFuncEvals = 0;
        obj.numScalarizationFuncEvals = 0;
        obj.useCPP = false;
        obj.includeZ = false;
        obj.functionID = id;
        obj.imageTitle = 'results\foo';
        obj.symFreeAnalytic = false;
        
        basewd = "<working directory including final />";
        if strcmp(basewd, "<working directory including final />")
            disp("FATAL (MappingFunction.m): placeholder working directory detected. Please update.")
        end

        % Objective function definition
        
        % ZDT1 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt1/index.php
        if id == 1
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = ZDT1_b(rD);
            [fun2, deriv2, f2] = ZDT1_a(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.imageTitle = 'results\ZDT1';
        % ZDT2 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt2/index.php
        elseif id == 2
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = ZDT1_b(rD);
            [fun2, deriv2, f2] = ZDT2_a(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.imageTitle = 'results\ZDT2';
        % ZDT3 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt3/index.php
        elseif id == 3
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = ZDT1_b(rD);
            [fun2, deriv2, f2] = ZDT3_a(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.imageTitle = 'results\ZDT3';
        % ZDT4 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt4/index.php
        elseif id == 4
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = ZDT1_b(rD);
            [fun2, deriv2, f2] = ZDT4_a(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.imageTitle = 'results\ZDT4';
        % ZDT6 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/zdt6/index.php
        elseif id == 5
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = ZDT6_a(rD);
            [fun2, deriv2, f2] = ZDT6_b(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.imageTitle = 'results\ZDT5';
        % DTLZ1 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/dtlz1/      
        elseif id == 6
            rD = 6;
            obj.rD = rD;

            [fun1, deriv1, f1] = DTLZ1a(rD);
            [fun2, deriv2, f2] = DTLZ1b(rD);
            [fun3, deriv3, f3] = DTLZ1c(rD);

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
        % DTLZ2 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/dtlz2/      
        elseif id == 7
            rD = 6;
            obj.rD = rD;

            [fun1, deriv1, f1] = DTLZ2a(rD);
            [fun2, deriv2, f2] = DTLZ2b(rD);
            [fun3, deriv3, f3] = DTLZ2c(rD);

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
        % DTLZ3 Function - https://sop.tik.ee.ethz.ch/download/supplementary/testproblems/dtlz3/
        elseif id == 8
            rD = 6;
            obj.rD = rD;

            [fun1, deriv1, f1] = DTLZ3a(rD);
            [fun2, deriv2, f2] = DTLZ3b(rD);
            [fun3, deriv3, f3] = DTLZ3c(rD);

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
        % 2d function formulated as linear combinations of sinusoids with
        % randomized coefficients
        elseif id == 9
            rng(2)
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = sumSins(6,rD, -rD, rD);
            [fun2, deriv2, f2] = sumSins(6,rD, -rD, rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
        % 3d function formulated as linear combinations of sinusoids with
        % randomized coefficients
        elseif id == 10
            rng(3)
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = sumSins(3,rD, -rD, rD);
            [fun2, deriv2, f2] = sumSins(3,rD, -rD, rD);
            [fun3, deriv3, f3] = sumSins(3,rD, -rD, rD);

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
        % Kursawe Function - http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.47.8050&rep=rep1&type=pdf
        elseif id == 11      % unconstrained z value
            rD = 3;
            obj.rD = rD;
            [fun1, deriv1, f1] = Schulz_Kursawe_1(rD);
            [fun2, deriv2, f2] = Schulz_Kursawe_2(rD);
%             [fun1, deriv1, f1] = Kursawe_1(rD);
%             [fun2, deriv2, f2] = Kursawe_2(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.f1 = f1;
            obj.f2 = f2;
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = false;
            
        elseif id == 24 % constrained z value
            rD = 3;
            obj.rD = rD;
            [fun1, deriv1, f1] = Schulz_Kursawe_1(rD, z);
            [fun2, deriv2, f2] = Schulz_Kursawe_2(rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.f1 = f1;
            obj.f2 = f2;
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            
        % 3d Parametric lamp example
        elseif id == 12
            addpath meshModels\lamp;
            [fun1, f1, deriv1, fun2, f2, deriv2, fun3, f3, deriv3] = parametricLamp_helper();

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
        % 3d Parametric Brake Hub Example - https://people.csail.mit.edu/aschulz/optCAD/a11-schulz.pdf
        elseif id == 13
            rD = 3;
            obj.rD = 3;
            obj.fun = zeros(1, 3); %fix this lol
            cpp = cppCaller(rD, 3, 2);
            obj.cppCalls = cpp;
            obj.useCPP = true;
            obj.imageTitle = 'results\BrakeHub';
        % 3d Parametric Wrench Example - https://people.csail.mit.edu/aschulz/optCAD/a11-schulz.pdf
        elseif id == 14
            rD = 3;
            obj.rD = 3;
            obj.fun = zeros(1, 3); %fix this lol
            cpp = cppCaller(rD, 3, 1);
            obj.cppCalls = cpp;
            obj.useCPP = true;
            obj.imageTitle = 'results\Wrench3D';
        % Gamut-wrapped ZDT1 Function
        elseif id == 15
            rD = 3;
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_ZDT1_b(rD, z);
            [fun2, deriv2, f2] = Schulz_ZDT1_a(rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_ZDT1'; 
        % Gamut-wrapped ZDT2 Function
        elseif id == 16
            rD = 3;
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_ZDT1_b(rD, z);
            [fun2, deriv2, f2] = Schulz_ZDT2_a(rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_ZDT2';
        % Gamut-wrapped ZDT3 Function
        elseif id == 17
            rD = 3;
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_ZDT1_b(rD, z);
            [fun2, deriv2, f2] = Schulz_ZDT3_a(rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_ZDT3';
        % Gamut-wrapped ZDT4 Function
        elseif id == 18
            rD = 15;
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_ZDT1_b(rD, z);
            [fun2, deriv2, f2] = Schulz_ZDT4_a(rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_ZDT4'; 
        elseif id == 19
            rng(2)
            rD = 2;   % previously 15
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_sumSins(6,rD, -rD, rD, z);
            [fun2, deriv2, f2] = Schulz_sumSins(6,rD, -rD, rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_sumSins_19';
        elseif id == 20
            rng(2)    % should result in same random coeffs
            rD = 3;   % previously 15
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_sumSins(6,rD, -rD, rD);    % dont pass z --> z unconstrained as design var
            [fun2, deriv2, f2] = Schulz_sumSins(6,rD, -rD, rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_sumSins_20';
        elseif id == 21       %% Schulz's rD+1 version of ID #9 above -- constrained z value
            rng(3)    % should result in same random coeffs
            rD = 3;   % previously 15
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_sumSins_fancyZ(6,rD, -rD, rD, z);
            [fun2, deriv2, f2] = Schulz_sumSins_fancyZ(6,rD, -rD, rD, z);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_sumSins_21';
            obj.id = 21;
        elseif id == 22       %% Schulz's unconstrained version of id 21 -- part 2! - unconstrained z value, with augmented perf space
            rng(3)    % should result in same random coeffs
            rD = 3;   % previously 15
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_sumSins_fancyZ(6,rD, -rD, rD);
            [fun2, deriv2, f2] = Schulz_sumSins_fancyZ(6,rD, -rD, rD);
            [fun3, deriv3, f3] = Schulz_appVarDim(rD);

            obj.fun = {fun1, fun2, fun3};
            obj.f = [f1, f2, f3];
            obj.deriv = [deriv1, deriv2, deriv3];
            obj.includeZ = true;
            obj.imageTitle = 'results\Schulz_sumSins_22'; 
            obj.id = 22;
        elseif id == 23                                         % unconstrained z values, non-augmented performance
            rng(3)    % should result in same random coeffs
            rD =3;   % previously 15
            obj.rD = rD;

            [fun1, deriv1, f1] = Schulz_sumSins_fancyZ(6,rD, -rD, rD);
            [fun2, deriv2, f2] = Schulz_sumSins_fancyZ(6,rD, -rD, rD);
%             [fun3, deriv3, f3] = Schulz_appVarDim(rD);

            obj.fun = {fun1, fun2};
            obj.f = [f1, f2];          
            
            obj.f1 = f1;
            obj.f2 = f2;
            
            obj.deriv = [deriv1, deriv2];
            obj.includeZ = false;
            obj.imageTitle = 'results\Schulz_sumSins_23';
            obj.id = 23;
        % id 24 is Kursawe above
        elseif id == 25     %% 2D L bracket - constrained z
            rD = 3; %for the 2d design case
            obj.rD = rD;
            obj.rA = 1;
            obj.id = 25;
            
            addpath meshBased/LBracket;
            addpath meshBased;
            addpath mathUtils;
            b = Schulz_LBracket(basewd, false, true, z); %3d, unconstrained
            
            [rfun, rf, rfuncsCell, rderiv, rhess] = b.createProblem(true); %write new files

            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            
            f1 = rfuncsCell{1};
            f2 = rfuncsCell{2};
            obj.f1 = rfuncsCell{1}; % need separate ones for quasi hessians
            obj.f2 = rfuncsCell{2};
            obj.includeZ = true;
            
        elseif id == 26     %% 2D L bracket - unconstrained z
            rD = 3; % for the 2d design case 
            obj.rD = rD;
            obj.rA = 1;
            obj.id = 26;
            
            addpath meshBased/LBracket;
            addpath meshBased;
            addpath mathUtils;
            b = Schulz_LBracket(basewd, false, false, z);%3d, unconstrained
            
            [rfun, rf, rfuncsCell, rderiv, rhess] = b.createProblem(true); %write new files

            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            
            f1 = rfuncsCell{1};
            f2 = rfuncsCell{2};
            obj.f1 = rfuncsCell{1}; % need separate ones for quasi hessians
            obj.f2 = rfuncsCell{2};
            obj.includeZ = false;
        elseif id == 27     %% 3D L bracket - constrained z
            rD = 4; 
            obj.rD = rD;
            obj.rA = 1;
            obj.id = 25; % CHANGE (there are some if statements that use this) 
            
            addpath meshBased/LBracket;
            addpath meshBased;
            addpath mathUtils;
            b = Schulz_LBracket(basewd, true, true, z); %3d, unconstrained
            
            [rfun, rf, rfuncsCell, rderiv, rhess] = b.createProblem(true); %write new files

            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            
            f1 = rfuncsCell{1};
            f2 = rfuncsCell{2};
            obj.f1 = rfuncsCell{1}; % need separate ones for quasi hessians
            obj.f2 = rfuncsCell{2};
            obj.includeZ = true;
            
            obj.deriv = [deriv1, deriv2];
        elseif id == 28     %% 3D L bracket - unconstrained z
            rD = 4;
            obj.rD = rD;
            obj.rA = 1;
            obj.id = 26; % CHANGE (there are some if statements that use this) 
            
            addpath meshBased/LBracket;
            addpath meshBased;
            addpath mathUtils;
            b = Schulz_LBracket(basewd, true, false, z); %3d, unconstrained
            
            [rfun, rf, rfuncsCell, rderiv, rhess] = b.createProblem(true); %write new files

            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            
            f1 = rfuncsCell{1};
            f2 = rfuncsCell{2};
            obj.f1 = rfuncsCell{1}; % need separate ones for quasi hessians
            obj.f2 = rfuncsCell{2};
            obj.includeZ = false;
            
        elseif id == 29    %% Lamp
            rD = 21;
            obj.rD = rD;
            perfMetricIDs = [1,2,3];
            addpath testFunctions/Schulz_Lamp;
            addpath mathUtils;
            b = Schulz_Lamp(basewd, perfMetricIDs, z); 

            [rfun, rf1, rf2, rf3, rderiv1, rderiv2, rderiv3] = b.createProblem(true);%write new files

            obj.fun = rfun;
            f1 = rf1;
            f2 = rf2;
            f3 = rf3;
            obj.f = [f1, f2, f3];
            obj.deriv = [rderiv1, rderiv2, rderiv3];
            obj.includeZ = true;
            
         elseif id == 30    %% Turbine
            rD = 3;
            obj.rD = rD;
            perfMetricIDs = [1,2];
            addpath testFunctions/Schulz_Turbine;
            addpath mathUtils;
            b = Schulz_Turbine(basewd, perfMetricIDs, z); 

            [rfun, rf, ~, rderiv, ~] = b.createProblem(true);%write new files

            obj.fun = rfun;
            f1 = rf{1};
            f2 = rf{2};
            obj.f = [f1, f2];
            obj.deriv = [rderiv{1}, rderiv{2}];
            obj.includeZ = true;
        elseif id == 31    %% BikeRocker
            rD = 3;
            obj.rD = rD;
            addpath testFunctions/Schulz_BikeRocker;
            addpath mathUtils;
            b = Schulz_BikeRocker(basewd, z); 

            [rfun, rf, ~, rderiv, ~] = b.createProblem(true);%write new files

            obj.fun = rfun;
            f1 = rf{1};
            f2 = rf{2};
            obj.f = [f1, f2];
            obj.deriv = [rderiv{1}, rderiv{2}];
            obj.includeZ = true;
            
        elseif id == 32    %% Bicopter
            b = bicopter(basewd);
            b.fixedz();
            b.set_app_var(z); % must set before createProblem()
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
            
        elseif id == 33    %% Gridshell
            metricsToUse = [1,2];
            numX = 12;
            numY = 12;
            addpath testFunctions/Schulz_GridShell;
            addpath mathUtils;
            b = Schulz_Gridshell_nonsymb(basewd, numX, numY, metricsToUse, z); 
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
            
        elseif id == 34    %% Gridshell
            metricsToUse = [1,2];
            numX = 15;
            numY = 15;
            addpath testFunctions/Schulz_GridShell;
            addpath mathUtils;
            b = Schulz_Gridshell_nonsymb(basewd, numX, numY, metricsToUse, z); 
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
        elseif id == 35    %% Gridshell
            metricsToUse = [1,2];
            numX = 20;
            numY = 20;
            addpath testFunctions/Schulz_GridShell;
            addpath mathUtils;
            b = Schulz_Gridshell_nonsymb(basewd, numX, numY, metricsToUse, z); 
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
        elseif id == 36    %% Gridshell
            metricsToUse = [1,2];
            numX = 25;
            numY = 25;
            addpath testFunctions/Schulz_GridShell;
            addpath mathUtils;
            b = Schulz_Gridshell_nonsymb(basewd, numX, numY, metricsToUse, z); 
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
            
        elseif id == 37    %% Gridshell
            metricsToUse = [1,2];
            numX = 27;
            numY = 27;
            addpath testFunctions/Schulz_GridShell;
            addpath mathUtils;
            b = Schulz_Gridshell_nonsymb(basewd, numX, numY, metricsToUse, z); 
            obj.app_var = z;
            
            [rfun, rf, ~, rderiv, rhess] = b.createProblem(true); % write new files
            rD = b.rD;
            obj.rD = rD;
            obj.fun = rfun;
            obj.f = rf;
            obj.deriv = rderiv;
            obj.hess = rhess;
            obj.includeZ = true;
            obj.symFreeAnalytic = true;
            
            
        else
            disp('Invalid function ID');
            pause;
        end
        
        rd = length(obj.fun);
        obj.rd = rd;

        obj.minBound = zeros(1,rD);
        obj.maxBound = ones(1,rD);

        % Generate symbolic design variables for derivative/hessian evaluation 
        if obj.includeZ == true
            % original 2 lines
            xs = sym('x', [1 obj.rD]);
            q = xs.';
        else
            xs = sym('x', [1 obj.rD]);
            zs = sym('z', [1 1]); %% to fix for more than 1 app var
            q = [xs.'; zs.'];
        end
            

        if (~obj.useCPP && ~obj.symFreeAnalytic)
            % Compute each hessian symbolically
            if rd == 2
                 hessian_1 = hessian(f1,q);
                 hessian_2 = hessian(f2,q); 
                 obj.hess = {hessian_1, hessian_2};
            elseif rd == 3
                 hessian_1 = hessian(f1,q);
                 hessian_2 = hessian(f2,q);         
                 hessian_3 = hessian(f3,q);
                 obj.hess = {hessian_1, hessian_2, hessian_3};
            else
                 obj.hess = {};
            end

            % Generate symbolic target variables for scalarization function 
            targets = sym('t', [1 obj.rd]);
            t = targets.';

            % Scalarization function: |f(x) - target(x)|^2
            if obj.includeZ % only need these for constrained z
                scalarFSymb(xs, targets) = sum((obj.f-targets).^2) / 2.0;
                gradFSymb = jacobian(scalarFSymb, q).'; % column gradf

                % Convert to Matlab function handle for easier use later
                obj.fhSymb = matlabFunction(scalarFSymb, gradFSymb,'Vars',{q, t});
            end
        end
    end
      
    %% Retrieve the image title for saved figures/data
    function imageTitle = getImageTitle(obj)
        imageTitle = obj.imageTitle;
    end

    %% Return the buffer size across one dimension
    % return numCells: 200x200x200 for 3-d cases, and 1000x1000 for 2-d
    % cases
    function numCells = getNumCells(obj)
        if obj.functionID == 29
            numCells = 50;
        elseif obj.functionID == 32
            numCells = 200;
        elseif obj.rd == 3
            numCells = 200;
        else
            numCells = 1000;
        end
    end
      
    %% Return the dimensionality of the problem
    % return rD: Dimensionality of design space
    % return rd: Dimensionality of objective space
    function [rD, rd] = getDimensions(obj)
        rD = obj.rD;
        rd = obj.rd;
    end
      
    %% Return the lower and upper bounds for design space
    %
    % return minBound: 1 x D vector representing the minimum allowed
    % values for the D-dimensions in design space
    % return minBound: 1 x D vector representing the maximum allowed
    % values for the D-dimensions in design space      
    function [minBound, maxBound] = getBounds(obj)
        minBound = obj.minBound;
        maxBound = obj.maxBound;
    end

    %% Return a function handle for local optimization
    function fh = getFunctionHandle(obj, target)
        if obj.useCPP
            fh = obj.cppCalls.getOptimizationHandle(target);
        elseif obj.symFreeAnalytic
            fh = @(x)obj.nonSymb_scalarF(x, target);
        else
            fh = @(q)obj.fhSymb(q', target');
        end
    end
    
    function [f, gradf] = nonSymb_scalarF(obj, pDes, target)
                        % hessian is given in the create_problem()
        pApp = obj.app_var;
        % Scalarization function: 1/2 * |f(x) - target|^2
        Dfull = obj.deriv(pDes.', pApp.');
        f = 0;
        gradf = 0;
        for i = 1:obj.rd
            F_i = obj.fun{i};
            alpha_i = F_i(pDes', pApp') - target(i);

            scalarizedFi = alpha_i.^2;
            f = f + scalarizedFi;

            scalarizedDFi = alpha_i * Dfull(i, :);
            gradf = gradf + scalarizedDFi;
        end
        if size(gradf, 2) ~= obj.rD
            gradf = gradf(:, 1:obj.rD);
        end
        f = f/2;
        gradf = gradf';
    end
      
    %% Validates a list of samples in design space are within the bounds of the objective functions
    %
    % @param points: an n x D list of n points
    % ---
    % @return V: an n x 1 binary vector where a 1 at index i denotes that
    % points(i,:) is a valid point in design space, and a 0 denotes that 
    % it is invalid.
    function V = validate(obj, points)
        V = ones(size(points,1), 1);
        [XMax, ~] = find(points > 1);
        [XMin, ~] = find(points < 0 );
        V(XMax) = 0;
        V(XMin) = 0;
    end
     
    %% Generate nSamples random samples
    function S = sample(obj, nSamples, z)
        P = obj.minBound + rand(nSamples, obj.rD).*(obj.maxBound-obj.minBound);
        S = [P obj.eval(P, z)];
    end
 
    %% Evaluate a list of points in design space
    % @param points: an n x D list of n points in design space
    % ---
    % @return result: n x d list of n evaluated points in objective space
    function [result] = eval(obj, points, z)
        obj.numFuncEvals = obj.numFuncEvals + size(points, 1);
        if obj.useCPP 
            result = obj.cppCalls.evalPoints(points);
        else         
            [X,~] = size(points);
            result = zeros(X, obj.rd);
            for i = 1:obj.rd
%             for i = 1:1
                F = obj.fun{i};
                if (obj.includeZ || obj.functionID == 26 || obj.functionID == 25) && obj.functionID ~= 29 && obj.functionID ~= 30 && obj.functionID ~= 31
                    result(:,i) = F(points', z')';
                else
                    result(:,i) = F(points')';
                end
            end
        end
    end
      
    %% Evaluate the hessian and gradient of the objective functions at a point P
    % @param P: 1 x D vector representing a point in design space
    % ---
    % @return G: d x D matrix where each row i is the value of the derivative of
    % the ith objective function evaluated at P.
    % @return H: D x D x d multidimensional array, where each third
    % dimension i is the hessian of the ith objective function evaluated
    % at P.
    function [G,H] = evalHessAndDeriv(obj, P)
        if obj.useCPP
            [~, G, H] = obj.cppCalls.evalPointsHessian(P);
            G = G';
        elseif obj.symFreeAnalytic
            pApp = obj.app_var;
            % Scalarization function: 1/2 * |f(x) - target|^2
            G = obj.deriv(P.', pApp.');
            G = G(:, 1:obj.rD);
            H = zeros(obj.rD,obj.rD,obj.rd);

            % try with preexisting matlab functions (instead of symbolic)
            for i=1:obj.rd
                Hfh = obj.hess{i}; %function handle
                Hfull = Hfh(P.', pApp.');
                rows = 1:obj.rD;
                cols = 1:obj.rD;     
                H(:, :, i) = Hfull(rows, cols);
            end
        else
            for i = 1:3
                eps = 0.00001;
                if P(i) - 0.5 < eps
                    P(i) = P(i) + eps;
                end
            end
            P = num2cell(P);

            G = eval(obj.deriv(P{:}))';
            h = obj.hess;
            H = zeros(obj.rD,obj.rD,obj.rd);
 
            for i = 1:obj.rd
                h_i = h{i};
                h_i_eval = double(h_i(P{:}));
                H(:,:,i) = h_i_eval;
            end
        end
    end
end

end
