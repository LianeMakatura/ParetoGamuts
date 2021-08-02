classdef AppMappingFunction < handle
%% Instance variables
properties     
    % Visualization
    imageTitle

    % Static constants
    minBound
    maxBound
    rD
    rA
    rd
    functionID
    varNames
    p               % problem class
    hasMeshes
    hasGroundTruth
    hasLinearConstraints
    useFD
    symFreeAnalytic
    onTheFlyJH
    
    % Function evaluation variables
    f
    funcsCell
    fun
    deriv
    hess
    fhSymb
    
    %statistics
    numFuncEvals
    numJacEvals
    numHessEvals
    numScalarizationFuncEvals
end
%% Functions
methods        
    %% Constructor.
    function obj = AppMappingFunction(id, basewd) % z only necessary sometimes
        if nargin==0 % needed for Object arrays
            return
        end
        
        obj.numFuncEvals = 0;
        obj.numJacEvals = 0;
        obj.numHessEvals = 0;
        obj.numScalarizationFuncEvals = 0;
        
        obj.functionID = id;
        obj.imageTitle = 'results\foo';
        onTheFlyJH = false;
        
        % matlab file writing
        createFiles = true;
        
        if id >=1 && id <= 5 % ZDT 1-4 & 6 test functions
            % these parameters are allowed to change for any ZDT problem
            rD = 2;
            rA = 1;
            
            % all have 2 performance metrics
            if id == 1
                p = ZDT1(rD, rA, basewd);
%                 p = ZDT1_noApp(rD, basewd);
            elseif id == 2
                p = ZDT2(rD, rA, basewd);
            elseif id == 3
                p = ZDT3(rD, rA, basewd);
            elseif id == 4
                p = ZDT4(rD, rA, basewd);
            elseif id == 5
                p = ZDT6(rD, rA, basewd);
            end
            
        elseif id == 7
            p = activeDTLZ1(basewd);
            warning('off','all');
        elseif id == 8
            p = activeKNO1(basewd);
            warning('off','all');
        elseif id == 9
            p = ZDT1_fhtest(2, basewd);
        elseif id == 10
            rD = 2;
            rA = 1;
            rd = 2;
            polyOrder = 6;
            p = Fourier(rD, rA, rd, polyOrder, basewd); %% currently only written for rA=1
        elseif id == 11     %% bicopter
            p = bicopter(basewd);
            
        elseif id == 12     %% bicopter_2_contexts
            p = bicopter_2_contexts(basewd);

        elseif id == 26     %% 2D L bracket - unconstrained z
            p = LBracket(basewd, false, false); %3d, unconstrained
            onTheFlyJH = true;
            
        elseif id == 27     %% 3D L bracket - unconstrained z
            p = LBracket(basewd, true, false); %3d, unconstrained
            onTheFlyJH = true;
                        
        elseif id == 28     %% 3D L bracket - stress metric
            p = LBracket_stress(basewd, true, false); %3d, unconstrained
            onTheFlyJH = true;
                        
        elseif id == 29     %% Lamp
            metricsToUse = [1,2,3]; %1=Stability, 2= Mass, 3 = Focal Point (choose 2 or all 3 of them)
            p = Lamp(basewd, metricsToUse);
            onTheFlyJH = true;
            
        elseif id == 30     %% Acoustics
            tic;
            p = Acoustics(basewd);
            elapsed = toc;
            fprintf('time to create Acoustics object: %s\n', durationString(elapsed));
            
        elseif id == 31
            metricsToUse = [1,2]; %1=Mass, 2=Ideal Power Output
            p = Turbine(basewd, metricsToUse);
            onTheFlyJH = true;
        elseif id == 32
            p = BikeRocker(basewd);
        elseif id == 33
            p = BikeStay(basewd);
            
        elseif id == 36
            metricsToUse = [1,2];
            numX = 9;
            numY = 9;
            p = Gridshell(basewd, numX, numY, metricsToUse);
        elseif id == 37
            metricsToUse = [1,2];
            numX = 10;
            numY = 10;
            p = Gridshell_nonsymb(basewd, numX, numY, metricsToUse);
        elseif id == 38
            metricsToUse = [1,2];
            numX = 13;
            numY = 13;
            p = Gridshell_nonsymb(basewd, numX, numY, metricsToUse);
        elseif id == 39
            metricsToUse = [1,2];
            numX = 18;
            numY = 18;
            p = Gridshell_nonsymb(basewd, numX, numY, metricsToUse);
        elseif id == 40
            metricsToUse = [1,2];
            numX = 23;
            numY = 23;
            p = Gridshell_nonsymb(basewd, numX, numY, metricsToUse);
        elseif id == 41
            metricsToUse = [1,2];
            numX = 25;
            numY = 25;
            p = Gridshell_nonsymb(basewd, numX, numY, metricsToUse);
        else
            disp('Invalid function ID');
            pause;
        end
        
        % ====== universal processing ======= 
        obj.rD = p.rD;
        obj.rA = p.rA;
        obj.rd = p.rd;
                
        rD = obj.rD; 
        rA = obj.rA;
        rd = obj.rd;
                       
        obj.p = p;
        obj.varNames = p.varNames;
        obj.hasGroundTruth = p.hasGroundTruth;
        obj.hasMeshes = p.hasMeshes;
        obj.hasLinearConstraints = p.hasLinearConstraints;
        obj.useFD = p.useFD;
        obj.symFreeAnalytic = p.symFreeAnalytic;
        
        obj.onTheFlyJH = onTheFlyJH;
        
        tic;
        [rfun, rf, rfuncsCell, rderiv, rhess] = p.createProblem(createFiles);%write new files
        elapsed = toc;
        fprintf('time to create problem: %s\n', durationString(elapsed));

        obj.fun = rfun;
        obj.f = rf;
        obj.funcsCell = rfuncsCell;
        obj.deriv = rderiv;
        obj.hess = rhess;

        obj.minBound = zeros(1,rD+rA);
        obj.maxBound = ones(1,rD+rA);

        % Generate symbolic design variables for derivative/hessian evaluation  
        xs = obj.p.xSymb;
        zs = obj.p.zSymb;

        % Compute each hessian symbolically %may use if we store it
        % ahead of time, but would need multiple variants (dxdx, dxdz)
%             q = [xs.'; zs.'];
%             hess = cell([rd, 1]);
%             for i=1:rd
%                 hess{i} = hessian(obj.funcsCell{i}, q);
%             end


        % Scalarization function: 1/2 * |f(x) - target(x)|^2
        if ~obj.useFD && ~obj.symFreeAnalytic
            % generate scalarization function
            targets = sym('t', [1 obj.rd]);   % symbolic target variables
            t = targets.';
            tic;
            scalarFSymb(xs, zs, targets) = sum((obj.f-targets).^2) / 2.0;
            elapsed = toc;
            fprintf('time to create target scalarization: %s\n', durationString(elapsed));

            % compute the analytic gradient to aid fmincon
            tic;
            gradFSymb = jacobian(scalarFSymb, xs).'; % column gradf, only wrt x bc z, targ fixed
            elapsed = toc;
            fprintf('time to create scalarization jacobian: %s\n', durationString(elapsed));
            tic;
            obj.fhSymb = matlabFunction(scalarFSymb, gradFSymb,'vars',{xs.', zs.', t}, 'outputs',{'f','gradf'}); % Convert to Matlab function handle for easier use later
            elapsed = toc;
            fprintf('time to create scalarization matlab function: %s\n', durationString(elapsed));
        end
    end
    

    %% Retrieve the image title for saved figures/data
    function imageTitle = getImageTitle(obj)
        imageTitle = obj.imageTitle;
    end

    %% Return the buffer size across one dimension
    function numCells = getNumCells(obj)
        if obj.rd + obj.rA <= 3
            numCells = 200;
        else
            numCells = 100;
        end
    end
      
    %% Return the dimensionality of the problem
    % return rD: Dimensionality of design space
    % return rd: Dimensionality of objective space
    function [rD, rd, rA] = getDimensions(obj)
        rD = obj.rD;
        rA = obj.rA;
        rd = obj.rd;
    end
      
    
    %% Return the problem's ground truth points
    function [pts_global, pts_gamut] = getGroundTruth(obj)
        [pts_global, pts_gamut] = obj.p.getGroundTruth();
    end
    
    %% write the meshes for this problem
    function [] = writeMeshes(obj, pts)
        if obj.hasMeshes
            disp("Writing meshes corresponding to Pareto optimal designs.")
            obj.p.writeMeshes(pts);
        else
            disp("This problem does not require any mesh generation.");
        end
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
    function fh = getFunctionHandle(obj, target, fixedAppValues)
        if obj.useFD
            fh = @(x)obj.FD_scalarF(x, fixedAppValues, target); %don't transpose, already correct
        elseif obj.symFreeAnalytic
            fh = @(x)obj.nonSymb_scalarF(x, fixedAppValues, target);
        else
            fh = @(x)obj.fhSymb(x', fixedAppValues', target');
        end
    end
    
    %% Evaluate the scalarization function for a problem with non-symbolic metric(s)
    % @param pDes: an 1 x D point in design space
    % @param pApp: an 1 x A point in application space
    % @param target: a 1 x d target point in performance space 
    % ---
    % @return result: scalarized evaluated point in objective space
    function f = FD_scalarF(obj, pDes, pApp, target)
        % use fake app value if no application variable (to allow standard
        % MOOP)
        if obj.rA == 0
            pApp = 0;
        end
        
        % Scalarization function: 1/2 * |f(x) - target(x)|^2
        f = 0;
        for i = 1:obj.rd
            F_i = obj.fun{i};
            scalarizedFi = (F_i(pDes', pApp') - target(i))^2;
            f = f + scalarizedFi;
        end
        f = f/2;
    end
      
    %% Evaluate the scalarization function for a problem with non-symbolic but analytic metrics
    function [f, gradf] = nonSymb_scalarF(obj, pDes, pApp, target)
        % use fake app value if no application variable (to allow standard
        % MOOP)
        if obj.rA == 0
            if obj.p.app_var
                pApp = obj.p.app_var;
            else
                pApp = 0;
            end
        end
        
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
        V = ones(size(points, 1), 1);
        
        % take care of box constraints
        [XMax, ~] = find(points > 1);
        [XMin, ~] = find(points < 0 );
        V(XMax) = 0;
        V(XMin) = 0;
        
        % take care of linear constraints
        if obj.hasLinearConstraints
            validIDs = find(V==1); % row indices of valid points
           
            [~, invalidIndices] = obj.p.checkLinearConstraints(points(validIDs, 1:obj.rD), points(validIDs, obj.rD+1:obj.rD+obj.rA));
            V(invalidIndices) = 0; %update full validity vector
        end
    end
     
    %% Generate nSamples random samples in rD+rA space
    function S = sample(obj, nSamples)
        if obj.hasLinearConstraints
            % collect samples that don't violate constraints
            P = zeros(nSamples, obj.rD+obj.rA);
            i=1;
            while i < nSamples+1
                augSamp = obj.minBound + 0.001* rand(1, obj.rD+obj.rA).*(obj.maxBound-obj.minBound) + 0.5;
                [numViolated, arg] = obj.p.checkLinearConstraints(augSamp(1, 1:obj.rD), augSamp(obj.rD+1:obj.rD+obj.rA));
                if numViolated == 0 %valid sample
                    P(i, :) = augSamp;
                    i = i+1;
                end
                if numViolated == -1 % only for bicopter example
                    P(i, :) = arg;
                    perf_val = obj.eval(P(i, 1:obj.rD), P(i, obj.rD+1 : obj.rD+obj.rA));
                    if perf_val(1) < 2
                        i = i+1;
                    end
                end
            end
        else
            P = obj.minBound + rand(nSamples, obj.rD+obj.rA).*(obj.maxBound-obj.minBound);
        end
        pDes = P(:, 1:obj.rD);

        % use fake app value if no application variable (to allow standard MOOP)
        if obj.rA == 0
            if obj.p.app_var
                pApp = obj.p.app_var;
            else
                pApp = zeros(size(P, 1), 1);
            end
        else
            pApp = P(:, obj.rD+1 : obj.rD+obj.rA);
        end
        
        S = [P obj.eval(pDes, pApp)];
    end
 
    %% Evaluate a list of points in design space
    % @param pDes: an n x D list of n points in design space
    % @param pApp: an n x A list of n points in application space
    % ---
    % @return result: n x d list of n evaluated points in objective space
    function [result] = eval(obj, pDes, pApp)
        obj.numFuncEvals = obj.numFuncEvals + size(pDes, 1);
%         tic;
        % use fake app value if no application variable (to allow standard MOOP)
        if obj.rA == 0
            if obj.p.app_var
                pApp = obj.p.app_var;
            else
                pApp = zeros(size(pDes, 1), 1);
            end
        end
        
        if obj.hasLinearConstraints  % this was actually only for the acoustics problem (30) originally, bc can't eval more than 1 point at a time due to matrix math. May not always need.
            numPts = size(pDes, 1);
            result = zeros(numPts, obj.rd);
            for i=1:obj.rd
                F = obj.fun{i};
                parfor j=1:numPts
                    result(j, i) = F((pDes(j, :))', (pApp(j, :))')';
                end
            end
        else         
            numPoints = size(pDes, 1);
            result = zeros(numPoints, obj.rd);
            for i = 1:obj.rd
                F = obj.fun{i};
                result(:,i) = F(pDes', pApp')';
            end
        end
%         elapsed = toc;
%         fprintf('time for function evaluation: %0.6f', elapsed);
    end
    
    %% Evaluate a list of points in design space
    % @param pDes: an n x D list of n points in design space
    % @param pApp: an n x A list of n points in application space
    % @param i:     index (1..rd) of performance metric to eval
    % ---
    % @return result: n x 1 list of n evaluated points in objective space
    function [result] = eval_ith_metric(obj, pDes, pApp, i)
        % currently only used in FD hessian, don't add this to the
        % evaluation count (hessian evals counted separately)
        % use fake app value if no application variable (to allow standard MOOP)
        if obj.rA == 0
            if obj.p.app_var
                pApp = obj.p.app_var;
            else
                pApp = zeros(size(pDes, 1), 1);
            end
        end
        
        if obj.hasLinearConstraints        
            % evaluate each point separately
            F = obj.fun{i};
            numPts = size(pDes, 1);
            result = zeros(numPts, 1);
            for i=1:numPts
                result(i, 1) = F((pDes(i, :))', (pApp(i, :))')';
            end
        else
            % evaluate all together
            F = obj.fun{i};
            result = F(pDes', pApp')';
        end
    end
    
    
    function [J] = evalJacobian(obj, wrt_vars, cpDes, cpApp, eps, fd_type)
        obj.numJacEvals = obj.numJacEvals + size(cpDes, 1);
         % use fake app value if no application variable (to allow standard MOOP)
        if obj.rA == 0
            if obj.p.app_var
                cpApp = obj.p.app_var;
            else
                cpApp = zeros(size(cpDes, 1), 1);
            end
        end
        
        if obj.useFD
            if ~exist('eps', 'var')
                disp('Must provide epsilon for Finite Differencing.')
            end            
            J = fd_jacobian(obj, wrt_vars, cpDes, cpApp, eps, fd_type);
        elseif obj.onTheFlyJH
            xSymb = obj.p.xSymb;
            zSymb = obj.p.zSymb;
            if strcmp(wrt_vars, 'des_only')
                DF(xSymb, zSymb) = jacobian(obj.funcsCell, xSymb);
            elseif strcmp(wrt_vars, 'app_only')
                DF(xSymb, zSymb) = jacobian(obj.funcsCell, zSymb);
            elseif strcmp(wrt_vars, 'both')
                DF(xSymb, zSymb) = jacobian(obj.funcsCell, [xSymb, zSymb]);
            end
            J = double(evalAtPt(DF, cpDes, cpApp));
        else
            % use preexisting matlab functions (instead of symbolic)
            Jfull = obj.deriv(cpDes.', cpApp.');
            if strcmp(wrt_vars, 'des_only')
                J = Jfull(:, 1:obj.rD); % only cols for x derivatives
            elseif strcmp(wrt_vars, 'app_only')
                J = Jfull(:, obj.rD+1:end); % only cols for z derivatives
            elseif strcmp(wrt_vars, 'both')
                J = Jfull; % use all
            end
        end
    end
    
    function [H] = evalSecondPartials(obj, wrt_vars1, wrt_vars2, cpDes, cpApp, eps, fd_type)
        obj.numHessEvals = obj.numHessEvals + size(cpDes, 1);
         % use fake app value if no application variable (to allow standard MOOP)
        if obj.rA == 0
            if obj.p.app_var
                cpApp = obj.p.app_var;
            else
                cpApp = zeros(size(cpDes, 1), 1);
            end
        end
        
        % construct the quasi-hessian d^2/dudv(f_i) = J_v(grad_u(f_i))^T
        if obj.useFD
            if ~exist('eps', 'var')
                disp('Must provide epsilon for Finite Differencing.')
            end    
            H = fd_hessian(obj, wrt_vars1, wrt_vars2, cpDes, cpApp, eps, fd_type);
        elseif obj.onTheFlyJH
            xSymb = obj.p.xSymb;
            zSymb = obj.p.zSymb;
            
            H = cell(obj.rd, 1);
            
            for i=1:obj.rd
                fi = obj.funcsCell{i};
                % take the first partials
                if strcmp(wrt_vars1, 'des_only')
                    grad_fi = gradient(fi, xSymb);
                elseif strcmp(wrt_vars1, 'app_only')
                    grad_fi = gradient(fi, zSymb);
                elseif strcmp(wrt_vars1, 'both')
                    grad_fi = gradient(fi, [xSymb, zSymb]);
                end
                
                % take the second partials
                if strcmp(wrt_vars2, 'des_only')
                    H_fi(xSymb, zSymb) = jacobian(grad_fi, xSymb);
                elseif strcmp(wrt_vars2, 'app_only')
                    H_fi(xSymb, zSymb) = jacobian(grad_fi, zSymb);
                elseif strcmp(wrt_vars2, 'both')
                    H_fi(xSymb, zSymb) = jacobian(grad_fi, [xSymb, zSymb]);
                end              
                H{i} = double(evalAtPt(H_fi, cpDes, cpApp));
            end
        else
            H = cell(obj.rd, 1);

            % try with preexisting matlab functions (instead of symbolic)
            for i=1:obj.rd
                Hfh = obj.hess{i}; %function handle
                Hfull = Hfh(cpDes.', cpApp.'); % for large hessians
                % take the rows for appropriate first partials
                if strcmp(wrt_vars1, 'des_only')
                    rows = 1:obj.rD;
                elseif strcmp(wrt_vars1, 'app_only')
                    rows = obj.rD+1:obj.rD+obj.rA;
                elseif strcmp(wrt_vars1, 'both')
                    rows = 1:obj.rD+obj.rA;
                end
                
                % take cols for appropriate second partials
                if strcmp(wrt_vars2, 'des_only')
                    cols = 1:obj.rD;
                elseif strcmp(wrt_vars2, 'app_only')
                    cols = obj.rD+1:obj.rD+obj.rA;
                elseif strcmp(wrt_vars2, 'both')
                    cols = 1:obj.rD+obj.rA;
                end              
                H{i} = Hfull(rows, cols);
            end
        end
    end
      
end

end
