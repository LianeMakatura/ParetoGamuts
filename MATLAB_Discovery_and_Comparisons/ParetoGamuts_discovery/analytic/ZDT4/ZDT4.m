classdef ZDT4
    %ZDT3 analytic test function
    % https://esa.github.io/pagmo2/docs/cpp/problems/zdt.html
    %   DESIGN VARS
    %
    %   APP VARS
    % 
    %   PERFORMANCE METRICS
    
    
    properties
        % ==== design vars
        rD
        xSymb
        
        % ==== application vars
        rA
        zSymb
        
        % ==== performance metrics
        rd 
        
        varNames
        wd
        hasGroundTruth
        hasMeshes
        hasLinearConstraints
        problemName
        visName
        useFD
        symFreeAnalytic
    end
    
    methods
        function obj = ZDT4(rD, rA, basewd)
            % Constructor 
            obj.problemName = "ZDT4";
            obj.visName = obj.problemName;
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = rD; 
            desNames = []; % for plot labels and filenames
            for j=1:obj.rD
                var_id = "x" + num2str(j);
                desNames = [desNames, var_id];
            end
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % No aliasing necessary; range adjusted in each perf metric
            
            
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = rA; 
            appNames = [];
            for j=1:obj.rA
                var_id = "z" + num2str(j);
                appNames = [appNames, var_id];
            end
            obj.zSymb = sym('z', [1 obj.rA]);
            
            % No aliasing necessary; range adjusted in each perf metric


            
            % ========== PERFORMANCE METRICS ============
            obj.rd = 2; 
            perfNames = ["f1", "f2"];
            
            
            % ============ GENERIC SETUP ============
            obj.varNames = [desNames, appNames, perfNames];
            
            % working directory; make sure it ends with a "/"
            probwd = "analytic/" + obj.problemName + "/";
            obj.wd = basewd + probwd;
            obj.hasGroundTruth = true;
            obj.hasMeshes = false;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
        end
        
        % =========== DEFINE PERFORMANCE METRICS ===========
        function f1 = performanceMetric1(obj)
            % PerformanceMetric1 - define the performance metric.
            % f(x, z) = x1
            % already normalized, x1 \in [0,1]
        
            f1 = obj.xSymb(1);             
        end
        
        function f2 = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            % x1 \in [0,1]
            % x_2.. x_rD & z_1 .. z_rA \in [-5, 5]
            
            x1 = obj.xSymb(1); 
            
            % sum over remaining (rescaled) design vars
            rescale = @(q)(10*q - 5); %from q\in [0,1] to q\in [-5, 5]
            summand = @(q)(q^2 - 10 * cos(4 * pi * q) );
            s = 0;
            for i=2:obj.rD
                s = s + summand(rescale(obj.xSymb(i)));
            end
            for i=1:obj.rA
                s = s + summand(rescale(obj.zSymb(i)));
            end
            
            m = obj.rD + obj.rA; 
            g = 1 + 10*(m-1) + s;
            
            f2 = g * (1 - sqrt(x1 / g));
            
            % normalize so all possible outputs are between 0 and 1
            scalingFactor = (m-1)*30;
            f2 = f2 / scalingFactor; 
        end
                
         % ========== GROUND TRUTH SOLUTION ===================
        function [gt_pts_global, gt_pts_gamut] = getGroundTruth(obj, numPts)
            % recall that x1 in [0,1], x2..xD&z1..zA in [-5, 5]
            if nargin < 2 
                numPts = 500;
            end
            gt_f1 = linspace(0,1,numPts);
            scalingFactor = (obj.rD + obj.rA - 1)*30; %% make sure this matches the one in perfMetric2
            
            % ground truth for lower envelope
            gt_f2 = (1 - sqrt(gt_f1))/scalingFactor; %ground truth ZDT 4
            gt_pts_global = {[gt_f1', gt_f2']};
            
            % get the analytic gamut solution --  only visualizable with 1
            % app var, no sum over app vars needed in g
            gt_z = linspace(0,1,numPts);
            [F1, Z] = meshgrid(gt_f1, gt_z);
            Zscaled = lerp(Z, -5, 5);
            g = 1 + 10 + (Zscaled.^2 - 10.*cos(4.*pi.*Zscaled)); % all other factors of 10 cancel with x_i=0
            F2 = g .* (1 - sqrt(F1 ./ g));
            F2 = F2 ./ (scalingFactor);
  
            gt_pts_gamut = {{F1, F2, Z}};
        end
        
        % =========== WRITE AND ALTER METRIC FILES ===========
        function [rfun, rf, rfuncsCell, rderiv, rhess] = createProblem(obj, createFiles)
            % create the usable matlab function files from symbolic equiv.
            % defined in this problem class
            % Called from the mapping function
            if ~exist(obj.wd + "generated/", 'dir')
                mkdir(obj.wd + "generated/");
            end
            addpath(obj.wd+"generated/");

            
            firstPerfIdx = obj.rD + obj.rA + 1;
            pName = strrep(obj.problemName, " ", "_");

            f1 = obj.performanceMetric1();
            f2 = obj.performanceMetric2();
            f = {f1, f2};
            
            rf = [];
            
            for i=1:obj.rd
                fname = pName + "_" + strrep(obj.varNames(firstPerfIdx + i-1)," ","_");
                fhandle = f{i}; %pick out the correct metric
                
                [fun, ffunc, hess] = obj.setupFunction(fhandle, fname, createFiles);
                
                rfun{i} = fun;
                rf = [rf, ffunc]; % used for scalarization function
                rfuncsCell{i} = ffunc;
                rhess{i} = hess;
            end
            
            % write out the full jacobian matrix as a file function too
            deriv = jacobian(rfuncsCell, [obj.xSymb, obj.zSymb]); 
            deriv_name = pName + "_deriv";
            matlabFunction(deriv, 'File', char(obj.wd +"generated/"+ deriv_name), 'vars',{obj.xSymb.', obj.zSymb.'});
            rderiv = str2func(deriv_name);
        end
        
        function [fun, f, hessh] = setupFunction(obj, fHandle, fname, createFile)
            if createFile
                matlabFunction(fHandle, 'File', char(obj.wd + "generated/"+fname), 'vars',{obj.xSymb.', obj.zSymb.'},... 
                                'Comments', sprintf('Automatically generated by %s.createProblem()', obj.problemName));
%                 obj.fixFileInputs(fname, obj.wd);
            end
            x = [obj.xSymb];
            z = [obj.zSymb];
            
            f_func = str2func(fname);
            fun = f_func;
            f(x, z) = f_func(x.', z.');
            
            % try with matlab function handles
            hess([x,z]) = hessian(f, [x,z]);
            hess_name = fname + "_hess";
            matlabFunction(hess, 'File', char(obj.wd + "generated/"+hess_name), 'vars',{x.', z.'});
            hessh = str2func(hess_name);
        end
        
        function [] = writeMeshes(obj, pts)
            %WRITEMESHES Takes design points (and application points, if necessary) and
            %writes out the STL files for the visualization tool.
            %   Detailed explanation goes here
            if ~obj.hasMeshes
                disp('No meshes associated with this problem.');
                return
            end
        end
    end
end



