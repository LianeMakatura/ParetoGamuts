classdef ZDT3
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
        function obj = ZDT3(rD, rA, basewd)
            % Constructor 
            obj.problemName = "ZDT3";
            obj.visName = obj.problemName;
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = rD; 
            desNames = []; % for plot labels and filenames
            for j=1:obj.rD
                var_id = "x" + num2str(j);
                desNames = [desNames, var_id];
            end
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % No aliasing necessary; range is already [0,1]
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = rA; 
            appNames = [];
            for j=1:obj.rA
                var_id = "z" + num2str(j);
                appNames = [appNames, var_id];
            end
            obj.zSymb = sym('z', [1 obj.rA]);
            
            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary

            
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

            x1 = obj.xSymb(1); 
            s = sum(obj.xSymb(2:end)) + sum(obj.zSymb);
            g = 1 + 9*s ./ (obj.rD + obj.rA - 1);
            
            x1_g = x1 ./ g;
            
            f2 = g * (1 - sqrt(x1_g)  - (x1_g .* sin(10 * pi * x1)));
            
            % normalize so all possible outputs are between 0 and 1
            f2 = f2 / 11 + 0.5; 
        end
        
        % ========== GROUND TRUTH SOLUTION ===================
        function [gt_pts_global, gt_pts_gamut] = getGroundTruth(obj, numPts)
            if nargin < 2 
                numPts = 10;
            end
            % Ground truth for the ZDT 3 problem is disconnected on these
            % intervals
            gt_x1 = linspace(0,.083001, numPts);
            gt_x2 = linspace(.1822, .2577, numPts);
            gt_x3 = linspace(.4093, .4538, numPts);
            gt_x4 = linspace(.6183, .6525, numPts);
            gt_x5 = linspace(.8233, .85183,numPts);

            % global front (lower envelope)
            gt_f2 = @(gt_f1)( (1 - sqrt(gt_f1) - gt_f1.*(sin(10*pi*gt_f1))) ./ 11 + 0.5); %ground truth ZDT 3, with normalization on f2
            gt_pts_global = {[gt_x1', gt_f2(gt_x1)'], ...
                      [gt_x2', gt_f2(gt_x2)'], ...
                      [gt_x3', gt_f2(gt_x3)'], ...
                      [gt_x4', gt_f2(gt_x4)'], ...
                      [gt_x5', gt_f2(gt_x5)'], ...
                      };
                  
            % get the analytic gamut solution
            gt_pts_gamut = ZDT3_gamutGroundTruth();
%             f2_gamut = @(f1_pts)obj.gt_f2_gamut(f1_pts, numPts);
% 
%             gt_pts_gamut = {f2_gamut(gt_x1), ...
%                 f2_gamut(gt_x2),...
%                 f2_gamut(gt_x3),...
%                 f2_gamut(gt_x4),...
%                 f2_gamut(gt_x5),...
%                 };
        end
        
        % helper function
        function gt_surf_pts = gt_f2_gamut(obj, f1pts, numPts)
            gt_z = linspace(0,1,numPts);
            [F1, Z] = meshgrid(f1pts, gt_z);
            
            % --  only visualizable with 1
            % app var, no sum over app vars needed in g
            g = (1 + 9 / (obj.rD + obj.rA - 1) .* Z);
            F2 = g .* (1 - sqrt(F1./g) - (F1 ./ g).*sin(10.*pi.*F1) );
            F2 = F2 ./ 11 + 0.5;
            gt_surf_pts = {F1, F2, Z};
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



