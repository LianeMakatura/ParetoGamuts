classdef ZDT1_fhtest
    %ZDT1 analytic test function
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
        function obj = ZDT1_fhtest(rD, basewd)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "ZDT1_fhtest";
            obj.visName = "ZDT1_fhtest";
        
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
            obj.rA = 0; 
            appNames = [];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % must have at least 1 symb var for z for the rest of the framework to function (even if you don't use it in the objectives)
            
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
            % s = sum(x_2..x_rD) + sum(z_1..z_rA) 
            % g(x) = 1 +  9s /(rD + rA - 1)
            % f(x, z) = g(x, z) * [1 - \sqrt(x1 / g(x, z))]

            x1 = obj.xSymb(1); 
            s = sum(obj.xSymb(2:end));
            g = 1 + 9*s / (obj.rD + obj.rA - 1);
            
            f2 = g * (1 - sqrt(x1 / g));
            
            % normalize so all possible outputs are between 0 and 1
            f2 = f2 / 10; 
        end
        
        % ========== GROUND TRUTH SOLUTION ===================
        function [gt_pts_global, gt_pts_gamut] = getGroundTruth(obj, numPts)
            if nargin < 2 
                numPts = 100;
            end
            gt_f1 = linspace(0,1,numPts);   %f1 values of optimal points
            
            % ground truth for lower envelope
            gt_f2 = (1 - sqrt(gt_f1))/10; %ground truth ZDT 1
            gt_pts_global = {[gt_f1', gt_f2']};
            
            % get the analytic gamut solution --  only visualizable with 1
            % app var, no sum over app vars needed in g
            gt_z = linspace(0,1,numPts);
            [F1, Z] = meshgrid(gt_f1, gt_z);
            g = (1 + 9 / (obj.rD + obj.rA - 1) .* Z);
            F2 = g .* (1 - sqrt(F1 ./ g)) / 10;

            gt_pts_gamut = {{F1, F2, Z}};
        end
        
        
        
        % =========== WRITE AND ALTER METRIC FILES ===========
        function [rfun, rf, rfuncsCell, rderiv, rhess] = createProblem(obj, createFiles)
            % create the usable matlab function files from symbolic equiv.
            % defined in this problem class
            % Called from the mapping function
            
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
            matlabFunction(deriv, 'File', char(obj.wd +deriv_name), 'vars',{obj.xSymb.', obj.zSymb.'});
            rderiv = str2func(deriv_name);
        end
        
        function [fun, f, hessh] = setupFunction(obj, fHandle, fname, createFile)
            if createFile
                matlabFunction(fHandle, 'File', char(obj.wd + fname));
                obj.fixFileInputs(fname, obj.wd);
            end
            x = [obj.xSymb];
            z = [obj.zSymb];
            
            f_func = str2func(fname);
            fun = f_func;
            f(x, z) = f_func(x.', z.');
            
            % try with matlab function handles
%             derivh = matlabFunction(deriv, 'vars',{x.', z.'}); % build up full jacobian before converting to matlab function
            hess([x,z]) = hessian(f, [x,z]);
            hess_name = fname + "_hess";
            matlabFunction(hess, 'File', char(obj.wd + hess_name), 'vars',{x.', z.'});
            hessh = str2func(hess_name);
        end
        
        
        function [] = fixFileInputs(obj, filename, wd)
            % FIXFILEINPUTS - change matlab function to only take x, z
            % vectors, then parse into xi, zi within the function itself
            % Internal function, only called from createProblem
            % shouldn't need to alter
    
            fullname = strcat(wd, filename, '.m');
            data = fileread(fullname);
            lines = strsplit(data, '\n');
            
            % write a new file with the altered signature / headers
            fileh = fopen(fullname, 'w');
            
            % fix the function signature to take x (and z if unconstrained)
            sig = strsplit(lines{1}, '(');
            newsig = strcat(sig{1}, '(x,z)');
            fprintf(fileh, '%s\n', newsig);
                        
            % output modified comment block
            for i=2:5
                fprintf(fileh, '%s\n', lines{i});
            end
            fprintf(fileh, '%% Includes scripted file modifications (%s.fixFileInputs) to match system specs\n\n', obj.problemName);
            
            % parse the variable blocks into individual sym. vars
            fprintf(fileh, '%% unpack the block vars into the variable names used by the function\n');
            for j=1:obj.rD
                var_id = strcat('x', num2str(j));
                fprintf(fileh, '%s = x(%d, :).'';\n', var_id, j);
            end
            for j=1:obj.rA
                var_id = strcat('z', num2str(j));
                fprintf(fileh, '%s = z(%d, :).'';\n', var_id, j);
            end
            
            % write out the rest of the auto-generated function
            fprintf(fileh, '\n%% Original content from Symbolic Math Toolbox\n');
            for i=6:length(lines)
                fprintf(fileh, '%s\n', lines{i});
            end
                
            fclose(fileh);
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

