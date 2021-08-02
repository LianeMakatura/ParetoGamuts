classdef ProblemTemplate
    %ZDT1 analytic test function
    %
    %   DESIGN VARS
    %
    %   APP VARS
    % 
    %   PERFORMANCE METRICS
    
    
    properties
        % ==== design vars
        rD
        xSymb
        %+alias vars for xSymb (eg, obj.length, width)
        
        % ==== application vars
        rA
        zSymb
        %+alias vars for zSymb (eg, obj.angle)

        
        % ==== performance metrics
        rd 
        
        % generic setup
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
        function obj = ProblemTemplate(basewd)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "ProblemTemplate"; %shouldn't have any spaces or special chars except underscore
            obj.hasGroundTruth = false;
            obj.hasMeshes = false;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
            
            % working directory; make sure it ends with a "/"
            probwd = "analytic/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = 2; 
            desNames = ["", ""]; % for plot labels and filenames
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = lerp(obj.xSymb(1), 2, 14);  %2 + x1*(14-2);
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = 1; 
            appNames = [""];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.

            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = 2; 
            perfNames = ["", ""];

            
            
            % ============ GENERIC SETUP ============
            obj.varNames = [desNames, appNames, perfNames];
            
        end
        
        % ==================================================
        %            DEFINE PERFORMANCE METRICS 
        %
        % the pareto front is assumed to lie between 0 and 1 in all
        % performance metrics. You will have to normalize them accordingly.
        %
        % Add extra perf metric functions if necessary. If you add any, 
        % you'll have to change the createProblem function to account 
        % for the extra metric(s)
        %
        % Follow ProblemTemplate.m if fully symbolic (relying on obj.x/zSymb)
        % If one or more metrics are not symbolic, follow example in
        % LBracket_stress.m (you will have to change more in eg createProblem)
        % If using matrix operations, be careful -- may not be able to
        % evaluate multiple points at once.
        % ==================================================

        function f1 = performanceMetric1(obj)
            % PerformanceMetric1 - define the performance metric.
            
            f1 = obj.xSymb(1); % use the aliases if the vars should be scaled
            
            % normalize so all possible outputs are between 0 and 1
        end
        
        function f2 = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            
            f2 = obj.zSymb(1); 
            
            % normalize so all possible outputs are between 0 and 1
        end
        
        
        
        
        
        % ==================================================
        %            DEFINING LINEAR CONSTRAINTS
        %
        % If any linear constraints on design/application variables.
        % Box constraints on each variable are handled autmoatically,
        % do not specify them here. Acoustics.m has an example.
        % ==================================================

        function [G, DxG, DzG] = constraints(obj)
            % constraints satisfied if <= 0 (active if == 0)
            % G defined symbolically, careful to use appropriately scaled
            % (aliased) variables!
            G = {
                -0.8*obj.xSymb(1) -3;
                };
            
            DxG = jacobian(G, obj.xSymb);
            DzG = jacobian(G, obj.zSymb);
        end
        
        function [A, b] = matrixLinearConstraints(obj)
            A = [ -0.8, 0]; % coefficients of x1..x5 (scaled appropriately) (same as DxG)
            b = [3];            % constant term(s) in the G equations (with scaled variables) -- create object p=ProblemTemplate(), then print p.constraints to obtain easily
        end
        
        % shouldn't have to modify
        % Report the active linear constraints on a single point
        function [numActive, DxG_linear, DzG_linear] = getActiveLinearConstraints(obj, xVals, zVals)
            if size(xVals, 1) > 1
                disp('Only written for single point right now');
                return
            end
            
            % use fake app value if no application variable (to allow standard MOOP)
            if obj.rA == 0
                zVals = zero(size(xVals, 1), 1);
            end
            
            epsilon = .001;
            [~, DxG, DzG] = obj.constraints();

            % check which ones are active
            [A,b] = obj.matrixLinearConstraints();
            
            % need (A*xVals' - b) < 0; active if >-eps
            activeMask = bsxfun(@gt, A*xVals' - b, -epsilon); %bit mask, 1 if constraint active
            activeIDs = find(activeMask);
            
            %return only active rows
            numActive = length(activeIDs);
            DxG_linear = DxG(activeIDs, :);
            DzG_linear = DzG(activeIDs, :);
        end

        % shouldn't have to modify
        function [numViolated, invalidIndices] = checkLinearConstraints(obj, xVals, zVals)
            % use fake app value if no application variable (to allow standard MOOP)
            if obj.rA == 0
                zVals = zero(size(xVals, 1), 1);
            end
            
            [A,b] = obj.matrixLinearConstraints();
            
            % need each column of A*xVals' to be <= b
            violated = bsxfun(@gt, A*xVals', b); %bit mask, 1 if constraint violated
            
            %columns with all 0 are valid points
            invalidIndices = find(any(violated==1))'; % all zeros 
            numViolated = length(invalidIndices);
        end
        
        
        
        
        
        
        % ========== GROUND TRUTH SOLUTION ===================
        % for extablished benchmark functions, this "ground truth" is only
        % for the fixed context situation; here, it should be the lower
        % envelope (if the application var is just one of the original
        % design vars) or some particular front (if the app var was originally
        % a constant, and we're varying it)
        % ignore if no ground truth.
        function gt_pts = getGroundTruth(obj, numPts)
            if ~obj.hasGroundTruth
                disp("Current problem has no ground truth.")
                gt_pts = [];
                return
            end
            
            % otherwise, return the ground truth points [example]
            if nargin < 2 
                numPts = 100;
            end
            gt_x = linspace(0,1,numPts);
            gt_y = gt_x^2;
            gt_pts = {[gt_x', gt_y']}; % return in cell array for distinct patches
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
        
        % shouldn't need to edit
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
        

        
        % shouldn't need to edit (unless you want to write meshes, then go for it)
        function [] = writeMeshes(obj, pts)
            %WRITEMESHES Takes design points (and application points, if necessary) and
            %writes out the STL files for the visualization tool.
            %   Detailed explanation goes here
            if ~obj.hasMeshes
                disp('No meshes associated with this problem.');
                return
            end
            
            numMeshes = size(pts, 1);
            for i=1:numMeshes
                p = pts(i, :);

                % generate the target geometry file name
                MESHout = "mesh.stl";
                
                % write out file
                outfile = obj.wd + "meshes/" + MESHout;
                
            end
        end
        
    end
end


