classdef bicopter < handle
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
        hasoptimizesampling
        problemName
        visName
        useFD
        symFreeAnalytic
        % only for bicopter
        objfunc
        fixedz
        app_var
    end
    
    methods
        function obj = bicopter(basewd)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "bicopter"; %shouldn't have any spaces or special chars except underscore
            obj.hasGroundTruth = false;
            obj.hasMeshes = false;
            obj.hasLinearConstraints = false; % we use this function to sample `good' initial points
            obj.hasoptimizesampling = false; % we use this function to sample `good' initial points
            obj.useFD = false;
            obj.symFreeAnalytic = true;
            obj.fixedz = false;
            
            % working directory; make sure it ends with a "/"
            probwd = "analytic/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            
            % ============ DESIGN VARIABLES ============ 
            total_time_steps = 16;
            obj.rD = total_time_steps*2; % 16 time steps * 2 propellers
            desNames = [];
            for i = 1:total_time_steps
                desNames = [desNames, strcat("u1_t", num2str(i)), strcat("u2_t", num2str(i))]; % for plot labels and filenames
            end
            
            obj.xSymb = sym('x', [obj.rD, 1]);
            % TODO maybe there is a proper way to achieve this
%             for i = 1:obj.rD
%                 eval(['obj.xsymb_', num2str(i), ' = lerp(obj.xSymb(', num2str(i), '), -5, 5);']);
%             end
                
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = lerp(obj.xSymb(1), 2, 14);  %2 + x1*(14-2);
            % TODO check do we really need this scale? seems that the
            % vars have not been used
            
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = 1; 
            appNames = ["length of the bicopter"];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.
            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = 2; 
            perfNames = ["distance to goal", "energy needed"];

            
            
            % ============ GENERIC SETUP ============
            obj.varNames = [desNames, appNames, perfNames];
            probwd = "analytic/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
        end

        function set_app_var(obj, app_var)
            obj.app_var = app_var;
        end

        function set_fixedz(obj)
            obj.fixedz = true;
            obj.zSymb = sym('z', [1 1]);
            obj.rA = 0;
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

        function [f1, g1, h1] = performanceMetric1(obj)
            % PerformanceMetric1 - define the performance metric.
            [f1, g1, h1] = bicopter_setup(obj.xSymb, obj.zSymb, 'position');
            % normalize so all possible outputs are between 0 and 1
        end
        
        function [f2, g2, h2] = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            [f2, g2, h2] = bicopter_setup(obj.xSymb, obj.zSymb, 'energy'); 
            
            % normalize so all possible outputs are between 0 and 1
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
            % only for bicopter
            if obj.hasoptimizesampling
                obj.objfunc = cell(2);
                [obj.objfunc{1}, obj.objfunc{2}, ~]= bicopter_setup(sym('x_tmp', [obj.rD,1]), sym('z_tmp', [obj.rA,1]), 'both');
            end
            
            [fun1, deriv1, hess1] = obj.performanceMetric1(); 
            f1=fun1; 
            
            [fun2, deriv2, hess2] = obj.performanceMetric2(); 
            f2=fun2;
            
            rfun = {fun1, fun2};
            
            rf = [];
            rfuncsCell = {};%{f1, f2};
            rderiv = @(x,z)obj.createJacobian({deriv1, deriv2}, x, z);
            rhess = {hess1, hess2};
        end
        
        function [F, J] = optimized_sampling_points(obj, x, z, F_, J_)
            F = F_(x,z);
            J = J_(x,z);
%             H_ = @(x, z)fast_hess(F_, J_, x, z);
%             H = H_(x,z);
        end
        
        function J = createJacobian(obj, derivs, x, z)
            J = [];
            for i=1:length(derivs)
                d_fi = derivs{i};
                J = [J; d_fi(x,z)];
            end
        end
        
        function [numViolated, optimized_output] = checkLinearConstraints(obj, pDes, pApp)
            numViolated = -1; % a flag...
            options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'SpecifyObjectiveGradient',true);%, 'MaxIterations', 30); %'StepTolerance',5e-2
            disp('time to optimize for sampling')
            lb = [0 * ones(size(pDes',1), 1)];
            ub = [1.0 * ones(size(pDes',1), 1)];
            [optimized_points, fe] = fmincon(@(x)obj.optimized_sampling_points(x, pApp, obj.objfunc{1}, obj.objfunc{2}), pDes',[],[] ,[],[],lb, ub, [], options);
            optimized_output = [optimized_points', pApp];
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


