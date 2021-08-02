classdef BikeStay
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
        drop_h
        tube_h
        stay_depth
        
        % ==== application vars
        rA
        zSymb
        %+alias vars for zSymb (eg, obj.angle)
        pivotOffsetX
        pivotOffsetY
        
        % ==== performance metrics
        rd 
        
        % ==== other variables
        pivotbase
        stay_b
        drop_w
        ur_drop
        finalPivot
        tube_length
        
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
        function obj = BikeStay(basewd)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "BikeStay"; %shouldn't have any spaces or special chars except underscore
            obj.visName = obj.problemName;
            obj.hasGroundTruth = false;
            obj.hasMeshes = true;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
            
            % working directory; make sure it ends with a "/"
            probwd = "meshBased/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = 3; 
            desNames = ["Dropout Height", "Tube Height", "Stay Depth"]; % for plot labels and filenames
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            obj.drop_h = lerp(obj.xSymb(1), 50, 80);
            obj.tube_h = lerp(obj.xSymb(2), 10, 20);
            obj.stay_depth = lerp(obj.xSymb(3), 6, 16);
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = 1; %% this problem does have 2 context vars, but honestly I've never tested that. It should work? but you may find bugs. I'd just use the one for now, unless we're totally done with time and energy to kill
            appNames = ["Pivot Offset X"];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.

            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary
%             obj.pivotOffsetX = lerp(obj.zSymb(1), 0, -40); % offsetting to the left, negative x direction
%             if obj.rA < 2
%                 obj.pivotOffsetY = 0;
%             else
%                 obj.pivotOffsetY = lerp(obj.zSymb(1), -9, 9);
%             end
            obj.pivotOffsetY = lerp(obj.zSymb(1), -9, 40); % offsetting vertically
            if obj.rA < 2
                obj.pivotOffsetX = 0;
            else
                obj.pivotOffsetX = lerp(obj.zSymb(1), 0, -40);
            end

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = 2; 
            perfNames = ["Mass", "Compliance"];

            % ========== OTHER CONSTANTS ============
            % Dropout height: x1
            % Dropout width: 30 fixed
            % Dropout thickness : x3
            % Tube Height: x2
            % Tube Width: x3
            % Tube length: ||finalPivot - ur_drop||_2
            obj.pivotbase = [0,0];
            obj.stay_b = [obj.pivotbase(1) - 180, obj.pivotbase(2) - 175];
            obj.drop_w = obj.drop_h;
            obj.ur_drop = [obj.stay_b(1) + obj.drop_w, obj.stay_b(2) + obj.drop_h];
            obj.finalPivot = [obj.pivotbase(1) + obj.pivotOffsetX, obj.pivotbase(2) + obj.pivotOffsetY];
%             obj.tube_length = sqrt(sum((obj.finalPivot - obj.ur_drop).^2));
            diff = (obj.finalPivot - obj.ur_drop);
            obj.tube_length = sqrt( diff(1)^2 + (2*diff(2))^2 );
            
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
            % MASS PROXY - define the performance metric.
            drop_mass = obj.drop_h * obj.drop_w * obj.stay_depth;
            tube_mass = obj.tube_h * obj.tube_length * obj.stay_depth;
            
            f1 = drop_mass + tube_mass;
            f1 = (f1 - 2.35e+04) / (1.16e+05 - 2.35e+04) / 1.7;
            % normalize so all possible outputs are between 0 and 1
        end
        
        function f2 = performanceMetric2(obj)
            % Stiffness PROXY 
            tube_stiff = obj.tube_length / (obj.tube_h * obj.stay_depth);
            drop_stiff = 1 / (obj.drop_h * obj.drop_w);
            f2 = tube_stiff + drop_stiff;
            f2 = (f2 - 41/64) / (55/12 - 41/64) / 1.7;
            % normalize so all possible outputs are between 0 and 1
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


