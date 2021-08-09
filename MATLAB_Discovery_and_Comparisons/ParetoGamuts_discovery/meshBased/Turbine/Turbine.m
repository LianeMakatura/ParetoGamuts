classdef Turbine
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
        radius
        height
        pitch
        
        % ==== application vars
        rA
        zSymb
        %+alias vars for zSymb (eg, obj.angle)
        vwind
        
        % ==== performance metrics
        rd 
        metricIDs
        
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
        function obj = Turbine(basewd, perfMetricIDs)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "Turbine"; %shouldn't have any spaces or special chars except underscore
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
            desNames = ["Turbine Radius", "Turbine Height", "Pitch Angle"]; % for plot labels and filenames
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = lerp(obj.xSymb(1), 2, 14);  %2 + x1*(14-2);
            obj.radius = lerp(obj.xSymb(1), 40, 100);  %[m]
            obj.height = lerp(obj.xSymb(2), 6, 15);    % [m]
            obj.pitch = lerp(obj.xSymb(3), 5, 20) * pi / 180;     % [rad] (4-20 degrees, then converted to radians)
            
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = 1; 
            appNames = ["Wind Speed"];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.

            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary
            obj.vwind = lerp(obj.zSymb(1), 4, 6); % [m/s] %(prev 2-20)

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = length(perfMetricIDs);
            perfNames = ["Mass", "Power", "Lift", "Stiffness", "RadiusPenalty"];
            perfNames = perfNames(perfMetricIDs); % choose the ones you want to use
            obj.metricIDs = perfMetricIDs; 

            
            
            % ============ GENERIC SETUP ============
            obj.varNames = [desNames, appNames, perfNames];
            ext = "";
            for i=1:obj.rd
                ext = ext + "_" + strrep(perfNames(i)," ","_");
            end
            obj.visName = obj.problemName + ext;
            
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
            % Mass - just use volume
%             blade_thickness = 4;
%             blade_arch = obj.radius/11; % curvature in blade profile 
% 
%             blade_curve_radius = (obj.radius^2 + 4*blade_arch^2)/(8*blade_arch);
%             
%             sketch_area = blade_thickness *  blade_curve_radius * ...
%                             2 * asin(obj.radius / (2*blade_curve_radius)); % area of a single slice
%                         
%             unit_extrusion_area = 1/2 * obj.pitch * obj.radius^2;
%             f1 = sketch_area * unit_extrusion_area * obj.height;
%             min_mass = 1.8263e+04; % min value when all params =0
%             max_mass = 1.0701e+07; % max value when all params =1

            f1 = obj.radius * obj.height * 1/obj.pitch;
            min_mass = 229.1831; % 0,0,1
            max_mass = 2.1486e+04; %1,1,0

            % normalize so all possible outputs are between 0 and 1
            f1 = (f1 - min_mass) / (max_mass - min_mass);
        end
        
        function f2 = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            % approximate power output
            betz = false;
            
            fluidDensity = 1; % can really disregard, assumed constant
            sweptArea = pi * obj.radius^2;
            
            if betz
                Cp = 16/27; % BetzCoeff; theoretical limit, for ideal only
                min_power = -2.0106e+06; % for experiment with wind in 4-6m/s
                max_power = -9.5318e+04;
            else
                Cp = obj.performanceCoefficient();
                min_power = -1096672.537507; % for experiment with wind in 4-6m/s
                max_power = 28929.459366;
            end

            
            P_out = Cp * 1/2 * fluidDensity * sweptArea * obj.vwind^3; 
            f2 = -P_out; % want to maximize power, so we minimize the negative

            % normalize so all possible outputs are between 0 and 1
            f2 = (f2 - min_power) / (max_power - min_power);
        end
        
        function cp = performanceCoefficient(obj)
            TSR = obj.tipSpeedRatio();
            degree = 6;
            p = cp_calc(degree, false); %8th order polynomial fit to empirical Cp data, don't plot
            
            evald = 0;
            for i=1:degree+1
                evald = evald + p(i)*TSR^(degree - (i-1));
            end
            
            cp = evald;
        end
        
        function tsr = tipSpeedRatio(obj)
            % function of tip-speed ratio (TSR, aka lambda)
            rpm = 8; %chosen to preserve tip ratio range
            omegam = 2 * pi * rpm / 60;
            tsr = obj.radius * omegam / obj.vwind;
        end
        
        function f3 = performanceMetric3(obj)
            % PerformanceMetric2 - define the performance metric.
            % drag proxy - encourage large radius / height, small pitch
            fluidDensity = 1; % can really disregard, assumed constant
            
%             unit_extrusion_area = 1/2 * obj.pitch * obj.radius^2;
%             frontalSurfaceArea = unit_extrusion_area * obj.height;
            frontalSurfaceArea = obj.radius * obj.height * 1/obj.pitch;
            
            lift = 1/2 * fluidDensity * frontalSurfaceArea * obj.vwind^2;
            f3 = lift; % want to maximize lift, so we minimize the negative
            
            % normalize so all possible outputs are between 0 and 1
            max_lift = 3.0940e+05; % min value when params are 1,1,0,1 (bc negated)
            min_lift = 5.5004e+03; % max value when params are 0,0,1,0 (bc negated)
            f3 = (f3 - min_lift) / (max_lift - min_lift);
        end
        
        function f4 = performanceMetric4(obj)
            % PerformanceMetric2 - define the performance metric.
            % blade stiffness proxy - aspect ratio. want a low one.
%             unit_extrusion_area = 1/2 * obj.pitch * obj.radius^2;
%             frontalSurfaceArea = unit_extrusion_area * obj.height;

            frontalSurfaceArea = obj.radius * obj.height * 1/obj.pitch;
            
            aspectRatio = obj.radius^2 / frontalSurfaceArea;
            f4 = aspectRatio;
            
            % normalize so all possible outputs are between 0 and 1
            min_stiffness = 0.3820; % min value when all params =1 (bc tall)
            max_stiffness = 14.3239; % max value when all params =0 (bc blades thin, flat (lots of resistance))
%             f4 = (f4 - min_stiffness) / (max_stiffness - min_stiffness);
            
            % unclear -- think more about this one. Penalizes surface area
            % instead of incentivizing it
        end
        
        function f5 = performanceMetric5(obj)
            % penalize large radius
            exponent = 3;
            f5 = obj.radius^exponent; %want to minimize radius for cost/transportation
            
            % normalize so all possible outputs are between 0 and 1
            min_penalty = 40^exponent; % min value when radius "0" (40)
            max_penalty = 100^exponent; % max value when radius "1" (100)
            f5 = (f5 - min_penalty) / (max_penalty - min_penalty);
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
        end
        
        function [A, b] = matrixLinearConstraints(obj)
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
            
            firstPerfIdx = obj.rD + obj.rA + 1;
            pName = strrep(obj.problemName, " ", "_");

            f1 = obj.performanceMetric1();
            f2 = obj.performanceMetric2();
            f3 = obj.performanceMetric3();
            f4 = obj.performanceMetric4();
            f5 = obj.performanceMetric5();
            f = {f1, f2, f3, f4, f5};
            
            rf = []; rderiv = [];
            
            for i=1:obj.rd
                fname = pName + "_" + strrep(obj.varNames(firstPerfIdx + i-1)," ","_");
                fhandle = f{obj.metricIDs(i)}; %pick out the correct metric
                
                [fun, ffunc, deriv] = obj.setupFunction(fhandle, fname, createFiles);
                
                rfun{i} = fun;
                rf = [rf, ffunc];
                rfuncsCell{i} = ffunc;
                rderiv = [rderiv, deriv];
            end
            rhess = {};
        end
        
        function [fun, f, deriv] = setupFunction(obj, fHandle, fname, createFile)
            if createFile
                matlabFunction(fHandle, 'File', char(obj.wd + fname));
                obj.fixFileInputs(fname, obj.wd);
            end
            x = [obj.xSymb];
            z = [obj.zSymb];
            
            f_func = str2func(fname);
            fun = f_func;
            f(x, z) = f_func(x.', z.');
            deriv([x,z]) = gradient(f, [x,z]);
        end
        
        % shouldn't need to edit
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


