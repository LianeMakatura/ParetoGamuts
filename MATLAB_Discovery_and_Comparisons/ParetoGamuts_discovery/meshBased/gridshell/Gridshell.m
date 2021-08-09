classdef Gridshell
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
        pts
        numX
        numY
        N
        
        % ==== application vars
        rA
        zSymb
        %+alias vars for zSymb (eg, obj.angle)
        morningSun
        eveningSun
        houseTheta
        regularizerWeight
        
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
        function obj = Gridshell(basewd, numX, numY, perfMetricIDs)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "Gridshell"; %shouldn't have any spaces or special chars except underscore
            obj.hasGroundTruth = false;
            obj.hasMeshes = true;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
            
            % working directory; make sure it ends with a "/"
            probwd = "meshBased/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            
            % ============ DESIGN VARIABLES ============ 
            obj.numX = numX;
            obj.numY = numY;
            obj.N = obj.numX * obj.numY;
            numFixedBoundaryPts = obj.N - (numX-2)*(numY-2);
            
            obj.rD = obj.N - numFixedBoundaryPts; 
            desNames = []; % for plot labels and filenames
            for j=1:obj.numY-2
                for i=1:obj.numX-2
                    desNames = [desNames, "z" + num2str(i) + "," + num2str(j)];
                end
            end
            obj.xSymb = sym('x', [1 obj.rD]);
            variablePts = reshape(lerp(obj.xSymb, 0, 5), numX-2, numY-2);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = lerp(obj.xSymb(1), 2, 14);  %2 + x1*(14-2);
            
            % ========= APPLICATION VARIABLES ==========
            obj.morningSun=[-0.5, -0.5, -0.6];
            obj.eveningSun = [-0.5, 0.5, -0.6];
            obj.regularizerWeight = 0.2;
            obj.rA = 1; 
            appNames = ["House Orientation"];
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.
            
            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary
            obj.houseTheta = lerp(obj.zSymb(1), 0, 2*pi);
            
            % ====== combine into the full points array 
            obj.pts = sym(zeros(obj.numX, obj.numY));
            obj.pts(2:obj.numX-1, 2:obj.numY-1) = variablePts;
            
            
            % ========== PERFORMANCE METRICS ============
            obj.rd = 2; 
            perfNames = ["Morning Power Output", "Evening Power Output"];
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
        % ==================================================

        function f1 = performanceMetric1(obj)
            % PerformanceMetric1 - define the performance metric.
            f1 = obj.powerOutput(obj.morningSun) + obj.regularizerWeight * obj.smoothness();
            
            % normalize so all possible outputs are between 0 and 1
        end
        
        function f2 = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            f2 = obj.powerOutput(obj.eveningSun) + obj.regularizerWeight * obj.smoothness();
            
            % normalize so all possible outputs are between 0 and 1
        end
        
        function p = powerOutput(obj, sun_dir)
            sun_dir = sun_dir / norm(sun_dir);
            % rotate normalized direction by -theta to account for house
            % orientation
            % transpose to col vec then back to row vec
            sun_dir = obj.rotateAboutZ(-obj.houseTheta, sun_dir')';
            
            % each quad split into 2 triangles
            %normals for tri A -- use upper left of vec1,2
            vec1 = obj.pts(:, 2:end) - obj.pts(:, 1:(end-1)); %(i, j+1)-(i,j)
            vec2 = obj.pts(2:end,:) - obj.pts(1:(end-1),:); %(i+1, j)-(i,j)
            incidenceA = sym(zeros(obj.numX-1, obj.numY-1));% matrix of normals
            for i=1:obj.numX-1
                for j=1:obj.numY-1
                    n = cross([0,1,vec1(i,j)], [1,0, vec2(i,j)]);
                    n = n / norm(n);
                    incidenceA(i,j) = 1 - dot(n, sun_dir)^2; % 1 - dot^2 bc want dot to be +-1, need to minimize f
                end
            end
            
            %normals for tri B -- use lower right of vec1,2
            vec1 = -vec1(2:end, :); %(i, j-1)-(i,j), chop first row
            vec2 = -vec2(:, 2:end); %(i-1, j)-(i,j), chop first col
            incidenceB = sym(zeros(obj.numX-1, obj.numY-1));% matrix of normals
            for i=1:obj.numX-1
                for j=1:obj.numY-1
                    n = cross([0,-1,vec1(i,j)], [-1,0, vec2(i,j)]);
                    n = n / norm(n);
                    incidenceB(i,j) = 1 - dot(n, sun_dir)^2; % 1 - dot^2 bc want dot to be +-1, need to minimize f
                end
            end
            
            p = sum(incidenceA(:)) + sum(incidenceB(:)); 
            p = p / (2*numel(incidenceA)); % each tri can contribute up to 1
        end
        
        function s = smoothness(obj)
            % SMOOTHNESS - wants second derivative to be small
            diff_y = (obj.pts(2:end,:) - obj.pts(1:(end-1),:))*obj.numY;
            diff_y2 = (diff_y(2:end,:) - diff_y(1:(end-1),:))*obj.numY;
            
            diff_x = (obj.pts(:, 2:end) - obj.pts(:, 1:(end-1)))*obj.numX;
            diff_x2 = (diff_x(:, 2:end) - diff_x(:, 1:(end-1)))*obj.numX;

            s = sum(diff_y2(:).^2)/obj.numY + sum(diff_x2(:).^2)/obj.numX;% scale by size of x so that it's invariant to N
            % normalize so all possible outputs are between 0 and 1
            s = s / 1000000;
        end
        
        % vec must be column vector
        function vecOut = rotateAboutZ(obj, theta, vec)
            rotZ = [cos(theta), -sin(theta), 0;
                    sin(theta),  cos(theta), 0;
                       0,           0,       1];
            vecOut = rotZ * vec;
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

            optimizeMatlabFuncs = false;
            firstPerfIdx = obj.rD + obj.rA + 1;
            pName = strrep(obj.problemName, " ", "_");

            f1 = obj.performanceMetric1();
            f2 = obj.performanceMetric2();
            f = {f1, f2};
            
            rf = [];
            
            for i=1:obj.rd
                fname = pName + "_" + strrep(obj.varNames(firstPerfIdx + i-1)," ","");
                fhandle = f{obj.metricIDs(i)}; %pick out the correct metric
                
                [fun, ffunc, hess] = obj.setupFunction(fhandle, fname, createFiles, optimizeMatlabFuncs);
                
                rfun{i} = fun;
                rf = [rf, ffunc]; % used for scalarization function
                rfuncsCell{i} = ffunc;
                rhess{i} = hess;
            end
            
            % write out the full jacobian matrix as a file function too
            deriv = jacobian(rfuncsCell, [obj.xSymb, obj.zSymb]); 
            deriv_name = pName + "_deriv";
            if createFiles
                matlabFunction(deriv, 'File', char(obj.wd +"generated/"+ deriv_name), 'vars',{obj.xSymb.', obj.zSymb.'}, 'Optimize', optimizeMatlabFuncs);
            end
            rderiv = str2func(deriv_name);
        end
        
        
        % shouldn't need to edit
        function [fun, f, hessh] = setupFunction(obj, fHandle, fname, createFile, optimizeMatlabFuncs)
            if createFile
                dispstat(sprintf('writing matlab function for %s', fname),'keepthis','keepprev','timestamp');
                matlabFunction(fHandle, 'File', char(obj.wd + "generated/"+fname), 'vars',{obj.xSymb.', obj.zSymb.'},... 
                                'Optimize', optimizeMatlabFuncs, ...
                                'Comments', sprintf('Automatically generated by %s.createProblem()', obj.problemName));
                dispstat(sprintf('done writing matlab function for %s', fname),'keepthis','keepprev','timestamp');
            end
            x = [obj.xSymb];
            z = [obj.zSymb];
            
            f_func = str2func(fname);
            fun = f_func;
            f(x, z) = f_func(x.', z.');
            
            % try with matlab function handles
            hess([x,z]) = hessian(f, [x,z]);
            hess_name = fname + "_hess";
            if createFile
                dispstat(sprintf('writing hess matlab function for %s', fname),'keepthis','keepprev','timestamp');
                matlabFunction(hess, 'File', char(obj.wd + "generated/"+hess_name), 'vars',{x.', z.'}, ...
                                'Optimize', optimizeMatlabFuncs);
                dispstat(sprintf('done writing hess matlab function for %s', fname),'keepthis','keepprev','timestamp');
            end
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


