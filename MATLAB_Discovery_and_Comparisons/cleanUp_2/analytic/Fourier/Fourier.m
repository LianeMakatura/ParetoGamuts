classdef Fourier
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
        order
        
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
        function obj = Fourier(rD, rA, rd, polyOrder, basewd)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "Fourier"; %shouldn't have any spaces or special chars except underscore
            obj.hasGroundTruth = false;
            obj.hasMeshes = false;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
            
            % working directory; make sure it ends with a "/"
            probwd = "analytic/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = rD; 
            desNames = []; % for plot labels and filenames
            for j=1:obj.rD
                var_id = "x" + num2str(j);
                desNames = [desNames, var_id];
            end
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = lerp(obj.xSymb(1), 2, 14);  %2 + x1*(14-2);
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = rA; 
            appNames = [];
            for j=1:obj.rA
                var_id = "z" + num2str(j);
                appNames = [appNames, var_id];
            end
            obj.zSymb = sym('z', [1 max(1, obj.rA)]); % need to have at least 1 symbolic var for the framework to be happy. Just don't use it if there shouldn't be one.

            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = rd;
            perfNames = [];
            for j=1:obj.rd
                var_id = "f" + num2str(j);
                perfNames = [perfNames, var_id];
            end
            obj.order = polyOrder;
            
            
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

        % just defined directly with the fourierFancyZ functions in
        % obj.createProblem
        
        
        
        

        
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

            f = cell(1, obj.rd);
            rng(3)    % should result in same random coeffs
            for i=1:obj.rd
                f{i} = fourierFancyZ(obj.order, obj.rD, obj.rA);
            end
            
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


