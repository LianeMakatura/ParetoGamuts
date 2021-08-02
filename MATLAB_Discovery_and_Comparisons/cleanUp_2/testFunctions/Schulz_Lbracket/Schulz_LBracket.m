classdef Schulz_LBracket
    %LBRACKET Summary of this class goes here
    %  for now, assume both arms of L bracket are symmetric
    %  
    %    v1 -- v2
    %     |     |
    %     |     |
    %     |    v4 ------ v6
    %     |              |
    %     v3 ---------- v5    same order for far side, verts v7-12
    %
    %   DESIGN VARS
    %   length = length of arms, eg from v3-v1 or v5-v3
    %   width = distance between front and back (eg between v3-v9)
    %   thickness = distance between eg v1 and v2
    %
    %   APP VARS
    %   angle = angle between outer arm plates; given in radians 
    %
    %   PERFORMANCE METRICS
    %   mass
    %   [tbd]

    
    properties
        % general
        ref_pt
        mesh
        constrainedZ
        zValFixed
        
        % design vars
        rD
        xSymb
        length
        width
        thickness
        
        % application vars
        rA
        zSymb
        angle
        
        % performance metrics
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
        function obj = Schulz_LBracket(basewd, threeD_bool, constrained_bool, zValFixed)
            %LBRACKET Construct an instance of this class
            %   @param constrained_bool - whether app vars are constrained
            %   @param zValFixed - value of app vars *if constrained*
            obj.constrainedZ = constrained_bool;
            obj.ref_pt = [0,0,0];
            
            obj.problemName = "Schulz_LBracket";
            desNames = ["Arm length", "Arm width", "Arm thickness"];
            appNames = ["Bend Angle"];
            perfNames = ["Mass", "Utopia Distance"];
            
            % DESIGN VARIABLES ============= 
            if threeD_bool
                obj.rD = 3;
                obj.visName = "Lbracket_3D";
            else
                obj.rD = 2;
                desNames = desNames(1:2);
                obj.visName = "Lbracket_2D";
            end
            obj.varNames = [desNames, appNames, perfNames];
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % create aliases for the variables so they're meaningful
            % ensure it's longer than it is thick (no inverted tris)
            obj.length = 0.4 + (obj.xSymb(1) * 0.6); % allowed from 0.4-1, 0 nonsense
            obj.width = 0.2 + (obj.xSymb(2) * 0.8); %allowed from 0.2-1, 0 is nonsense
            if threeD_bool
                obj.thickness = 0.1 + (obj.xSymb(3)*0.2); %allowed from 0.1-0.3
            else
                obj.thickness = 0.1;
            end
            
            % APPLICATION VARIABLES =========
            eps = 1e-3;
            if ~obj.constrainedZ %unconstrained, use symbolic z
                obj.rA = 1;
                obj.zSymb = sym('z', [1 obj.rA]);

                % create alias 
                obj.angle = obj.zSymb(1);
                
                % normalize so zSymb between 0-1 gives angle pi/2 to pi/4
                % eps offsets divie by zero errors that happen at pi/4,
                % pi/2
                obj.angle = pi/4 + eps + (obj.angle * (pi/4 - 2*eps));
            else
                obj.rA = 0;
                obj.zSymb = [];
                obj.angle = pi/4 + eps + (zValFixed * (pi/4-2*eps));
                obj.zValFixed = zValFixed;
            end
            
            % create the mesh object
            obj.mesh = Mesh(obj.calcSymbolicVertices_arbAngle(), ...
                            obj.createTopology());
            
            probwd = "testFunctions/" + obj.problemName + "/";
            obj.wd = basewd + probwd; 
            
            obj.hasGroundTruth = false;
            obj.hasMeshes = true;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
        end
        
        function vertsSymb = calcSymbolicVertices_rightAngle(obj)
            % Calculates the vertex positions wrt the symbolic design parameters
            % @return vertsSymb -- 12x3 matrix of symb verts, row=vertex
            
            v9 = obj.ref_pt;                        % far lower left
            v3 = obj.ref_pt + [obj.width, 0, 0];        % near lower left

            v1 = v3 + [0, 0, obj.length];           % near upper left
            v7 = v9 + [0, 0, obj.length];           % far upper left 

            v5 = v3 + [0, obj.length, 0];           % near lower right
            v11 = v9 + [0, obj.length, 0];          % far lower right

            v6 = v5 + [0, 0, obj.thickness];        % near center right
            v12 = v11 + [0, 0, obj.thickness];      % far center right

            v2 = v1 + [0, obj.thickness, 0];        % near upper center
            v8 = v7 + [0, obj.thickness, 0];        % far upper center

            v4 = v3 + [0, obj.thickness, obj.thickness]; % near center center
            v10 = v9 + [0, obj.thickness, obj.thickness]; % far center center

            vertsSymb = [v1; v2; v3; v4; v5; v6; v7; v8; v9; v10; v11; v12];
        end
        
        function vertsSymb = calcSymbolicVertices_arbAngle(obj)
            rotMat = rot([1, 0, 0], obj.angle); % rotation matrix about x axis

            offset = [obj.width, 0, 0];
            % define invariant side
            v3 = obj.ref_pt + offset;
            v5 = v3 + [0, obj.length, 0];
            v6 = v5 + [0, 0, obj.thickness];

            % define pts to be rotated
            v1 = v5;
            v2 = v1 - [0, 0, obj.thickness];

            % rotate the "vertical" beam as desired
            v1 = (rotMat * v1')';
            v2 = (rotMat * v2')';

            % compute location of v4; intersection of inner bracket faces
            % we have 1 known pt on each, and direction parallel to outer face
            % only in-plane (yz) so pass last 2 coords
            m_rot = slope2D(v1(2:3), v3(2:3));
            m_horz = slope2D(v5(2:3), v3(2:3)); % should always be 0
            [y,z] = intersectLines2D(v2(2:3), m_rot, v6(2:3), m_horz);

            v4 = [obj.width, y, z];

            % extrude into far plane
            v7 = v1 - offset;
            v8 = v2 - offset;
            v9 = v3 - offset;
            v10 = v4 - offset;
            v11 = v5 - offset;
            v12 = v6 - offset;

            vertsSymb = [v1; v2; v3; v4; v5; v6; v7; v8; v9; v10; v11; v12];
        end
        
        function faces = createTopology(obj)
            faces = [1, 4, 2;    % near side
                    4, 1, 3;
                    4, 3, 5;
                    4, 5, 6;
                    7, 8, 10;   % far side
                    7, 10, 9;
                    9, 10, 11;
                    10, 12, 11;
                    1, 7, 3;    % back
                    3, 7, 9;
                    10, 8, 2;   % front (inner bracket face)
                    10, 2, 4;
                    11, 12, 6;  % front (edge)
                    11, 6, 5;
                    9, 11, 5;   % bottom
                    9, 5, 3;    
                    8, 7, 1;    % top (edge)
                    8, 1, 2;
                    12, 10, 4;  % top (inner bracket face)
                    12, 4, 6;
                    ];
        end
        
        function m = mass(obj)
            %MASS return the mass of the design
            %   We actually just use volume for now (equiv. for
            %   homogeneous)
            m = obj.mesh.getVolume();
            
            % normalize so all outputs between 0, 1 (needed by old code)
            % already in range, because can't exceed unit cube
        end
        
        function rcom = relativeCenterOfMass(obj)
            % measure distance between current (appx) center of mass and
            % the utopia point
            utopia = [1, 1, 1];
            diff = utopia - obj.mesh.getAvgVertPosition();
            rcom = sum(diff.^2);
            
            % normalize so all outputs between 0, 1 (needed by old code)
            % extremes are 0 (perfect match, must be positive) and  3,
            % since design space is limited so (0,0,0) worst possible
            rcom = rcom / 3;
        end
        
        
        % ========== GROUND TRUTH SOLUTION ===================
        % for extablished benchmark functions, this "ground truth" is only
        % for the fixed context situation; here, it should be the lower
        % envelope (if the application var is just one of the original
        % design vars) or some particular front (if the app var was originally
        % a constant, and we're varying it)
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
        
        
        % ========== WRITE AND MODIFY THE PROBLEM FILES ===================
        function [rfun, rf, rfuncsCell, rderiv, rhess] = createProblem(obj, createFiles)
            % create the usable matlab function files from symbolic equiv.
            % defined in this problem class
            % only write the files once each time the functions or underlying
            % parametrizations change 

            f1_name = 'f1_mass_';
            f2_name = 'f2_rcom_';
            if obj.constrainedZ
                x = [obj.xSymb];
                z = sym('q', [1, 1]); % fake sym variable, needed for function sig
                
                disp(obj.zValFixed)
                disp(mod(obj.zValFixed, 1))
                if mod(obj.zValFixed, 1)==0 %integer
                    zvalstr = num2str(obj.zValFixed);
                else
                    zvals = strsplit(num2str(obj.zValFixed), '.');
                    zvalstr = strcat(zvals{1}, '_', zvals{2});
                end
                suffix = strcat('constrained_z_', zvalstr);     % these will change for every app var value
            else
                x = [obj.xSymb];
                z = [obj.zSymb];
                suffix = 'unconstrained';
            end
            
            f1_name = strcat(f1_name, suffix);
            f2_name = strcat(f2_name, suffix);
            if createFiles           
                f1 = obj.mass();
                f2 = obj.relativeCenterOfMass();
                
                matlabFunction(f1, 'File', strcat(obj.wd, f1_name));
                matlabFunction(f2, 'File', strcat(obj.wd, f2_name));
                
                obj.fixFileInputs(f1_name, obj.wd);
                obj.fixFileInputs(f2_name, obj.wd);
            end
            
            f1_func = str2func(f1_name);
            fun1 = f1_func;
            f1(x, z) = f1_func(x.', z.');
            deriv1([x,z]) = gradient(f1, [x,z]);
            
            f2_func = str2func(f2_name);
            fun2 = f2_func;
            f2(x, z) = f2_func(x.', z.');
            deriv2([x,z]) = gradient(f2, [x,z]);
            
            % create structures to return
            rfun = {fun1, fun2};
            rf = [f1, f2];
            rfuncsCell = {f1, f2};
            rderiv = [deriv1, deriv2];
            
            rhess = {}; %use on the fly anyway
        end
        
        function [] = fixFileInputs(obj, filename, wd)
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
            fprintf(fileh, '%% Includes scripted file modifications (LBracket.fixFileInputs) to match system specs\n\n');
            
            % parse the variable blocks into individual sym. vars
            fprintf(fileh, '%% unpack the block vars into the variable names used by the function\n');
            for j=1:obj.rD
                var_id = strcat('x', num2str(j));
                fprintf(fileh, '%s = x(%d, :).'';\n', var_id, j);
            end
            if ~obj.constrainedZ
                for j=1:obj.rA
                    var_id = strcat('z', num2str(j));
                    fprintf(fileh, '%s = z(%d, :).'';\n', var_id, j);
                end
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
        %writes out the OBJ files for the visualization tool.
        %   Detailed explanation goes here
            numMeshes = size(pts, 1);
            for i=1:numMeshes
                p = pts(i, :);

                % generate the target geometry file name
                fname_fmt = '';
                for j=1:obj.rD+obj.rA
                    fname_fmt = [fname_fmt '%0.6f '];
                end
                fname_fmt = strcat(fname_fmt, '.obj');
                if obj.rD == 2
                    filename = sprintf(fname_fmt, p(1), p(2), p(3));
                else
                    filename = sprintf(fname_fmt, p(1), p(2), p(3), p(4));
                end
                destFolder = obj.wd + "meshes/";

                MESHout = strcat(destFolder, filename);
                
                % write out file
                obj.writeOBJ(p(1:obj.rD), p(obj.rD+1:obj.rD+obj.rA), MESHout);
            end
        
        end
        
        function writeOBJ(obj, x, z, obj_pathAndName)
            % write out a specific instance of the LBracket
            % **ONLY applies to symbolic instances of LBracket
            % x is (cell) array of design variables [x1, x2, x3 (maybe)]
            % z is (cell) array of app variables [z1]
            
            % use cells to feed them in, convert if necessary
            if iscell(x)
                x = cell2mat(x);
                z = cell2mat(z);
            end
            
            if ~(length(x) == obj.rD)
                disp('Not enough design inputs!')
            end
            if ~(length(z) == obj.rA)
                disp('Not enough app inputs!')
            end
            
            vert_func = matlabFunction(obj.mesh.vertices); % make symb vertices into callable matlab func
            
            if obj.rD == 2
                vertices = vert_func(x(1), x(2), z);
            else 
                vertices = vert_func(x(1), x(2), x(3), z);
            end
            
            m = Mesh(vertices, obj.mesh.faces); % should avoid this
            write_obj(m, obj_pathAndName, 'w');
        end

    end %methods
end % class

