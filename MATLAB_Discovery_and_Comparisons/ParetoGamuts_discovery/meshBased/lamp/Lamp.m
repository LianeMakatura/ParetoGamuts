classdef Lamp
    % CAD lamp example
    %
    %   DESIGN VARS (21)
    %       [base_x, base_y, base_z] - vector to end of the "base" beam (central support to ground)
    %
    %       [arm1_x, arm1_y, arm1_z] - vector to end of each "arm" beam(connected to base)
    %       [arm2_x, arm2_y, arm2_z]   defined relative to the end effector
    %       [arm3_x, arm3_y, arm3_z]   of the base
    %
    %       [hand1_x, hand1_y, hand1_z] - vector to end of each "hand" beam (connected to arm)
    %       [hand2_x, hand2_y, hand2_z]   defined relative to the end effector
    %       [hand3_x, hand3_y, hand3_z]   of the ancestor arm
    %
    %   APP VARS
    %       focalPoint_z - coordinate specifying the location of the focal
    %                       point, relevant to the third performance metric
    % 
    %   PERFORMANCE METRICS
    
    
    properties
        % ==== design vars
        rD
        xSymb
        %+alias vars for xSymb (eg, obj.length, width)
        base
        a1
        a2
        a3
        h1
        h2
        h3
        
        % fixed design variables (not changing)
        radii
        bMax
        aMax
        hMax
        wMax
        endVolume
        beamRanges
        
        % ==== application vars
        rA
        zSymb
        %+alias vars for zSymb (eg, obj.angle)
        focalPoint_z
        focalPoint

        % ==== performance metrics
        rd 
        metricIDs
        
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
        function obj = Lamp(basewd, perfMetricIDs)
            %PROBLEMTEMPLATE - constructor for the class. Named same as
            % class.
            %   Detailed explanation goes here
            obj.problemName = "Lamp"; %shouldn't have any spaces or special chars except underscore
            
            % ============ DESIGN VARIABLES ============ 
            obj.rD = 21; 
            desNames = ["base_x", "base_y", "base_z", ...
                        "arm1_x", "arm1_y", "arm1_z", ...
                        "arm2_x", "arm2_y", "arm2_z", ...
                        "arm3_x", "arm3_y", "arm3_z", ...
                        "hand1_x", "hand1_y", "hand1_z", ...
                        "hand2_x", "hand2_y", "hand2_z", ...
                        "hand3_x", "hand3_y", "hand3_z", ...
                        ]; % for plot labels and filenames
            obj.xSymb = sym('x', [1 obj.rD]);
            
            % Create aliases for the variables so they're meaningful
            % also, scale as necessary so the range xi=[0,1] yields correct
            %       range for design vars
            % eg.   obj.length = 0.4 + (obj.xSymb(1) * 0.6); % allowed from 0.4-1
            
            % set the maximum dimensions for each component type
            obj.bMax = 4; %2;
            obj.aMax = 2;
            obj.hMax = 1;
            obj.wMax = 0.05;
            obj.endVolume = 0.3.^3;

            % vector definition (start point is predecessor endpt)
            obj.beamRanges = struct('base', struct('x', [0,1], 'y', [0,1], 'z', [1,obj.bMax]),...
                                    'a1', struct('x', [0,1], 'y', [0,1], 'z', [1,2]), ...
                                    'a2', struct('x', [-1,0], 'y', [0,1], 'z', [1,2]), ...
                                    'a3', struct('x', [0,1], 'y', [-1, 0], 'z', [1,2]), ...
                                    'h1', struct('x', [0.4, 1], 'y', [0.4, 1], 'z', [-2,-1]), ...
                                    'h2', struct('x', [-1, -0.4], 'y', [0.4, 1], 'z', [-2,-1]), ...
                                    'h3', struct('x', [0.4, 1], 'y', [-1, -0.4], 'z', [-2,-1]) ...
                                );
            obj = obj.rescaleAndAliasDesignVars();
            
            
            % ========= APPLICATION VARIABLES ==========
            obj.rA = 1; 
            appNames = ["focal_point_z"]; %will save x and y for later
            obj.zSymb = sym('z', [1 obj.rA]);
            
            % Create aliases for the variables so they're meaningful, 
            % and scale as necessary
            scale = obj.bMax/2+obj.aMax; %obj.bMax+obj.aMax;
            obj.focalPoint_z = (1 + obj.zSymb(1)*4) * scale; %between (3-5)*scale; original paper had 4*scale  
            obj.focalPoint = [2*scale, ...
                              2*scale, ...
                              obj.focalPoint_z];

            
            % ========== PERFORMANCE METRICS ============
            obj.rd = length(perfMetricIDs);
            perfNames = ["Stability", "Mass", "Focal Point"];
            perfNames = perfNames(perfMetricIDs); % choose the ones you want to use
            obj.metricIDs = perfMetricIDs;
            
            % ============ GENERIC SETUP ============
            obj.varNames = [desNames, appNames, perfNames];
            ext = "";
            for i=1:obj.rd
                ext = ext + "_" + strrep(perfNames(i)," ","_");
            end
            obj.visName = obj.problemName + ext;

            
            % working directory; make sure it ends with a "/"
            probwd = "meshBased/" + obj.problemName + "/";
            obj.wd = basewd + probwd;   
            obj.hasGroundTruth = false;
            obj.hasMeshes = true;
            obj.hasLinearConstraints = false;
            obj.useFD = false;
            obj.symFreeAnalytic = false;
        end
        
        % =========== SETUP DESIGN + APP VARIABLES ================
        function obj = rescaleAndAliasDesignVars(obj)
            r = obj.beamRanges;
            obj.base = [obj.rescaleParam(r.base.x(1), r.base.x(2), obj.bMax, 1), ...
                        obj.rescaleParam(r.base.y(1), r.base.y(2), obj.bMax, 2) , ...
                        obj.rescaleParam(r.base.z(1), r.base.z(2), obj.bMax, 3) ];
            
            obj.a1 = [obj.rescaleParam(r.a1.x(1), r.a1.x(2), obj.aMax, 4), ...
                      obj.rescaleParam(r.a1.y(1), r.a1.y(2), obj.aMax, 5), ...
                      obj.rescaleParam(r.a1.z(1), r.a1.z(2), obj.aMax, 6) ];
            
            obj.a2 = [obj.rescaleParam(r.a2.x(1), r.a2.x(2), obj.aMax, 7), ...
                      obj.rescaleParam(r.a2.y(1), r.a2.y(2), obj.aMax, 8) , ...
                      obj.rescaleParam(r.a2.z(1), r.a2.z(2), obj.aMax, 9) ];
            
            obj.a3 = [obj.rescaleParam(r.a3.x(1), r.a3.x(2), obj.aMax, 10), ...
                      obj.rescaleParam(r.a3.y(1), r.a3.y(2), obj.aMax, 11) , ...
                      obj.rescaleParam(r.a3.z(1), r.a3.z(2), obj.aMax, 12) ];
            
            obj.h1 = [obj.rescaleParam(r.h1.x(1), r.h1.x(2), obj.hMax, 13), ...
                      obj.rescaleParam(r.h1.y(1), r.h1.y(2), obj.hMax, 14) ,...
                      obj.rescaleParam(r.h1.z(1), r.h1.z(2), obj.hMax, 15) ];
            
            obj.h2 = [obj.rescaleParam(r.h2.x(1), r.h2.x(2), obj.hMax, 16), ...
                      obj.rescaleParam(r.h2.y(1), r.h2.y(2), obj.hMax, 17) , ...
                      obj.rescaleParam(r.h2.z(1), r.h2.z(2), obj.hMax, 18) ];
            
            obj.h3 = [obj.rescaleParam(r.h3.x(1), r.h3.x(2), obj.hMax, 19), ...
                      obj.rescaleParam(r.h3.y(1), r.h3.y(2), obj.hMax, 20) , ...
                      obj.rescaleParam(r.h3.z(1), r.h3.z(2), obj.hMax, 21) ];

            obj.radii = ones(1,4)*obj.wMax*2;
        end
        
        % Remap the symbolic design vars (in 0-1) to the proper variable
        % value range for the problem
        function scaledParam = rescaleParam(obj, minVal, maxVal, scaleFactor, xID)
            possibleRange = maxVal - minVal;
            rangeToUse = obj.xSymb(xID) * possibleRange; % x in [0,1] is a percentage
            scaledParam = (minVal + rangeToUse) * scaleFactor; 
        end
        
        
        % =========== DEFINE PERFORMANCE METRICS ===========
        function f1 = performanceMetric1(obj)
            % PerformanceMetric1 - define the performance metric.
            %       Name may be altered, only called from inside the class.
            %       If altered, change createProblem assignments to f1, f2
            %   STABILITY - distance between projected center of mass (COM), and
            %   the center of the base. Assume center of mass is origin, so
            %   just (COM_x^2 + COM_y^2)
            
            % for each element, compute mass (just volume, bc density
            % constant), then treat as point mass at COM of the element.
            % This COM is the geometric center since constant density,
            % radially symmetric.
            masses = [norm(obj.base)*obj.radii(1).^2,...
                    norm(obj.a1)*obj.radii(2).^2, ...
                    norm(obj.a2)*obj.radii(3).^2, ...
                    norm(obj.a3)*obj.radii(4).^2, ...
                    norm(obj.h1)*obj.radii(2).^2, ...
                    norm(obj.h2)*obj.radii(3).^2, ...
                    norm(obj.h3)*obj.radii(4).^2
                    ];
            
            centerOfMass = obj.base*0.5* masses(1) + ...
               (obj.base + obj.a1*0.5)* masses(2) +...
               (obj.base + obj.a2*0.5)* masses(3) + ...
               (obj.base + obj.a3*0.5)* masses(4) + ...
               (obj.base + obj.a1 + obj.h1*0.5)* masses(5)+ ...
               (obj.base + obj.a2 + obj.h2*0.5)* masses(6) + ...
               (obj.base + obj.a3 + obj.h3*0.5)* masses(7) + ...
               (obj.base + obj.a1 + obj.h1)*obj.endVolume + ...
               (obj.base + obj.a2 + obj.h2)*obj.endVolume + ...
               (obj.base + obj.a3 + obj.h3)*obj.endVolume;
           
            % normalize by the mass to get the location only
            totalMass = sum(masses) + 3 * obj.endVolume;
            centerOfMass = centerOfMass / totalMass; 
           
            f1 = (centerOfMass(1).^2 + centerOfMass(2).^2);   
           

            % normalize so all possible outputs are between 0 and 1
            f1 = (f1 / 15) / 3; % observed via jscad
            
        end
        
        
        function f2 = performanceMetric2(obj)
            % PerformanceMetric2 - define the performance metric.
            %       Name may be altered, only called from inside the class.
            %       If altered, change createProblem assignments to f1, f2
            % MASS - use volume, because density assumed constant
            
            % may need to change
            minMass = ( norm([0,0,1])*( obj.aMax*3 + obj.bMax) + norm([0.4,0.4,1] * obj.hMax*3) )*2*(obj.wMax);
            maxMass  = norm([1,1,obj.bMax])*(obj.hMax*3 + obj.aMax*3 + obj.bMax)*(2*obj.wMax);

            mass = norm(obj.base)*obj.radii(1)+ ...
                    norm(obj.a1)*obj.radii(2) + ...
                    norm(obj.a2)*obj.radii(3) + ...
                    norm(obj.a3)*obj.radii(4) + ...
                    norm(obj.h1)*obj.radii(2) + ...
                    norm(obj.h2)*obj.radii(3) + ...
                    norm(obj.h3)*obj.radii(4);

            % normalize so all possible outputs are between 0 and 1
            f2 = (mass -minMass)/(maxMass - minMass);
        end
        
        function f3 = performanceMetric3(obj)
           % average distance to the focal point
            endPosCost = norm(obj.base+obj.a1+obj.h1 -obj.focalPoint) ...
                        + norm(obj.base+obj.a2+obj.h2 -obj.focalPoint) ...
                        + norm(obj.base+obj.a3+obj.h3 -obj.focalPoint);
            endPosCost = endPosCost / 3;
                    
            % estimated in jscad viewer 
            f3 = endPosCost/20;%(endPosCost/3- 10)/ 8;
        end
        
        
        % ========== GROUND TRUTH SOLUTION ===================
        % for established benchmark functions, this "ground truth" is only
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
        
        % =========== WRITE AND ALTER METRIC FILES ===========
        function [rfun, rf, rfuncsCell, rderiv, rhess] = createProblem(obj, createFiles)
            % create the usable matlab function files from symbolic equiv.
            % defined in this problem class
            % Called from the mapping function
            
            firstPerfIdx = obj.rD + obj.rA + 1;
            pName = strrep(obj.problemName, " ", "_");

            f1 = obj.performanceMetric1(); %stability
            f2 = obj.performanceMetric2(); %mass
            f3 = obj.performanceMetric3(); %focal
            f = {f1, f2, f3};
            
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
            rhess = [];
        end
        
        function [fun, f, deriv] = setupFunction(obj, fHandle, fname, createFile)
            if createFile
                matlabFunction(fHandle, 'File', char(obj.wd + fname), 'vars',{obj.xSymb.', obj.zSymb.'},... 
                                'Comments', sprintf('Automatically generated by %s.createProblem()', obj.problemName));
            end
            x = [obj.xSymb];
            z = [obj.zSymb];
            
            f_func = str2func(fname);
            fun = f_func;
            f(x, z) = f_func(x.', z.');
            deriv([x,z]) = gradient(f, [x,z]);
        end
       
        
        function filename = generateFileName(obj, p)
            % generate the target geometry file name
            fname_fmt = '';
            for j=1:obj.rD+obj.rA
                fname_fmt = [fname_fmt '%0.6f '];
            end
            fname_fmt = strcat(fname_fmt, '.stl');
            
            filename = sprintf(fname_fmt, p(1:obj.rD + obj.rA));
        end
        
        function [] = writeMeshes(obj, pts, use_jscad)
            %WRITEMESHES Takes design points (and application points, if necessary) and
            %writes out the STL files for the visualization tool.
            %   Detailed explanation goes here
            
            if ~exist('use_jscad', 'var')
                use_jscad = true;
            end
            
            % read in the JSCAD file
            if use_jscad
                import matlab.net.*
                import matlab.net.http.*
                jscad_text = fileread(obj.wd + "lamp_template.jscad");
            end
            
            numMeshes = size(pts, 1);
            for i=1:numMeshes
                p = pts(i, :);
                
                MESHout = "mesh.stl";
                outfile = obj.wd + "meshes/" + MESHout;

                if use_jscad
                    jscad_result = obj.writeLampJSCAD(p, jscad_text);

                    tic;
                    body = matlab.net.http.MessageBody(jscad_result);
                    method = matlab.net.http.RequestMethod.POST;
                    header = [];
                    request = matlab.net.http.RequestMessage(method,header,body);

%                     uri = URI('http://internal.hanneshergeth.com:8034');
                    uri = URI('http://localhost:8034');
                    resp = send(request,uri);
                    toc;

                    if resp.StatusCode == 200 %code is OK
                        stl_data = resp.Body.Data;
                        meshfile = fopen(outfile, 'w');
                        fwrite(meshfile, stl_data);
                        fclose(meshfile);
                    else
                        disp("JSCAD status code not ok");
                        return
                    end
                else % use openSCAD
                    SCADout = obj.wd + "SCADout.scad";
                    obj.writeLampSCAD(p, SCADout);

                    % call openSCAD to generate design + write out file
                    scadpath = "/usr/local/bin/openscad";               
                    system(sprintf("%s -o %s %s", scadpath, outfile, SCADout));
                end
                
                % rename the file to fit
                finalname = obj.generateFileName(p);
                finalout = obj.wd + "meshes/" + finalname;
                system(sprintf('mv "%s" "%s"', outfile, finalout));
            end
        end

        % only requires the 21 design points 
        function writeLampSCAD(obj, pt, filename)
            % write out the lamp.scad file specifying the parametrized design
            fileh = fopen(filename, 'w');
            fprintf(fileh, 'use <commonLampModules.scad>\n\n');

            fprintf(fileh, '// ====== set the parameter values ======\n');
            fprintf(fileh, '$fn = 10;\n\n');  % sets the mesh resolution. should stay below 100 for performance. even 20 takes a long time.

            fprintf(fileh, 'radii = [%0.6f, %0.6f, %0.6f, %0.6f];\n', obj.radii);
            fprintf(fileh, 'stand_h = 0.25;\n');
            fprintf(fileh, 'stand_r = 1;\n');
            fprintf(fileh, 'shade_end_radius = 0.4;\n');
            fprintf(fileh, 'shade_height = 0.8;\n');
            fprintf(fileh, 'shade_thickness=0.1;\n\n');

            fprintf(fileh, '// vectors defined wrt the ancestor endpoint\n');
            fprintf(fileh, 'base = [%0.6f, %0.6f, %0.6f];\n', pt(1:3));
            fprintf(fileh, 'a1 = [%0.6f, %0.6f, %0.6f];\n', pt(4:6));
            fprintf(fileh, 'a2 = [%0.6f, %0.6f, %0.6f];\n', pt(7:9));
            fprintf(fileh, 'a3 = [%0.6f, %0.6f, %0.6f];\n', pt(10:12));
            fprintf(fileh, 'h1 = [%0.6f, %0.6f, %0.6f];\n', pt(13:15));
            fprintf(fileh, 'h2 = [%0.6f, %0.6f, %0.6f];\n', pt(16:18));
            fprintf(fileh, 'h3 = [%0.6f, %0.6f, %0.6f];\n\n', pt(19:21));

            fprintf(fileh, 'createLampGeometry(base, a1, a2, a3, h1, h2, h3,\n');
            fprintf(fileh, '\t\t\tradii, stand_h, stand_r,\n');
            fprintf(fileh, '\t\t\tshade_end_radius, shade_height, shade_thickness\n'); 
            fprintf(fileh, '\t\t\t);\n\n');

%             fprintf(fileh, 'echo("done creating geometry");\n');
            fclose(fileh);
        end
        
        % only requires the 21 design points 
        function outtext = writeLampJSCAD(obj, pt, jscad_text)
            % alter the template jscad text with the parametrized design
            % values; don't need to write a file, just return the string
            resolution = 10;  % sets the mesh resolution. should stay below 100 for performance. even 20 takes a long time.
            
            q = "";
            
            q = sprintf("%s\tg_fn = %d;\n\n", q, resolution);  % sets the mesh resolution. should stay below 100 for performance. even 20 takes a long time.

            data_format = "0.6f";
            q = sprintf("%s\tvar dp = %s;\n", q, sprintArray(pt(1:obj.rD), data_format));
            q = sprintf("%s\tvar da = %s;\n", q, sprintArray(pt(obj.rD+1 : obj.rD+obj.rA), data_format));
            
            % insert this string into the template jscad file
            outtext = strrep(jscad_text, "//INSERT_PARAMETER_VALUES", q);
        end
        
    end
end


