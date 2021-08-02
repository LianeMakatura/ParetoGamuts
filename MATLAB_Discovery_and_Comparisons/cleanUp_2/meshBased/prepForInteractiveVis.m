function [] = prepForInteractiveVis(mFunc, pts, labels, mesh_provider, mesh_radius)
%PREPFORINTERACTIVEVIS Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('mesh_radius', 'var') && strcmp(mesh_provider, "jscad")
        sprintf("Error: need to provide a mesh radius if using jscad to generate meshes");
        return;
    end
        
    rD = mFunc.rD; 
    rA = mFunc.rA;
    rd = mFunc.rd;
    
    problemName = mFunc.p.problemName;
    visFolderName = mFunc.p.visName; %in most cases, same as problem name.
    problem_dir = mFunc.p.wd; 
    
    vis_topdir = "<location of ParetoGamutJSVisualizer> including final /"; % where the JS vis files is located
    JSONdestination_dir = vis_topdir + "vis/data/" + visFolderName + "/"; %where to write the resulting JSONs
    if ~exist(JSONdestination_dir, 'dir')
        system(sprintf('mkdir %s', JSONdestination_dir));
    end

        
    % =================== JSCAD template ========================
    % need a clean way to name these files, also named mesh template in JS
    % library right now
    if strcmp(mesh_provider, "jscad")
        curr_template = problemName + "_template.jscad";
        jscad_filepath = problem_dir + curr_template;
        system(sprintf('cp %s %s', jscad_filepath, JSONdestination_dir));
        system(sprintf('mv %s %s', JSONdestination_dir + curr_template, JSONdestination_dir + "mesh_template.jscad"));
    end
    
    % =================== Problem description ===================
    outfile = JSONdestination_dir + "problem.json";
%     writeProblemFile(mFunc, outfile, mesh_provider, mesh_radius);
    
    
    % =================== Front description ===================
    j = jsonencode([pts, labels]);
    fh = fopen(JSONdestination_dir + "fronts.json", 'w');
    fprintf(fh, j);
    fclose(fh);
    
    
    % =================== Mesh description ===================
    if strcmp(mesh_provider, "files")
        % write out meshes corresponding to pts
        mesh_dir = problem_dir + "meshes/";
        if ~exist(mesh_dir, 'dir')
            system(sprintf('mkdir %s', mesh_dir));
        end

        mFunc.writeMeshes(pts(:, 1:rD+rA));

        % convert the meshes to a JSON
        scriptname = vis_topdir + "scripts/convert_meshes.py";

        system(sprintf('python %s %s %s', ...
                        scriptname, mesh_dir, JSONdestination_dir));
    end

end


% Writes a JSON for the javascript visualization describing the problem
% specification
function writeProblemFile(mFunc, outfile, mesh_provider, mesh_radius)
    if ~exist('mesh_provider', 'var')
        mesh_provider = "jscad"; 
%         mesh_provider = "files";
    end
    if ~exist('mesh_radius', 'var')
        mesh_radius = 10;
    end

    % minmax vals set to [0,1] for now; change for each var if necessary
    % later
    minval = 0; 
    maxval = 1; 

    rD = mFunc.rD; 
    rA = mFunc.rA; 
    rd = mFunc.rd;
    names = mFunc.varNames;
    
    entryfmt = '\t\t{"name": "%s", "min": %0.6f, "max": %0.6f},\n';
    lastentryfmt = '\t\t{"name": "%s", "min": %0.6f, "max": %0.6f}\n';

    fileh = fopen(outfile, 'w');
    fprintf(fileh, '{\n');
    
    % ====== write out the design variable information
    fprintf(fileh, '\t"design_variables": [\n');
    for i=1:rD-1
        fprintf(fileh, entryfmt, ...
                        names(i), minval, maxval);
    end
    fprintf(fileh, lastentryfmt, names(rD), minval, maxval);
    fprintf(fileh, '\t],\n');
    
    % ====== write out the application variable information
    fprintf(fileh, '\t"application_variables": [\n');
    for i=rD+1:rD+rA-1
        fprintf(fileh, entryfmt, ...
                        names(i), minval, maxval);
    end
    fprintf(fileh, lastentryfmt, names(rD+rA), minval, maxval);
    fprintf(fileh, '\t],\n');
    
    % ====== write out the performance variable information
    fprintf(fileh, '\t"performance_metrics": [\n');
    for i=rD+rA+1:rD+rA+rd-1
        fprintf(fileh, entryfmt, ...
                        names(i), minval, maxval);
    end
    fprintf(fileh, lastentryfmt, names(rD+rA+rd), minval, maxval);
    fprintf(fileh, '\t],\n');
        
    % ====== additional specifications
    fprintf(fileh, '\t"mesh_provider": "%s",\n', mesh_provider);
    fprintf(fileh, '\t"max_mesh_radius": %f,\n', mesh_radius);
    fprintf(fileh, '\t"design_chart_options": [{"type": "chart2d"}, {"type": "radar"}],\n');
    fprintf(fileh, '\t"performance_chart_options": [{"type": "chart2d"}]\n');
        
    fprintf(fileh, '}\n');

    fclose(fileh);
end

