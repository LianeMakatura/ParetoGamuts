function [] = write_obj(mesh, obj_name, write_mode)
%WRITE_OBJ Summary of this function goes here
%   Detailed explanation goes here
    fileH = fopen(obj_name, write_mode);
    
    fprintf(fileH, '# Simple Wavefront file\n');
    
    numVerts = length(mesh.vertices);
    for v=1:numVerts
        vertex = mesh.getVertex(v);
        fprintf(fileH, 'v %.6f %.6f %.6f\n', ...
                vertex(1), vertex(2), vertex(3));
    end
    
    numFaces = length(mesh.faces);
    for f=1:numFaces
        face = mesh.getFace(f);
        numVerts = length(face);
        f_str = "f";
        for fv=1:numVerts
            f_str = strcat(f_str, " ", num2str(face(fv)));
        end
        f_str = strcat(f_str, "\n");
        fprintf(fileH, f_str);
    end
    
    fclose(fileH);
end

