function C = colorByDistance(surf, queryPts, triStruct) 
    if ~exist('triStruct', 'var')
        triStruct = false;
    end
    
    if triStruct
        f = surf.Faces;
        v = surf.Vertices;
    else
        % test surface to patch 
        [f,v,~] = surf2patch(surf, 'triangles');
    end
    
    %create point2trimesh input struct
    args.Faces = f;
    args.Vertices = v;
    args.QueryPoints = queryPts;
    
    [distances, ~] = point2trimesh( args );
    
    distances = abs(distances);
    C = distances;
end