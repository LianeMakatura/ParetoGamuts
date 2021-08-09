function structural_results = max_stress(meshstruct, fixed_faces, load_faces, load)
    structuralmodel = createpde('structural','static-solid');
    
    % 
    nodes = meshstruct.vertices';
    elements = meshstruct.faces';
    geometryFromMesh(structuralmodel,nodes,elements);

    pdegplot(structuralmodel,'FaceLabels','on','CellLabels','on','FaceAlpha',0.5)

    %Specify the Young's modulus and Poisson's ratio for each metal.
    structuralProperties(structuralmodel,'Cell',1,'YoungsModulus',110E9, ...
                                                  'PoissonsRatio',0.28);
                                              
    %Specify that faces 1 and 4 are fixed boundaries.
    structuralBC(structuralmodel,'Face',fixed_faces,'Constraint','fixed');
    
    %Specify the surface traction for faces 2 and 5.
    structuralBoundaryLoad(structuralmodel,'Face',load_faces,'SurfaceTraction',load);
    
    %Generate a mesh and solve the problem.
    tic;
    
    generateMesh(structuralmodel);
    toc;
    %solver finds the values of displacement, stress, strain, and von Mises stress 
    %at the nodal locations. To access these values, use structuralresults.Displacement, 
    %structuralresults.Stress, and so on. The displacement, stress, and strain values at
    %the nodal locations are returned as structure arrays with fields representing their components.
    structural_results = solve(structuralmodel);

    %Plot the deformed shape with the z-component of normal stress.
    figure;
    pdeplot3D(structuralmodel,'ColorMapData',structuralresults.VonMisesStress, ...
                              'Deformation',structuralresults.Displacement)
end
