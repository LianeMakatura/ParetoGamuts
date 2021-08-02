function gt_gamut = ZDT3_gamutGroundTruth()
    rD = 2; 
    rA = 1;

    % get the full set of candidate optimal points
    f1pts = linspace(0, 1, 500);
    numZpts = 200;
    gt_z = linspace(0,1,numZpts);
    gamut_map = gt_f2_gamut(f1pts, gt_z, rD, rA);
%     surf(gamut_map{1}, gamut_map{2}, gamut_map{3})

    % weed out the dominated points on this patch at each context
    endPatchIdx = @(X) (find(X==0, 1) - 1); % one before first nan index in context
    beginPatchIdx = @(X) (find(X==1, 1)); %first non-nan index in context
    vertices = cell(5, 1); % know there are 5 patches
    for c=1:numZpts
        F = [gamut_map{1}(c, :)', gamut_map{2}(c, :)']; % get all points corresponding to this context

        % must check for dominated points in each context
        % pass performance points in, to get non-dominate points
        paretoIndices = logical(getParetoIndices(F)); %gives logical array        
        for i=1:5
            % get starting index and remove all non-dominated points until then
            patchStartIdx = beginPatchIdx(paretoIndices); 
            paretoIndices = paretoIndices(patchStartIdx:end);
            F = F(patchStartIdx:end, :);
            
            patchEndIdx = endPatchIdx(paretoIndices); % get finish index
            
            % update the gamut & remove used points
            vertices{i} = [vertices{i}; F(1:patchEndIdx, :), repmat(gt_z(c), patchEndIdx, 1)];
            paretoIndices = paretoIndices(patchEndIdx+1:end);
            F = F(patchEndIdx+1:end, :);
        end
    end
    
%     figure;hold on;
    gt_gamut = cell(5, 1);
    for i=1:5
        T = delaunay(vertices{i}(:, 1), vertices{i}(:, 2));
        gt_gamut{i} = struct('Vertices', vertices{i}, ...
                             'Faces', T, ...
                             'triStruct', true...
                         );
%         scatter3(gt_gamut{i}(:, 1), gt_gamut{i}(:, 2), gt_gamut{i}(:, 3));
    end
%     hold off;
%     
    figure; hold on;
    for i=1:5
        t = trisurf(gt_gamut{i}.Faces, gt_gamut{i}.Vertices(:, 1), gt_gamut{i}.Vertices(:, 2), gt_gamut{i}.Vertices(:, 3));
        t.EdgeColor = 'none';
    end
    hold off;
end



function gt_surf_pts = gt_f2_gamut(f1pts, gt_z, rD, rA)
    [F1, Z] = meshgrid(f1pts, gt_z);

    % --  only visualizable with 1
    % app var, no sum over app vars needed in g
    g = (1 + 9 / (rD + rA - 1) .* Z);
    F2 = g .* (1 - sqrt(F1./g) - (F1 ./ g).*sin(10.*pi.*F1) );
    F2 = F2 ./ 11 + 0.5;
    gt_surf_pts = {F1, F2, Z};
end