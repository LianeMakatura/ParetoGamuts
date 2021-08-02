function [] = lowerEnvelopeCompare(buffer, gamutPts, gt_xlim, gt_ylim)
    % get the ground truth (global only)
    [pts_global, pts_gamut] = buffer.mFunc.getGroundTruth();
    
    
    
    % compute the lower envelope of our pareto gamut
    % ignore z axes, and get Pareto indices 
    flattenedGamut = gamutPts(:, buffer.rD+buffer.rA+1:buffer.rD+buffer.rA+buffer.rd); % only extract performance points, no app vars
    nonDominated = getParetoIndices(flattenedGamut); %gives logical array
    paretoIndices = find(nonDominated); % keep only the non-zero (non-dominated) indices
    fLowerEnv = flattenedGamut(paretoIndices, :);
    
    %% compute the error between them 
    numFrontPieces = length(pts_global);
    C = [];
    
    queryPts = fLowerEnv; % m x rd
    for i=1:numFrontPieces
        piece = pts_global{i};
        
        if buffer.rd == 2
            pts = piece; % n x 2
            [closestPt,dist,~] = distance2curve(pts,queryPts,'linear');
            C = [C, dist];
            
        elseif fLowerEnv.rd == 3
            s = surf(piece{1},piece{2},piece{3});
            C = [C, colorByDistance(s, queryPts)];
        end
    end

    % compute the closest distances over all surface pieces (each column result of different surface patch)
    if size(C, 2) > 1
       C =  min(C,[],2);
    end

    %% ====== Plot the result
    figure;
    h(1) = subplot(2,1,1);
    hold on;
    numFrontPieces = length(pts_global);
    for i=1:numFrontPieces
        piece = pts_global{i};
        if buffer.rd == 2
            plot(piece(:, 1),piece(:, 2), 'black','LineWidth',4);
        end
    end
    xlim(gt_xlim)
    ylim(gt_ylim)
    xlabel('f1');
    ylabel('f2');
    set(gca, 'FontSize', 18);
    
    % color the scatter plot by found points.
    markerSize = 40;
    h(2) = subplot(2,1,2);

    if buffer.rd == 2
        hold on
        for i=1:numFrontPieces
            piece = pts_global{i};
            if buffer.rd == 2
                plot(piece(:, 1),piece(:, 2), 'black','LineWidth',1);
            end
        end
        
        
        scatter(queryPts(:, 1), queryPts(:, 2), markerSize, C, 'filled');
        colormap(jet)
        caxis([0 max(C)]);
        xlabel('f1');
        ylabel('f2');
        set(gca, 'FontSize', 18);
        xlim(gt_xlim);
        ylim(gt_ylim);
        cbh = colorbar(h(2));
        hold off
        
%         % Reposition to figure's left edge, centered vertically
%         cbh.Position(1) = .95-cbh.Position(3);
%         cbh.Position(2) = 0.5-cbh.Position(4)/2;
        % decrease horizontal extent of subplots to 92% of their current width
        set(h(1), {'Position'}, mat2cell(vertcat(h(1).Position) .* [1 1 .865, 1], ones(size(h(1))),4))
        
    elseif buffer.rd == 3
        scatter3(queryPts(:, 1), queryPts(:, 2), queryPts(:, 3), markerSize, C, 'filled');
        colormap(jet)
        colorbar
        caxis([0 max(C)]);
        xlabel('f1');
        ylabel('f2');
        zlabel('f3');
        set(gca, 'FontSize', 18);
    end
    hold off;
    
end