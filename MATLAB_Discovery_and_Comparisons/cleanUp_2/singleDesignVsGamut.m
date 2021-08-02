
% @param buffer
% @param gamutPts
% @param designPt - NxrD point in design space

% assumes a single context variable
function singleDesignVsGamut(buffer, gamutPts, gamutLabels, numContexts, numSamples, designPt)
    % get a particular point, map it into all contexts
    if numSamples == size(gamutPts, 1)
        % visualize all samples
        disp(numSamples);
        designPt = gamutPts(:, 1:buffer.rD);
    elseif ~exist('designPt', 'var')
        idx = ceil(rand(numSamples, 1) * size(gamutPts, 1));
        designPt = gamutPts(idx, 1:buffer.rD);
    end
    
    contextPts = linspace(0, 1, numContexts)';
    
    samples = cell(numSamples, 1);
    for i=1:numSamples
        ptDesMultiC = repmat(designPt(i, :), [numContexts, 1]);
        ptPerfMultiC = buffer.mFunc.eval(ptDesMultiC, contextPts);
        
        samples{i} = struct('des', ptDesMultiC, ...
                            'perf', ptPerfMultiC);
    end
    
    % find the closest point on each Pareto front
    sac_curves = cell(numSamples, 1);
    sac_closestAugPerfPts = cell(numSamples, 1);
    for i=1:numSamples
        s = samples{i};
        [sac_curves{i}, sac_closestAugPerfPts{i}] = findSacrificeCurves(buffer, gamutPts, s.perf, contextPts);
        percentOptimalCurve{i} = nadir_error_curve(s.perf, sac_closestAugPerfPts{i});
    end
    
    writeDesignExpToJSON(samples, sac_curves, sac_closestAugPerfPts);
    
    % plot the error curve to the closest point
    figure; 
    title('Sacrifice Curves (Dist from eval to closest pt on gamut)')
    xlabel('Context Value');
    ylabel('Distance to closest point on gamut');
    hold on;
    for i=1:numSamples
        s = sac_curves{i};
        plot(contextPts, s, 'LineWidth', 3);
    end
    set(gca, 'FontSize', 18)
    hold off;
    
    % plot the percent optimal curve
    figure; 
    title('% Optimal Performance Achieved')
    xlabel('Context Value');
    ylabel('% Optimal Performance');
    hold on;
    l_entries = cell(numSamples, 1);
    for i=1:numSamples
        s = percentOptimalCurve{i};
        plot(contextPts, s, 'LineWidth', 3);
        l = sprintf("r=%0.2f, h=%0.2f, a=%0.2f", samples{i}.des(1, :));
        l_entries{i} = l;
    end
    set(gca, 'FontSize', 18);
    ylim([0, 1.05]);    
    legend(l_entries);
    hold off;
    
    
    % plot the gamut with the single design curve over it
    idxToVis = struct('designIndex1', 1, ...
                  'designIndex2', 2, ...
                  'appIndex', 1, ...
                  'perfIndex1', 1, ...
                  'perfIndex2', 2);
    [des_figh, perf_figh] = buffer.visualizeAppBuff(idxToVis, false, gamutPts, gamutLabels);
    
    % plot the single design curve, and the points used for its measurement
    figure(perf_figh);
    hold on;
    for i=1:numSamples
        s = samples{i};
        plot3(s.perf(:, idxToVis.perfIndex1), ...
             s.perf(:, idxToVis.perfIndex2), ...
             contextPts(:, idxToVis.appIndex), ...
             '-o', ...
             'LineWidth', 3);
    end
    set(gca, 'FontSize', 18)
    title('Pareto Gamut')
    hold off;
    
    figure(des_figh);
    hold on;
    for i=1:numSamples
        s = samples{i};
        plot3(s.des(:, idxToVis.designIndex1), ...
             s.des(:, idxToVis.designIndex2), ...
             contextPts(:, idxToVis.appIndex), ...
             '-o', ...
             'LineWidth', 3);
    end
    set(gca, 'FontSize', 18)
    hold off;
    
end



% only pass in the performance points, no context
function percentOptimal = nadir_error_curve(fx, fx_opt)
    rd = size(fx_opt, 2);
    nadirpt = ones(1, rd);
        
    percentOptimal = zeros(length(fx_opt), 1);
    for i=1:length(fx_opt)
        xoptVol = prod(nadirpt - fx_opt(i, :));
        xVol = prod(nadirpt - fx(i, :));
        percentOptimal(i) = min(xVol / xoptVol, 1);
    end
end





function writeDesignExpToJSON(samples, sac_curves, sac_closestAugPerfPts)
    numsamples = length(samples);
    numContexts = size(samples{1}.perf, 1); % all have same # contexts
    
    % for each sample (single design point)
    jdata = cell(numsamples, 1);
    for i=1:numsamples
%         % construct the perf evals as nested cell array
%         perfEvals = cell(numContexts, 1);
%         for j=1:numContexts
%             perfEvals{j} = samples{i}.perf(j, :); %eval at jth context val
%         end

%         js = struct('designPt', samples{i}.des(1, :), ...
%                     'perfEvalsPerContext', samples{i}.perf,...
%                     'errorPerContext', sac_curves{i} ...
%                 );
        js = struct('designPt', samples{i}.des(1, :), ...
            'perfEvalsPerContext', samples{i}.perf,...
            'closestPtPerContext', sac_closestAugPerfPts{i} ...
        );    
        
        jdata{i} = js;
    end
    
    j = jsonencode(jdata);
%     fh = fopen(JSONdestination_dir + "designExp.json", 'w');
    fh = fopen("designExp.json", 'w');
    fprintf(fh, j);
    fclose(fh);
end




% @param perfMultiC - numContext x rd array mapping fixed design in every context
% @param contextPts - context values over range
function [sacCurve, closestPerfPts] = findSacrificeCurves(buffer, gamutPts, perfMultiC, contextPts)
    numContextpts = size(contextPts, 1);
    sacCurve = zeros(numContextpts, 1);
%     closestPerfPts = zeros(numContextpts, buffer.rd + buffer.rA);
    closestPerfPts = zeros(numContextpts, buffer.rd);

    cellSize = 1 / numContextpts;
    
%     figure;
%     hold on;
    for c=1:numContextpts
        % assume 1 context dimension 
        center = (c - 0.5) * cellSize;
        % find all points near that context
        FCidx = find( (abs(gamutPts(:, buffer.rD+1:buffer.rD+buffer.rA) - center) - cellSize*0.5 ) < 0 );
        FCperfpoints = gamutPts(FCidx, buffer.rD+buffer.rA+1:end);
        
        if isempty(FCperfpoints)
            sacCurve(c) = NaN;
            closestPerfPts(c, :) = NaN([1, buffer.rd+buffer.rA]);
            continue
        end
        
        % find the point closest to our query pt's performance at this
        % context
        queryPerfPt = perfMultiC(c, :);
        
        % define a triangle mesh with the points, then find the minimum
        % distance to that
        if buffer.rd == 2
             [~, I] = sort(FCperfpoints(:, 1)); % sort in ascending order to ensure valid curve; works bc pareto optimality
             FCperfpoints = FCperfpoints(I, :);
             
             [closestPt,dist,~] = distance2curve(FCperfpoints,queryPerfPt,'linear');
%              plot(FCperfpoints(:, 1), FCperfpoints(:, 2));
%              keyboard;
             sacCurve(c) = dist;
%              closestPerfPts(c, :) = [closestPt, c];
             closestPerfPts(c, :) = closestPt;
        else
            %compute Euclidean distances to all points:
            dist2 = sum((FCperfpoints - queryPerfPt).^ 2, 2);
            
            %find the smallest distance and use that as an index into B:
            closestDist = min(dist2);
            closestGamutPt = FCperfpoints(dist2 == closestDist,:);
            sacCurve(c) = closestDist;
%             closestPerfPts(c, :) = [closestGamutPt(1, :), c]; %small fudge factor in app var
            closestPerfPts(c, :) = closestGamutPt(1, :); %small fudge factor in app var

        end        
    end
%     hold off
%     figure
%     title('closest gamut point')
%     plot(closestPerfPts(:, 1), (closestPerfPts(:, 2)))
end

