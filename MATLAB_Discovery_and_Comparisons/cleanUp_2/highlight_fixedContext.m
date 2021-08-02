% [buffer, des_figh, perf_figh, pts, labels, filledPtsArray, idxToVis] = appDiscover(31);
% 
% highlightFixedContext(0.04, buffer, des_figh, perf_figh, filledPtsArray, idxToVis);
% highlightFixedContext(0.5, buffer, des_figh, perf_figh, filledPtsArray, idxToVis);
% highlightFixedContext(0.99, buffer, des_figh, perf_figh, filledPtsArray, idxToVis);

% highlightFixedContext(0.04, t_buffer, t_des_figh, t_perf_figh, t_filledPtsArray, t_idxToVis);
% highlightFixedContext(0.5, t_buffer, t_des_figh, t_perf_figh, t_filledPtsArray, t_idxToVis);
% highlightFixedContext(0.99, t_buffer, t_des_figh, t_perf_figh, t_filledPtsArray, t_idxToVis);
% pause;
% 
% highlightFixedContext(0.02, buffer_r, des_figh_r, perf_figh_r, filledPtsArray_r, idxToVis_r);
% highlightFixedContext(0.47, buffer_r, des_figh_r, perf_figh_r, filledPtsArray_r, idxToVis_r);
% highlightFixedContext(0.97, buffer_r, des_figh_r, perf_figh_r, filledPtsArray_r, idxToVis_r);
%
highlightFixedContext(0.02, buffer_s, des_figh_s, perf_figh_s, filledPtsArray_s, idxToVis_s);
highlightFixedContext(0.47, buffer_s, des_figh_s, perf_figh_s, filledPtsArray_s, idxToVis_s);
highlightFixedContext(0.97, buffer_s, des_figh_s, perf_figh_s, filledPtsArray_s, idxToVis_s);


function [FCpoints] = highlightFixedContext(contextValue, buffer, desFigH, perfFigH, filledPtsArray, idxToVis)
    % ===== pick indices to visualize
    d1=idxToVis.designIndex1;
    d2=idxToVis.designIndex2; 
    a1=idxToVis.appIndex;
    p1=idxToVis.perfIndex1;
    p2=idxToVis.perfIndex2;
    
    des1_idx = d1;
    des2_idx = d2;
    app_idx = buffer.rD + a1;
    perf1_idx = buffer.rD+buffer.rA+p1;
    perf2_idx = buffer.rD+buffer.rA+p2;
    
    markerSize = 10;
    
    
    % find the subscripted context cell index
    % figure out the appropriate contextual cell for this point
    contextCell = buffer.getLinearCellIdx(contextValue, buffer.cylCellDims(buffer.rd-1+1:end));

    % get Pareto optimal points out of this set
    appDims = buffer.cylCellDims{buffer.rd-1+1:buffer.rd+buffer.rA};
    cellSize = 1./appDims;
    
    FCpoints = [];
    for patchID = 1:length(filledPtsArray)
        allPts = filledPtsArray{patchID};
        if isempty(allPts)
            continue;
        end
        
        % get center of context, and find all points within cell
        if length(appDims) == 1
            center = (contextCell - 0.5) * cellSize;
            % find all points near that context
            FCidx = find( (abs(allPts(:, buffer.rD+1:buffer.rD+buffer.rA) - center) - cellSize*0.8 ) < 0 );
            FCpoints = [FCpoints; allPts(FCidx, :)];
        else
           [cellIdx{:}] = ind2sub(appDims, contextCell{:});
           center = (cellIdx - 0.5) * cellSize;

           for pInd=1:size(allPts, 1)
               keep = max(abs(allPts(pInd, buffer.rD+1:buffer.rD+buffer.rA) - center) - cellSize*0.5 ) < 0;
               if keep
                   FCpoints = [FCpoints; allPts(pInd, :)];
               end
           end
        end
    end
    
    if ~isempty(FCpoints)    
        figure(desFigH);
        hold on;
        plot3(FCpoints(:, des1_idx), FCpoints(:, des2_idx), FCpoints(:, app_idx), ...
                    '-o',...
                    'Color','k',...
                    'MarkerSize',markerSize,...
                    'MarkerFaceColor','k');
        hold off;
        
        figure(perfFigH);
        hold on;
        plot3(FCpoints(:, perf1_idx), FCpoints(:, perf2_idx), FCpoints(:, app_idx), ...
                    'o',...
                    'Color','k',...
                    'MarkerSize',markerSize,...
                    'MarkerFaceColor','k');
        hold off;
    end
end
