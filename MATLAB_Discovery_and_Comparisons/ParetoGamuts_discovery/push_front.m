function [appxDesignPts, appxPerfPts] = push_front(funcID, buffer_z, z, target_z)
%PUSH_FRONT Summary of this function goes here
%   Detailed explanation goes here

    % % extract (pareto optimal) design and performance points from buffer
    idx = buffer_z.getParetoInd();
    design_pts = [];
    perf_pts = [];
    for i=1:size(idx, 1)
        if ~buffer_z.emptyCells(idx(i))
            if funcID == 19
                design_pts = [design_pts; [buffer_z.Buff(idx(i)).minD, z]];
                perf_pts = [perf_pts; [buffer_z.Buff(idx(i)).minF, z]];
            elseif funcID == 21
                design_pts = [design_pts; buffer_z.Buff(idx(i)).minD];
                perf_pts = [perf_pts; [buffer_z.Buff(idx(i)).minF, z]];
            end
        end
    end

    %% sample from the buffer points (just keep them all for now)

    %% create new (rD+1)-dimensional design space & fixed params (needed for exploration)
    mFunc = MappingFunction(funcID+1);  % the unconstrained z version of whatever sample I picked

    imageTitle = mFunc.getImageTitle();
    [rD,rd] = mFunc.getDimensions();
    numCells = mFunc.getNumCells(); 
    [minBound, maxBound] = mFunc.getBounds();
    fixedParams = struct('rD', rD, 'rd', rd, 'imageTitle', imageTitle, 'minBound', minBound, ...
                         'maxBound', maxBound);

    mFunc_targ_z = MappingFunction(funcID, target_z); % for easier evaluation

    %% project to desired z value -- samples *must* be between 0 and 1 for this to work
    appxDesignPts = [];
    appxPerfPts = [];
    num_samples = size(design_pts,  1);
    skipped = 0;
    for i=1:num_samples
        centerPoint = design_pts(i, :);
        explorationDirection = getExplorationDirection(centerPoint(1:rD),mFunc, fixedParams);
        if size(explorationDirection, 2) < 2    % not enough directions
            disp('skipped pt')
            skipped = skipped + 1;
            continue;
        end

        % % ==== NOTE: only valid for 2+1 D right now =====
        % find intersection between exploration manifold and target z plane
        n_targz = [0,0,1];
        n_expDirs = cross(explorationDirection(:, 1)', explorationDirection(:,2)');
        d = cross(n_targz, n_expDirs); % direction of intersection line

        % find shortest path along manifold to get to z=c; follow until z=c
        % TODO: fix for general case
        new_des_pt = plane_intersection(n_expDirs, centerPoint, n_targz, [0,0,target_z], centerPoint);
        new_des_pt = new_des_pt';
        %new_des_pt = [centerPoint(1:2), z];

        % % === End dimension specific code ====

        % save point on appx Pareto set
        appxDesignPts = [appxDesignPts; new_des_pt];

        % map to perf space and save to appx Pareto front
        if funcID == 19
            new_perf_pt = mFunc_targ_z.eval(new_des_pt(1:rD-1), new_des_pt(rD));   %using unconstrainted z
        elseif funcID == 21
            new_perf_pt = mFunc_targ_z.eval(new_des_pt, new_des_pt(rD));
        end
        appxPerfPts = [appxPerfPts; new_perf_pt];
    end

end

