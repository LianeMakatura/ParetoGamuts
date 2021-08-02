function [design_pts, perf_pts, stats] = ground_truth(funcID, imageTitle, zmin, zmax, step, constrainted_z)
    % construct multiple fronts
    tic
%     zlim_plot = [zmin, zmax];
    stats = struct();
    design_pts = [];
    perf_pts = [];
    hypervolume = cell(1, length(zmin:step:zmax));
    n_evals = 0;
    n_scalar_evals = 0;
    n_pareto_font_points = {};
    j = 1;
    for z=zmin:step:zmax
        disp(strcat('Creating the front for z=', num2str(z)))
        [buffer_z, fig_obj_z, ~, hv_fixed_context, n_evals_fixed_context, n_scalar_evals_fixed_context] = single_z_front(funcID, z, constrainted_z);
%         save(strcat('Lamp_Schulz_50_15run_', num2str(fix((z-zmin)/step)), '.mat'));
        % % extract (pareto optimal) design and performance points from buffer
        idx = buffer_z.getParetoInd();
        pointsD = [];
        pointsF = [];
        for i=1:size(idx, 1)
            if ~buffer_z.emptyCells(idx(i))
                pointsD = [pointsD; buffer_z.Buff(idx(i)).minD];
                pointsF = [pointsF; [buffer_z.Buff(idx(i)).minF, z]];
            end
        end
        design_pts = [design_pts; pointsD]; % each row of _pts is (x,y,z)
        perf_pts = [perf_pts; pointsF];
        hypervolume{j} = hv_fixed_context;
        n_evals = n_evals + n_evals_fixed_context;
        n_scalar_evals = n_scalar_evals + n_scalar_evals_fixed_context;
        n_pareto_font_points{j} = size(idx, 1);
        j = j + 1;
    end
    timeElaspes = toc;
    stats.timeElaspes = timeElaspes;
    stats.n_evals = n_evals + n_scalar_evals;
    stats.n_pareto_font_points = n_pareto_font_points;
    stats.hypervolume = hypervolume;
    
%     zlim(zlim_plot);
%     hold off;
   
%     % plot sample points from the buffer (these should match the pareto front)
%     figure
%     subplot(1, 2, 1)
%     scatter3(design_pts(:, 1), design_pts(:, 2), design_pts(:,3));
%     title('Sampled Points: Design Space')
%     subplot(1,2,2)
%     scatter3(perf_pts(:,1), perf_pts(:,2), perf_pts(:,3));
%     title('Sampled Points: Performance Space')
%     
end