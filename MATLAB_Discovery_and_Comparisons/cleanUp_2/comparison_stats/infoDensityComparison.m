
%% Turbine
problem_name = 'Turbine';
Schulz_id = 30;
Ours_id = 31;
rd = 2;
max_eval_constraints=1e10;



%% 
maxRuns = 1;
visBool = false;

numFronts = 10;

zmin = 0.0;
zmax = 1.0;
contextStepSize = (zmax-zmin)/(numFronts);

zmin = zmin + 0.5*contextStepSize;
zmax = zmax - 0.5*contextStepSize; % center the samples

z_range = [zmin, zmax, contextStepSize]; %min, max, stepsize
zVals = z_range(1):z_range(3):z_range(2);
if numFronts ~= length(zVals)
    disp('numFronts do not match')
    return
end



% [perf_Schulz, time_Schulz, n_evals_Schulz, hypervolume_Schulz] = comparison_Schulz(Schulz_id, z_range, rd, visBool, max_eval_constraints);
[perf_Ours, time_ours, n_evals_ours, hypervolume_ours] = comparison_ours(Ours_id, z_range, rd, visBool);


figure
hold on;
markerSize = 50;

plottingColors = perf_Ours(:, 3); % color by context value
                        caxis([0,1]);
                        colormap jet;
scatter3(perf_Ours(:,1), perf_Ours(:,2), perf_Ours(:,3), markerSize, plottingColors, 'filled');

markerSize = 80;
scatter3(perf_Schulz(:,1), perf_Schulz(:,2), perf_Schulz(:,3), markerSize, 'k', 'filled');
xlim([0,1]);
ylim([0,1]);
xlabel('Mass');
ylabel('Power');
zlabel('Wind Speed')
hold off;

