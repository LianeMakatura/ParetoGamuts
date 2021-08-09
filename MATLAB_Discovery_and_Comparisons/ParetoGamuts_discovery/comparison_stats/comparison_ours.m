% z_range = [0.0025, 0.5025, 0.1];
z_range = [0.1, 1.0, 0.2];
rd = 2;
problemID = 32;
num_runs = 3;
visualize = false;
time_all = zeros(num_runs, 1);
n_evals_all = zeros(num_runs, 1);
hypervolume_all = {};
for n_ = 1:num_runs
    [perf, time_ours, n_evals_ours, hypervolume_ours] = comparison_ours_single(problemID, z_range, rd, visualize);
    time_all(n_) = time_ours;
    n_evals_all(n_) = n_evals_ours;
    hypervolume_all{n_} = hypervolume_ours;
end

function [perf, time_ours, n_evals_ours, hypervolume_ours] = comparison_ours_single(problemID, z_range, rd, visualize)
    [time_ours, n_evals_ours, buffer, ~, ~, pts, ~, ~, ~] = appDiscover(problemID, false);
    perf = [pts(:,buffer.rD+buffer.rA+1:buffer.rD+buffer.rA+rd), pts(:,buffer.rD+buffer.rA)];
    hypervolume_ours = [];
    nadirPt = ones(rd, 1);
    for z = z_range(1):z_range(3):z_range(2)
        idx = find(abs(perf(:, end) - z) < 1/200);
        pointF = [];
        for i = 1:rd
            pointF = [pointF, perf(idx, i)];
        end
        hypervolume_ours = [hypervolume_ours, lebesgue_measure(pointF', nadirPt)];
    end
    if visualize
        figure
        scatter3(perf(:,1), perf(:,2), perf(:,3), 'r');
        xlim([0,1]);
        ylim([0,1]);
    end
    title('Ours: Performance Space')
end