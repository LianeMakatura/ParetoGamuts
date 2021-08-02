% Compare Pareto Gamut to single-context ground truth solutions

% appFunc = 1; scFunc = 15;   % ZDT1
% appFunc = 2; scFunc = 16;     % ZDT2
% appFunc = 3; scFunc = 17;     % ZDT3
% appFunc = 4; scFunc = 18;     % ZDT4 - doesn't work yet, and Schulz takes forever, ask about the settings they needed
% appFunc = 10; scFunc = 21;      % Fourier
% zmin = 0.0;
% zmax = 1.0;
appFunc = 10; scFunc = 32; % bicopter
zmin = 0.0025; 
zmax = 0.5025;
vis = false;
time_all = zeros(10, 1);
n_evals_all = zeros(10, 1);
hypervolume_all = {};
for n_ = 1 : 10
    [timeElaspes, n_evals, hypervolume] = groundTruthComparison(appFunc, scFunc, zmin, zmax, vis);
    time_all(n_) = timeElaspes;
    n_evals_all(n_) = n_evals;
    hypervolume_all{n_} = hypervolume;
end


%%
function [timeElaspes, n_evals, hypervolume] = groundTruthComparison(appFuncID, singleContextFuncID, zmin, zmax, vis)
    [gt_des, gt_perf, timeElaspes, n_evals, hypervolume] = ground_truth(singleContextFuncID, 'title', zmin, zmax, 0.1,  singleContextFuncID ~= 29 && singleContextFuncID ~= 30&& singleContextFuncID ~= 32);
    if vis
        figure(); 
        set(gca,'FontSize', 14);
        scatter3(gt_perf(:, 1), gt_perf(:, 2), gt_perf(:, 3));
    %     scatter3(zeros(50, 1), zeros(50, 1), linspace(0, 1, 50)', 80, 'k', 'filled'); % reference line for optimization

        figure(); 
        set(gca,'FontSize', 14);
        scatter3(gt_des(:, 1), gt_des(:, 2), gt_des(:, 3));
    end
    
end





%% Visualize the single context functions
% visualizeSingleContextFuncs('eval_ZDT1_f2');

% %% Function visualization
% % pass in the meshgrid points for the design points x1, x2
% % and grid of points for the context value (a repmat)
% 
% function f2 = eval_ZDT1_f2(X1,X2,Z)
%     s = X1 + X2 + Z;
%     g = 1 + 4.5*s;
% 
%     f2 = g * (1 - sqrt(X1 ./ g));
%     f2 = f2 / 10;       % normalize so all possible outputs are between 0 and 1
% end
% 
% 
% function visualizeSingleContextFuncs(f)
%     x1 = 0:0.01:1;
%     x2 = 0:0.01:1;
%     [X, Y] = meshgrid(x1, x2);
%     z = 0:0.1:1;
% 
%     figure; hold on;
%     for i=1:length(z)
%         zstar = z(i);
%         Z = repmat(zstar, size(X));
%         evald_f = feval(f, X,Y,Z);
%         
%         surf(X,Y,evald_f);
%         xlabel('x1');
%         ylabel('x2');
%         zlabel('f2');
%         drawnow;
%         keyboard;
%     end
%     hold off;
% end


