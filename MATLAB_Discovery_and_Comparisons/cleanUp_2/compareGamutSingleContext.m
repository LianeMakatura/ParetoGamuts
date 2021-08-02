% Compare Pareto Gamut to single-context ground truth solutions

% appFunc = 1; scFunc = 15;   % ZDT1
% appFunc = 2; scFunc = 16;     % ZDT2
% appFunc = 3; scFunc = 17;     % ZDT3
% appFunc = 4; scFunc = 18;     % ZDT4
% appFunc = 10; scFunc = 21;      % Fourier
zmin = 0.1;
% zmax = 0.5;
zmax = 1.0;
vis = false;
num_run = 1;
appFunc = 10;
% appFunc = 10; scFunc = 32; % bicopter
% appFunc = 10; scFunc = 31;
% appFunc = 10; scFunc = 33; % GridShell
% zmin = 0.0025; 
% zmax = 0.5025;
% scFunc_list = [33, 34, 35];
scFunc_list = [33];
for scFunc_i = 1
    stats_idx = {};
    for n_ = 1 : num_run
        stats_single_run = groundTruthComparison(appFunc, scFunc_list(scFunc_i), zmin, zmax, vis);
        stats_idx{n_} = stats_single_run;
    end
    save(strcat('stats_Schulz_', num2str(scFunc_i), '.mat'), 'stats_idx');
end


%%
function stats = groundTruthComparison(appFuncID, singleContextFuncID, zmin, zmax, vis)
    [gt_des, gt_perf, stats] = ground_truth(singleContextFuncID, 'title', zmin, zmax, 0.2,  singleContextFuncID ~= 29 && singleContextFuncID ~= 30&&singleContextFuncID ~= 31&& singleContextFuncID ~= 32&& singleContextFuncID ~= 33&& singleContextFuncID ~= 34&& singleContextFuncID ~= 35);
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


