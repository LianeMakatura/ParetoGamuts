
function plotNEvalsComparison(maxRuns, frontsPerTest, problemName)
    addpath('raacampbell-shadedErrorBar-263537a/')

    SchulzCol = 1;
    NSGACol = 2;
    OurCol = 3;

    numTests = length(frontsPerTest);
    avg_evals_NSGA = zeros(1, numTests);
    avg_evals_Schulz = zeros(1, numTests);
    avg_evals_Ours = zeros(1, numTests);
    
    evals_NSGA = zeros(maxRuns, numTests);
    evals_Schulz = zeros(maxRuns, numTests);
    evals_Ours = zeros(maxRuns, numTests);
    
    for i=1:length(frontsPerTest)
        numFCFronts = frontsPerTest(i);
        filename = sprintf('./%s_numFronts_%d.mat', problemName, numFCFronts);
        stats = load(filename);
        
        avg_evals_NSGA(1, i) = stats.avg_evals(NSGACol);
        evals_NSGA(:, i) = stats.all_nevals(:, NSGACol);

        avg_evals_Schulz(1, i) = stats.avg_evals(SchulzCol);
        evals_Schulz(:, i) = stats.all_nevals(:, SchulzCol);

    end
    
    % load our results
    filename = sprintf('./%s_ours.mat', problemName);
    stats = load(filename);
    problemEvalsOurs = stats.avg_evals(OurCol);
    avg_evals_Ours(1, :) = repmat(problemEvalsOurs, 1, numTests);
    evals_Ours(:, :) = repmat(stats.all_nevals(:, OurCol), 1, numTests);
    
    % Plot actual figure
    figure('DefaultAxesFontSize',14);
    hold on;

%     plot(frontsPerTest, avg_evals_Schulz, 'LineWidth', 2);
%     plot(frontsPerTest, avg_evals_NSGA, 'LineWidth', 2);

%     minmax = @std;
    minmax = @(x) [max(x) - mean(x); mean(x) - min(x)];
    shadedErrorBar(frontsPerTest, evals_Schulz, {@mean,minmax}, 'lineprops', {'LineWidth', 2})
    

    %% temporary
    notDoneYet = 5:5:201;
    remainingNSGA = [];
    for i=1:length(notDoneYet)
        numF = notDoneYet(i);
        NSGA_dir = "../../NSGA_comparisons/";
        NSGA_filename = NSGA_dir + problemName + "_numFronts" + numF + "_NSGA.mat";
        NSGA_log = load(NSGA_filename);
        n_evals_NSGA = NSGA_log.total_num_evals;
        remainingNSGA = [remainingNSGA, n_evals_NSGA];
        
        evals_NSGA(:, i) = NSGA_log.total_num_evals;
    end
%     plot([frontsPerTest, notDoneYet], [avg_evals_NSGA, remainingNSGA], 'LineWidth', 2);
%     plot([frontsPerTest, notDoneYet], [avg_evals_Ours, repmat(problemEvalsOurs, 1, length(notDoneYet))], 'r--', 'LineWidth', 2)
    % ========= end temporary 
    
    %%
    shadedErrorBar(frontsPerTest, evals_NSGA, {@mean,minmax}, 'lineprops', {'LineWidth', 2})

    scatter(200, problemEvalsOurs, 80, 'r', 'filled');
%     plot(frontsPerTest, avg_evals_Ours, 'r--', 'LineWidth', 2);
    
    shadedErrorBar(frontsPerTest, evals_Ours, {@mean,minmax}, 'lineprops', {'r--', 'LineWidth', 2})

    
    xlabel('Number (N) of Fixed-Context Fronts Computed')
    ylabel('Number of Function Evals Required')
    legend('Schulz: N Fixed-Contexts', 'NSGA: N Fixed-Contexts', 'Ours: Pareto Gamut (w/ 200 context cells)', 'Ours: Reference Line')
    hold off;
end