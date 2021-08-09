function plotTimeComparison(maxRuns, frontsPerTest, problemName)
    addpath('raacampbell-shadedErrorBar-263537a/')

    SchulzCol = 1;
    NSGACol = 2;
    OurCol = 3;

    numTests = length(frontsPerTest);
    avg_times_NSGA = zeros(1, numTests);
    avg_times_Schulz = zeros(1, numTests);
    avg_times_Ours = zeros(1, numTests);
    
    times_NSGA = zeros(maxRuns, numTests);
    times_Schulz = zeros(maxRuns, numTests);
    times_Ours = zeros(maxRuns, numTests);
    
    for i=1:length(frontsPerTest)
        numFCFronts = frontsPerTest(i);
        filename = sprintf('./%s_numFronts_%d.mat', problemName, numFCFronts);
        stats = load(filename);
        
        avg_times_NSGA(1, i) = stats.avg_time(NSGACol);
        times_NSGA(:, i) = stats.all_time(:, NSGACol);
        avg_times_Schulz(1, i) = stats.avg_time(SchulzCol);
        times_Schulz(:, i) = stats.all_time(:, SchulzCol);
    end
    
    % load our results
    filename = sprintf('./%s_ours.mat', problemName);
    stats = load(filename);
    problemTimeOurs = stats.avg_time(OurCol);
    avg_times_Ours(1, :) = repmat(problemTimeOurs, 1, numTests);
    
    % Plot actual figure
    figure;
    hold on;

    minmax = @(x) [max(x) - mean(x); mean(x) - min(x)];
%     plot(frontsPerTest, avg_times_Schulz, 'LineWidth', 2);
%     plot(frontsPerTest, avg_times_NSGA, 'LineWidth', 2);
    shadedErrorBar(frontsPerTest, times_Schulz, {@mean,minmax}, 'lineprops', {'LineWidth', 2})
    shadedErrorBar(frontsPerTest, times_NSGA, {@mean,minmax}, 'lineprops', {'LineWidth', 2})
%     shadedErrorBar(frontsPerTest, times_Ours, {@mean,minmax}, 'lineprops', {'LineWidth', 2})
    plot(frontsPerTest, avg_times_Ours, 'r--', 'LineWidth', 2);


%     %% temporary
%     notDoneYet = 90:5:201;
%     remainingNSGA = [];
%     for i=1:length(notDoneYet)
%         numF = notDoneYet(i);
%         NSGA_dir = "../../NSGA_comparisons/";
%         NSGA_filename = NSGA_dir + problemName + "_numFronts" + numF + "_NSGA.mat";
%         NSGA_log = load(NSGA_filename);
%         n_evals_NSGA = NSGA_log.totalTime;
%         remainingNSGA = [remainingNSGA, n_evals_NSGA];
%     end
%     plot([frontsPerTest, notDoneYet], [avg_times_NSGA, remainingNSGA], 'LineWidth', 2);
%     plot([frontsPerTest, notDoneYet], [avg_times_Ours, repmat(problemTimeOurs, 1, length(notDoneYet))], 'r--', 'LineWidth', 2)
    % ========= end temporary 

    scatter(200, problemTimeOurs, 'r', 'filled');
    xlabel('Number of Fixed-Context Fronts Computed')
    ylabel('Time (s)')
    legend('Schulz: N Fixed-Contexts', 'NSGA: N Fixed-Contexts', 'Ours: Pareto Gamut (w/ 200 context cells)')
    hold off;
end


% y=randn(256,80)*5; 
% x=(1:size(y,2));
% yP = cos( linspace(-2*pi,2*pi,length(x)) )*10;
% y = bsxfun(@plus,y,yP);
% 
% 
% shadedErrorBar(x, y, {@mean,@std}, 'lineprops', '-r')
% 
% hold on
% 
% y=mean(y)+16;
% errbar = [2*ones(1,length(x)) ; 4*ones(1,length(x))];
% 
% shadedErrorBar(x, y, errbar, 'lineprops', '-g')