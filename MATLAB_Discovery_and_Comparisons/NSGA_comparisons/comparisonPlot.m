times = [];
evals = [];

problemName = "Turbine";
problemTimeOurs = 52.1;
problemEvalsOurs = 124713;

frontsPerTest = 5:5:195;
numtests = length(frontsPerTest);

for i=frontsPerTest
    filename = "./" + problemName + "_numFronts" + i + "_run0.mat";
    stats=load(filename);
    
    times = [times, stats.totalTime];
    evals = [evals, stats.total_num_evals];
end

figure;
hold on;
plot(frontsPerTest, times, 'LineWidth', 2);
plot(frontsPerTest, repmat(problemTimeOurs, numtests), 'r--', 'LineWidth', 2);
scatter(200, problemTimeOurs, 'r', 'filled');
xlabel('Number of Fixed-Context Fronts Computed')
ylabel('Time (s)')
legend('Fixed-Context NSGA', 'Pareto Gamut (w/ 200 context cells)')
hold off;

figure;
hold on;
plot(frontsPerTest, evals, 'LineWidth', 2);
plot(frontsPerTest, repmat(problemEvalsOurs, numtests), 'r--', 'LineWidth', 2);
scatter(200, problemEvalsOurs, 'r', 'filled');
xlabel('Number of Fixed-Context Fronts Computed')
ylabel('Number of Function Evals Required (s)')
% set(gca, 'YScale', 'log')
legend('Fixed-Context NSGA', 'Pareto Gamut (w/ 200 context cells)')
hold off;