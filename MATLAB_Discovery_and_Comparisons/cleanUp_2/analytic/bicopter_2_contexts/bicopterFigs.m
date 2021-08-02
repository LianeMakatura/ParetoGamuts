pts_file = 'bicopter_C2_Dec9_noLabels.json'; %'copter2context.json';
hv_file = 'bicopter_C2_Dec9_HVhistory.json'; %'hypervolumes.json';
c1_file = 'bicopter_C2_Dec9_selected_points_unswapped.json'; %'selected_points_unswapped.json';

% Problem information
rD = 32;
rA = 2;
rd = 2;
perf1_idx = rD+rA+1;
perf2_idx = rD+rA+2;
app1_idx = rD+1;
app2_idx = rD+2;


% Load the files
pts = readFromJSON(pts_file);
hv = readFromJSON(hv_file);
hv= hv(end, :);
c1_pts = readFromJSON(c1_file);

%% ======= Hypervolume visualization
% create color maps
c = parula(256);

% custom dark
n = 100;               %// number of colors
R = linspace(0,0.2422,n);
G = linspace(0,0.1504,n);  
B = linspace(0, 0.6603,n);
c2 = [R(:), G(:), B(:)];  %// create colormap

%original
figure
imshow(reshape(hv, 100, 100));
colormap(gca, c);
hcb = colorbar(gca);
caxis([min(hv), 1]);
xlabel('Length');
ylabel('Density');

colorTitleHandle = get(hcb,'Title');
titleString = 'FC HV';
set(colorTitleHandle ,'String',titleString);
set(gca,'FontSize',18)



% % ======= pareto gamut visualization
appToVis = app2_idx;
color_idx = app1_idx;
contextColorRange = [0,1];
monochrome = false;


figure;
plottingColors = pts(:, color_idx); % color by context value
% caxis(contextColorRange);
if monochrome
    colormap gray;
else
    colormap jet;
end
scatter3(pts(:, perf1_idx), pts(:,perf2_idx), pts(:,appToVis),...
    25, plottingColors);
grid off;
xlabel('Distance To Goal');
ylabel('Energy');
zlabel('Length');
hcb=colorbar(gca);
colorTitleHandle = get(hcb,'Title');
titleString = 'Density';
set(colorTitleHandle ,'String',titleString);
set(gca,'FontSize',14)

% ======= visualize select density
figure;
% plottingColors = [106/255, 189/255, 69/255];
plottingColors = c1_pts(:, appToVis);
if monochrome
    colormap gray;
else
    colormap jet;
end
scatter3(c1_pts(:, perf1_idx), c1_pts(:,perf2_idx), c1_pts(:,appToVis),...
    25, plottingColors);
grid off;
xlabel('Distance To Goal');
ylabel('Energy');
zlabel('Length');
hcb=colorbar(gca);
colorTitleHandle = get(hcb,'Title');
titleString = 'Density';
set(colorTitleHandle ,'String',titleString);
set(gca,'FontSize',14)

