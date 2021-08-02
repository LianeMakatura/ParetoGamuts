%% Pareto Gamut
pts_origin = pts;
val_idx = find(pts_origin(:, 34) <= 1);
pts = pts_origin(:, :);
pts(:, 34) = pts(:, 34) / 1.5;
idx_1 = find(abs(pts(:, 33) - 0.2575) < 1 / 200);
idx_2 = find(abs(pts(:, 33) - 0.5325) < 1 / 200);
idx_3 = find(abs(pts(:, 33) - 0.8325) < 1 / 200);
figure
color_base = jet(fix((max(pts(:,33)) - min(pts(:,33))) / 0.005) + 1);
color = color_base(int64((pts(:,33) - min(pts(:,33))) / 0.005 + 1), :);
scatter3(pts(:, 34), pts(:, 35), pts(:,33), 80, color)
hold on
grid off
scatter3(pts(idx_1, 34), pts(idx_1, 35), pts(idx_1,33), 200, [0,0,0], 'filled')
scatter3(pts(idx_2, 34), pts(idx_2, 35), pts(idx_2,33), 200, [0,0,0], 'filled')
scatter3(pts(idx_3, 34), pts(idx_3, 35), pts(idx_3,33), 200, [0,0,0], 'filled')
xlabel('Distance to goal','FontSize', 22);
ylabel('Energy','FontSize', 22);
zlabel('Length','FontSize', 22);
set(gcf,'Position',[680,373,722,598]);
set(gca, 'FontSize',18);

%% Design Space
figure
color_base = jet(fix((max(pts(:,33)) - min(pts(:,33))) / 0.005) + 1);
color = color_base(int64((pts(:,33) - min(pts(:,33))) / 0.005 + 1), :);
% color_base = jet(fix((max(pts(:,33)) - min(pts(:,33))) / 0.005) + 1);
design_1 = mean((1-pts(:, 1:2:32)), 2);
design_2 = mean((1-pts(:, 2:2:32)), 2);
% design_1 = (design_1-0.6)*20 + 1;
% design_2 = (design_2-0.6)*20 + 1;
scatter3(design_1, design_2, pts(:,33), 80, color)
hold on
grid off
scatter3(design_1(idx_1), design_2(idx_1), pts(idx_1,33), 160, [0,0,0], 'filled')
scatter3(design_1(idx_2), design_2(idx_2), pts(idx_2,33), 160, [0,0,0], 'filled')
scatter3(design_1(idx_3), design_2(idx_3), pts(idx_3,33), 160, [0,0,0], 'filled')
set(gca,'FontSize',18);
% zlim([0,1.0]);
xlabel('Actuator 1','FontSize', 22);
ylabel('Actuator 2','FontSize', 22);
zlabel('Length','FontSize', 22);
set(gcf,'Position',[680,373,722,598]);

%% Design slice
figure
% hold on
% idx = idx_3;
% color = jet(length(idx));
% pts_norm = (1-pts - (1-0.38)) * 15 + 1.15;
pts_norm = (1 - pts - (1-0.38)) * 7 + 0.55;

idx_1_1 = find(abs(pts(idx_1, 34) - 4e-4) < 1 / 200);
idx_2_2 = find(abs(pts(idx_2, 34) - 4e-4) < 1 / 200);
idx_3_3 = find(abs(pts(idx_3, 34) - 4e-4) < 1 / 200);
pts_1 = pts(idx_1, :);
pts_2 = pts(idx_2, :);
pts_3 = pts(idx_3, :);
idx_1_1_max = find(pts_1(idx_1_1, 35) == max(pts_1(idx_1_1, 35)));
% idx_2_2_mean = find(abs(pts_2(idx_2_2, 35) - (max(pts_2(idx_2_2, 35)) * 0.99 + min(pts_2(idx_2_2, 35)) * 0.01)) < 1 / 10);
% idx_2_2_mean = idx_2_2_mean(1);
idx_2_2_mean = find(pts_2(idx_2_2, 35) == max(pts_2(idx_2_2, 35)));
idx_3_3_min = find(pts_3(idx_3_3, 35) == max(pts_3(idx_3_3, 35)));
pts_1_1 = pts_1(idx_1_1, :);
pts_2_2 = pts_2(idx_2_2, :);
pts_3_3 = pts_3(idx_3_3, :);
xlim([-0.5,15.5]);
ylim([0,1]);
hold on
plot(16, 16, 'b--o', 'Color', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth',2);
plot(16, 17, 'c*-', 'Color', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth',2);
grid on;
plot([0:15], pts_norm(idx_1(idx_1_1(idx_1_1_max)), 1:2:32), 'b--o', 'Color', color(idx_1(end),:), 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth',5);
plot([0:15], pts_norm(idx_1(idx_1_1(idx_1_1_max)), 2:2:32), 'c*-', 'Color', color(idx_1(end),:), 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth',5);
plot([0:15], pts_norm(idx_2(idx_2_2(idx_2_2_mean)), 1:2:32), 'b--o', 'Color', [128 / 255, 160/255, 53/255], 'MarkerEdgeColor', [128 / 255, 160/255, 53/255], 'MarkerSize', 15, 'LineWidth',5);
plot([0:15], pts_norm(idx_2(idx_2_2(idx_2_2_mean)), 2:2:32), 'c*-', 'Color', [128 / 255, 160/255, 53/255], 'MarkerEdgeColor', [128 / 255, 160/255, 53/255], 'MarkerSize', 15, 'LineWidth',5);
plot([0:15], pts_norm(idx_3(idx_3_3(idx_3_3_min)), 1:2:32), 'b--o', 'Color', color(idx_3(1),:), 'MarkerEdgeColor', color(idx_3(1),:), 'MarkerSize', 15, 'LineWidth',5);
plot([0:15], pts_norm(idx_3(idx_3_3(idx_3_3_min)), 2:2:32), 'c*-', 'Color', color(idx_3(1),:), 'MarkerEdgeColor', color(idx_3(1),:), 'MarkerSize', 15, 'LineWidth',5);
ylim(([1-0.455,1-0.39] -1 + 0.38) * 15 + 1.1);
legend({}, 'String', {'Actuator 1', 'Actuator 2'}, 'TextColor', 'black', 'FontSize', 20);
set(gca,'FontSize',18);
set(gca,'xtick',[0:15])
ylabel('Actuating Force','FontSize', 20);
xlabel('Time Step','FontSize', 20);

% legend('b--o', 'Actuator 1');
% legend('c*-', 'Actuator 2');

% scatter([0:15], pts_norm(idx_1(1), 2:2:32), 80, color(i,:),'filled');
% 
% idx = idx_3;
% for i = 1:length(idx)
%     plot([min(ac_1_norm(idx(i), :)), max(ac_1_norm(idx(i), :))], [design_2(idx(i)), design_2(idx(i))]);
%     plot([design_1(idx(i)), design_1(idx(i))], [min(ac_2_norm(idx(i), :)), max(ac_1_norm(idx(i), :))]);
% end

%% 
figure
% hold on
idx = idx_3;
% color = jet(length(idx));
% pts_norm = (1-pts - (1-0.38)) * 15 + 1.1;
pts_norm = (1 - pts - (1-0.38)) * 7 + 0.55;
% pts_norm = pts;
for i = 1:length(idx)
    ylim(([1-0.455,1-0.39] -1 + 0.38) * 7 + 0.55);
    zlim(([1-0.455,1-0.39] -1 + 0.38) * 7 + 0.55);
    xlim([0,15]);
    scatter3(0, pts_norm(idx(i), 1), pts_norm(idx(i), 2), 80, color(i, :), 'filled');
    hold on;
    grid on;
    set(gca,'xtick',[0:15])
end

for i = 1:50:101
    scatter3([0:15], pts_norm(idx(i), 1:2:32), pts_norm(idx(i), 2:2:32), 80, color(i, :), 'filled');
    plot3([0:15], pts_norm(idx(i), 1:2:32), pts_norm(idx(i), 2:2:32), '-o', 'color', color(i, :));
end
set(gca,'FontSize',12);
ylabel('Actuator 1','FontSize', 14);
zlabel('Actuator 2','FontSize', 14);
xlabel('Time Step','FontSize', 14);
view(-14.399999938687216,31.200000005911065);

%% Pareto slice
figure 
hold on
idx = idx_3;
% color = jet(length(idx));
xlim([0,1]);
ylim([0,0.4]);
grid on
% line([1e-1, 1e-1], [0, 0.4], 'Color', 'black', 'LineWidth',5);
for i = 1:length(idx)
    scatter(pts(idx(i), 34), pts(idx(i), 35), 80, color(idx(i), :));%[128 / 255, 160/255, 53/255]);%
end
set(gca,'FontSize',12);
xlabel('Distance to goal','FontSize', 14);
ylabel('Energy','FontSize', 14);
