idx_1 = find(pts(:,33)==0.0175);
idx_2 = find(pts(:,33)==0.2575);
idx_3 = find(pts(:,33)==0.5325);
color_base = jet(fix((pts(end,33) - pts(1,33)) / 0.005) + 1);
color = color_base(int64((pts(:,33) - pts(1,33)) / 0.005 + 1), :);
%% Pareto Gamut
figure
hold on
ylim([0.4,1]);
zlim([0,0.6]);

color_base = jet(fix((pts(end,33) - pts(1,33)) / 0.005) + 1);
color = color_base(int64((pts(:,33) - pts(1,33)) / 0.005 + 1), :);
sc = scatter3(pts(:, 34), pts(:, 35), pts(:,33), 80, color);
grid off
% xlabel('Distance to goal','FontSize', 22);
% ylabel('Energy','FontSize', 22);
% zlabel('Length','FontSize', 22);
set(gca, 'FontSize',18);
set(gcf,'Position',[680,373,722,598]);
view(42.551475222565813,34.648459890860877);
delete(sc);

%%
writerObj = VideoWriter('bicopter_video.mp4','MPEG-4');
writerObj.Quality = 100;
% writerObj.VideoCompressionMethod = 'H.264';
framerate = 4;
writerObj.FrameRate = framerate;
Frames = [];

figure('Position', [1,72,1792,935])
set(gcf,'color','w');
dt = 0.25;
pts_norm = (1-pts - (1-0.38)) * 15 + 1.15;
x = [0:15];
idx_1_actuator_1 =  pts_norm(idx_1(1), 1:2:32);
idx_1_actuator_2 =  pts_norm(idx_1(1), 2:2:32);
idx_2_actuator_1 =  pts_norm(idx_2(51), 1:2:32);
idx_2_actuator_2 =  pts_norm(idx_2(51), 2:2:32);
idx_3_actuator_1 =  pts_norm(idx_3(101), 1:2:32);
idx_3_actuator_2 =  pts_norm(idx_3(101), 2:2:32);

u1 =  pts(idx_1(1), 1:32);
u1 = u1 * 10 - 5;
r1 =  0.0175*0.5 + 0.5;

u2 =  pts(idx_2(51), 1:32);
u2 = u2 * 10 - 5;
r2 = 0.2575*0.5 + 0.5;

u3 =  pts(idx_3(101), 1:32);
u3 = u3 * 10 - 5;
r3 = 0.5325*0.5 + 0.5;
dt = 0.25;

all_qs_1 = dynamics(u1, dt, r1);
all_qs_2 = dynamics(u2, dt, r2);
all_qs_3 = dynamics(u3, dt, r3);

n = size(u1, 2) / 2;

sub1=subplot(2,3,1);
idx_1_curve_1 = animatedline('Marker', 'o', 'Color', color(idx_1(1),:), 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth',2);
idx_1_curve_2 = animatedline('Marker', '*', 'Color', color(idx_1(1),:), 'MarkerEdgeColor', 'b', 'MarkerSize', 15, 'LineWidth',2);
xlim([0,15])
ylim(([1-0.455,1-0.39] -1 + 0.38) * 15 + 1.1);
ylabel('Actuating Force','FontSize', 20);
xlabel('Time Step','FontSize', 20);
set(gca,'FontSize',18);
grid on;
legend({}, 'String', {'Actuator 1', 'Actuator 2'}, 'TextColor', 'black', 'FontSize', 20);

subplot(2,3,2);
idx_2_curve_1 = animatedline('Marker', 'o', 'Color', color(idx_2(51),:), 'MarkerEdgeColor', color(idx_2(51),:), 'MarkerSize', 15, 'LineWidth',2);
idx_2_curve_2 = animatedline('Marker', '*', 'Color', color(idx_2(51),:), 'MarkerEdgeColor', color(idx_2(51),:), 'MarkerSize', 15, 'LineWidth',2);
xlim([0,15])
ylim(([1-0.455,1-0.39] -1 + 0.38) * 15 + 1.1);
ylabel('Actuating Force','FontSize', 20);
xlabel('Time Step','FontSize', 20);
set(gca,'FontSize',18);
grid on;
legend({}, 'String', {'Actuator 1', 'Actuator 2'}, 'TextColor', 'black', 'FontSize', 20);

subplot(2,3,3);
idx_3_curve_1 = animatedline('Marker', 'o', 'Color', color(idx_3(101),:), 'MarkerEdgeColor', color(idx_3(101),:), 'MarkerSize', 15, 'LineWidth',2);
idx_3_curve_2 = animatedline('Marker', '*', 'Color', color(idx_3(101),:), 'MarkerEdgeColor', color(idx_3(101),:), 'MarkerSize', 15, 'LineWidth',2);
xlim([0,15])
ylim(([1-0.455,1-0.39] -1 + 0.38) * 15 + 1.1);
ylabel('Actuating Force','FontSize', 20);
xlabel('Time Step','FontSize', 20);
set(gca,'FontSize',18);
grid on;
legend({}, 'String', {'Actuator 1', 'Actuator 2'}, 'TextColor', 'black', 'FontSize', 20);

subplot4 = subplot(2,3,4);
s1_4 = scatter(0, 0, 180, 'o', 'b', 'filled');
hold on
ts1_4 = text(-0.8, 0.2, 'start','Color', 'b', 'FontSize',20);
s2_4 = scatter(1, 0, 180, '*', 'r');
ts2_4 = text(1.2, 0.2, 'destination','Color', 'r', 'FontSize',20);
ylim(subplot4, [-1.5, 2.0]);
xlim(subplot4, [-1.5, 2.0]);
set(gca,'FontSize',18);
grid on

subplot5 = subplot(2,3,5);
s1_5 = scatter(0, 0, 180, 'o', 'b', 'filled');
hold on
ts1_5 = text(-0.8, 0.2, 'start','Color', 'b', 'FontSize',20);
s2_5 = scatter(1, 0, 180, '*', 'r');
ts2_5 = text(1.2, 0.2, 'destination','Color', 'r', 'FontSize',20);
ylim(subplot5, [-1.5, 2.0]);
xlim(subplot5, [-1.5, 2.0]);
set(gca,'FontSize',18);
grid on

subplot6 = subplot(2,3,6);
s1_6 = scatter(0, 0, 180, 'o', 'b', 'filled');
hold on
ts1_6 = text(-0.8, 0.2, 'start','Color', 'b', 'FontSize',20);
s2_6 = scatter(1, 0, 180, '*', 'r');
ts2_6 = text(1.2, 0.2, 'destination','Color', 'r', 'FontSize',20);
ylim(subplot6, [-1.5, 2.0]);
xlim(subplot6, [-1.5, 2.0]);
set(gca,'FontSize',18);
grid on

pause(2)

for i = 1:2*framerate
    Frames = [Frames, getframe(gcf)];
end
    
delete([s1_4, ts1_4, s2_4, ts2_4, s1_5, ts1_5, s2_5, ts2_5, s1_6, ts1_6, s2_6, ts2_6]);

for i = 1:0.5*framerate
    Frames = [Frames, getframe(gcf)];
end

for k = 1:length(x)+1
    
    subplot(2,3,4);
    x1 = all_qs_1(1, k) - r1 * cos(all_qs_1(3, k));
    x2 = all_qs_1(1, k) + r1 * cos(all_qs_1(3, k));
    y1 = all_qs_1(2, k) - r1 * sin(all_qs_1(3, k));
    y2 = all_qs_1(2, k) + r1 * sin(all_qs_1(3, k));
    
    [h1_4, h2_4, h3_4, h4_4, h5_4] = copter_plot(x1,y1,x2,y2, r1, all_qs_1(3, k), [0.8824, 0.5843, 0.0863]);
    
    subplot(2,3,5);
    x1 = all_qs_2(1, k) - r2 * cos(all_qs_2(3, k));
    x2 = all_qs_2(1, k) + r2 * cos(all_qs_2(3, k));
    y1 = all_qs_2(2, k) - r2 * sin(all_qs_2(3, k));
    y2 = all_qs_2(2, k) + r2 * sin(all_qs_2(3, k));
    
    [h1_5, h2_5, h3_5, h4_5, h5_5] = copter_plot(x1,y1,x2,y2, r2, all_qs_2(3, k), [0.8824, 0.5843, 0.0863]);
    
    subplot(2,3,6);
    x1 = all_qs_3(1, k) - r3 * cos(all_qs_3(3, k));
    x2 = all_qs_3(1, k) + r3 * cos(all_qs_3(3, k));
    y1 = all_qs_3(2, k) - r3 * sin(all_qs_3(3, k));
    y2 = all_qs_3(2, k) + r3 * sin(all_qs_3(3, k));
    
    [h1_6, h2_6, h3_6, h4_6, h5_6] = copter_plot(x1,y1,x2,y2, r3, all_qs_3(3, k), [0.8824, 0.5843, 0.0863]);
    
    if k > 1
        addpoints(idx_1_curve_1, x(k-1), idx_1_actuator_1(k-1));
        addpoints(idx_1_curve_2, x(k-1), idx_1_actuator_2(k-1));

        addpoints(idx_2_curve_1, x(k-1), idx_2_actuator_1(k-1));
        addpoints(idx_2_curve_2, x(k-1), idx_2_actuator_2(k-1));
        addpoints(idx_3_curve_1, x(k-1), idx_3_actuator_1(k-1));
        addpoints(idx_3_curve_2, x(k-1), idx_3_actuator_2(k-1));
        drawnow
    end
%     saveas(gcf, strcat('fuck.png'))
    Frames = [Frames, getframe(gcf)];
    pause(dt)
    if k ~= length(x)+1
        delete([h1_4, h2_4, h3_4, h4_4, h5_4, h1_5, h2_5, h3_5, h4_5, h5_5, h1_6, h2_6, h3_6, h4_6, h5_6]);
    end

end
open(writerObj); 
writeVideo(writerObj, Frames);
close(writerObj);

%%

hold on
xlim([-3.0, 3.0])
ylim([-3.0, 3.0])
n = size(u1, 2) / 2;
for i = 1:1:n+1
    cla;
    x1 = all_qs(1, i) - r * cos(all_qs(3, i));
    x2 = all_qs(1, i) + r * cos(all_qs(3, i));
    y1 = all_qs(2, i) - r * sin(all_qs(3, i));
    y2 = all_qs(2, i) + r * sin(all_qs(3, i));

    [h1, h2, h3, h4, h5] = copter_plot(x1,y1,x2,y2, r, all_qs(3, i), color(11500, :));
    pause(0.25)
end

function [h1, h2, h3, h4, h5] = copter_plot(x1, y1, x2, y2, r,  theta, color)
    h1 = plot([x1; x2], [y1; y2], 'LineWidth',10, 'Color', [0.5,0.5,0.5]);
    x_ac_1 = x1 + r * 0.2 * cos(theta);
    y_ac_1 = y1 + r * 0.2 * sin(theta);
    x_ac_2 = x2 - r * 0.2 * cos(theta);
    y_ac_2 = y2 - r * 0.2 * sin(theta);
    
    bar = 0.3;
    x_ac_1_c = x_ac_1 + bar * sin(-1*theta);
    y_ac_1_c = y_ac_1 + bar * cos(-1*theta);
    x_ac_2_c = x_ac_2 + bar * sin(-1*theta);
    y_ac_2_c = y_ac_2 + bar * cos(-1*theta);
    ellipse_a = r * 0.2 * 2;
    ellipse_b = 0.2 * ellipse_a * 2;
    
    h2 = plot([x_ac_1;x_ac_1_c], [y_ac_1;y_ac_1_c], 'LineWidth',3, 'Color', [0.5,0.5,0.5]);
    h3 = plot([x_ac_2;x_ac_2_c], [y_ac_2;y_ac_2_c], 'LineWidth',3, 'Color', [0.5,0.5,0.5]);
    h4 = plot_ellipse(ellipse_a, ellipse_b, x_ac_1_c, y_ac_1_c, -1*theta, color);
    h5 = plot_ellipse(ellipse_a, ellipse_b, x_ac_2_c, y_ac_2_c, -1*theta, color);
    
end

function all_qs = dynamics(u, dt, r)
    u = reshape(u, 2, []);
    q_init = [0,0,0,0,0,0]';
    m = r*1.0;
    g = -9.81;
    
    inertia_func = @(m, r)(1/12 * m * (2 * r)^2);
    I = inertia_func(m,r);
    
    last_q = q_init;

    [~, n] = size(u);

    f = 0;
    f1 = 0;
    f2 = 0;

    all_qs = zeros(6, n + 1);

    for i = 1:1:n
        x = last_q(1);
        y = last_q(2);
        theta = last_q(3);
        old_vel = last_q(4:end);
        delta_x = old_vel * dt;
        all_qs(:, i) = last_q;
        u1 = u(1, i);
        u2 = u(2, i);
        vx = -(u1 + u2) * sin(theta) / m * dt;
        vy = (u1 + u2) * cos(theta) - m * g * dt;
        vtheta = r * (u1 - u2) / I * dt;
        last_q = last_q + [delta_x; vx; vy; vtheta];
    end
    all_qs(:,end) = last_q;
end
