function [f, f1, f2] = quadcopter_dynamics(q_init, u, dt, m, r, running_weight, ...
    running_func, final_weight, final_func,...
    inertia_func, visualize)
g = -9.81;
last_q = q_init;
I = inertia_func(m, r);

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
    
    if running_weight
        run_val = running_weight * running_func(u(:, i));
        f = f + run_val;
        f2 = f2 + run_val;
    end
    
    if visualize
        all_qs(:, i) = last_q;
    end
    
    u1 = u(1, i);
    u2 = u(2, i);
    
    vx = -(u1 + u2) * sin(theta) / m * dt;
    vy = (u1 + u2) * cos(theta) - m * g * dt;
    vtheta = r * (u1 - u2) / I * dt;
    
    last_q = last_q + [delta_x; vx; vy; vtheta];
    
    
    
end

f1 = 0;
if final_weight
    final_val = final_weight * final_func(last_q);
    f = f + final_val;
    f1 = final_val;
end

f1 = f1 / 1e3 / 5; % /2.5;
f2 = f2 / 10; % / 5;

if visualize
    figure
    all_qs(:, end) = last_q;
    pause();
    
    hold on
    xlim([-3.0, 3.0])
    ylim([-3.0, 3.0])
    for i = 1:1:n+1
        cla;
        x1 = all_qs(1, i) - r * cos(all_qs(3, i));
        x2 = all_qs(1, i) + r * cos(all_qs(3, i));
        y1 = all_qs(2, i) - r * sin(all_qs(3, i));
        y2 = all_qs(2, i) + r * sin(all_qs(3, i));
        
        h = plot([x1; x2], [y1; y2]);
        
        saveas(h, strcat('visualize_',num2str(i, '%02d'),'.jpg'))
        pause(0.5)
    end
end

end