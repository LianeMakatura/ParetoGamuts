function [fh, gradh, hessh] = bicopter_setup(design_variables, app_variables, performance_type)
    % design_variables: [32, 1], x * y * 16
    % app_variables: context (length), [1, 1]
    % some hard-code hyper-params
    position_goal = [1.; 0.; 0.; 0.; 0.; 0.];
    mass_density = 1.0;
    dt = 0.25;
    q_init = zeros(6, 1);
   
    % Scale to proper range
    design_variables_scale = lerp(design_variables, -5.0, 5.0);
    app_variables_scale = lerp(app_variables, 0.5, 1.0); % 10.5 (m=r, false)
    func1 = @(q_final)(norm(q_final - position_goal)^2); %* 2.5e-6 (m=r, false); %* 2.5e-4(m=0.5); %
    func2 = @(u)(0.5 * norm(u)^2) * 0.7; % / 20;% * 5e-2;
    inertia_func = @(m, r)(1/12 * m * (2 * r)^2);
    
    total_time_steps = fix(size(design_variables_scale, 1) / 2);
    u = reshape(design_variables_scale, 2, total_time_steps);
    m = app_variables_scale * mass_density; %; 
    r = app_variables_scale;
    
    % f1: position performance
    % f2: energy performance
    [f, f1, f2] = quadcopter_dynamics(q_init, u, dt, m, r, 0.1, func2, 0.5, func1, inertia_func, false);
    
    % TODO: is simplify worth it?
    if strcmp(performance_type, 'position')
        f_symp = simplify(f1); 
    elseif strcmp(performance_type, 'energy')
        f_symp = simplify(f2);
    elseif strcmp(performance_type, 'both')
        f_symp = simplify(f);
    else
        error('bicopter_setup only supports position & energy');
    end
    
    all_vars = [design_variables; app_variables];
    all_vars_cell = {design_variables; app_variables};
    
    fh = matlabFunction(f_symp, 'Vars', all_vars_cell);
%     df = jacobian(f_symp, design_variables);
%     gradh = matlabFunction(df, 'Vars', all_vars_cell); % 1 x 32
    
    df_xz = jacobian(f_symp, all_vars);
    grad_xz = matlabFunction(df_xz, 'Vars', all_vars_cell); % 1 x 33
    gradh = grad_xz;
    % Test1: matlab built-in hessian(); GOODBYE
%     ddf = hessian(f_symp, all_vars); % the time cost is unbearable
    % Test2: FD; GOODBYE
    % TODO validate the result
%     ddf = fd_for_hess(fh, design_variables, app_variables);
    % Test3: Andy's implementation
    % DONOT do this
%     ddf = fast_hess(fh, gradh, design_variables, app_variables);
%     hessh = matlabFunction(ddf, 'Vars', all_vars_cell);
    hessh = @(x, z)fast_hess(fh, grad_xz, x, z);
end

function [hess] = fd_for_hess(metric, design_variables, app_variables)
    % hard code the eps and fd type
    eps = 1e-3; % central fd type
    
    % Note here vars is a row-vector
    [~, rD] = size(design_variables);
    [~, rA] = size(app_variables);
    
    numVars = rD + rA;
    changes = eye(numVars);
    des_e = changes(:, 1:rD);           % change one design var at a time
    app_e = changes(:, rD+1:rD+rA);     % change one app var at a time
    
    % Work for both i and j
    num_i = numVars;
    num_j = numVars;
    des_ei_all = des_e;
    app_ei_all = app_e; 
    des_ej_all = des_e;
    app_ej_all = app_e; 
    
    hess = [];
    for i=1:num_i
        ei_des = des_ei_all(i, :);
        ei_app = app_ei_all(i, :);
        hess_row = [];
        for j=1:num_j
            i,j
            ej_des = des_ej_all(j, :);
            ej_app = app_ej_all(j, :);
            h = eps;
            t1 = metric(design_variables + h*ei_des + h*ej_des, app_variables + h*ei_app + h*ej_app);
            t2 = metric(design_variables + h*ei_des - h*ej_des, app_variables + h*ei_app - h*ej_app);
            t3 = metric(design_variables - h*ei_des + h*ej_des, app_variables - h*ei_app + h*ej_app);
            t4 = metric(design_variables - h*ei_des - h*ej_des, app_variables - h*ei_app - h*ej_app);
            dxidxj = (t1 - t2 - t3 + t4) / (4*h^2);
            hess_row = [hess_row dxidxj];
        end
        hess = [hess; hess_row];
    end
end
