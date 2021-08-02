function [hess] = fd_hessian(mFunc, wrt_vars1, wrt_vars2, pDes, pApp, eps, fd_type)
% computes a finite difference hessian
%
% @param mFunc -
% @param wrt_vars1 - 'des_only', 'app_only', 'both'
% @param wrt_vars2 - 'des_only', 'app_only', 'both'
% @param pDes - 
% @param pApp - 
% @param eps - 
% @param fd_type - 'forward', 'backward' or 'central'
%
% @return hess - cell array of d hessians (one for each metric) 
%                     each hessian is numVars1 x numVars2; numVarsx is
%               D if wrt_varsx='des_only', 
%               A if wrt_varsx = 'app_only' and 
%               D+A if wrt_varsx = 'both

    rd = mFunc.rd;
    hess = cell([rd, 1]);
    
    for i=1:rd
        hess{i} = hessian_singlef(mFunc, i, pDes, pApp, eps, wrt_vars1, wrt_vars2, fd_type);
    end
end

function [hess] = hessian_singlef(mFunc, funcID, cpDes, cpApp, eps, wrt_vars1, wrt_vars2, fd_type)
    [num_i, des_ei_all, app_ei_all] = getMasks(mFunc, wrt_vars1);
    [num_j, des_ej_all, app_ej_all] = getMasks(mFunc, wrt_vars2);

    hess = zeros(num_i, num_j);
    for i=1:num_i
        ei_des = des_ei_all(i, :);
        ei_app = app_ei_all(i, :);
        for j=1:num_j
            ej_des = des_ej_all(j, :);
            ej_app = app_ej_all(j, :);
            
            if strcmp(fd_type, 'central')
                dxidxj = centralSecondOrder(mFunc, funcID, cpDes, cpApp, ...
                                            ei_des, ei_app, ...
                                            ej_des, ej_app, eps);
            elseif strcmp(fd_type, 'forward')
                disp('Foward differences for FDHessian not implemented yet.')
                return
            elseif strcmp(fd_type, 'backward')
                disp('Backward differences for FDHessian not implemented yet.')
                return
            else
                disp('FD type not recognized.')
                return
            end
            
            hess(i, j) = dxidxj;
        end
    end
end

function [numVars, des_e, app_e] = getMasks(mFunc, wrt_vars)
    rD = mFunc.rD;
    rA = mFunc.rA;

    if strcmp(wrt_vars, 'des_only')
        numVars = rD;
        changes = eye(numVars);
        des_e = changes;              % change one design var at a time
        app_e = zeros([numVars, rA]);       % application always same (no perturbation)
        
    elseif strcmp(wrt_vars, 'both')
        numVars = rD + rA;
        changes = eye(numVars);
        des_e = changes(:, 1:rD);              % change one design var at a time
        app_e = changes(:, rD+1:rD+rA);       % change one app var at a time
    
    elseif strcmp(wrt_vars, 'app_only')
        numVars = rA;
        changes = eye(numVars);
        des_e = zeros([numVars, rD]);              % change one design var at a time
        app_e = changes;       % application always same (no perturbation)
    
    else
        disp('wrt_vars not recognized');
        return
    end
end


function d2_dxidxj = centralSecondOrder(mFunc, i, cpDes, cpApp, ei_des, ei_app, ej_des, ej_app, h) 
    t1 = mFunc.eval_ith_metric(cpDes + h*ei_des + h*ej_des, ...
                               cpApp + h*ei_app + h*ej_app, i);
    t2 = mFunc.eval_ith_metric(cpDes + h*ei_des - h*ej_des, ...
                               cpApp + h*ei_app - h*ej_app, i);
    t3 = mFunc.eval_ith_metric(cpDes - h*ei_des + h*ej_des, ...
                               cpApp - h*ei_app + h*ej_app, i);
    t4 = mFunc.eval_ith_metric(cpDes - h*ei_des - h*ej_des, ...
                               cpApp - h*ei_app - h*ej_app, i);
    d2_dxidxj = (t1 - t2 - t3 + t4) / (4*h^2);
end

