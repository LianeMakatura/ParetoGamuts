function [J] = fd_jacobian(mFunc, wrt_vars, cpDes, cpApp, eps, fd_type)
% computes a finite difference jacobian
%
% @param mFunc - 
% @param wrt_vars - 
% @param cpDes - 
% @param cpApp - 
% @param eps - 
% @param fd_type -
%
% @return jac - the dx(numVars) jacobian, where numVars is 
%               D if wrt_vars='des_only', 
%               A if wrt_vars = 'app_only' and 
%               D+A if wrt_vars = 'both'

    if size(cpDes, 1) > 1 || size(cpApp, 1) > 1
        disp('Error using fd_jacobian; can only request 1 point at a time.');
        return 
    end
    rD = mFunc.rD;
    rA = mFunc.rA;
        
    if strcmp(wrt_vars, 'des_only')
        numVarsVarying = rD;
        changes = eye(numVarsVarying);
        desPtbnMask = changes;              % change one design var at a time
        appPtbnMask = zeros([numVarsVarying, rA]);       % application always same (no perturbation)
        
    elseif strcmp(wrt_vars, 'app_only')
        numVarsVarying = rA;
        changes = eye(numVarsVarying);
        desPtbnMask = zeros([numVarsVarying, rD]);              % design always same
        appPtbnMask = changes;       % change one app var at a time
        
    elseif strcmp(wrt_vars, 'both')
        numVarsVarying = rD+rA;
        changes = eye(numVarsVarying);
        desPtbnMask = changes(:, 1:rD);              % change one design var at a time
        appPtbnMask = changes(:, rD+1:rD+rA);       % change one app var at a time

    else
        disp('wrt_vars not recognized');
        return
    end
    
    % compute the jacobian
    J = fd_helper(mFunc, cpDes, cpApp, desPtbnMask, appPtbnMask, ...
                            eps, numVarsVarying, fd_type);
end


function jac = fd_helper(mFunc, cpDes, cpApp, desPtbnMask, appPtbnMask, ...
                            eps, numVarsVarying, fd_type)
    if strcmp(fd_type, 'forward')
        % ith row is orig pt w/ ptbn in ith design/app var
        perturbedDes = repmat(cpDes, [numVarsVarying, 1]) + (desPtbnMask * eps); 
        perturbedApp = repmat(cpApp, [numVarsVarying, 1]) + (appPtbnMask * eps);
        
        % compute the jacobian by finite differencing
        fperturbed = mFunc.eval(perturbedDes, perturbedApp)';  % eval all metrics at perturbed points, and transpose        
        fx0 = mFunc.eval(cpDes, cpApp)'; %eval original point, transpose so dx1
        orig = repmat(fx0, [1, numVarsVarying]); %dxD
        jac = fperturbed - orig;
        
    elseif strcmp(fd_type, 'backward')
        % ith row is orig pt w/ ptbn in ith design/app var
        perturbedDes = repmat(cpDes, [numVarsVarying, 1]) - (desPtbnMask * eps); 
        perturbedApp = repmat(cpApp, [numVarsVarying, 1]) - (appPtbnMask * eps);
        
        % compute the jacobian by finite differencing
        fperturbed = mFunc.eval(perturbedDes, perturbedApp)';  % eval all metrics at perturbed points, and transpose
        fx0 = mFunc.eval(cpDes, cpApp)'; %evaluate the original point, transpose so dx1
        orig = repmat(fx0, [1, numVarsVarying]);
        jac = orig - fperturbed;
        
    elseif strcmp (fd_type, 'central')
        halfeps = eps/2;
        
        % ith row is orig pt w/ ptbn in ith design/app var
        perturbedDesPlus = repmat(cpDes, [numVarsVarying, 1]) + (desPtbnMask * halfeps); 
        perturbedAppPlus = repmat(cpApp, [numVarsVarying, 1]) + (appPtbnMask * halfeps);
        
        perturbedDesMinus = repmat(cpDes, [numVarsVarying, 1]) - (desPtbnMask * halfeps); 
        perturbedAppMinus = repmat(cpApp, [numVarsVarying, 1]) - (appPtbnMask * halfeps);
        
        % compute the jacobian by finite differencing
        fperturbedPlus = mFunc.eval(perturbedDesPlus, perturbedAppPlus)';  % evaluate all metrics at perturbed points, and transpose
        fperturbedMinus = mFunc.eval(perturbedDesMinus, perturbedAppMinus)';  % evaluate all metrics at perturbed points, and transpose
        jac = fperturbedPlus - fperturbedMinus;
    else
        disp('FD type not recognized.')
        return
    end
    
    % normalize by the step size in all cases
    jac = jac / eps;
end
