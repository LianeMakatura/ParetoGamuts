function [des_extents, perf_extents, principal_dirs, principal_extents] ...
                = getExtents(centerPt, alphaStar, betaStar, ...
                            alphaPrimes, betaPrimes, ...
                            explorationDirections, mFunc_uc, rD, rd, rA)
% get the extents that we trust 

    num_dirs = size(alphaPrimes, 2); % should be rd - 1

%     %% bound the extent of the design space we trust
%     % ----- all alpha_i >= 0 -------
%     for dir=1:num_dirs
%         alpha_extents = -alphaStar ./ alphaPrimes(:, dir); % extent at which each alpha crosses 0
%         des_extents(1, dir) = min(abs(alpha_extents));
%     end

    %% bound the extent of the design space we trust (maybe we need beta too?)
    % ----- all alpha_i >= 0 -------
    alphaBar = getSimplexCorners(rd); %corners, target points
    alpha0 = repmat(alphaStar, 1, rd);
    targ = alphaBar - alpha0;
    
    % solve for the coefficients c1..cd to reach target
    % know it's rank deficient so just use the first rd-1 equations for some
    % solution
    % and constrain z' to be 0 (only makes sense in single context)
    zprimes = explorationDirections(rD+1:end, :); 
    C = [alphaPrimes(1:rd-1, :); zprimes ]\ [targ(1:rd-1, :); zeros(rA, size(targ, 2))]; 
    
    if max(max(abs(alphaPrimes*C - targ))) > 1e-5
        disp('error in extent calculation')
    end
    
    % new design directions (linear comb. to corners; extent is 0-1)
    corner_exp_dirs = explorationDirections * C;
    disp('design points corresponding to extreme alpha points: ')
    disp(corner_exp_dirs)
    
    des_extents = [-1,1];


    %% bound the extent of the linearized performance space we trust
    funcsCell = {mFunc_uc.f1, mFunc_uc.f2};
    H = zeros(rD+rA);
    for i=1:rd
        % construct the quasi-hessian d^2/dx^2(f_i) = J_x(grad_x(f_i))^T
        % then evaluate, scale by alpha, and add to H_x
        Hxx_fi = jacobian(gradient(funcsCell{i}));
        Hxx_fi_evald = double(evalAtPt(Hxx_fi, centerPt(1:rD), centerPt(rD+1:rD+rA)));
        H = H + Hxx_fi_evald;
    end

    tol = 0.01;
    for dir=1:num_dirs
        d = explorationDirections(:,dir);
        t = sqrt( tol / (d' * H * d) );  
        t = min(5, t); % prevent giant expansion in linear regions
        perf_extents(1, dir) = t;
    end
    
    %% find min/max eigenvector of H, for principal directions
    [VL, S, VR] = svd(H);
    
    principal_dirs(:, 1) = VL(:, 1);
    principal_dirs(:, 2) = VL(:, 1);
    
    for dir=1:2
        d = principal_dirs(:,dir);
        t = sqrt( tol / (d' * H * d) );  
        t = min(5, t); % prevent giant expansion in linear regions
        principal_extents(1, dir) = t;
    end
    
end

function corners = getSimplexCorners(rd)
    corners = eye(rd);
end



% function [des_extents, perf_extents, principal_dirs, principal_extents] ...
%                 = getExtents(centerPt, alphaStar, betaStar, ...
%                             alphaPrimes, betaPrimes, ...
%                             explorationDirections, mFunc_uc, rD, rd, rA)    
% % get the extents that have positive alpha in simplex
% 
%     num_dirs = size(alphaPrimes, 2); % should be rd + rA - 1
% 
%     %% bound the extent of the design space we trust (maybe we need beta too?)
%     % ----- all alpha_i >= 0 -------
%     alphaBar = getSimplexCorners(rd); %corners, target points
%     alpha0 = repmat(alphaStar, 1, rd);
%     targ = alphaBar - alpha0;
%     
%     % solve for the coefficients c1..cd to reach target
%     % know it's singular so just use the first rd-1 equations for some
%     % solution
%     C = alphaPrimes(1:rd-1, :) \ targ(1:rd-1, :); 
%     
%     if max(max(abs(alphaPrimes*C - targ))) > 1e-5
%         disp('error in extent calculation')
%     end
%     
%     % new design directions (linear comb. to corners; extent is 0-1)
%     corner_exp_dirs = explorationDirections * C;
%     
%     des_extents = [-1,1];
% 
%     
%     %% bound the extent of the linearized performance space we trust
%     funcsCell = {mFunc_uc.f1, mFunc_uc.f2};
%     for i=1:rd
%         % construct the quasi-hessian d^2/dx^2(f_i) = J_x(grad_x(f_i))^T
%         % then evaluate, scale by alpha, and add to H_x
%         H_fi = hessian(funcsCell{i});%jacobian(gradient(funcsCell{i}));
%         H_fi_evald = double(evalAtPt(H_fi, centerPt(1:rD), centerPt(rD+1:rD+rA)));
%         H{i} = H_fi_evald;
%     end
% 
%     tol = 0.01;
%     max_extent = inf*ones(1, rd); % max extent in each direction
%     
%     %% find min/max eigenvector of H, for principal directions
%     for i=1:rd
%         subH = explorationDirections' * H{i} * explorationDirections; % restrict to exploration subspace
%         [VL, S, VR] = svd(subH); 
% 
%         principal_dirs(:, 1) = explorationDirections * VL(:, 1); % map back to global
%         principal_dirs(:, 2) = explorationDirections * VL(:, end);
% 
%         for dir=1:2
%             d = principal_dirs(:,dir);
%             t = sqrt( abs( tol / (d' * H{i} * d) ) );  % abs might be bad...
%             t = min(5, t); % prevent giant expansion in linear regions
%             principal_extents(1, dir) = t;
%         end
%     end
% end
% 
% 
% 
