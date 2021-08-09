function [directions, numActiveC] = appGetExplorationDirections_xz(mFunc, ...
                            xVals, zVals, alpha, beta, ...
                            rD, rd, rA, userParams)
    % note that beta from the optimize function includes *all* beta, not
    % just the active ones. Non-active are 0 as expected. Active are
    % non-zero.
                        

    % determine the active constraints
    % KKT only informs movement wrt x
    % z constraints are used in the full jacobian, DG, which comes in from
    % complementary slackness
    [activeBoxC_x, activeBoxC_z, upper_x, upper_z] = getActiveBoxConstraints(xVals, zVals);
    if mFunc.hasLinearConstraints
        [numActiveLinearC, DxG_linear, DzG_linear] = mFunc.p.getActiveLinearConstraints(xVals, zVals);
    else
        numActiveLinearC = 0; DxG_linear = []; DzG_linear = [];
    end
    if rA == 0
        DzG_linear = [];
    end
    numActiveBoxC_x = length(activeBoxC_x);
    numActiveBoxC_z = length(activeBoxC_z);
    numActiveC = numActiveBoxC_x + numActiveBoxC_z + numActiveLinearC;
    constrained = false;
    if numActiveC > 0
        constrained = true;
%         disp('CONSTRAINED = TRUE');
%         disp([xVals, zVals]);
    end    
              
    % ==== compute DxF^T at (x*, z*) ====
    DxF_evald = mFunc.evalJacobian('des_only', xVals, zVals, userParams.fd_eps, userParams.fd_type);
    
    % ==== Constrained case: compute DxG ====
    if constrained
        DxG_box = boxConstraintJacobian_x(activeBoxC_x, upper_x, numActiveBoxC_x, numActiveBoxC_z, rD);
        DxG = [DxG_box; DxG_linear];
        Dx = [DxF_evald', DxG'];
        % for now, don't need to evaluate at a point, because it's just the
        % box/linear constraints on x-- constant wrt x,z
    else
        Dx = DxF_evald';
    end

    % ==== compute quasi-Hessian terms ====
    H_x = zeros(rD, rD); % for H_x = sum( alpha_i * d^2/dx^2(f_i) )
    H_z = zeros(rD, rA); % for H_z = sum( alpha_i * d/dz(d/dx(f_i)) )

    % construct the quasi-hessian d^2/dx^2(f_i) = J_x(grad_x(f_i))^T
    % construct the "hessian" d/dz(d/dx (F)) = J_z(grad_x(f_i))^T
    Hxx = mFunc.evalSecondPartials('des_only', 'des_only', xVals, zVals, userParams.fd_eps, userParams.fd_type);
    Hxz = mFunc.evalSecondPartials('des_only', 'app_only', xVals, zVals, userParams.fd_eps, userParams.fd_type);
    
    for i=1:rd
        % then evaluate, scale by alpha, and add to Hx / Hz
        H_x = H_x + alpha(i)*Hxx{i};
        H_z = H_z + alpha(i)*Hxz{i};
    end
    
%     if constrained
%         % For box and linear constraints, hessians are also 0 -- no need to compute
%         % IMPT: eval G second partials not implemented in mFunc, just here for
%         % reference if extended to non-linear constraints
% 
%         % construct the quasi-hessian d^2/dx^2(g_j) = J_x(grad_x(g_j))^T
%         % construct the "hessian" d/dz(d/dx g_j) = J_z(grad_x(g_j))^T
%         % gives cell array containing the quasi-Hessian for each active
%         % constraint, evaluated at the center point (xVals,zVals)
%         Hxx = mFunc.evalGSecondPartials('des_only', 'des_only', xVals, zVals, userParams.fd_eps, userParams.fd_type);
%         Hxz = mFunc.evalGSecondPartials('des_only', 'app_only', xVals, zVals, userParams.fd_eps, userParams.fd_type);
% 
%         for i=1:numActiveC
%             % then evaluate, scale by beta, and add to H_x / H_z
%             H_x = H_x + beta(i)*Hxx{i};
%             H_z = H_z + beta(i)*Hxz{i};
%         end
%     end

    % ==== find directions: stack into DxHxHz, find null space ====
    % recall: #cols in DxHxHz = (rd + numActiveC_x + numActiveC_z + rD + rA)
    alpha_constraint = [ones(1, rd), zeros(1, numActiveC + rD + rA)]; % enforce sum [alpha] = 1;
    
    if constrained
        % add in complementary slackness condition
        DzG_box = boxConstraintJacobian_z(activeBoxC_z, upper_z, numActiveBoxC_x, numActiveBoxC_z, rA);
        DzG = [DzG_box; DzG_linear];
        DG = [DxG, DzG];     % full Jacobian of constraints
        comp_slack_constraint = [zeros(numActiveC, rd + numActiveC), DG];
        
        DxHxHz = [alpha_constraint;
                  comp_slack_constraint;
                  Dx, H_x, H_z];
    else
        DxHxHz = [alpha_constraint;
                    Dx, H_x, H_z];
    end
    
    directions = null(DxHxHz);
    
    % ==== Possible: if z constraints are active, enforce ==== 
    % can still expand in one direction... maybe we don't deal with this
    % here, but upstream.
    
%     directions(rd +numActiveC+1:rd+numActiveC+rD+rA, :) = ...
%         orth( directions(rd +numActiveC+1:rd+numActiveC+rD+rA, :) ); %orthonormalize for stability and even
%     sampling -- want to normalize the directions, not the alpha / beta /
%     directions... also, null returns orthonormal vectors so this was
%     redundant

end




%% Report the active box constraints on a point
%
% @param xVals: design point to check constraints for
% @param zVals: application values to check constraints for
% ---
% @return activeConstraints: array of indices in P with active constraints
function [activeC_x, activeC_z, upper_x, upper_z] = getActiveBoxConstraints(xVals, zVals)
    epsilon = .001;
    rD = length(xVals);
    rA = length(zVals);
    
    minBound_x = zeros(1,rD);
    maxBound_x = ones(1,rD);
    
    minBound_z = zeros(1, rA);
    maxBound_z = ones(1, rA);
    
    UC_x = maxBound_x - xVals < epsilon;
    LC_x = xVals - minBound_x < epsilon;
    C_x = UC_x + LC_x;
    
    UC_z = maxBound_z - zVals < epsilon;
    LC_z = zVals - minBound_z < epsilon;
    C_z = UC_z + LC_z;

    activeC_x = find(C_x);
    activeC_z = find(C_z);
    
    % for active values, whether upper or lower bound active
    upper_x = UC_x(activeC_x);
    upper_z = UC_z(activeC_z);
end


%% Construct the constraint jacobian wrt design variables (D_{x}G)
% currently only handles the box constraints for design vars
% 1*x'_{i}=0 for all constrained vars x_i
% differentiating constraints on z has no effect (no x, KKT says nothing)
function DxG = boxConstraintJacobian_x(activeC_x, upper_x, numActiveConstraints_x, numActiveConstraints_z, rD)
    Kp = numActiveConstraints_x + numActiveConstraints_z; % K prime in derivation
    DxG = zeros(Kp, rD);

    % assumes that all x constraints come first, all z constraints later
    for l = 1:numActiveConstraints_x
        constrainedVarIndex = activeC_x(l);
        upperBoundActive = upper_x(l);
        constraint = zeros(1,rD);
        if upperBoundActive %(g(x) = x - c <= 0) --> d/dx(g(x)) = 1 
            constraint(constrainedVarIndex) = 1;
        else %lower bound active; (g(x) = c-x <= 0) --> d/dx(g(x)) = -1 
            constraint(constrainedVarIndex) = -1;
        end
        DxG(l, :) = constraint;
    end
end

%% Construct the constraint jacobian wrt application variables (D_{z}G)
function DzG = boxConstraintJacobian_z(activeC_z, upper_z, numActiveC_x, numActiveC_z, rA)
    Kp = numActiveC_x + numActiveC_z; % K prime in derivation
    DzG = zeros(Kp, rA);

    % assumes that all x constraints come first, all z constraints later
    % offset z indices by numActive_x
    for l = 1:numActiveC_z
        constrainedVarIndex = activeC_z(l);
        upperBoundActive = upper_z(l);

        constraint = zeros(1,rA);
        if upperBoundActive %(g(x) = z - c <= 0) --> d/dz(g(z)) = 1 
            constraint(constrainedVarIndex) = 1;
        else %lower bound active; (g(x) = c-z <= 0) --> d/dz(g(z)) = -1 
            constraint(constrainedVarIndex) = -1;
        end
        constraint(constrainedVarIndex) = 1;
        DzG(l + numActiveC_x, :) = constraint;
    end
end

