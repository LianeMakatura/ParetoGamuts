function [directions, numActiveC] = getExplorationDirections_xz(funcsCell, ...
                            xSymb, xVals, zSymb, zVals, ...
                            alpha, beta, ...
                            rD, rd, rA)
    % note that beta from the optimize function includes *all* beta, not
    % just the active ones. Non-active are 0 as expected. Active are
    % non-zero.
                        

    % determine the active constraints
    % KKT only informs movement wrt x
    % z constraints are used in the full jacobian, DG, which comes in from
    % complementary slackness
    [activeC_x, activeC_z, upper_x, upper_z] = getActiveConstraints(xVals, zVals);
    numActiveC_x = length(activeC_x);
    numActiveC_z = length(activeC_z);
    numActiveC = numActiveC_x + numActiveC_z;
    constrained = false;
    if numActiveC > 0
        constrained = true;
        disp('CONSTRAINED = TRUE');
        disp([xVals, zVals]);
    end    
              
    % ==== compute DxF^T at (x*, z*) ====
    DxF(xSymb, zSymb) = jacobian(funcsCell, xSymb);
    DxF_evald = double(evalAtPt(DxF, xVals, zVals));
    
    % ==== Constrained case: compute DxG ====
    if constrained
        DxG = constraintJacobian_x(activeC_x, upper_x, numActiveC_x, numActiveC_z, rD);
        Dx = [DxF_evald', DxG'];
        % for now, don't need to evaluate at a point, because it's just the
        % box constraints on x
    else
        Dx = DxF_evald';
    end

    % computing the alpha, beta in optimize now
%     [alpha, beta] = getKKTDualVars(Dx, constrained, rd, numActiveC_x);


    % ==== compute quasi-Hessian terms ====
    H_x = zeros(rD, rD); % for H_x = sum( alpha_i * d^2/dx^2(f_i) )
    H_z = zeros(rD, rA); % for H_z = sum( alpha_i * d/dz(d/dx(f_i)) )

    for i=1:rd
        % construct the quasi-hessian d^2/dx^2(f_i) = J_x(grad_x(f_i))^T
        % then evaluate, scale by alpha, and add to H_x
        grad_x_fi = gradient(funcsCell{i}, xSymb);
        
        Hxx_fi(xSymb, zSymb) = jacobian(grad_x_fi, xSymb);
        Hxx_fi_evald = double(evalAtPt(Hxx_fi, xVals, zVals));
        H_x = H_x + alpha(i)*Hxx_fi_evald;

        % construct the "hessian" d/dz(d/dx (F)) = J_z(grad_x(f_i))^T
        % then evaluate, scale by alpha, and add to H_z
        Hxz_fi(xSymb, zSymb) = jacobian(grad_x_fi, zSymb);
        Hxz_fi_evald = double(evalAtPt(Hxz_fi, xVals, zVals));
        H_z = H_z + alpha(i)*Hxz_fi_evald;
    end
    
    if constrained
        % don't need to do this with box constraints only, because the values
        % are always 0; adding these matrices doesn't change anything.
        % also this code will break, since there are no symbolic
        % constraints.
%         for i=1:numActiveC_x
%             % construct the quasi-hessian d^2/dx^2(g_j) = J_x(grad_x(g_j))^T
%             % then evaluate, scale by beta, and add to H_x
%             grad_x_gi = gradient(constrCell{i}, xSymb);
% 
%             Hxx_gi(xSymb, zSymb) = jacobian(grad_x_gi, xSymb);
%             Hxx_gi_evald = double(evalAtPt(Hxx_gi, xVals, zVals));
%             H_x = H_x + beta(i)*Hxx_gi_evald;
% 
%             % construct the "hessian" d/dz(d/dx g_j) = J_z(grad_x(g_j))^T
%             % then evaluate, scale by beta, and add to H_z
%             Hxz_gi(xSymb, zSymb) = jacobian(grad_x_gi, zSymb);
%             Hxz_gi_evald = double(evalAtPt(Hxz_gi, xVals, zVals));
%             H_z = H_z + beta(i)*Hxz_gi_evald;
%         end
    end

    % ==== find directions: stack into DxHxHz, find null space ====
    % recall: #cols in DxHxHz = (rd + numActiveC_x + numActiveC_z + rD + rA)
    alpha_constraint = [ones(1, rd), zeros(1, numActiveC_x + numActiveC_z + rD + rA)]; % enforce sum [alpha] = 1;
    
    if constrained
        % add in complementary slackness condition
        DzG = constraintJacobian_z(activeC_z, upper_z, numActiveC_x, numActiveC_z, rA);
        DG = [DxG, DzG];     % full Jacobian of constraints
        comp_slack_constraint = [zeros(numActiveC, rd + numActiveC_x + numActiveC_z), DG];
        
        DxHxHz = [alpha_constraint;
                  comp_slack_constraint;
                  Dx, H_x, H_z];
    else
        DxHxHz = [alpha_constraint;
                    Dx, H_x, H_z];
    end
    
    directions = null(DxHxHz);
    if size(directions, 2) > 2
        disp('too many')
%         keyboard;
    end
        
end


%% Compute the KKT dual variables 
function [alpha, beta]= getKKTDualVars(Dx, constrained, rd, numActiveC_x)
    % ==== compute the KKT dual variables, alpha*, beta*, at (x*, z*) ====
    if constrained
        alpha_beta = null(Dx);
        alpha_beta = alpha_beta(:, 1); %only take one column, in case there are multiple. 
        if alpha_beta(1) < 0
            alpha_beta = -alpha_beta; % flip the vector (legal bc anything parallel still in subspace)
        end
        alpha = alpha_beta(1:rd); % first rd entries
        alpha_norm = sum(alpha);
        alpha = alpha / alpha_norm;
        beta = alpha_beta(rd+1:rd+numActiveC_x) / alpha_norm; % last K' entries
    else
        % in null(J^T) by KKT conditions
        alpha = null(Dx);
        if alpha(1) < 0
            alpha = -alpha; % flip the vector (legal bc anything parallel still in subspace)
        end
        alpha = alpha / sum(alpha);
    end
    if ~isempty(find(alpha < 0))
        fprintf('alpha vector not positive!');
        keyboard;
    elseif abs(sum(alpha) - 1) > 1e-8
        fprintf('alpha doesnt sum to 1');
        keyboard;
    end
end


%% Report the active constraints on a point
%
% @param xVals: design point to check constraints for
% @param zVals: application values to check constraints for
% ---
% @return activeConstraints: array of indices in P with active constraints
function [activeC_x, activeC_z, upper_x, upper_z] = getActiveConstraints(xVals, zVals)
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
function DxG = constraintJacobian_x(activeC_x, upper_x, numActiveConstraints_x, numActiveConstraints_z, rD)
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
function DzG = constraintJacobian_z(activeC_z, upper_z, numActiveC_x, numActiveC_z, rA)
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

