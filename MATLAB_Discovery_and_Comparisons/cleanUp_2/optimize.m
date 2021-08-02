%% Locally optimize a set of points; uses Matlab's fmincon function
% @param samples: n x D+d list of n points
% @param targets: n x d list of n optimization targets for the samples
% ---
% @return: n x D+d list of n points after performing a local
% optimization on each using fmincon.
function [projectedSamples, alphas, betas] = optimize(samples, targets, dirs, mFunc, fixedParams, z) 
    rD = fixedParams.rD;
    rd = fixedParams.rd;

    minBound = fixedParams.minBound;
    maxBound = fixedParams.maxBound;
    projectedPoints = zeros(size(samples, 1), rD);
    alphas = zeros(size(samples, 1), rd);
    betas = zeros(size(samples, 1), rD);

    for  i=1:size(samples,1)
        P = samples(i,:);
        P(1:rD) = min(P(1:rD), maxBound);
        P(1:rD) = max(P(1:rD), minBound);

        p = P(1, 1:rD);  

        target = targets(i,:);
        dir = dirs(i, :);

        options = optimoptions('fmincon', ...
            'GradObj','on', ...
            'HessianApproximation', {'lbfgs',10}, ...
            'Algorithm','interior-point', ...
            'Display','off');

        fh = mFunc.getFunctionHandle(target);
        
        % Linear constraints if any
        if fixedParams.constrain_z
            Aeq = zeros(1, rD);
            Aeq(rD) = 1;
            beq = [z];
        else
            Aeq = [];
            beq = [];
        end
        
        [xopt,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(fh,p,[],[],Aeq,beq, minBound, maxBound, [],options);
        
        mFunc.numScalarizationFuncEvals = mFunc.numScalarizationFuncEvals + output.funcCount;
        
        fopt = mFunc.eval(xopt, z);
        while min(fopt - target) < 0 
            % update the target; wasn't utopian, so xopt not optimal
%             target = target + 0.2*dir;
            targworse = find((fopt - target) < 0);
            target(targworse) = fopt(targworse) - 0.1; % set the target values along offending axes to something better than the values fopt reached 

            target = min(target, 1); % prevent the target from going out of bounds, else fopt will always be better than it on this axis

            disp('Optimize: target wasn''t utopian. Trying again')
            
            fh = mFunc.getFunctionHandle(target);
            [xopt, ~, ~, output, ~, ~, ~] = fmincon(fh,p,[],[],Aeq,beq, minBound, maxBound, [],options);
            mFunc.numScalarizationFuncEvals = mFunc.numScalarizationFuncEvals + output.funcCount;
            fopt = mFunc.eval(xopt, z);
        end

        projectedPoints(i,:) = xopt;
        
        % get alpha values at xopt
        diff = fopt - target;
        norm_a = sum(diff);
        alphas(i, :) = diff / norm_a;   % KKT dual weights, known to be positive
        
        % obtain the beta values for active BOX constraints (given by lambda)
        tol = 1e-3;
        for j=1:rD  % loop over all possible box constrained vals
            maxActive = abs(lambda.upper(j)) > tol;
            minActive = abs(lambda.lower(j)) > tol;
            if maxActive && minActive 
                disp('claiming upper and lower constraint of same param both active');
            elseif maxActive
                betas(i, j) = lambda.upper(j);
            elseif minActive
                betas(i, j) = lambda.lower(j);
            end
        end
        % normalize by same factor so alpha*f' + beta*g' = 0
        betas(i, :) = betas(i, :) / norm_a;    
    end

    projectedSamples = [projectedPoints, mFunc.eval(projectedPoints, z)];
end