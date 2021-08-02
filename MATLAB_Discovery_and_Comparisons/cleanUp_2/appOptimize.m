%% Locally optimize a set of points; uses Matlab's fmincon function
% @param samples: n x D+A+d list of n points
% @param targets: n x d list of n optimization targets for the samples
% ---
% @return: n x D+A+d list of n points after performing a local
% optimization on each using fmincon.
function [projectedSamples, alphas, betas] = appOptimize(samples, targets, dirs, mFunc, fixedParams, userParams) 
    eps = 1e-3;
    
    rD = fixedParams.rD;
    rA = fixedParams.rA;
    rd = fixedParams.rd;

    numSamples = size(samples, 1);
    
    minBound = fixedParams.minBound;
    maxBound = fixedParams.maxBound;
    projectedPoints = zeros(numSamples, rD);
    projectedPtsApp = zeros(numSamples, rA);
    
    alphas = zeros(numSamples, rd);    % weights on perf metrics
    betas = zeros(numSamples, rD);     % constraints (only box in design; may need box in app too)

    for  i=1:numSamples
        P = samples(i,:);
        
        P(1:rD+rA) = min(P(1:rD+rA), maxBound); % ensure design, app in range
        P(1:rD+rA) = max(P(1:rD+rA), minBound);
        
        pDes = P(1:rD);                 % only the design values
        pApp = P(rD+1:rD+rA);           % only the application values
        pPerf = P(rD+rA+1:rD+rA+rd);    % only perf

        target = targets(i,:);  % only in perf space (rd), with same app values
        dir = dirs(i, :);       % only in perf space (rd), with same app values

        if mFunc.useFD %|| mFunc.symFreeAnalytic
            gradOption = 'off';
        else
            gradOption = 'on';
        end
        
        
        options = optimoptions('fmincon', ...
            'GradObj',gradOption, ...
            'CheckGradients', false,...
            'MaxIterations', 3000,...
            'MaxFunctionEvaluations', 3000,...
            'HessianApproximation', {'lbfgs',10}, ...
            'Algorithm','interior-point', ...
            'Display','notify'); % notify

        fh = mFunc.getFunctionHandle(target, pApp); % scalarized opt function
        
        % inequality constraints
        if mFunc.hasLinearConstraints
            [A, b] = mFunc.p.matrixLinearConstraints();
        else
            A = []; b=[];
        end
        
        % equality constraints
        Aeq = []; % don't need constraint on z anymore because we won't let z change
        beq = [];
        
        [xopt,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(fh,pDes,A,b,Aeq,beq, ...
                    minBound(1:rD), maxBound(1:rD), [],options);
        mFunc.numScalarizationFuncEvals = mFunc.numScalarizationFuncEvals + output.funcCount;
        fopt = mFunc.eval(xopt, pApp);
        
%         figure;
%         hold on;
%             scatter(pPerf(1), pPerf(2), 'r'); % starting point
%             scatter(target(1), target(2), 'b'); % target
%             scatter(fopt(1), fopt(2), 'g'); % obtained
%             xlim([0,1]); ylim([0,1]); 
%             keyboard;
        
        while norm(fopt - target) < eps || min(fopt - target) < 0
            if min(fopt - target) < 0
                targworse = find((fopt - target) < 0);
                target(targworse) = fopt(targworse) - userParams.directionDelta; %0.1; % set the target values along offending axes to something better than the values fopt reached 
            else
                target = target + userParams.directionDelta*dir;
            end
            
            target = min(target, 1); % prevent the target from going out of bounds, else fopt will always be better than it on this axis

            disp('Optimize: target wasn''t utopian. Trying again')
            fh = mFunc.getFunctionHandle(target, pApp);
            xopt_tmp = xopt;
            [xopt,fval,exitflag,output,lambda,grad,hessian] = ...
                fmincon(fh,xopt_tmp,A,b,Aeq,beq, minBound(1:rD), maxBound(1:rD), [],options);
            mFunc.numScalarizationFuncEvals = mFunc.numScalarizationFuncEvals + output.funcCount;
            fopt = mFunc.eval(xopt, pApp);
            
%             scatter(target(1), target(2), 'b'); % target
%             scatter(fopt(1), fopt(2), 'g'); % obtained
%             xlim([0,1]); ylim([0,1]); 
%             keyboard;
        end
        
%         hold off;

        projectedPoints(i,:) = xopt;
        projectedPtsApp(i, :) = pApp;
        
        % get alpha values at xopt
        diff = fopt - target;
        norm_a = sum(diff);
        alphas(i, :) = diff / norm_a;   % KKT dual weights, known to be positive
        
        % obtain the beta values for active BOX constraints (given by lambda)
        tol = 5e-3;
        for j=1:rD  % loop over all possible box constrained vals
            maxActive = abs(lambda.upper(j)) > tol;
            minActive = abs(lambda.lower(j)) > tol;
            if maxActive && minActive 
                disp('claiming upper and lower constraint of same param both active');
                if abs(lambda.upper(j)) > abs(lambda.lower(j))
                    maxActive = true; 
                    minActive = false;
                else
                    minActive = true;
                    maxActive = false;
                end
            end
            if maxActive
                betas(i, j) = lambda.upper(j);
            elseif minActive
                betas(i, j) = lambda.lower(j);
            end
        end
        % normalize by same factor so alpha*f' + beta*g' = 0
        betas(i, :) = betas(i, :) / norm_a;    
    end

    projectedSamples = [projectedPoints, projectedPtsApp, mFunc.eval(projectedPoints, projectedPtsApp)];
end