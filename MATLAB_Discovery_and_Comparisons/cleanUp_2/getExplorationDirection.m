%% Obtain the direction to locally expand a manifold in design space
%
% @p: The point around which to locally expand the manifold.
% @mFunc: MappingFunction object used to obtain Gradient/Hessian
%   information about the objective functions
% @params: Object containing various static parameters (currently only used
%   to obtain the dimensionality of the design/objective space)
% ---
% @return direction: D x n matrix with each column as a
%   direction in design space across which to expand the manifold, where
%   D is the dimensionality of the design space.
  function direction = getExplorationDirection(p, mFunc, params)
    rD = params.rD;
    rd = params.rd;
    constrainZ = params.constrain_z;
   
    [J,h] = mFunc.evalHessAndDeriv(p);         
    % Get indices for active constraints
    C = getActiveConstraints(p, constrainZ);
    numActiveConstraints = length(C);

    % CASE(1) General case when there are active constraints on the 
    % design variables.
    if numActiveConstraints > 0
        DG = zeros(rD, numActiveConstraints);
        DG_new = zeros(numActiveConstraints, rD);

        for l = 1:numActiveConstraints
            constraintIndex = C(l);
            constraint = zeros(rD,1);
            constraint(constraintIndex) = 1;
            DG(:, l) = constraint;
            DG_new(l, :) = constraint;
        end

        DG = DG';

        % Use SVD decomposition to ensure that the jacobian is rank
        % deficient by zeroing out the smallest singular value(s)
        [U, S, V] = svd(J); 
        S(rd,rd) = 0;
        J = U*S*V';
        alpha = null(J');

        J = [J; -DG];

        omega = null(null(J)'); % w = Im(DF') 

        H = zeros(rD);
        for i = 1:rd
            hessian_i = h(:,:,i);
            H = H + alpha(i,1)*hessian_i; % H = sum[alpha(i)*hessian(i)]
        end

        H = double(H); % Convert matrix of symbolic expressions to doubles
        
        % doing this because of a numerical precision bug
        [U, S, V] = svd(H);
        S_inv = zeros(size(S));
        pos = find(S > 0.0001);
        S_inv(pos) = 1./S(pos);
        H_inv = U*S_inv*V';

        cprime = H_inv*omega;

        [sX, sY] = size(cprime);
        finalOmega = zeros(sX,sY);

        for l = 1:length(cprime(1,:))
            vec = cprime(:,l);
            for r = 1:numActiveConstraints
                constraint = DG(r,:);
                proj = dot(vec,constraint);
                vec = vec - proj*constraint';
            end
            
            finalOmega(:,l) = vec;
        end   
        
        
        finalOmega = reduce(finalOmega);
        direction = orth(null(null(finalOmega')'));

    % CASE(2) Special case when there are no active constraints. We can
    % therefore assume that the matrix of the objective derivative
    % functions is degenerate.
    else
        % Use SVD decomposition to ensure that the jacobian is rank
        % deficient by zeroing out the smallest singular value(s)
        [U, S, V] = svd(J);
        S(rd,rd) = 0;
        J = U*S*V';
        omega = null(null(J)'); % w = Im(DF')
        alpha = null(J'); % scaling factors alpha for the hessians

        % Get the hessians for the objective functions

        % Matlab is weird about storing arrays of functions, so we
        % have to loop through a cell-array and compute each
        % hessian separately.
        H = zeros(rD);
        for i = 1:rd
            hessian_i = h(:,:,i);
            H = H + alpha(i)*hessian_i; % H = sum[alpha(i)*hessian(i)]
        end

        H = double(H); % Convert matrix of symbolic expressions to doubles 
        H_inv = pinv(H); % Get H^-1
        direction = H_inv*omega; % H_inverse * w 
        direction = orth([direction null(H)]);
    end
 end

%% Report the active constraints on a point
%
% @param P: point P to check constraints for
% ---
% @return activeConstraints: array of indices in P with active constraints
function activeConstraints = getActiveConstraints(P, constrainZ)
    epsilon = .001;
    rD = length(P);
    
    minBound = zeros(1,rD);
    maxBound = ones(1,rD);
    
    UC = maxBound - P < epsilon;
    LC = P - minBound < epsilon;
    C = UC + LC;
    
    if constrainZ
        C(rD) = 1; %% only works for single app var
    end

    activeConstraints = find(C);
end


%% Reduce the dimensionality of a matrix using SVD decomposition 
%
% @param A: Matrix A to be reduced
% ---
% @return A_reduced: Matrix A with insignificant singular values zeroed
% out.
function A_reduced =  reduce(A)
    [U, S, V] = svd(A);
    M = S(1,1);
    S(S/M < 0.0001) = 0;
    A_reduced = U*S*V';
end