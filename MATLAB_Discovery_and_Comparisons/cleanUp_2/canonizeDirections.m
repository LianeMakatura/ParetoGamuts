% Take in the [x', z'] directions and alter them
function [canonical_exp_dirs, degenerate] = canonizeDirections(exp_dirs_orig, rd, rD, rA)  
    eps = 1e-6;
    num_dirs = size(exp_dirs_orig, 2);
    
    % swap order of app / design vars, transpose whole matrix
    exp_dirs = [exp_dirs_orig(rD+1:rD+rA, :); exp_dirs_orig(1:rD, :)]';
    
    % order the rows in decreasing order of context vars, first column
    exp_dirs = sortrows(exp_dirs, 'descend');
    
    % do gauss elimination
    [exp_dirs, ~, success] = solveGauss(exp_dirs, zeros(num_dirs, 1));
    if ~success || any(abs(sum(exp_dirs, 2)) < eps) % at least one direction linearly dependent
        degenerate = true;
    else
        degenerate = false;
    end
    
    % normalize the non-context rows
    exp_dirs(rA+1:end, :) = exp_dirs(rA+1:end, :)./sqrt(sum(exp_dirs(rA+1:end, :).^2,2));
    
    %transpose and put context vars at end again
    canonical_exp_dirs = [exp_dirs(:, rA+1:rA+rD), exp_dirs(:, 1:rA)]';
end

function [A, b, success] = solveGauss(A,b)
    success = true;

    s = size(A, 1);
    for j = 1:(s-1)
        for i = s:-1:j+1
            if A(j,j) < eps
                success = false;
                break;
            end
            m = A(i,j)/A(j,j);
            A(i,:) = A(i,:) - m*A(j,:);
            b(i) = b(i) - m*b(j);
        end
        if success == false
            break
        end
    end  
end