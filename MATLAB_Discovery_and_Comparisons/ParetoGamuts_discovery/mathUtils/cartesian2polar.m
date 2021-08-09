function [phi, r] = cartesian2polar(x)
    %CARTESIAN2POLAR Summary of this function goes here
    %   Detailed explanation goes here
    tol = 1e-8;
    use_acos = true;
    
    n = length(x);
    
    phi = zeros(1, n-1);
    for k=1:n-1
        phi(k) = phi_k(x, k, tol, use_acos);
    end
    
    r = norm(x, 2);
end

function [phik] = phi_k(x, k, tol, use_acos)
    n = length(x);
    abovek_allzeros = isempty( find( abs( x(k+1:n) )>tol ) ); 
    
    if abovek_allzeros 
        % all of x_{k+1} to x_n are 0; ambiguous, decide phi_k by x_k
        if abs(x(k)) < tol  % x_k = 0, let phi_k = 0
            phik = 0;
            return;
        elseif x(k) > 0      % x_k > 0, let phi_k = 0
            phik = 0;
            return;
        elseif x(k) < 0     % x_k < 0, let phi_k = pi (180 degrees)
            phik = pi; 
            return;
        end
        
    elseif k < n-1                   % no ambiguity
        num = x(k);
        if use_acos
            denom = norm(x(k:n), 2);
            phik = acos(num / denom);
            return;
        else
            denom = norm(x(k+1:n), 2);
            phik = acot(num / denom); 
            return;
        end
    elseif k == n-1
        if use_acos
            num = x(k);
            denom = norm(x(k:n), 2);
            phik = acos(num / denom);
            if num < 0
                phik = 2*pi - phik;
            end
            return;
        else
            num = x(k) + norm(x(k:n), 2);
            denom = x(n);
            phik = 2*acot(num / denom);
            return;
        end
    end
end