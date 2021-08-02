function [x] = polar2cartesian(phi, r)
    if ~exist('r', 'var')
        r = 1;
    end
    
    n = length(phi) + 1; % cartesian coords require 1 more dim than polar
    x = zeros(1, n);
    
    sinphi = sin(phi); %precompute sin(phi_i) for all i
    for i=1:n-1
            x(i) = r*cos(phi(i))*prod(sinphi(1:i-1));
    end
    x(n) = r*prod(sinphi);
end