function [ out ] = Schulz_Kursawe_1_helper( x, z )
% x=x(1:3,:);

    [X,Y] = size(x);
    x(:,:) = 10*x(:,:) - 5; % rescale design vars to [0,1]
    rD = size(x,1);

    if exist('z', 'var') %% we should fix z (overwrite x(rD))
        x(rD, :) = z;
    end

    summation = zeros(1,Y);
    
    for i=1:X-1
        xi = x(i,:);
        xj = x(i+1,:);
        rad = sqrt(xi.^2 + xj.^2);
        p = exp(-0.2*rad);
        
        summation = summation + p;
    end
   
    out = -10.0 * summation;
    out = out./ 16.0;
    out = out + (5.0/4.0);
end

